import argparse
import hail as hl
import os
import pandas as pd

def get_covariates(stats_path, covar_path, dataset=os.getenv("WORKSPACE_CDR"), 
                   ancestry_pred_path="gs://fc-aou-datasets-controlled/v6/wgs/vcf/aux/ancestry/ancestry_preds.tsv"):
    ancestry_pred = hl.import_table(ancestry_pred_path,
                                    key="research_id", 
                                    impute=True, 
                                    types={"research_id":"tstr","pca_features":hl.tarray(hl.tfloat)},
                                    min_partitions=50)

    if hl.hadoop_exists(os.path.join(covar_path, '_SUCCESS')):
        sample_covariates = hl.read_table(covar_path)
    else:
        ht_samp_flat = hl.import_table(stats_path,
                                       key="s", 
                                       impute=True, 
                                       types={"s":"tstr"},
                                       min_partitions=50)
#        ht_samp_flat = ht_samp_flat.annotate(isFemale = hl.if_else(ht_samp_flat.sex_at_birth == 'Female', 1, 
#                                                        hl.if_else(ht_samp_flat.sex_at_birth == 'Male', 0, hl.missing(hl.tint32))))
        ht_samp_flat = ht_samp_flat.annotate(isFemale = hl.if_else(ht_samp_flat.sex_at_birth == 'F', 1, 
                                                        hl.if_else(ht_samp_flat.sex_at_birth == 'M', 0, hl.missing(hl.tint32))))

        person_sql = f"""
        SELECT  person.person_id,
                person.birth_datetime,
                p_location_concept.concept_name as loc,
                p_site_concept.concept_name as site
            FROM
                `{dataset}.person` person 
            LEFT JOIN
                `{dataset}.concept` p_location_concept 
                    on person.location_id = p_location_concept.CONCEPT_ID 
            LEFT JOIN
                `{dataset}.concept` p_site_concept 
                    on person.care_site_id = p_site_concept.CONCEPT_ID
            WHERE
                person.PERSON_ID IN (
                    select
                        person_id  
                    from
                        `{dataset}.cb_search_person` cb_search_person  
                    where
                        cb_search_person.person_id in (
                            select
                                person_id 
                            from
                                `{dataset}.cb_search_person` p 
                            where
                                has_whole_genome_variant = 1 
                        ) 
                    )"""

        wgs_demog = pd.read_gbq(person_sql, dialect="standard")

        sample_covariates = ht_samp_flat.select(isFemale = ht_samp_flat.isFemale, 
                                                mtdna_mean_coverage = ht_samp_flat.mean_coverage, 
                                                nucdna_mean_coverage = ht_samp_flat.nuc_mean_coverage)
        sample_covariates = sample_covariates.annotate(mtcn = 2 * sample_covariates.mtdna_mean_coverage / sample_covariates.nucdna_mean_coverage)
        wgs_demog['approx_age'] = 2021 - pd.DatetimeIndex(wgs_demog.birth_datetime).year
        age_table = wgs_demog[['person_id', 'approx_age']]
        age_ht = hl.Table.from_pandas(age_table)
        age_ht = age_ht.annotate(person_id = hl.str(age_ht.person_id)).key_by('person_id')
        sample_covariates = sample_covariates.annotate(approx_age = age_ht[sample_covariates.s].approx_age)
        sample_covariates = sample_covariates.annotate(age_isFemale = sample_covariates.approx_age * sample_covariates.isFemale,
                                                    age2 = sample_covariates.approx_age**2)
        sample_covariates = sample_covariates.annotate(age2_isFemale = sample_covariates.age2 * sample_covariates.isFemale)

        sample_covariates = sample_covariates.annotate(**{f'PC{str(idx+1)}': ancestry_pred[sample_covariates.s].pca_features[idx] for idx in range(0,10)})
        sample_covariates = sample_covariates.annotate(pop = ancestry_pred[sample_covariates.s].ancestry_pred)
        sample_covariates = sample_covariates.repartition(100).checkpoint(covar_path, overwrite=True)
    return sample_covariates

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('stats_path', help='GCS URI to the total merged stats file that is output by aou_collate_tables.py')
    parser.add_argument('gwas_dir', help='GCS URI to the location where the covariates.ht file should be written.')
    parser.add_argument('--export_tsv_path', help='Path to tsv file where the covariates should be exported to the local file system.')
    parser.add_argument('--ancestry_pred_path', default="gs://fc-aou-datasets-controlled/v6/wgs/vcf/aux/ancestry/ancestry_preds.tsv")
    args = parser.parse_args()

    covar_path = os.path.join(args.gwas_dir, 'covariates.ht')
    dataset = os.getenv('WORKSPACE_CDR')
    ht = get_covariates(args.stats_path, covar_path, dataset=dataset, ancestry_pred_path=args.ancestry_pred_path)
    ht = ht.select(sex = ht.isFemale,
                   age = ht.approx_age,
                   pop = ht.pop)
    if args.export_tsv_path:
    #    f'file:///home/jupyter/covariates_export.tsv'
        ht.export(f'file://{os.path.abspath(args.export_tsv_path)}')