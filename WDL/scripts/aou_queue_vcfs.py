''' This script collects individual sample VCFs for re-merging. We are using it to re-merge
a subset of the batches using new merging code that will pull additional fields from the 
individual sample VCFs into the merged VCFs.
'''

import pandas as pd
import argparse
import subprocess, sys
import re
import os
import gzip
import multiprocessing

from tqdm import tqdm

try:
    import gcsfs
except ImportError:
    subprocess.call([sys.executable, '-m', 'pip', 'install', 'gcsfs'])


def generate_regex(merged_vcf_path):
    this_search_pre = re.search('.+/(?=call-MergeMitoMultiSampleOutputsInternal)', merged_vcf_path)[0] 
    if len(this_search_pre) < 1:
        raise ValueError('ERROR: the regex search did not work.')
    this_search_str = this_search_pre + \
        'call-MitochondriaPipeline_v2_5/*/MitochondriaPipeline/*/call-LiftOverAfterSelf/**/out/*.self.ref.split.selfToRef.final.vcf'
    return this_search_str


def get_vcf_paths(this_search_str, merged_path, fs):
    run_out = subprocess.run(['gsutil', '-u', os.getenv("GOOGLE_PROJECT"), 'ls', this_search_str], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    vcf_list_pre = [x for x in run_out.stdout.split('\n') if len(x) > 0]
    sample_names_pre = [re.sub('.self.ref.split.selfToRef.final.vcf$', '', os.path.basename(x)) for x in vcf_list_pre]

    # edge case -- address cases where multiple attempts were successful by taking the latest one
    # example: shard-31, 7ca1cbb7-7c25-4683-825c-4a3522bba71c
    attempt_count = []
    for x in vcf_list_pre:
        search = re.search('(?<=call-LiftOverAfterSelf/)attempt-([0-9]{1,})/(?=out/)', x)
        this_att = int(search[1]) if search else 1
        attempt_count.append(this_att)
    df_attempts = pd.DataFrame({'sample': sample_names_pre, 'vcf': vcf_list_pre, 'attempt_number': attempt_count})
    df_attempts_filt = df_attempts.iloc[df_attempts.groupby(['sample'])['attempt_number'].idxmax()].sort_index(axis=0)
    vcf_list = list(df_attempts_filt['vcf'])
    sample_names = [re.sub('.self.ref.split.selfToRef.final.vcf$', '', os.path.basename(x)) for x in vcf_list]

    # confirm that all samples were identified
    merged_names = get_sample_names_from_vcf(fs, merged_path)
    if sorted(sample_names) != sorted(merged_names):
        raise ValueError(f'ERROR: for {merged_path}, sample names from file paths are not identical to sample names from the merged VCF.')
    
    return pd.DataFrame({'sample': sample_names, 'vcf': vcf_list})


def get_sample_names_from_vcf(fs, merged_path):
    this_line = ''
    with fs.open(merged_path, 'rb') as f:
        g = gzip.GzipFile(fileobj=f)
        for line in g:
            line = line.decode('utf-8')
            if re.search('^#CHROM', line):
                this_line = line
                break
    this_line = re.sub(r'\n','',this_line)
    return this_line.split('\t')[9:]


def generate_output_paths(df, output_dir):
    # Generates gs:// output paths for sample and VCF flat files
    df['sample_list_file'] = df["batch"].map(lambda x: f'{output_dir}{x}/sample_list.txt')
    df['vcf_list_file'] = df["batch"].map(lambda x: f'{output_dir}{x}/vcf_list.txt')
    df['write_success'] = df["batch"].map(lambda x: f'{output_dir}{x}/_SUCCESS')
    return df


def produce_lists(df, overwrite):
    for _, row in tqdm(df.iterrows()):
        list_writer_core(overwrite=overwrite, row=row)


def list_writer_core(overwrite, row):
    fs = gcsfs.GCSFileSystem(project=os.getenv('GOOGLE_PROJECT'))
    if (overwrite) or (not fs.exists(row['write_success'])):
        df_per_batch = get_vcf_paths(row['search_str'], row['vcf'], fs=fs)
        df_per_batch['sample'].to_csv(row['sample_list_file'], sep='\t', index=False, header=False)
        df_per_batch['vcf'].to_csv(row['vcf_list_file'], sep='\t', index=False, header=False)
        with fs.open(row['write_success'], 'w') as f:
            f.write('0\n')


parser = argparse.ArgumentParser()
parser.add_argument('--batch-paths', required=True, type=str, 
                    help='Path to tab_batch_file_paths.tsv, which is an output of "collate_tables.py."')
parser.add_argument('--table-output', required=True, type=str, 
                    help='Local path to output flat file containing a modified batch paths file now containing paths to per-batch sample and individual VCF flat files.')
parser.add_argument('--flat-file-output', required=True, type=str, 
                    help='gs:// path to a folder that will contain per-batch sample and individual VCF flat files.')
parser.add_argument('--overwrite', action='store_true',
                    help='If enabled, will overwrite per-batch files.')
parser.add_argument('--serial', action='store_true',
                    help='If enabled, will avoid parallel processing.')
args = parser.parse_args()


def _internal_list_writer(row):
    _, row_this = row
    return list_writer_core(overwrite=args.overwrite, row=row_this)


if __name__ == "__main__":
    df = pd.read_csv(args.batch_paths, sep='\t')
    print(f'Obtained file with {str(df.shape[0])} batches.')

    df['search_str'] = df.vcf.map(generate_regex)
    df = generate_output_paths(df, args.flat_file_output)

    print(f'Generating per-batch files...')
    if args.serial:
        _ = produce_lists(df=df, overwrite=args.overwrite)
    else:
        multiprocessing.set_start_method('spawn')
        num_cores = multiprocessing.cpu_count()
        p = multiprocessing.Pool(processes=num_cores)
        print(f'Using {str(num_cores)} cores...')
        results = tqdm(p.imap(_internal_list_writer, df.iterrows()), total=df.shape[0])
        tuple(results)

    print(f'Outputting local table with paths...')
    df.to_csv(args.table_output, sep='\t', index=False)
