version 1.0

workflow MergeVCFs {

  meta {
    description: "Takes in an array of VCFs and outputs 1 merged VCF."
    allowNestedInputs: true
  }

  input {
    File vcf_list
    File sample_name_list

    File MergePerBatch

    #Docker and version arguments
    String genomes_cloud_docker = "docker.io/rahulg603/genomes_cloud_bcftools"

    #Optional runtime arguments
    Int? preemptible_tries
  }

  Array[File] vcf = read_lines(vcf_list)
  Array[String] sample_name = read_lines(sample_name_list)

  call MergeVCFsInternal {
    input:
      sample_name = sample_name,
      variant_vcf = vcf,
      MergePerBatch = MergePerBatch,
      preemptible_tries = preemptible_tries,
      genomes_cloud_docker = genomes_cloud_docker
  }

  output {
    # merged files
    File merged_calls = MergeVCFsInternal.merged_calls
  }
}

task MergeVCFsInternal {
  # this task merges arrayed inputs for variant calls
  input {
    Array[String] sample_name
    Array[File] variant_vcf

    File MergePerBatch
    Int? preemptible_tries
    String genomes_cloud_docker
  }

  Int disk_size = ceil(size(variant_vcf, "GB") * 2) + 20

  command <<<
    set -e

    # now produce the relevant inputs for the per-batch MT script
    R --vanilla <<CODE
      sample_ids <- read.csv("~{write_lines(sample_name)}", sep='\t', stringsAsFactors=F, header=F)[[1]]
      paths_vcf <- read.csv("~{write_lines(variant_vcf)}", sep='\t', stringsAsFactors=F, header=F)[[1]]

      write.table(data.frame(s=sample_ids, path=paths_vcf), sep='\t', row.names=F, file='vcf_paths.tsv', quote=F)
    CODE

    mkdir tmp

    # Merge variant calls
    python3.7 ~{MergePerBatch} \
    --run-variants \
    --input-tsv vcf_paths.tsv \
    --temp-dir tmp/ \
    --output-flat-file batch_merged_mt_calls.vcf.bgz
  >>>

  output {
    File merged_calls = "batch_merged_mt_calls.vcf.bgz"
  }

  runtime {
    disks: "local-disk " + disk_size + " HDD"
    memory: "2 GB"
    docker: genomes_cloud_docker
    preemptible: select_first([preemptible_tries, 5])
  }
}
