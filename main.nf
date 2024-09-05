#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

// Print help message, supply typical command line usage for the pipeline

log.info """
    |            #################################################
    |            #    _  _                             _         #
    |            #   | \\| |  ___  __ __  ___   _ __   (_)  __    #
    |            #   | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<   #
    |            #   |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/   #
    |            #                                               #
    |            #################################################
    |
    | compress_fq: Compress fastq files from a directory and check integrity.
    |
    |""".stripMargin()

if (params.help) {
  log.info paramsHelp("nextflow run nexomis/compress_fq --in_dir /path/to/fastq/dir --out_dir /path/to/out/dir")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

process CHECK {
  container 'ghcr.io/nexomis/check_fastq:1.0'

  publishDir "${params.out_dir}/.check_fastq/", mode: 'link', pattern: "*.yml"

  input:
  tuple val(meta), path(spring_files, arity: 1..2)
  tuple val(meta2), path(original_files, arity: 1..2)

  output:
  tuple val(meta), path("${meta.id}.yml", arity: 1)

  script:
  """
  #!/usr/bin/bash

  check_fastq.py \\
  -f1 ${spring_files[0]} ${spring_files.size() > 1 ? '-r1 ' + spring_files[1] : '' } \\
  -f2 ${original_files[0]} ${original_files.size() > 1 ? '-r2 ' + original_files[1] : '' } \\
  ${params.hash_algo != "" ? "--hash " + params.hash_algo : ''} \\
  ${params.n_bytes != 0 ? "--n_bytes " + params.n_bytes : ''} \\
  > ${meta.id}.yml

  # Read the first line of the file
  first_line=\$(head -n 1 "${meta.id}.yml")

  # Check if the first line matches the expected text
  if [ "\$first_line" = "sequences_equal: True" ]; then
    echo "First line matches expected text."
  else
    echo "Error: First line does not match the expected text."
    exit 1
  fi

  """

  stub:
  """
  #!/usr/bin/bash

  echo "sequences_equal: True" > ${sample_name}.yml
  echo "${params.hash_algo != "" ? "--hash " + params.hash_algo : ''}" >> ${sample_name}.yml
  echo "${params.n_bytes != 0 ? "--n_bytes " + params.n_bytes : ''}" >> ${sample_name}.yml

  # Read the first line of the file
  first_line=\$(head -n 1 "${sample_name}.yml")

  # Check if the first line matches the expected text
  if [ "\$first_line" = "sequences_equal: True" ]; then
    echo "First line matches expected text."
  else
    echo "Error: First line does not match the expected text."
    exit 1
  fi
  """

}

include { PARSE_SEQ_DIR } from './modules/subworkflows/parse_seq_dir/main.nf'
include { SPRING_COMPRESS } from './modules/process/spring/compress/main.nf'
include { SPRING_DECOMPRESS } from './modules/process/spring/decompress/main.nf'

workflow {

  Channel.fromPath(params.in_dir, type: 'dir', checkIfExists: true)
  | PARSE_SEQ_DIR

  reads = PARSE_SEQ_DIR.out.fastq

  SPRING_COMPRESS(reads)
  SPRING_DECOMPRESS(SPRING_COMPRESS.out)
  
  SPRING_DECOMPRESS.out
  | map {[it[0].id, it]}
  | join(reads.map {[it[0].id, it]})
  | set { joinReads }
  
  CHECK(joinReads.map {it[1]}, joinReads.map {it[2]})

}
