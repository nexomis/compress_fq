#!/usr/bin/env nextflow

include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

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
  log.info paramsHelp("nextflow run nexomis/compress_fq --input_dir /path/to/fastq/dir --output_dir /path/to/out/dir")
  exit 0
}
validateParameters()
log.info paramsSummaryLog(workflow)

// Validate input parameters
if (params.input_dir == null) {
  error """No input directory provided. Use --input_dir.""".stripMargin()
}

if (params.output_dir == null) {
  error """No output directory provided. Use --output_dir.""".stripMargin()
}

process COMPRESS {
  container 'ghcr.io/nexomis/spring:1.1.1'

  publishDir "${params.output_dir}/", mode: 'link', pattern: "${sample_name}.spring"

  label 'cpu_x4'
  label 'mem_16G'

  input:
  tuple val(sample_name), path(files, arity: 1..2)

  output:
  tuple val("${sample_name}"), path("${sample_name}.spring*fq.gz", arity: 1..2), emit: fastq
  tuple val("${sample_name}"), path("${sample_name}.spring", arity: 1), emit: spring

  script:
  """
  #!/usr/bin/bash

  spring -q ${params.quality_mode} -t ${task.cpus} -c -i ${files} -o ${sample_name}.spring \\
    ${params.drop_order ? '-r' : ''} \\
    ${params.drop_ids ? '--no-ids' : ''} \\
    ${files[0].toString().endsWith('.gz') || files[0].toString().endsWith('.gzip') || files[0].toString().endsWith('.z') ? '-g' : ''}

  spring -d -g -t ${task.cpus} -i ${sample_name}.spring -o ${sample_name}.spring${files.size() > 1 ? '.1' : '' }.fq.gz ${files.size() > 1 ? sample_name + '.spring.2.fq.gz' : '' }

  """

  stub:
  """
  #!/usr/bin/bash

  touch ${sample_name}.spring ${sample_name}.spring${files.size() > 1 ? '.1' : '' }.fq ${files.size() > 1 ? sample_name + '.spring.2.fq' : '' }

  """
}

process CHECK {
  container 'ghcr.io/nexomis/check_fastq:1.0'

  label 'cpu_x1'
  label 'mem_2G'

  publishDir "${params.output_dir}/.check_fastq/", mode: 'link', pattern: "${sample_name}.yml"

  input:
  tuple val(sample_name), path(spring_files, arity: 1..2), path(original_files, arity: 1..2)

  output:
  tuple val("${sample_name}"), path("${sample_name}.yml", arity: 1)

  script:
  """
  #!/usr/bin/bash

  check_fastq.py \\
  -f1 ${spring_files[0]} ${spring_files.size() > 1 ? '-r1 ' + spring_files[1] : '' } \\
  -f2 ${original_files[0]} ${original_files.size() > 1 ? '-r2 ' + original_files[1] : '' } \\
  ${params.hash_algo != "" ? "--hash " + params.hash_algo : ''} \\
  ${params.n_bytes != 0 ? "--n_bytes " + params.n_bytes : ''} \\
  > ${sample_name}.yml

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

workflow {

  Channel.fromPath(params.input_dir, type: 'dir')
  | PARSE_SEQ_DIR
  | set {reads}

  reads.paired
  | concat(reads.single)
  | set {allReads}

  COMPRESS(allReads).fastq
  | map { return tuple(it[0], it[1].sort { a, b -> a.name <=> b.name }) }
  | join(allReads, by: 0)
  | CHECK
}
