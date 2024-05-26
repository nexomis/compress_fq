#!/usr/bin/env nextflow

if ( params.help ) {
  help = """
    |________________________________________________________________________________
    |                _  _                             _      
    |               | \\| |  ___  __ __  ___   _ __   (_)  ___
    |               | .` | / -_) \\ \\ / / _ \\ | '  \\  | | (_-<
    |               |_|\\_| \\___| /_\\_\\ \\___/ |_|_|_| |_| /__/
    |                                      
    |_#################______________________________________________________________
    | ## compress_fq ##
    |_#################______________________________________________________________
    |
    | A workflow for the compression of fastq files with a step to check the sequence
    |  identity after decompression even with loosy compression regardeing, read 
    |  order, read identities and qualities.  
    |
    |_#################______________________________________________________________
    | ## Input/Ouput ##
    |_#################______________________________________________________________
    |
    |  --input_dir      Input directory with uncompressed fastq files.
    |                    Potentially in gzip, auto id paired and single reads
    |
    |  --output_dir     Output directory with compressed files (.spring)
    |
    |_###########################____________________________________________________
    | ## Compression Arguments ##
    |_###########################____________________________________________________
    |
    |  --quality_mode   "lossless" (default)
    |                   "qvz qv_ratio" (QVZ lossy compression, parameter qv_ratio 
    |                     roughly corresponds to bits used per quality value)
    |                   "ill_bin" (Illumina 8-level binning)
    |                   "binary thr" high low (binary (2-level) thresholding, quality
    |                     binned to high if >= thr and to low if < thr)
    |
    |  --drop_order     true/false, whether to drop reads ordering
    |
    |  --drop_ids       true/false, whether to drop reads ids
    |
    |_#####################__________________________________________________________
    | ## Check Arguments ##
    |_#####################__________________________________________________________
    |
    |  --hash_algo      Algorithm to build an hash digest of reads before  checking 
    |                    identity based on occurences. To spare memory. Must be in 
    |                    python hashlib. Default is None.
    |                   
    |  --n_bytes        Only in combination with hash_algo. Number of bytes to use to
    |                    count occurences. Default is to take all bytes from the hash
    |________________________________________________________________________________
    """.stripMargin()
  // Print the help with the stripped margin and exit
  println(help)
  exit(0)
}

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
  tuple val("${sample_name}"), path("${sample_name}.spring*fq", arity: 1..2), emit: fastq
  tuple val("${sample_name}"), path("${sample_name}.spring", arity: 1), emit: spring

  script:
  """
  #!/usr/bin/bash

  spring -q ${params.quality_mode} -t ${task.cpus} -c -i ${files} -o ${sample_name}.spring \\
    ${params.drop_order ? '-r' : ''} \\
    ${params.drop_ids ? '--no-ids' : ''} \\
    ${files[0].toString().endsWith('.gz') || files[0].toString().endsWith('.gzip') || files[0].toString().endsWith('.z') ? '-g' : ''}

  spring -d -t ${task.cpus} -i ${sample_name}.spring -o ${sample_name}.spring${files.size() > 1 ? '.1' : '' }.fq ${files.size() > 1 ? sample_name + '.spring.2.fq' : '' }

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
