#!/usr/bin/env nextflow

if ( params.help ) {
  help = """
    |compress_fq: A workflow for the compression of fasq files.
    |
    |Required arguments:
    |  --input_dir      Reference sample for the annotation
    |  --output_dir     Output directory
    |  --quality_mode   "lossless" (default)
    |                   "qvz qv_ratio" (QVZ lossy compression,
    |                     parameter qv_ratio roughly corresponds to
    |                     bits used per quality value)
    |                   "ill_bin" (Illumina 8-level binning)
    |                   "binary thr" high low (binary (2-level)
    |                     thresholding, quality binned to high if >=
    |                     thr and to low if < thr)
    |  --drop_order     true/false, whether to drop reads ordering
    |  --drop_ids       true/false, whether to drop reads ids
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

process COMPRESS_PAIRED {
  container 'ghcr.io/nexomis/spring:1.1.1'

  publishDir "${params.output_dir}/", mode: 'link', pattern: "${sample_name}.spring"

  label 'cpu_x4'
  label 'mem_16G'

  input:
  tuple val(sample_name), path(files, arity: 2)

  output:
  tuple val("${sample_name}"), path("${sample_name}.spring.*.fq", arity: 2)

  script:
  """
  #!/usr/bin/bash

  spring -q ${params.quality_mode} -t ${task.cpus} -c -i ${files} -o ${sample_name}.spring \\
    ${params.drop_order ? '-r' : ''} \\
    ${params.drop_ids ? '--no-ids' : ''} \\
    ${files[0].toString().endsWith('.gz') || files[0].toString().endsWith('.gzip') || files[0].toString().endsWith('.z') ? '-g' : ''}

  spring -d -t ${task.cpus} -i ${sample_name}.spring -o ${sample_name}.spring.1.fq ${sample_name}.spring.2.fq

  """

  stub:
  """
  #!/usr/bin/bash

  touch ${sample_name}.spring ${sample_name}.spring.1.fq ${sample_name}.spring.2.fq

  """
}


process CHECK_PAIRED {
  container 'ghcr.io/nexomis/check_fastq:1.0'

  label 'cpu_x1'
  label 'mem_2G'

  publishDir "${params.output_dir}/", mode: 'link', pattern: ".check_${sample_name}.yml"

  input:
  tuple val(sample_name), path(spring_files, arity: 2), path(original_files, arity: 2)

  script:
  """
  #!/usr/bin/bash

  check_fastq.py -f1 ${spring_files[0]} -r1 ${spring_files[1]} -f2 ${original_files[0]} -r2 ${original_files[1]} \\
  ${params.hash_algo != "" ? "--hash " + params.hash_algo : ''} \\
  ${params.n_bytes != 0 ? "--n_bytes " + params.n_bytes : ''} \\
  > .check_${sample_name}.yml

  # Read the first line of the file
  first_line=\$(head -n 1 ".check_${sample_name}.yml")

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

  cat  ${spring_files[0]} ${spring_files[1]}
  cat  ${original_files[0]} ${original_files[1]}

  echo "sequences_equal: True" > .check_${sample_name}.yml
  echo "${params.hash_algo != "" ? "--hash " + params.hash_algo : ''}" >> .check_${sample_name}.yml
  echo "${params.n_bytes != 0 ? "--n_bytes " + params.n_bytes : ''}" >> .check_${sample_name}.yml

  # Read the first line of the file
  first_line=\$(head -n 1 ".check_${sample_name}.yml")

  # Check if the first line matches the expected text
  if [ "\$first_line" = "sequences_equal: True" ]; then
    echo "First line matches expected text."
  else
    echo "Error: First line does not match the expected text."
    exit 1
  fi
  """

}


workflow{
  Channel.fromFilePairs([
    params.input_dir + "/*[._]{1,2}*.fastq",
    params.input_dir + "/*[._]{1,2}*.fq",
    params.input_dir + "/*[._]{1,2}*.fastq.gz",
    params.input_dir + "/*[._]{1,2}*.fq.gz",
    params.input_dir + "/*[._][Rr]{1,2}*.fastq",
    params.input_dir + "/*[._][Rr]{1,2}*.fq",
    params.input_dir + "/*[._][Rr]{1,2}*.fastq.gz",
    params.input_dir + "/*[._][Rr]{1,2}*.fq.gz"
  ])
  | set {pairedInput}

  Channel.fromPath("${params.input_dir}/*.{fastq,fq}{.gz,}")
  | set { allFiles }

  pairedInput
  | map { it -> it[1] } 
  | flatten()
  | toList()
  | set {pairedFlat}

  allFiles.combine(pairedFlat)
  | filter { !it[1..-1].contains(it[0])}
  | map {  
      def basename = it[0].getName()
      def sampleName = basename.find(/(.+?)\.(fastq|fq)\.?gz?$/) { it[1] }
      return tuple(sampleName, tuple(it[0]))
    }
  | set {singleInput}

  COMPRESS_PAIRED(pairedInput)
  | join(pairedInput, by: 0)
  | CHECK_PAIRED

}
