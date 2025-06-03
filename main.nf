#!/usr/bin/env nextflow
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'
include { 
  PARSE_SEQ_DIR 
} from './modules/subworkflows/parse_seq_dir/main.nf'
include { SLIMFASTQ_COMPRESS } from './modules/process/slimfastq/compress/main.nf'

nextflow.preview.output = true

// Print help message, supply typical command line usage for the pipeline

workflow {
  main:
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

  Channel.fromPath(params.in_dir, type: 'dir', checkIfExists: true)
  | map { path -> 
       [ [
         depth: params.depth,
         parsingArgs: [
           tailPattern: params.tail_pattern,
           readPattern: params.read_pattern,
           lanePattern: params.lane_pattern
         ]
       ], path ] 
   }
  | PARSE_SEQ_DIR

  reads = PARSE_SEQ_DIR.out.fastq

  SLIMFASTQ_COMPRESS(reads)

  publish:
  slimfastq = SLIMFASTQ_COMPRESS.out
}

output {
  slimfastq {
    path "."
  }
}
