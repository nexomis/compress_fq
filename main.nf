#!/usr/bin/env nextflow
nextflow.preview.output = true
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

include { PARSE_SEQ_DIR } from './modules/subworkflows/parse_seq_dir/main.nf'
include { SLIMFASTQ_COMPRESS } from './modules/process/slimfastq/compress/main.nf'

workflow {

  Channel.fromPath(params.in_dir, type: 'dir', checkIfExists: true)
  | PARSE_SEQ_DIR

  reads = PARSE_SEQ_DIR.out.fastq

  SLIMFASTQ_COMPRESS(reads)

  publish:
  SLIMFASTQ_COMPRESS.out >> 'slimfastq'

}

def parseOutDir(outDir) {
    def idx = outDir.lastIndexOf('/')
    if (idx == -1) {
        return [path: '', prefix: outDir]
    } else {
        def prefix = outDir.substring(idx + 1)
        def path = outDir.substring(0, idx + 1)
        return [path: path, prefix: prefix]
    }
}

output {

  'slimfastq' {
    path ""
  }

}
