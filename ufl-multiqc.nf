#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.multiqc_params = [:]

include { MULTIQC_RUN } from './modules/multiqc/multiqc_run' addParams([*:params, "run_id" : params.run_id])

workflow MULTIQC {
    MULTIQC_RUN(
        params.run_dir
    )
}