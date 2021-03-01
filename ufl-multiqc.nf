#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

//params.run_id
//params.run_dir

include { MULTIQC_RUN } from './modules/multiqc/multiqc' addParams([*:params, "run_id" : params.run_id])

workflow MULTIQC {
    MULTIQC_RUN(
        params.run_dir
    )
}