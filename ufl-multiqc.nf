#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.multiqc_params = [:]

include { MULTIQC_RUN } from './modules/multiqc/multiqc' addParams([*:params, "run_id" : params.run_id])
include { STOP_EC2    } from './modules/ec2/stop'

workflow MULTIQC {
    MULTIQC_RUN(
        params.run_dir
    )
    STOP_EC2(
        MULTIQC_RUN.out.qc
    )
}