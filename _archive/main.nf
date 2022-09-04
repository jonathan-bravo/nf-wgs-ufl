#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { GERMLINE     } from './workflows/germline'
include { CNV_CONTROLS } from './workflows/cnv_controls'


workflow {
    if ( params.cnv_controls ) {
        CNV_CONTROLS()
    }
    else {
        GERMLINE()
    }
