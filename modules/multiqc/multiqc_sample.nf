#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process MULTIQC_SAMPLE {
    
    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/MultiQC", mode: 'copy'
    label 'multiqc'
    label 'small_process'

    input:
    tuple val(sample_id), path(files)

    output:
    file "${sample_id}.html"
    path "${sample_id}_data"

    script:
    """
    multiqc -n ${sample_id} .
    """
}