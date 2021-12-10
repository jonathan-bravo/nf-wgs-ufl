#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CREATE_CNV_CONTROLS {

    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'
    label 'cn_mops'
    label 'medium_process'

    input:
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_cnv_control.RData")

    script:
    """
    CNV_Controls.R
    """
}
