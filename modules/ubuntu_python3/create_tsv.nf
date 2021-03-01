#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process CREATE_TSV {

    tag "${sample_id}_${panel}"
    publishDir "${params.glue_dir}/${sample_id}/${panel}", mode: 'copy'
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), val(panel), file("${sample_id}_${panel}_OPL.vcf")
    tuple val(sample_id), val(panel), file("${sample_id}_${panel}.data"), file("${sample_id}_${panel}.info.tsv")

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.final.tsv")

    script:
    """
    createFinalTSV.sh ${sample_id} ${panel}
    """
} 