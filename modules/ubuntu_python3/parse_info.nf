#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process PARSE_INFO {

    tag "${sample_id}_${panel}"
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), val(panel), file("${sample_id}_${panel}_OPL.vcf")

    output:
    tuple val(sample_id), val(panel), file("${sample_id}_${panel}.data"), file("${sample_id}_${panel}.info.tsv"), emit: data

    script:
    """
    generateTmpFiles.sh ${sample_id} ${panel}
    parseInfo.py ${sample_id} ${panel}
    """
}