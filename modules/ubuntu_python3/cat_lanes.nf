#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CAT_LANES {

    tag "${sample_id}"
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), path(reads1)
    tuple val(sample_id), path(reads2)

    output:
    tuple val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz"), emit: read_pairs

    script:
    """
    cat ${reads1[0]} ${reads1[1]} > ${sample_id}_R1.fastq.gz
    cat ${reads2[0]} ${reads2[1]} > ${sample_id}_R2.fastq.gz
    """
}