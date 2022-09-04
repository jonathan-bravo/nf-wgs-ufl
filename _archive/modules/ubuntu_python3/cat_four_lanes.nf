#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CAT_FOUR_LANES {

    tag "${sample_id}"
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), path(reads1), path(reads2), path(reads3), path(reads4)

    output:
    tuple val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz"), emit: read_pairs

    script:
    """
    cat ${reads1[0]} ${reads2[0]} ${reads3[0]} ${reads4[0]} > ${sample_id}_R1.fastq.gz
    cat ${reads1[1]} ${reads2[1]} ${reads3[1]} ${reads4[1]} > ${sample_id}_R2.fastq.gz
    """
}