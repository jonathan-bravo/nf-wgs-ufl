#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_VIEW {

    tag "${sample_id}"
    label 'samtools'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}.sam")

    output:
    tuple val(sample_id), file("${sample_id}.bam"), emit: bam

    script:
    """
    samtools view -@ ${task.cpus} -b -o ${sample_id}.bam ${sample_id}.sam
    """
}