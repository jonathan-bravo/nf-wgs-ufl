#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_INDEX {

    tag "${sample_id}"
    label 'samtools'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}-sort.bam")

    output:
    tuple val(sample_id), file("${sample_id}-sort.bam.bai"), emit: index_sort_bam

    script:
    """
    samtools index -@ ${task.cpus} ${sample_id}-sort.bam
    """
}