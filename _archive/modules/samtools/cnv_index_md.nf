#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CNV_SAMTOOLS_INDEX_MD {

    tag "${sample_id}"
    label 'samtools'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}_md.bam")

    output:
    tuple val(sample_id), file("${sample_id}_md.bam.bai"), emit: index_md_bam

    script:
    """
    samtools index -@ ${task.cpus} ${sample_id}_md.bam
    """
}