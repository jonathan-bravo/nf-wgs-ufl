#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process SAMTOOLS_INDEX_MD {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/alignment", mode: 'copy'
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