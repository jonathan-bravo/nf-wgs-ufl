#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process SAMTOOLS_SORT {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/alignment", mode: 'copy'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}.bam")

    output:
    tuple val(sample_id), file("${sample_id}-sort.bam"), emit: sort_bam

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}-sort.bam ${sample_id}.bam
    """
}