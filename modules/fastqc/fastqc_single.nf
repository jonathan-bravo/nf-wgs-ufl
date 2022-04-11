#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_SINGLE {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/", mode: 'copy'
    label 'fastqc'
    label 'small_process'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("fastqc_${sample_id}_logs"), emit: qc

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc \
    -t ${task.cpus} \
    -o fastqc_${sample_id}_logs \
    -f fastq \
    ${reads[0]} \
    ${reads[1]}
    """
}