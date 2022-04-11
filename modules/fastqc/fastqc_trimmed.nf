#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC_TRIMMED {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/", mode: 'copy'
    label 'fastqc'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_R1-p_trimmed.fastq.gz"), file("${sample_id}_R2-p_trimmed.fastq.gz")

    output:
    tuple val(sample_id), path("fastqc_${sample_id}_trimmed_logs"), emit: qc_trimmed

    script:
    """
    mkdir fastqc_${sample_id}_trimmed_logs
    fastqc \
    -t ${task.cpus} \
    -o fastqc_${sample_id}_trimmed_logs \
    -f fastq \
    ${sample_id}_R1-p_trimmed.fastq.gz \
    ${sample_id}_R2-p_trimmed.fastq.gz
    """
}