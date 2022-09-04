#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CNV_TRIM_READS {

    tag "${sample_id}"
    label 'trimmomatic'
    label 'medium_process'

    input:
    tuple val(sample_id), path(reads)
    path trim_adapters

    output:
    tuple val(sample_id), file("${sample_id}_R1-p_trimmed.fastq.gz"), file("${sample_id}_R2-p_trimmed.fastq.gz"), emit: trimmed_paired_reads

    script:
    """
    TrimmomaticPE -threads ${task.cpus} \
    ${reads[0]} \
    ${reads[1]} \
    ${sample_id}_R1-p_trimmed.fastq.gz \
    ${sample_id}_R1-u_trimmed.fastq.gz \
    ${sample_id}_R2-p_trimmed.fastq.gz \
    ${sample_id}_R2-u_trimmed.fastq.gz \
    ILLUMINACLIP:${trim_adapters}:2:30:10:2:keepBothReads \
    SLIDINGWINDOW:4:20 CROP:149 MINLEN:36 2> ${sample_id}_trim_out.log
    """
}