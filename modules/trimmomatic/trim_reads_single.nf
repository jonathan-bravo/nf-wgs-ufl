#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRIM_READS_SINGLE {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Trimmomatic", mode: 'copy'
    label 'trimmomatic'
    label 'medium_process'

    input:
    tuple val(sample_id), path(reads)
    path trim_adapters

    output:
    tuple val(sample_id), file("${sample_id}_R1-p_trimmed.fastq.gz"), file("${sample_id}_R2-p_trimmed.fastq.gz"), emit: trimmed_paired_reads
    tuple file("${sample_id}_R1-u_trimmed.fastq.gz"), file("${sample_id}_R2-u_trimmed.fastq.gz"), emit: trimmed_unpaired_reads
    tuple val(sample_id), file("${sample_id}_trim_out.log"), emit: trim_log

    script:
    """
    TrimmomaticPE -threads ${task.cpus} \
    ${sample_id}_R1.fastq.gz \
    ${sample_id}_R2.fastq.gz \
    ${sample_id}_R1-p_trimmed.fastq.gz \
    ${sample_id}_R1-u_trimmed.fastq.gz \
    ${sample_id}_R2-p_trimmed.fastq.gz \
    ${sample_id}_R2-u_trimmed.fastq.gz \
    ILLUMINACLIP:${trim_adapters}:2:30:10:2:keepBothReads \
    SLIDINGWINDOW:4:20 CROP:149 MINLEN:36 2> ${sample_id}_trim_out.log
    """
}