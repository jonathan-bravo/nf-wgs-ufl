#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process TRIM_READS_SINGLE {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Trimmomatic", mode: 'copy'
    label 'small_process'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("${sample_id}_forward-paired.fastq.gz"), file("${sample_id}_reverse-paired.fastq.gz"), emit: trimmed_paired_reads
    tuple val(sample_id), file("${sample_id}_forward-unpaired.fastq.gz"), file("${sample_id}_reverse-unpaired.fastq.gz"), emit: trimmed_unpaired_reads
    tuple val(sample_id), file("${sample_id}_trim_out.log")

    script:
    """
    TrimmomaticPE -threads ${task.cpus} \
    ${reads[0]} \
    ${reads[1]} \
    ${sample_id}_forward-paired.fastq.gz \
    ${sample_id}_forward-unpaired.fastq.gz \
    ${sample_id}_reverse-paired.fastq.gz \
    ${sample_id}_reverse-unpaired.fastq.gz \
    ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
    LEADING:3 TRAILING:3 MINLEN:36 2> ${sample_id}_trim_out.log
    """
}