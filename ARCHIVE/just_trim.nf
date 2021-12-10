#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.run_id        = ""
params.bucket        = "s3://hakmonkey-genetics-lab"
params.outdir        = "${params.bucket}/Pipeline_Output/TRIM/${params.run_id}"
params.ref_dir       = "${params.bucket}/Pipeline/Reference"
params.trim_adapters = "${params.ref_dir}/trim/NEBNext.fa"

params.reads1 = "${params.bucket}/Fastqs/${params.run_id}*_L001_{R1,R2}_001.fastq.gz"
params.reads2 = "${params.bucket}/Fastqs/${params.run_id}*_L002_{R1,R2}_001.fastq.gz"
params.reads3 = "${params.bucket}/Fastqs/${params.run_id}*_L003_{R1,R2}_001.fastq.gz"
params.reads4 = "${params.bucket}/Fastqs/${params.run_id}*_L004_{R1,R2}_001.fastq.gz"

reads1_ch = Channel.fromFilePairs(params.reads1)
reads2_ch = Channel.fromFilePairs(params.reads2)
reads3_ch = Channel.fromFilePairs(params.reads3)
reads4_ch = Channel.fromFilePairs(params.reads4)

process CAT_FOUR_LANES {

    tag "${sample_id}"
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), path(reads1)
    tuple val(sample_id), path(reads2)
    tuple val(sample_id), path(reads3)
    tuple val(sample_id), path(reads4)

    output:
    tuple val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz"), emit: read_pairs

    script:
    """
    cat ${reads1[0]} ${reads2[0]} ${reads3[0]} ${reads4[0]} > ${sample_id}_R1.fastq.gz
    cat ${reads1[1]} ${reads2[1]} ${reads3[1]} ${reads4[1]} > ${sample_id}_R2.fastq.gz
    """
}

process TRIM_READS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Trimmomatic", mode: 'copy'
    label 'trimmomatic'
    label 'medium_process'

    input:
    tuple val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz")
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

workflow {
    CAT_FOUR_LANES(
        reads1_ch,
        reads2_ch,
        reads3_ch,
        reads4_ch,
    )

    TRIM_READS(
        CAT_FOUR_LANES.out.read_pairs,
        params.trim_adapters
    )
}