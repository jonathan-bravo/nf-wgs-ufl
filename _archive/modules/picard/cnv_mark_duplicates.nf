#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CNV_PICARD_MARK_DUPLICATES {

    tag "${sample_id}"
    label 'picard'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_md.bam"), emit: md_bam
    tuple val(sample_id), file("${sample_id}_md_metrics.txt"), emit: md_metrics

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx62g \
    /picard.jar MarkDuplicates \
    --TAGGING_POLICY All \
    -I ${sample_id}-sort.bam \
    -O ${sample_id}_md.bam \
    -M ${sample_id}_md_metrics.txt
    """
}