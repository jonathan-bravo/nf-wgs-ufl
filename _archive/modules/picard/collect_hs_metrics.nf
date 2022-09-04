#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PICARD_COLLECT_HS_METRICS {

    tag "${sample_id}"
    publishDir "${params.run_dir}_Exome/${sample_id}/hs_metrics", mode: 'copy'
    label 'picard'
    label 'high_mem'

    input:
    path reference
    path ref_fai
    path ref_gzi
    path target
    path bait
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_gatk_collect_hs_metrics.txt"), emit: hs_metrics

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx62g \
    /picard.jar CollectHsMetrics \
    -I ${sample_id}_md.bam \
    -O ${sample_id}_gatk_collect_hs_metrics.txt \
    -R ${reference} \
    --BAIT_INTERVALS ${bait} \
    --TARGET_INTERVALS ${target}
    """
}