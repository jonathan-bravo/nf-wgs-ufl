#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process PICARD_COLLECT_HS_METRICS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}_Exome/${sample_id}/hs_metrics", mode: 'copy'
    label 'high_mem'

    input:
    path reference
    path ref_fai
    path target
    path bait
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    file "${sample_id}_gatk_collect_hs_metrics.txt"

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx32g /picard.jar CollectHsMetrics \
    -I ${sample_id}-sort.bam \
    -O ${sample_id}_gatk_collect_hs_metrics.txt \
    -R ${reference} \
    --BAIT_INTERVALS ${bait} \
    --TARGET_INTERVALS ${target}
    """
}