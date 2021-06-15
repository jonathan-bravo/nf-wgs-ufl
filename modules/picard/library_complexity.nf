#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process PICARD_ESTIMATE_LIBRARY_COMPLEXITY {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/library_complexity", mode: 'copy'
    label 'picard'
    label 'high_mem'

    input:
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_est_lib_complex_metrics.txt"), emit: lib_complex_metrics

    script:
    """
    java -jar -XX:ParallelGCThreads=${task.cpus} -Xmx62g \
    /picard.jar EstimateLibraryComplexity \
    -I ${sample_id}_md.bam \
    -O ${sample_id}_est_lib_complex_metrics.txt
    """
}