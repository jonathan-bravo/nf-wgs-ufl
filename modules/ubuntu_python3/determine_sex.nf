#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process DETERMINE_SEX {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_m_or_f.txt"), emit: sex

    script:
    """
    sexTest.py -b ${sample_id}-sort.bam -s ${sample_id} -t ${task.cpus}
    """
}