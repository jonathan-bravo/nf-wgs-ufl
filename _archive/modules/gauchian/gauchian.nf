#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process GAUCHIAN {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/Gauchian_GBA", mode: 'copy'
    label 'gauchian'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    path ref_gzi
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("gauchian_out/${sample_id}_gauchian.json"), path("gauchian_out/${sample_id}_gauchian.tsv"), emit: gba

    script:
    """
    echo ${bam} > manifest.txt

    gauchian \
    --manifest manifest.txt \
    --genome 37 \
    --prefix ${sample_id}_gauchian \
    --outDir gauchian_out \
    --threads ${task.cpus}
    """
}