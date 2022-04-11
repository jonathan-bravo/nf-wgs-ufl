#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CYRIUS {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/Cyrius_CYP2D6", mode: 'copy'
    label 'cyrius'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    path ref_gzi
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), path("cyrius_out/${sample_id}_cyrius.json"), path("cyrius_out/${sample_id}_cyrius.tsv"), emit: cyp2d6

    script:
    """
    echo ${bam} > manifest.txt

    python3 /Cyrius-1.1.1/star_caller.py \
    --manifest manifest.txt \
    --genome 37 \
    --prefix ${sample_id}_cyrius \
    --outDir cyrius_out \
    --threads ${task.cpus}
    """
}