#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process ANNOTATE_EH {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/ExpansionHunter", mode: 'copy'
    label 'ubuntu_python3'
    label 'small_process'

    input:
    path variant_catalog
    tuple val(sample_id), file("${sample_id}_filtered_eh.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_eh_ann.vcf"), emit: eh_ann

    script:
    """
    annotate_eh.py -v ${sample_id}_filtered_eh.vcf -c ${variant_catalog}
    """
}