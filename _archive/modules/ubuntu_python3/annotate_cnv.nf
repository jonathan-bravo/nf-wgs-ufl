#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process ANNOTATE_CNV {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/cn_MOPS", mode: 'copy'
    label 'ubuntu_python3'
    label 'medium_process'

    input:
    path hs37d5_genes
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_cnv_ann.vcf"), emit: cnv_ann

    script:
    """
    annotate_cnv.py -v ${sample_id}_filtered_cnv.vcf -b ${hs37d5_genes}
    """
}