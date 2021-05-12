#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process ANNOTATE_CNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'ubuntu_python3'
    label 'high_mem'

    input:
    path hg19_genes
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf.gz"), file("${sample_id}_filtered_cnv.vcf.gz.tbi")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_cnv_ann.vcf"), emit: cnv_ann

    script:
    """
    annotate_cnv.py -v ${sample_id}_filtered_cnv.vcf.gz -b ${hg19_genes}
    """
}