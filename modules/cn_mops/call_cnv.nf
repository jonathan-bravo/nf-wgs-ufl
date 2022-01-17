#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CALL_CNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/cn_MOPS", mode: 'copy'
    label 'cn_mops'
    label 'medium_process'

    input:
    path cnv_control
    path header
    tuple val(sample_id), file("${sample_id}_md.bam"), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf"), emit: cnv_vcf
    tuple val(sample_id), file("${sample_id}_cnv.vcf"), emit: cnv_gvcf
    file("${sample_id}_cnv.pdf")

    script:
    """
    callCNV.R ${sample_id}
    csvToVCF.sh ${sample_id} ${header}
    toGRanges.sh ${sample_id}
    CNVPlot.R ${sample_id}
    """
}
