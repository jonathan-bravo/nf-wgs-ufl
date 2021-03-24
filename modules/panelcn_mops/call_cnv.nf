#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_CNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'panelcn_mops'
    label 'medium_process'

    input:
    path male_control
    path female_control
    path header
    tuple val(sample_id), file("${sample_id}_m_or_f.txt")
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_cnv.vcf"), emit: cnv
    file("${sample_id}_cnv.pdf")

    script:
    """
    callCNV.R ${sample_id} ${male_control} ${female_control}
    csvToVCF.sh ${sample_id} ${header}
    toGRanges.sh ${sample_id}
    CNVPlot.R ${sample_id}
    """
}
