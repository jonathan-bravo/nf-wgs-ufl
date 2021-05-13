#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process INDEX_CNV {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'bcftools_tabix'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf.gz"), file("${sample_id}_filtered_cnv.vcf.gz.tbi"), emit: cnv_index

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_filtered_cnv.vcf
    bcftools index --threads ${task.cpus} --tbi ${sample_id}_filtered_cnv.vcf.gz
    """
}