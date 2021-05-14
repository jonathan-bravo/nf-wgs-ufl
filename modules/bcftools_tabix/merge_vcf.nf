#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_VCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'bcftools_tabix'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_snpsift.vcf.gz"), file("${sample_id}_filtered_cnv_ann.vcf"), file("${sample_id}_filtered_eh.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_concat.vcf.gz"), file("${sample_id}_concat.vcf.gz.tbi"), emit: vcf

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_filtered_cnv_ann.vcf
    bgzip -@ ${task.cpus} ${sample_id}_filtered_eh.vcf

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_snpsift.vcf.gz

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_filtered_cnv_ann.vcf.gz

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_filtered_eh.vcf.gz

    bcftools concat --threads ${task.cpus} -a \
    -o ${sample_id}_concat.vcf \
    ${sample_id}_filtered_cnv_ann.vcf.gz \
    ${sample_id}_snpsift.vcf.gz \
    ${sample_id}_filtered_eh.vcf.gz

    bgzip -@ ${task.cpus} ${sample_id}_concat.vcf

    bcftools index --threads ${task.cpus} --tbi ${sample_id}_concat.vcf.gz
    """
}