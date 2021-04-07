#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_GVCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'bcftools_tabix'
    label 'small_process'

    input:
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf.gz")
    tuple val(sample_id), file("${sample_id}_cnv.vcf")
    tuple val(sample_id), file("${sample_id}_eh.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_concat.gvcf.gz"), file("${sample_id}_concat.gvcf.gz.tbi"), emit: vcf

    script:
    """
    bgzip -@ ${task.cpus} ${sample_id}_cnv.vcf
    bgzip -@ ${task.cpus} ${sample_id}_eh.vcf

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf.gz

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_cnv.vcf.gz

    bcftools index --threads ${task.cpus} --tbi \
    ${sample_id}_eh.vcf.gz

    bcftools concat --threads ${task.cpus} -a \
    -o ${sample_id}_concat.gvcf \
    ${sample_id}_cnv.vcf.gz \
    ${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf.gz \
    ${sample_id}_eh.vcf.gz

    bgzip -@ ${task.cpus} ${sample_id}_concat.gvcf

    bcftools index --threads ${task.cpus} --tbi ${sample_id}_concat.gvcf.gz
    """
}