#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_VCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'bcftools_tabix'
    label 'small_process'

    input:
    path("${sample_id}_strelka2/results/variants/genome.S1.vcf.gz")
    tuple val(sample_id), file("${sample_id}_cnv.vcf")
    tuple val(sample_id), file("${sample_id}_eh.vcf")

    output:
    tuple val(sample_id), file("${sample_id}_concat.vcf.gz"), file("${sample_id}_concat.vcf.gz.csi"), emit: vcf

    script:
    """
    # change SAMPLE1 to ${sample_id}-sort in Strelka2 VCF file

    gunzip ${sample_id}_strelka2/results/variants/genome.S1.vcf.gz

    sed -i s/SAMPLE1/${sample_id}-sort/g ${sample_id}_strelka2/results/variants/genome.S1.vcf

    bgzip -@ ${task.cpus} ${sample_id}_strelka2/results/variants/genome.S1.vcf
    bgzip -@ ${task.cpus} ${sample_id}_cnv.vcf
    bgzip -@ ${task.cpus} ${sample_id}_eh.vcf

    bcftools index --threads ${task.cpus} \
    ${sample_id}_strelka2/results/variants/genome.S1.vcf.gz

    bcftools index --threads ${task.cpus} \
    ${sample_id}_cnv.vcf.gz

    bcftools index --threads ${task.cpus} \
    ${sample_id}_eh_vcf.gz

    bcftools concat --threads ${task.cpus} -a \
    -o ${sample_id}_concat.vcf \
    ${sample_id}_cnv.vcf.gz \
    ${sample_id}_strelka2/results/variants/genome.S1.vcf.gz \
    ${sample_id}_eh.vcf.gz

    bgzip -@ ${task.cpus} ${sample_id}_concat.vcf

    bcftools index --threads ${task.cpus} ${sample_id}_concat.vcf.gz
    """
}