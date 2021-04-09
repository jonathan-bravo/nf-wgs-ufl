#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process ANNOTATE_VCF {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'snpeff_tabix'
    label 'high_mem'

    input:
    path dbNSFP
    path dbNSFP_tbi
    path dbNSFP_data_types
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/${sample_id}_variants.vcf.gz")

    output:
    tuple val(sample_id), file("${sample_id}_snpsift.vcf.gz"), emit: sift_vcf
    tuple val(sample_id), file("${sample_id}_snpeff_stats.csv"), emit: snpeff_stats

    script:
    """
    tabix ${sample_id}_strelka2/results/variants/${sample_id}_variants.vcf.gz

    java -jar /snpEff/snpEff.jar \
    -csvStats ${sample_id}_snpeff_stats.csv \
    -v -canon hg19 \
    ${sample_id}_strelka2/results/variants/${sample_id}_variants.vcf.gz \
    > ${sample_id}_snpeff.vcf

    bgzip -@ ${task.cpus} ${sample_id}_snpeff.vcf

    java -jar /snpEff/SnpSift.jar \
    dbnsfp -v -db ${dbNSFP} \
    ${sample_id}_snpeff.vcf.gz \
    > ${sample_id}_snpsift.vcf
    
    bgzip -@ ${task.cpus} ${sample_id}_snpsift.vcf
    """
}