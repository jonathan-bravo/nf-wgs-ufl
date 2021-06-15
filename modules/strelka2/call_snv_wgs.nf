#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_SNV_WGS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Strelka2", mode: 'copy'
    label 'strelka2'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/${sample_id}_variants.vcf.gz"), emit: snv_vcf
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf"), emit: snv_gvcf
    
    script:
    """
    mkdir ${sample_id}_strelka2

    /strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${sample_id}_md.bam \
    --referenceFasta ${reference} \
    --runDir ${sample_id}_strelka2

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${sample_id}_strelka2/runWorkflow.py

    ${sample_id}_strelka2/runWorkflow.py -m local -j ${task.cpus}

    mv ${sample_id}_strelka2/results/variants/variants.vcf.gz \
    ${sample_id}_strelka2/results/variants/${sample_id}_variants.vcf.gz
    
    mv ${sample_id}_strelka2/results/variants/genome.S1.vcf.gz \
    ${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf.gz

    gunzip ${sample_id}_strelka2/results/variants/${sample_id}_genome.S1.vcf.gz
    """
}