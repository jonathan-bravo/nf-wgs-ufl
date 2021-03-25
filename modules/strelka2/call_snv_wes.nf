#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_SNV_WES {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'strelka2'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/variants.vcf.gz"), emit: snv_vcf
    tuple val(sample_id), path("${sample_id}_strelka2/results/variants/genome.S1.vcf.gz"), emit: snv_gvcf
    
    script:
    """
    mkdir ${sample_id}_strelka2

    /strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${sample_id}-sort.bam \
    --referenceFasta ${reference} \
    --runDir ${sample_id}_strelka2 \
    --exome

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${sample_id}_strelka2/runWorkflow.py

    ${sample_id}_strelka2/runWorkflow.py -m local -j ${task.cpus}
    """
}