#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process MANTA_WGS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/Manta", mode: 'copy'
    label 'manta'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    path ref_gzi
    tuple val(sample_id), file("${sample_id}_md.bam"), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), path("${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz"), path("${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz.tbi"), emit: manta_vcf

    script:
    """
    mkdir ${sample_id}_manta

    /manta-1.6.0.centos6_x86_64/bin/configManta.py \
    --bam ${sample_id}_md.bam \
    --referenceFasta ${reference} \
    --runDir ${sample_id}_manta

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${sample_id}_manta/runWorkflow.py

    ${sample_id}_manta/runWorkflow.py -j ${task.cpus}

    mv ${sample_id}_manta/results/variants/diploidSV.vcf.gz \
    ${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz

    mv ${sample_id}_manta/results/variants/diploidSV.vcf.gz.tbi \
    ${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz.tbi
    """
}