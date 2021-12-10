#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.run_id     = ""
params.sample_bam = "s3://hakmonkey-genetics-lab/Pipeline_Output/${params.run_id}/${params.run_id}-*/alignment/*_md.bam"
params.sample_bai = "s3://hakmonkey-genetics-lab/Pipeline_Output/${params.run_id}/${params.run_id}-*/alignment/*_md.bam.bai"
params.reference  = "s3://hakmonkey-genetics-lab/Pipeline/Reference/hs37d5/hs37d5.fa.gz"
params.ref_fai    = "${params.reference}.fai"
params.ref_gzi    = "${params.reference}.gzi"

bam_ch = Channel.fromFilePairs([params.sample_bam, params.sample_bai])


process MANTA_WGS {

    tag "${sample_id}"
    publishDir "s3://hakmonkey-genetics-lab/Pipeline_Output/${params.run_id}/${sample_id}/Manta", mode: 'copy'
    label 'manta'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    path ref_gzi
    tuple val(sample_id), path(bam_bai)

    output:
    tuple val(sample_id), path("${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz"), path("${sample_id}_manta/results/variants/${sample_id}_diploidSV.vcf.gz.tbi"), emit: manta_vcf

    script:
    """
    mkdir ${sample_id}_manta

    /manta-1.6.0.centos6_x86_64/bin/configManta.py \
    --bam ${bam_bai[0]} \
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

workflow {

    MANTA_WGS(
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )
}