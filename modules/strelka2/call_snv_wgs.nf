#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_SNV_WGS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'medium_process'

    input:
    path reference
    path ref_fai
    tuple val(sample_id), file("${sample_id}-sort.bam")
    tuple val(sample_id), file("${sample_id}-sort.bam.bai")

    output:
    path "${sample_id}_strelka2/results/variants/variants.vcf.gz", emit: snv
    
    script:
    """
    mkdir ${sample_id}_strelka2

    /strelka-2.9.10.centos6_x86_64/bin/configureStrelkaGermlineWorkflow.py \
    --bam ${sample_id}-sort.bam \
    --referenceFasta ${reference} \
    --runDir ${sample_id}_strelka2

    sed -i s/"isEmail = isLocalSmtp()"/"isEmail = False"/g ${sample_id}_strelka2/runWorkflow.py

    ${sample_id}_strelka2/runWorkflow.py -m local -j ${task.cpus}
    """
}