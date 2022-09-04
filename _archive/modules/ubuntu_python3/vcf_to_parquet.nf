#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process VCF_TO_PARQUET {

    tag "${sample_id}"
    publishDir "s3://hakmonkey-datalake/germline_sample_datalake/", mode: 'copy'
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_concat.vcf.gz"), file("${sample_id}_concat.vcf.gz.tbi")

    output:
    tuple val(sample_id), file("${sample_id}.parquet"), emit: parquet
    

    script:
    """
    vcf_to_parquet.py -v ${sample_id}_concat.vcf.gz -s ${sample_id}
    """
}
