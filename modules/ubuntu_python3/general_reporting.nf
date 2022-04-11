#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process GENERAL_REPORTING {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/General_Report", mode: 'copy'
    label 'ubuntu_python3'
    label 'reporting'

    input:
    path reporting_resources
    tuple val(sample_id), path(vcf), path(tbi)

    output:
    tuple val(sample_id), file("${sample_id}_report.json"), file("${sample_id}_low-qc_report.json"), file("${sample_id}_report.xlsx"), file("${sample_id}_low-qc_report.xlsx")

    script:
    """
    reporting.py -v ${vcf} -t ${task.cpus} -s ${sample_id}

    json_to_csv.py -j ${sample_id}_report.json
    json_to_csv.py -j ${sample_id}_low-qc_report.json
    """
}