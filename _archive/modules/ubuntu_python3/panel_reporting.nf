#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

process PANEL_REPORTING {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/Panel_Reports", mode: 'copy'
    label 'ubuntu_python3'
    label 'reporting'

    input:
    path reporting_resources
    tuple val(sample_id), path(vcf), path(panel), path(tbi)

    output:
    tuple val(sample_id), file("${sample_id}_${panel}_report.json"), file("${sample_id}_${panel}_low-qc_report.json"), file("${sample_id}_${panel}_report.xlsx"), file("${sample_id}_${panel}_low-qc_report.xlsx")

    script:
    """
    reporting.py -v ${vcf[0]} -t ${task.cpus} -s ${sample_id} -p ${panel}

    json_to_csv.py -j ${sample_id}_${panel}_report.json
    json_to_csv.py -j ${sample_id}_${panel}_low-qc_report.json
    """
}