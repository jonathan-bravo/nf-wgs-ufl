#!/usr/bin/env nextflow

params.run_id
params.run_dir

log.info """\

         M U L T I Q C - U F L    P I P E L I N E
         ========================================
         run directory path : ${params.run_dir}
         run id             : ${params.run_id}
         
         """
         .stripIndent()

process multiqcRun {
    
    tag "${params.run_id}"
    publishDir "${params.run_dir}MultiQC"
    label 'small_process'

    input:
    path run_dir from params.run_dir

    output:
    file "multiqc_report.html" into multiqc_report
    file "multiqc_data"

    script:
    """
    multiqc -n ${params.run_id} .
    """
}