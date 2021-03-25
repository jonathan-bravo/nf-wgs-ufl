#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process CNV_CONTIGS {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/variants", mode: 'copy'
    label 'ubuntu_python3'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_filtered_cnv.vcf")

    output:
    file("${sample_id}_cnv_contigs_report.tsv")

    script:
    """
    cnv_contigs.py -c ${sample_id}_filtered_cnv.vcf -s ${sample_id}
    """
}