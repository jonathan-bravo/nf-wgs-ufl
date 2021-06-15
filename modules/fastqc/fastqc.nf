#!/usr/bin/env nextflow

nextflow.enable.dsl=2

#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

process FASTQC {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/", mode: 'copy'
    label 'fastqc'
    label 'small_process'

    input:
    tuple val(sample_id), file("${sample_id}_R1.fastq.gz"), file("${sample_id}_R2.fastq.gz")

    output:
    tuple val(sample_id), path("fastqc_${sample_id}_logs"), emit: qc

    script:
    """
    mkdir fastqc_${sample_id}_logs
    fastqc \
    -t ${task.cpus} \
    -o fastqc_${sample_id}_logs \
    -f fastq \
    ${sample_id}_R1.fastq.gz \
    ${sample_id}_R2.fastq.gz
    """
}