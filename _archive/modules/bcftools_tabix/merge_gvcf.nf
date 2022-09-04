#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process MERGE_GVCF {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/variants", mode: 'copy'
    label 'bcftools_tabix'
    label 'small_process'

    input:
    tuple val(sample_id), path(files)

    output:
    tuple val(sample_id), file("${sample_id}_concat.gvcf.gz"), file("${sample_id}_concat.gvcf.gz.tbi"), emit: gvcf

    script:
    """
    bgzip -@ ${task.cpus} ${files[0]}
    bgzip -@ ${task.cpus} ${files[1]}
    bgzip -@ ${task.cpus} ${files[2]}

    bcftools index --threads ${task.cpus} --tbi ${files[0]}.gz
    bcftools index --threads ${task.cpus} --tbi ${files[1]}.gz
    bcftools index --threads ${task.cpus} --tbi ${files[2]}.gz

    bcftools concat --threads ${task.cpus} -a \
    -o ${sample_id}_concat.gvcf \
    ${files[0]}.gz \
    ${files[1]}.gz \
    ${files[2]}.gz

    bgzip -@ ${task.cpus} ${sample_id}_concat.gvcf

    bcftools index --threads ${task.cpus} --tbi ${sample_id}_concat.gvcf.gz
    """
}