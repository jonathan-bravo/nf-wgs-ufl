#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_EH {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/ExpansionHunter", mode: 'copy'
    label 'expansion_hunter'
    label 'small_process'

    input:
    path reference
    path ref_fai
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_eh.vcf"), emit: eh_vcf
    tuple val(sample_id), file("${sample_id}_eh.vcf"), emit: eh_gvcf
    file "${sample_id}_eh_realigned.bam"
    file "${sample_id}_eh.json"

    shell:
    '''
    /ExpansionHunter-v4.0.2-linux_x86_64/bin/ExpansionHunter \
    --reads !{sample_id}_md.bam \
    --reference !{reference} \
    --variant-catalog /ExpansionHunter-v4.0.2-linux_x86_64/variant_catalog/hg19/variant_catalog.json \
    --output-prefix !{sample_id}_eh

    sed -i s/!{sample_id}-sort/SAMPLE1/g !{sample_id}_eh.vcf

    grep '^#' !{sample_id}_eh.vcf > !{sample_id}_filtered_eh.vcf
    grep '^chr' !{sample_id}_eh.vcf | \
    awk '{if ($5 !=".")print}' - >> !{sample_id}_filtered_eh.vcf
    '''
}