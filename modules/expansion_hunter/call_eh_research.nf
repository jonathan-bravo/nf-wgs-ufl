#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_EH_RESEARCH {

    tag "${sample_id}"
    publishDir "${params.outdir}/${params.run_id}/${sample_id}/ExpansionHunter_Research", mode: 'copy'
    label 'expansion_hunter'
    label 'small_process'

    input:
    path research_variant_catalog
    path reference
    path ref_fai
    tuple val(sample_id), file("${sample_id}_md.bam")
    tuple val(sample_id), file("${sample_id}_md.bam.bai")

    output:
    tuple val(sample_id), file("${sample_id}_filtered_eh.vcf"), emit: eh_vcf_research
    tuple val(sample_id), file("${sample_id}_eh.vcf"), emit: eh_gvcf_research
    file "${sample_id}_eh_research_realigned.bam"
    file "${sample_id}_eh_research.json"

    shell:
    '''
    /ExpansionHunter-v4.0.2-linux_x86_64/bin/ExpansionHunter \
    --reads !{sample_id}_md.bam \
    --reference !{reference} \
    --variant-catalog !{research_variant_catalog} \
    --output-prefix !{sample_id}_eh_research

    sed -i s/!{sample_id}_md/SAMPLE1/g !{sample_id}_eh.vcf

    grep '^#' !{sample_id}_eh_research.vcf > !{sample_id}_filtered_eh_research.vcf
    grep '^chr' !{sample_id}_eh_research.vcf | \
    awk '{if ($5 !=".")print}' - >> !{sample_id}_filtered_eh_research.vcf
    '''
}