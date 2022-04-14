#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process CALL_EHD {

    tag "${sample_id}"
    publishDir "${params.run_dir}/${sample_id}/ExpansionHunterDenovo", mode: 'copy'
    label 'expansion_hunter_denovo'
    label 'medium_process'

    input:
    path ehd_controls
    path reference
    path ref_fai
    path ref_gzi
    tuple val(sample_id), path(bam), path(bai)

    output:
    tuple val(sample_id), file("${sample_id}_ehd.vcf"), emit: ehd_vcf
    file "${sample_id}_ehd.str_profile.json"
    file "${sample_id}_ehd.locus.tsv"
    file "${sample_id}_ehd.motif.tsv"

    script:
    """
    /ExpansionHunterDenovo/bin/ExpansionHunterDenovo profile \
    --reads ${bam} \
    --reference ${reference} \
    --output-prefix ${sample_id}_ehd \
    --min-anchor-mapq 50 \
    --max-irr-mapq 40

    echo -e "${sample_id}\tcase\t${sample_id}_ehd.str_profile.json" >> manifest.txt

    /ExpansionHunterDenovo/bin/ExpansionHunterDenovo merge \
    --reference ${reference} \
    --manifest manifest.txt \
    --output-prefix sample_dataset

    /ExpansionHunterDenovo/scripts/casecontrol.py locus \
    --manifest manifest.txt \
    --multisample-profile sample_dataset.multisample_profile.json \
    --output sample_dataset.casecontrol_locus.tsv

    locus_finder.py -s ${sample_id} -t sample_dataset.casecontrol_locus.tsv
    """
}