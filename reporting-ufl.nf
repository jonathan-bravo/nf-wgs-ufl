#!/usr/bin/env nextflow

params.run_dir
params.panels_dir

pairs_ch = Channel.from()

process applyPanel {

    tag "${sample_id}"
	publishDir "${params.run_dir}/${sample_id}/variants", mode: 'copy'
	label 'small_process'
    
    input:
    tuple sample_id, panel from pairs_ch
    path panel_dir from params.panels_dir
    path sample_path from params.run_dir

    output:
    tuple sample_id, file("${sample_id}_${panel}.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel_dir}/!{panel} | tr '\n' '|')

    zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_!{panel}.vcf

	zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz > !{sample_id}_eh_!{panel}.vcf

    zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_!{panel}.vcf

	zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_eh_!{panel}.vcf
    '''
}

/* 
	3. Visualize filtered output
	4. Generate report template based on panel(s) selected and data found
	   - This will require that I parse the hits provided depending on how much of the report
	   is being automatically generated

	Comment from Lee:
	So in the share drive Clinical Genomics > Lab Med & Path > Epilepsy exomes each of those folders has a report. The template looked reasonable to me.

	Word docs or PDFs
*/


/*
process visualizeVCF {
	
	tag "${sample_id}"
	
	input:
	// filtered vcf file from applyPanel process
	
	output:
	// page that is the visualization of the filtered data
	
	script:
	"""
	vcftools --gvcf ${sample_id}.vcf.gz --out ${sample_id} --plink
	gemini -v ${sample_id}.vcf.gz -p ${sample_id}.ped -t snpEff --cores ${task.cpus} ${sample_id}.db
	
	"""
}

process generateReportTemplate {
	
	tag "${sample_id}"
	
	input: 
	//Output from parseSampleSheet
	//filtered vcf file from applyPanel process (maybe)
	
	output:
	//report template that contains genes, commands run [tools] (maybe),
	//and most problimatic phenotypes/ mutations (maybe)
	
	script:
	"""
	"""
}
*/