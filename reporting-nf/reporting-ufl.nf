#!/usr/bin/env nextflow

pairs_ch = Channel.from([['NQ-20-10-BC702503-143_S28','ataxia']])

params.panels_dir = "s3://hakmonkey-genetics-lab/Pipeline/Reference/panels"





params.vcf_dirs = "s3://hakmonkey-genetics-lab/Pipeline_Output/_processed/**/variants"

// !{vcf_dir}/!{sample_id}/variants/

process applyPanel {

    tag "${sample_id}"
    echo true
    
    input:
    tuple sample_id, panel from pairs_ch
    path panel_dir from params.panels_dir
    path vcf_paths from params.vcf_dirs

    // !{vcf_paths}/!{sample_id}_concat_snpsift.vcf.gz

    output:
    tuple sample_id, file("${sample_id}_${panel}.vcf")

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel_dir}/!{panel} | tr '\n' '|')

    zgrep '#' !{vcf_paths}/!{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_!{panel}.vcf

    zcat !{vcf_paths}/!{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_!{panel}.vcf
    '''
}


/* 
	1. Parse sample sheet
	   - This is done in the python script, which will automatically populate
	   the pairs_ch in this nextflow file and will then reset it back to
	   `pairs_ch = Channel.from()`
	2. Apply correct filters
	3. Visualize filtered output
	4. Generate report template based on panel(s) selected and data found
	   - This will require that I parse the hits provided depending on how much of the report
	   is being automatically generated

	Comment from Lee:
	So in the share drive Clinical Genomics > Lab Med & Path > Epilepsy exomes each of those folders has a report. The template looked reasonable to me.

	Word docs or PDFs
*/

/*
// This is the s3 bucket that contains the wgs pipeline output
params.vcf_dirs = "s3://hakmonkey-genetics-lab/Pipeline_Output"
// This is the s3 bucket that holds all the panels
params.panels_dir = "s3://hakmonkey-genetics-lab/Pipeline/Reference/panels"
*/
/*
	These are the panels we currently have:
	- ataxia
	- dementia
	- dystonia
	- epilepsy
	- hsp
	- neuromuscular
	- neuropathy
	- parkinsons
	- sma
*/

// This is a list of lists that containes the sample_id and the panel name
// This will be used to determine which folder to look into and which panel to
// apply
/*
pairs_ch = Channel.from()

process filterVCF {
	
	tag "${sample_id}"
	publishDir "${params.outdir}/${sample_id}/filteredVariants", mode: 'copy'
	
	input:
	// output from parse sample sheet
	tuple sample_id, test from pairs_ch
	path panel_dir from params.panels_dir
	path vcf from 
	
	output:
	// filtered sample into *2* channels
	
	script:
	"""
	"""
}

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