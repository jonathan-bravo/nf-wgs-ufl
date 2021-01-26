#!/usr/bin/env nextflow

params.run_dir
params.panels_dir
params.bucket      = ""
params.ref_dir     = "${params.bucket}/Pipeline/Reference"
params.reference   = "${params.ref_dir}/hg19/hg19.fa"
params.ref_fai     = "${params.reference}.fai"


pairs_ch = Channel.from()

log.info """\

         R E P O R T I N G - U F L    P I P E L I N E
         ============================================
         run directory    : ${params.run_dir}
         panels directory : ${params.panels_dir}
		 
         """
         .stripIndent()

process applyPanel {

    tag "${sample_id}-${panel}"
	publishDir "${params.run_dir}/${sample_id}/${panel}", mode: 'copy'
	label 'small_process'
    
    input:
    tuple sample_id, panel from pairs_ch
    path panel_dir from params.panels_dir
    path sample_path from params.run_dir

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.vcf"), file("${sample_id}_eh_${panel}.vcf") into paneled_ch

    shell:
    '''
    GENES=$(tr -d '\r' <!{panel_dir}/!{panel} | tr '\n' '|')

    zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_!{panel}.vcf

	zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz > !{sample_id}_eh_!{panel}.vcf

    zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_!{panel}.vcf

	zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_eh_!{panel}.vcf
    '''
}

process onePerLine {

	tag "${sample_id}-${panel}"
	publishDir "${params.run_dir}/${sample_id}/${sample_id}-${panel}", mode: 'copy'
	label 'small_process'

	input:
	tuple sample_id, panel, file("${sample_id}_${panel}.vcf"), file("${sample_id}_eh_${panel}.vcf") from paneled_ch

	output:
	tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf"), file("${sample_id}_eh_${panel}_OPL.vcf") into (opl_ch1, opl_ch2)

	script:
	"""
	cat ${sample_id}_${panel}.vcf | /snpEff/scripts/vcfEffOnePerLine.pl > ${sample_id}_${panel}_OPL.vcf

	cat ${sample_id}_eh_${panel}.vcf | /snpEff/scripts/vcfEffOnePerLine.pl > ${sample_id}_eh_${panel}_OPL.vcf
	"""
}

// Going to do a qc filtered version and a non qc filtered version

process qualityFilter {

	tag "${sample_id}-${panel}"
	publishDir "${params.run_dir}/${sample_id}/${sample_id}-${panel}", mode: 'copy'
	label 'small_process'

	input:
	tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf"), file("${sample_id}_eh_${panel}_OPL.vcf") from opl_ch1
	path reference from params.reference
	path ref_fai from params.ref_fai

	output:
	tuple sample_id, panel, file("${sample_id}_${panel}_hq.vcf"), file("${sample_id}_eh_${panel}_hq.vcf") into hq_ch

	script:
	"""
	bcftools view -R ${reference} -i'FILTER="PASS"' --threads ${task.cpus} ${sample_id}_${panel}_OPL.vcf -o ${sample_id}_${panel}_hq.vcf

	bcftools view -R ${reference} -i'FILTER="PASS"' --threads ${task.cpus} ${sample_id}_eh_${panel}_OPL.vcf -o ${sample_id}_eh_${panel}_hq.vcf
	"""
}


process visualizePanel {

	tag "${sample_id}-${panel}"
	publishDir "${params.run_dir}/${sample_id}/${sample_id}-${panel}", mode: 'copy'
	label 'small_process'

	input:
	tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf"), file("${sample_id}_eh_${panel}_OPL.vcf") from opl_ch2
	tuple sample_id, panel, file("${sample_id}_${panel}_hq.vcf"), file("${sample_id}_eh_${panel}_hq.vcf") into hq_ch

	output:

	script:
	"""
	#!/usr/bin/env python3

	import pandas as pd
	"""
}

// Also want a component for visualizing `.bam` file using IGV

/* 
	3. Visualize filtered output
	4. Generate report template based on panel(s) selected and data found
	   - This will require that I parse the hits provided depending on how much of the report
	   is being automatically generated

	Comment from Lee:
	So in the share drive Clinical Genomics > Lab Med & Path > Epilepsy exomes each of those folders has a report. The template looked reasonable to me.

	Word docs or PDFs
*/

// Want to add a component that stores all variants from "panel filtered" data
// into a database that can be searched by gene or by sample


/*

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