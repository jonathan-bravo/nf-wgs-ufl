#!/usr/bin/env nextflow

params.run_dir    = ""
params.panels_dir = ""
params.bucket     = ""
params.ref_dir    = "${params.bucket}/Pipeline/Reference"
params.reference  = "${params.ref_dir}/hg19/hg19.fa"
params.ref_fai    = "${params.reference}.fai"
params.glue_dir   = "${params.bucket}/Pipeline_Output/_SampleTSV"


pairs_ch = Channel.from()

log.info """\

         R E P O R T I N G - U F L    P I P E L I N E
         ============================================
         run directory    : ${params.run_dir}
         panels directory : ${params.panels_dir}
         
         """
         .stripIndent()

process applyPanel {

    tag "${sample_id}_${panel}"
    label 'small_process'
    
    input:
    tuple sample_id, panel from pairs_ch
    path panel_dir from params.panels_dir
    path sample_path from params.run_dir

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf") into (panel_ch1, panel_ch2)

    script:
    """
    applyPanel.sh ${sample_id} ${panel} ${panel_dir} ${sample_path}
    """
}

process generateTmpFiles {

    tag "${sample_id}_${panel}"
    label 'small_process'

    input:
    tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf") from panel_ch1

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.data") into data_ch
    tuple sample_id, panel, file("${sample_id}_${panel}.info.txt") into info_ch

    script:
    """
    generateTmpFiles.sh ${sample_id} ${panel}
    """
}

process parseInfo {

    tag "${sample_id}_${panel}"
    label 'small_process'
    echo true

    input:
    tuple sample_id, panel, file("${sample_id}_${panel}.info.txt") from info_ch

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.info.tsv") into info_tsv_ch

    script:
    """
    parseInfo.py ${sample_id} ${panel}
    """
}

process createFinalTSV {

    tag "${sample_id}_${panel}"
    publishDir "${params.glue_dir}/${sample_id}/${panel}", mode: 'copy'
    label 'small_process'
    echo true

    input:
    tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf") from panel_ch2
    tuple sample_id, panel, file("${sample_id}_${panel}.data") from data_ch
    tuple sample_id, panel, file("${sample_id}_${panel}.info.tsv") from info_tsv_ch

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.final.tsv") into final_tsv_ch

    script:
    """
    createFinalTSV.sh ${sample_id} ${panel}
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
