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

    tag "${sample_id}_${panel}"
	label 'small_process'
    
    input:
    tuple sample_id, panel from pairs_ch
    path panel_dir from params.panels_dir
    path sample_path from params.run_dir

    output:
    tuple sample_id, panel, file("${sample_id}_${panel}.vcf"), file("${sample_id}_eh_${panel}.vcf") into paneled_ch

    shell:
    '''
    GENES=$(tr -d '\\r' <!{panel_dir}/!{panel} | tr '\n' '|')

    zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz > !{sample_id}_!{panel}.vcf

	zgrep '#' !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz > !{sample_id}_eh_!{panel}.vcf

    zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_concat_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_!{panel}.vcf

	zcat !{sample_path}/!{sample_id}/variants/!{sample_id}_eh_snpsift.vcf.gz | awk -v g="${GENES}" '$0 ~ g' - >> !{sample_id}_eh_!{panel}.vcf
    '''
}

process onePerLine {

	tag "${sample_id}_${panel}"
	publishDir "${params.run_dir}/${sample_id}/${panel}", mode: 'copy'
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
/*
process qualityFilter {

	tag "${sample_id}_${panel}"
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
*/

process generateTmpFiles {

	tag "${sample_id}_${panel}"
	label 'small_process'

	input:
	tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf"), file("${sample_id}_eh_${panel}_OPL.vcf") from opl_ch1
	//tuple sample_id, panel, file("${sample_id}_${panel}_hq.vcf"), file("${sample_id}_eh_${panel}_hq.vcf") from hq_ch

	output:
	tuple sample_id, panel, file("${sample_id}_${panel}.data") into data_ch
	tuple sample_id, panel, file("${sample_id}_${panel}.info.txt") into info_ch

	shell:
	'''
	zgrep '^#CHROM' !{sample_id}_!{panel}_OPL.vcf > !{sample_id}_!{panel}.headers
	sed -i 's/\tINFO//' !{sample_id}_!{panel}.headers
	sed -i 's/#//' !{sample_id}_!{panel}.headers

	cat !{sample_id}_!{panel}.headers | tr '\n' '\t' > !{sample_id}_!{panel}.data
	echo >> !{sample_id}_!{panel}.data

	zgrep '^chr' !{sample_id}_!{panel}_OPL.vcf | awk -F'\t' '{
    	print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10"\t";
	}' >> !{sample_id}_!{panel}.data

	sed -rn 's/^##INFO=<ID=([0-9a-zA-Z.-_\\;]+),.*/\\1/gp' !{sample_id}_!{panel}_OPL.vcf > !{sample_id}_!{panel}.info.headers

	cat !{sample_id}_!{panel}.info.headers | tr '\n' '\t' > !{sample_id}_!{panel}.info.txt
	echo >> !{sample_id}_!{panel}.info.txt

	zgrep '^chr' !{sample_id}_!{panel}_OPL.vcf | awk -F'\t' '{print $8}' | sed -e s'/\\;/\t/g' >> !{sample_id}_!{panel}.info.txt
	'''
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
	#!/usr/bin/env python3
header = []
f = open('${sample_id}_${panel}.info.tsv', 'w')
with open('${sample_id}_${panel}.info.txt', 'r') as fp:
    for cnt, line in enumerate(fp):
        if cnt == 0:
            header.append(line.split('\\t'))
            header[0].pop()
            header[0] = header[0]
            header = sorted(header[0])
            for i in header:
                f.write('{}\\t'.format(i))
            f.write('\\n')
        else:
            sorted_line = sorted(line.split('\\t'))
            sorted_line.pop()
            for i, element in enumerate(sorted_line):
                if "\\n" in element:
                    sorted_line[i] = element[:-2]
            for i in range(len(header) - len(sorted_line)):
                sorted_line.append('.')
            for col_num, val in enumerate(header):
                if not val in sorted_line[col_num] and sorted_line[col_num] not in ".":
                    sorted_line[col_num], sorted_line[col_num + 1] = sorted_line[col_num + 1], sorted_line[col_num]
            for i in sorted_line:
                f.write('{}\\t'.format(i))
            f.write('\\n')
f.close()
	"""
}

process createFinalTSV {

	tag "${sample_id}_${panel}"
	publishDir "${params.run_dir}/${sample_id}/${panel}", mode: 'copy'
	label 'small_process'

	input:
	tuple sample_id, panel, file("${sample_id}_${panel}_OPL.vcf"), file("${sample_id}_eh_${panel}_OPL.vcf") from opl_ch2
	tuple sample_id, panel, file("${sample_id}_${panel}.data") from data_ch
	tuple sample_id, panel, file("${sample_id}_${panel}.info.tsv") from info_tsv_ch

	output:
	tuple sample_id, panel, file("${sample_id}_${panel}.tsv")

	shell:
	'''
	ann=$(zgrep '^##INFO=<ID=ANN' !{sample_id}_!{panel}_OPL.vcf | cut -c75-316 | sed -e 's/\\s//g' | sed -e 's/|/\t/g')

	sed -r -i 's/([0-9a-zA-Z.-_]+=)//g' !{sample_id}_!{panel}.info.tsv

	sed -i "s|ANN\t|$ann|" !{sample_id}_!{panel}.info.tsv

	sed -i 's/|/\t/g' !{sample_id}_!{panel}.info.tsv

	awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' !{sample_id}_!{panel}.info.tsv > !{sample_id}_!{panel}.info.final.tsv

	paste !{sample_id}_!{panel}.data !{sample_id}_!{panel}.info.final.tsv > !{sample_id}_!{panel}.tsv
	'''
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