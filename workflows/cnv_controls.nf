#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.outdir = "${params.bucket}/CNV_Controls/RData/"
params.reads_path = "${params.bucket}/CNV_Controls/Fastqs/*_{1,2}.fastq.gz"

reads_ch = Channel.fromFilePairs(params.reads_path)

include { CNV_TRIM_READS             } from '../modules/trimmomatic/cnv_trim_reads'
include { SAMTOOLS_SORT              } from '../modules/samtools/sort'                 
include { SAMTOOLS_INDEX             } from '../modules/samtools/index'                
include { CNV_SAMTOOLS_INDEX_MD      } from '../modules/samtools/cnv_index_md'             
include { CNV_PICARD_MARK_DUPLICATES } from '../modules/picard/cnv_mark_duplicates'        
include { ALIGN_TRIMMED_READS        } from '../modules/bwa/align_trimmed_reads'
include { SAMTOOLS_VIEW              } from '../modules/samtools/view'
include { CREATE_CNV_CONTROLS        } from '../modules/cn_mops/create_cnv_controls' addParams(outdir : params.outdir)

workflow CNV_CONTROLS {

    CNV_TRIM_READS(
        reads_ch,
        params.trim_adapters
    )

    ALIGN_TRIMMED_READS(
        params.reference,
        params.bwa_amb,
        params.bwa_ann,
        params.bwa_bwt,
        params.bwa_pac,
        params.bwa_sa,
        CNV_TRIM_READS.out.trimmed_paired_reads
    )

    SAMTOOLS_VIEW(
        ALIGN_TRIMMED_READS.out.sam
    )

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam
    )

    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.sort_bam
    )

    CNV_PICARD_MARK_DUPLICATES(
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
    )

    CNV_SAMTOOLS_INDEX_MD(
        CNV_PICARD_MARK_DUPLICATES.out.md_bam
    )

    CREATE_CNV_CONTROLS(
        CNV_PICARD_MARK_DUPLICATES.out.md_bam,
        CNV_SAMTOOLS_INDEX_MD.out.index_md_bam
    )
}   
