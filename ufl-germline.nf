#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.single_lane    = ""
params.exome          = ""
params.match          = ""
params.bucket         = ""
params.run_id         = ""
params.ref_dir        = "${params.bucket}/Pipeline/Reference"

params.reference      = "${params.ref_dir}/hg19/hg19.fa"
params.bwa_amb        = "${params.reference}.amb"
params.bwa_ann        = "${params.reference}.ann"
params.bwa_bwt        = "${params.reference}.bwt"
params.bwa_pac        = "${params.reference}.pac"
params.bwa_sa         = "${params.reference}.sa"
params.ref_fai        = "${params.reference}.fai"

params.cnv_control    = "${params.ref_dir}/cnv/sliding_windows_sim_control.RData"
params.cnv_vcf_header = "${params.ref_dir}/cnv/cnv_vcf_header.tsv"

params.bait           = "${params.ref_dir}/exome_targets/bait.interval_list"
params.target         = "${params.ref_dir}/exome_targets/target.interval_list"

params.dbnsfp         = "${params.ref_dir}/dbnsfp/dbNSFP4.1a.txt.gz"
params.dbnsfp_tbi     = "${params.dbnsfp}.tbi"
params.dbnsfp_dt      = "${params.dbnsfp}.data_types"

params.outdir         = "${params.bucket}/Pipeline_Output"



include { CAT_LANES                  } from './modules/cat/cat_lanes'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC                     } from './modules/fastqc/fastqc'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_SINGLE              } from './modules/fastqc/fastqc_single'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS                 } from './modules/trimmomatic/trim_reads'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS_SINGLE          } from './modules/trimmomatic/trim_reads_single' addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ALIGN_TRIMMED_READS        } from './modules/bwa/align_trimmed_reads'       addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_VIEW              } from './modules/samtools/view'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_SORT              } from './modules/samtools/sort'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_INDEX             } from './modules/samtools/index'                addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_WGS_METRICS } from './modules/picard/collect_wgs_metrics'    addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_HS_METRICS  } from './modules/picard/collect_hs_metrics'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WGS               } from './modules/strelka2/call_snv_wgs'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WES               } from './modules/strelka2/call_snv_wes'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_CNV                   } from './modules/panelcn.mops/call_cnv'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_EH                    } from './modules/expansion_hunter/call_eh'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_VCF                  } from './modules/bcftools_tabix/merge_vcf'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_VCF               } from './modules/snpeff_tabix/annotate_vcf'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])



if (params.single_lane == "YES"){
    if (params.exome == "YES"){
        params.reads_path = "${params.bucket}/Exome_Fastqs/${params.run_id}*${params.match}"
    }
    else{
        params.reads_path = "${params.bucket}/Fastqs/${params.run_id}*${params.match}"
    }

    reads_ch  = Channel.fromFilePairs(params.reads_path)
}
else {
    if (params.exome == "YES"){
        params.reads1 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_{L001,L002}_R1_001.fastq.gz"
        params.reads2 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_{L001,L002}_R2_001.fastq.gz"
    }
    else {
        params.reads1 = "${params.bucket}/Fastqs/${params.run_id}*_{L001,L002}_R1_001.fastq.gz"
        params.reads2 = "${params.bucket}/Fastqs/${params.run_id}*_{L001,L002}_R2_001.fastq.gz"
    }

    reads1_ch = Channel.fromFilePairs(params.reads1)
    reads2_ch = Channel.fromFilePairs(params.reads2)
}



workflow GERMLINE {

    if (params.single_lane == "NO"){
        CAT_LANES(
            reads1_ch,
            reads2_ch
        )

        FASTQC(
            CAT_LANES.out.read_pairs
        )

        TRIM_READS(
            CAT_LANES.out.read_pairs
        )
    }
    else {
        FASTQC_SINGLE(
            reads_ch
        )

        TRIM_READS_SINGLE(
            reads_ch
        )
    }

    ALIGN_TRIMMED_READS(
        params.reference,
        params.bwa_amb,
        params.bwa_ann,
        params.bwa_bwt,
        params.bwa_pac,
        params.bwa_sa,
        TRIM_READS.out.trimmed_paired_reads
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

    if (params.exome == "YES"){
        PICARD_COLLECT_HS_METRICS(
            params.reference,
            params.ref_fai,
            params.target,
            params.bait,
            SAMTOOLS_SORT.out.sort_bam,
            SAMTOOLS_INDEX.out.index_sort_bam
        )

        CALL_SNV_WES(
            params.reference,
            params.ref_fai,
            SAMTOOLS_SORT.out.sort_bam,
            SAMTOOLS_INDEX.out.index_sort_bam
        )
    }
    else {
        PICARD_COLLECT_WGS_METRICS(
            params.reference,
            params.ref_fai,
            SAMTOOLS_SORT.out.sort_bam,
            SAMTOOLS_INDEX.out.index_sort_bam
        )
        CALL_SNV_WGS(
            params.reference,
            params.ref_fai,
            SAMTOOLS_SORT.out.sort_bam,
            SAMTOOLS_INDEX.out.index_sort_bam
        )
    }

    CALL_CNV(
        params.cnv_control,
        params.cnv_vcf_header,
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
    )

    CALL_EH(
        params.reference,
        params.ref_fai,
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
    )

    if (params.exome == "YES"){
        MERGE_VCF(
            CALL_SNV_WES.out.snv,
            CALL_CNV.out.cnv,
            CALL_EH.out.eh
        )
    }
    else {
        MERGE_VCF(
            CALL_SNV_WGS.out.snv,
            CALL_CNV.out.cnv,
            CALL_EH.out.eh
        )
    }
    
    ANNOTATE_VCF(
        params.dbnsfp,
        params.dbnsfp_tbi,
        params.dbnsfp_dt,
        MERGE_VCF.out.vcf
    )
}
