#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.germline_params = [:]

include { FASTQC                     } from './modules/fastqc/fastqc'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_SINGLE              } from './modules/fastqc/fastqc_single'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS                 } from './modules/trimmomatic/trim_reads'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS_SINGLE          } from './modules/trimmomatic/trim_reads_single' addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_SORT              } from './modules/samtools/sort'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_INDEX             } from './modules/samtools/index'                addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_WGS_METRICS } from './modules/picard/collect_wgs_metrics'    addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_HS_METRICS  } from './modules/picard/collect_hs_metrics'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WGS               } from './modules/strelka2/call_snv_wgs'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WES               } from './modules/strelka2/call_snv_wes'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_CNV                   } from './modules/panelcn_mops/call_cnv'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_EH                    } from './modules/expansion_hunter/call_eh'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_VCF                  } from './modules/bcftools_tabix/merge_vcf'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_VCF               } from './modules/snpeff_tabix/annotate_vcf'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { DETERMINE_SEX              } from './modules/ubuntu_python3/determine_sex'  addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CNV_CONTIGS                } from './modules/ubuntu_python3/cnv_contigs'    addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CAT_LANES                  } from './modules/ubuntu_python3/cat_lanes'
include { ALIGN_TRIMMED_READS        } from './modules/bwa/align_trimmed_reads'
include { SAMTOOLS_VIEW              } from './modules/samtools/view'




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

    DETERMINE_SEX(
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
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
        params.male_cnv_control,
        params.female_cnv_control,
        params.cnv_vcf_header,
        DETERMINE_SEX.out.sex
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
    )

    CNV_CONTIGS(
        CALL_CNV.out.cnv
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
