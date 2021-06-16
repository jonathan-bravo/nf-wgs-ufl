#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.germline_params = [:]

include { FASTQC                             } from './modules/fastqc/fastqc'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_TRIMMED                     } from './modules/fastqc/fastqc_trimmed'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_SINGLE                      } from './modules/fastqc/fastqc_single'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS                         } from './modules/trimmomatic/trim_reads'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS_SINGLE                  } from './modules/trimmomatic/trim_reads_single' addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_SORT                      } from './modules/samtools/sort'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_INDEX                     } from './modules/samtools/index'                addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_INDEX_MD                  } from './modules/samtools/index_md'             addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_WGS_METRICS         } from './modules/picard/collect_wgs_metrics'    addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_HS_METRICS          } from './modules/picard/collect_hs_metrics'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_MARK_DUPLICATES             } from './modules/picard/mark_duplicates'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_ESTIMATE_LIBRARY_COMPLEXITY } from './modules/picard/library_complexity'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WGS                       } from './modules/strelka2/call_snv_wgs'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WES                       } from './modules/strelka2/call_snv_wes'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_CNV                           } from './modules/cn_mops/call_cnv'              addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { INDEX_CNV                          } from './modules/bcftools_tabix/index_cnv'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_CNV                       } from './modules/ubuntu_python3/annotate_cnv'   addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_EH                            } from './modules/expansion_hunter/call_eh'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_VCF                          } from './modules/bcftools_tabix/merge_vcf'      addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_GVCF                         } from './modules/bcftools_tabix/merge_gvcf'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_VCF                       } from './modules/snpeff_tabix/annotate_vcf'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MULTIQC_SAMPLE                     } from './modules/multiqc/multiqc_sample'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CAT_LANES                          } from './modules/ubuntu_python3/cat_lanes'
include { ALIGN_TRIMMED_READS                } from './modules/bwa/align_trimmed_reads'
include { SAMTOOLS_VIEW                      } from './modules/samtools/view'




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
            CAT_LANES.out.read_pairs,
            params.trim_adapters
        )

        FASTQC_TRIMMED(
            TRIM_READS.out.trimmed_paired_reads
        )

        ALIGN_TRIMMED_READS(
            params.reference,
            params.bwa_amb,
            params.bwa_ann,
            params.bwa_bwt,
            params.bwa_pac,
            params.bwa_sa,
            TRIM_READS.out.trimmed_paired_reads
        )
    }
    else {
        FASTQC_SINGLE(
            reads_ch
        )

        TRIM_READS_SINGLE(
            reads_ch,
            params.trim_adapters
        )
        FASTQC_TRIMMED(
            TRIM_READS_SINGLE.out.trimmed_paired_reads
        )

        ALIGN_TRIMMED_READS(
            params.reference,
            params.bwa_amb,
            params.bwa_ann,
            params.bwa_bwt,
            params.bwa_pac,
            params.bwa_sa,
            TRIM_READS_SINGLE.out.trimmed_paired_reads
        )
    }

    SAMTOOLS_VIEW(
        ALIGN_TRIMMED_READS.out.sam
    )

    SAMTOOLS_SORT(
        SAMTOOLS_VIEW.out.bam
    )

    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.sort_bam
    )

    PICARD_MARK_DUPLICATES(
        SAMTOOLS_SORT.out.sort_bam,
        SAMTOOLS_INDEX.out.index_sort_bam
    )

    SAMTOOLS_INDEX_MD(
        PICARD_MARK_DUPLICATES.out.md_bam
    )

    PICARD_ESTIMATE_LIBRARY_COMPLEXITY(
        PICARD_MARK_DUPLICATES.out.md_bam,
        SAMTOOLS_INDEX_MD.out.index_md_bam
    )

    if (params.exome == "YES"){
        PICARD_COLLECT_HS_METRICS(
            params.reference,
            params.ref_fai,
            params.target,
            params.bait,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        CALL_SNV_WES(
            params.reference,
            params.ref_fai,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        ANNOTATE_VCF(
            params.dbnsfp,
            params.dbnsfp_tbi,
            params.dbnsfp_dt,
            CALL_SNV_WES.out.snv_vcf
        )
    }
    else {
        PICARD_COLLECT_WGS_METRICS(
            params.reference,
            params.ref_fai,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        CALL_SNV_WGS(
            params.reference,
            params.ref_fai,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        ANNOTATE_VCF(
            params.dbnsfp,
            params.dbnsfp_tbi,
            params.dbnsfp_dt,
            CALL_SNV_WGS.out.snv_vcf
        )
    }

    CALL_CNV(
        params.cnv_control,
        params.cnv_vcf_header,
        PICARD_MARK_DUPLICATES.out.md_bam,
        SAMTOOLS_INDEX_MD.out.index_md_bam
    )

    INDEX_CNV(
        CALL_CNV.out.cnv_vcf
    )

    ANNOTATE_CNV(
        params.hg19_genes,
        INDEX_CNV.out.cnv_index
    )

    CALL_EH(
        params.reference,
        params.ref_fai,
        PICARD_MARK_DUPLICATES.out.md_bam,
        SAMTOOLS_INDEX_MD.out.index_md_bam
    )

    vcf_concat_ch = ANNOTATE_VCF.out.sift_vcf
    vcf_concat_ch
        .mix(ANNOTATE_CNV.out.cnv_ann, CALL_EH.out.eh_vcf)
        .groupTuple()
        .set { vcf_concat_ch }

    MERGE_VCF(
        vcf_concat_ch
    )

    if (params.exome == "YES") {
        gvcf_concat_ch = CALL_SNV_WES.out.snv_gvcf
    }
    else {
        gvcf_concat_ch = CALL_SNV_WGS.out.snv_gvcf
    }

    gvcf_concat_ch
        .mix(CALL_CNV.out.cnv_gvcf, CALL_EH.out.eh_gvcf)
        .groupTuple()
        .set { gvcf_concat_ch }

    MERGE_GVCF(
        gvcf_concat_ch
    )

    qc_out_ch = ANNOTATE_VCF.out.snpeff_stats
        
    if (params.exome == "YES") {
        qc_out_ch.mix(PICARD_COLLECT_HS_METRICS.out.hs_metrics)
    }
    else {
        qc_out_ch.mix(PICARD_COLLECT_WGS_METRICS.out.wgs_metrics)
    }
    if (params.single_lane == "NO") {
         qc_out_ch.mix(TRIM_READS.out.trim_log, FASTQC.out.qc)
    }
    else {
        qc_out_ch.mix(TRIM_READS_SINGLE.out.trim_log, FASTQC_SINGLE.out.qc)
    }

    qc_out_ch
        .mix(
            FASTQC_TRIMMED.out.qc_trimmed,
            PICARD_MARK_DUPLICATES.out.md_metrics,
            PICARD_ESTIMATE_LIBRARY_COMPLEXITY.out.lib_complex_metrics
        )

    qc_out_ch
        .groupTuple()
        .set { qc_out_ch }

    MULTIQC_SAMPLE(
        qc_out_ch
    )
}   
