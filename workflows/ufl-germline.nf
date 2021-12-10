#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.germline_params = [:]

include { FASTQC                             } from '../modules/fastqc/fastqc'                     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_TRIMMED                     } from '../modules/fastqc/fastqc_trimmed'             addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { FASTQC_SINGLE                      } from '../modules/fastqc/fastqc_single'              addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS                         } from '../modules/trimmomatic/trim_reads'            addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { TRIM_READS_SINGLE                  } from '../modules/trimmomatic/trim_reads_single'     addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { SAMTOOLS_INDEX_MD                  } from '../modules/samtools/index_md'                 addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_WGS_METRICS         } from '../modules/picard/collect_wgs_metrics'        addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_COLLECT_HS_METRICS          } from '../modules/picard/collect_hs_metrics'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { PICARD_MARK_DUPLICATES             } from '../modules/picard/mark_duplicates'            addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WGS                       } from '../modules/strelka2/call_snv_wgs'             addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_SNV_WES                       } from '../modules/strelka2/call_snv_wes'             addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_CNV                           } from '../modules/cn_mops/call_cnv'                  addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { INDEX_CNV                          } from '../modules/bcftools_tabix/index_cnv'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_CNV                       } from '../modules/ubuntu_python3/annotate_cnv'       addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_EH                            } from '../modules/expansion_hunter/call_eh'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CALL_EH_RESEARCH                   } from '../modules/expansion_hunter/call_eh_research' addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_VCF                          } from '../modules/bcftools_tabix/merge_vcf'          addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MERGE_GVCF                         } from '../modules/bcftools_tabix/merge_gvcf'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { ANNOTATE_VCF                       } from '../modules/snpeff_tabix/annotate_vcf'         addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { MULTIQC_SAMPLE                     } from '../modules/multiqc/multiqc_sample'            addParams([*:params, "outdir" : params.outdir, "run_id" : params.run_id])
include { CAT_TWO_LANES                      } from '../modules/ubuntu_python3/cat_two_lanes'
include { CAT_FOUR_LANES                     } from '../modules/ubuntu_python3/cat_four_lanes'
include { ALIGN_TRIMMED_READS                } from '../modules/bwa/align_trimmed_reads'
include { SAMTOOLS_VIEW                      } from '../modules/samtools/view'
include { SAMTOOLS_SORT                      } from '../modules/samtools/sort'
include { SAMTOOLS_INDEX                     } from '../modules/samtools/index'


if (params.lanes == "ONE") {
    if (params.exome == "YES"){
        params.reads_path = "${params.bucket}/Exome_Fastqs/${params.run_id}*${params.match}"
    }
    else{
        params.reads_path = "${params.bucket}/Fastqs/${params.run_id}*${params.match}"
    }

    reads_ch  = Channel.fromFilePairs(params.reads_path)
}
else if (params.lanes == "TWO") {
    if (params.exome == "YES"){
        params.reads1 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L001_{R1,R2}_001.fastq.gz"
        params.reads2 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L002_{R1,R2}_001.fastq.gz"
    }
    else {
        params.reads1 = "${params.bucket}/Fastqs/${params.run_id}*_L001_{R1,R2}_001.fastq.gz"
        params.reads2 = "${params.bucket}/Fastqs/${params.run_id}*_L002_{R1,R2}_001.fastq.gz"
    }

    reads1_ch = Channel.fromFilePairs(params.reads1)
    reads2_ch = Channel.fromFilePairs(params.reads2)
}
else if (params.lanes == "FOUR") {
    if (params.exome == "YES"){
        params.reads1 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L001_{R1,R2}_001.fastq.gz"
        params.reads2 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L002_{R1,R2}_001.fastq.gz"
        params.reads3 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L003_{R1,R2}_001.fastq.gz"
        params.reads4 = "${params.bucket}/Exome_Fastqs/${params.run_id}*_L004_{R1,R2}_001.fastq.gz"
    }
    else {
        params.reads1 = "${params.bucket}/Fastqs/${params.run_id}*_L001_{R1,R2}_001.fastq.gz"
        params.reads2 = "${params.bucket}/Fastqs/${params.run_id}*_L002_{R1,R2}_001.fastq.gz"
        params.reads3 = "${params.bucket}/Fastqs/${params.run_id}*_L003_{R1,R2}_001.fastq.gz"
        params.reads4 = "${params.bucket}/Fastqs/${params.run_id}*_L004_{R1,R2}_001.fastq.gz"
    }

    reads1_ch = Channel.fromFilePairs(params.reads1)
    reads2_ch = Channel.fromFilePairs(params.reads2)
    reads3_ch = Channel.fromFilePairs(params.reads3)
    reads4_ch = Channel.fromFilePairs(params.reads4)
}



workflow GERMLINE {
    if (params.lanes == "ONE") {
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
    else if (params.lanes == "TWO") {
        CAT_TWO_LANES(
            reads1_ch,
            reads2_ch
        )

        FASTQC(
            CAT_TWO_LANES.out.read_pairs
        )

        TRIM_READS(
            CAT_TWO_LANES.out.read_pairs,
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
    else if (params.lanes == "FOUR") {
        CAT_FOUR_LANES(
            reads1_ch,
            reads2_ch,
            reads3_ch,
            reads4_ch,
        )

        FASTQC(
            CAT_FOUR_LANES.out.read_pairs
        )

        TRIM_READS(
            CAT_FOUR_LANES.out.read_pairs,
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

    if (params.exome == "YES"){
        PICARD_COLLECT_HS_METRICS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            params.target,
            params.bait,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        CALL_SNV_WES(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        ANNOTATE_VCF(
            CALL_SNV_WES.out.snv_vcf
        )
    }
    else {
        PICARD_COLLECT_WGS_METRICS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        CALL_SNV_WGS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        ANNOTATE_VCF(
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
        params.hs37d5_genes,
        INDEX_CNV.out.cnv_index
    )

    CALL_EH(
        params.variant_catalog,
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        PICARD_MARK_DUPLICATES.out.md_bam,
        SAMTOOLS_INDEX_MD.out.index_md_bam
    )

    CALL_EH_RESEARCH(
        params.research_variant_catalog,
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        PICARD_MARK_DUPLICATES.out.md_bam,
        SAMTOOLS_INDEX_MD.out.index_md_bam
    )

    vcf_concat_ch = ANNOTATE_VCF.out.snpeff_vcf
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
        picard_qc_ch = PICARD_COLLECT_HS_METRICS.out.hs_metrics
    }
    else {
        picard_qc_ch = PICARD_COLLECT_WGS_METRICS.out.wgs_metrics
    }
    if (params.lanes == "ONE") {
        trim_qc_ch = TRIM_READS_SINGLE.out.trim_log
        fastqc_qc_ch = FASTQC_SINGLE.out.qc
    }
    else {
        trim_qc_ch = TRIM_READS.out.trim_log
        fastqc_qc_ch = FASTQC.out.qc
    }

    qc_out_ch
        .mix(
            picard_qc_ch,
            trim_qc_ch,
            fastqc_qc_ch,
            FASTQC_TRIMMED.out.qc_trimmed,
            PICARD_MARK_DUPLICATES.out.md_metrics
        )
        .groupTuple()
        .set { qc_out_ch }

    MULTIQC_SAMPLE(
        qc_out_ch
    )
}   
