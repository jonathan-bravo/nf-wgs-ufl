#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC_TRIMMED             } from '../modules/fastqc/fastqc_trimmed'         addParams([*:params, "run_dir" : params.run_dir])
include { TRIM_READS                 } from '../modules/trimmomatic/trim_reads'        addParams([*:params, "run_dir" : params.run_dir])
include { TRIM_READS_SINGLE          } from '../modules/trimmomatic/trim_reads_single' addParams([*:params, "run_dir" : params.run_dir])
include { SAMTOOLS_INDEX_MD          } from '../modules/samtools/index_md'             addParams([*:params, "run_dir" : params.run_dir])
include { PICARD_COLLECT_WGS_METRICS } from '../modules/picard/collect_wgs_metrics'    addParams([*:params, "run_dir" : params.run_dir])
include { PICARD_COLLECT_HS_METRICS  } from '../modules/picard/collect_hs_metrics'     addParams([*:params, "run_dir" : params.run_dir])
include { PICARD_MARK_DUPLICATES     } from '../modules/picard/mark_duplicates'        addParams([*:params, "run_dir" : params.run_dir])
include { CAT_TWO_LANES              } from '../modules/ubuntu_python3/cat_two_lanes'
include { CAT_FOUR_LANES             } from '../modules/ubuntu_python3/cat_four_lanes'
include { ALIGN_TRIMMED_READS        } from '../modules/bwa/align_trimmed_reads'
include { SAMTOOLS_VIEW              } from '../modules/samtools/view'
include { SAMTOOLS_SORT              } from '../modules/samtools/sort'
include { SAMTOOLS_INDEX             } from '../modules/samtools/index'

if ( params.one ) {

    reads_base = Channel.fromFilePairs( params.reads1 )
}
else if ( params.two ) {
    
    reads1_ch = Channel.fromFilePairs( params.reads1 )
    reads2_ch = Channel.fromFilePairs( params.reads2 )

    reads1_ch
        .combine( reads2_ch, by: 0 )
        .set { reads_base }
}
else if ( params.four ) {

    reads1_ch = Channel.fromFilePairs( params.reads1 )
    reads2_ch = Channel.fromFilePairs( params.reads2 )
    reads3_ch = Channel.fromFilePairs( params.reads3 )
    reads4_ch = Channel.fromFilePairs( params.reads4 )

    reads1_ch
        .combine( reads2_ch, by: 0 )
        .combine( reads3_ch, by: 0 )
        .combine( reads4_ch, by: 0 )
        .set { reads_base }
}


workflow ALIGN {

    take:

    reads

    main:

    if ( reads == null ) {

        if ( params.one ) {

            reads_ch = reads_base
        }
        else if ( params.two ) {

            CAT_TWO_LANES(
                reads_base,
            )

            reads_ch = CAT_TWO_LANES.out.read_pairs
        }
        else if ( params.four ) {

            CAT_FOUR_LANES(
                reads_ch
            )

            reads_ch = CAT_FOUR_LANES.out.read_pairs
        }
    }
    else {

        reads_ch = reads
    }

    if ( params.one ) {

        TRIM_READS_SINGLE(
            reads_ch,
            params.trim_adapters
        )

        trimmed_reads_ch = TRIM_READS_SINGLE.out.trimmed_paired_reads
        trim_log_ch = TRIM_READS_SINGLE.out.trim_log
    }
    else {

        TRIM_READS(
            reads_ch,
            params.trim_adapters
        )

        trimmed_reads_ch = TRIM_READS.out.trimmed_paired_reads
        trim_log_ch = TRIM_READS.out.trim_log
    }

    FASTQC_TRIMMED(
        trimmed_reads_ch
    )

    ALIGN_TRIMMED_READS(
        params.reference,
        params.bwa_amb,
        params.bwa_ann,
        params.bwa_bwt,
        params.bwa_pac,
        params.bwa_sa,
        trimmed_reads_ch
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

    SAMTOOLS_SORT.out.sort_bam
        .combine( SAMTOOLS_INDEX.out.index_sort_bam, by: 0 )
        .set { bam_ch }

    PICARD_MARK_DUPLICATES(
        bam_ch
    )

    SAMTOOLS_INDEX_MD(
        PICARD_MARK_DUPLICATES.out.md_bam
    )

    PICARD_MARK_DUPLICATES.out.md_bam
        .combine( SAMTOOLS_INDEX_MD.out.index_md_bam, by: 0 )
        .set { bam_md_ch }

    if ( params.exome ){

        PICARD_COLLECT_HS_METRICS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            params.target,
            params.bait,
            PICARD_MARK_DUPLICATES.out.md_bam,
            SAMTOOLS_INDEX_MD.out.index_md_bam
        )

        picard_qc_ch = PICARD_COLLECT_HS_METRICS.out.hs_metrics
    }
    else {

        PICARD_COLLECT_WGS_METRICS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            bam_md_ch
        )

        picard_qc_ch = PICARD_COLLECT_WGS_METRICS.out.wgs_metrics
    }

    emit:

    bam          = bam_md_ch                             // channel: [ val(sample_id), bam, bai ]
    trim_log     = trim_log_ch                           // channel: [ val(sample_id), trim_log ]
    fastqc       = FASTQC_TRIMMED.out.qc_trimmed         // channel: [ val(sample_id), trimmed_logs ]
    picard_md_qc = PICARD_MARK_DUPLICATES.out.md_metrics // channel: [ val(sample_id), md_metrics ]
    picard_qc    = picard_qc_ch                          // channel: [ val(sample_id), wgs_metrics ]
}