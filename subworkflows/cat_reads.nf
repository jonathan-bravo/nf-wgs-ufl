#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { FASTQC_SINGLE  } from '../modules/fastqc/fastqc_single' addParams([*:params, "run_dir" : params.run_dir])
include { FASTQC         } from '../modules/fastqc/fastqc'        addParams([*:params, "run_dir" : params.run_dir])
include { CAT_TWO_LANES  } from '../modules/ubuntu_python3/cat_two_lanes'
include { CAT_FOUR_LANES } from '../modules/ubuntu_python3/cat_four_lanes'

workflow CAT_READS {

    main:

    if ( params.one ) {

        reads_ch = Channel.fromFilePairs(params.reads1)

        FASTQC_SINGLE(
            reads_ch
        )

        cat_reads_ch = reads_ch
        pre_fastqc_ch = FASTQC_SINGLE.out.qc
    }
    else if ( params.two ) {

        reads1_ch = Channel.fromFilePairs( params.reads1 )
        reads2_ch = Channel.fromFilePairs( params.reads2 )

        reads1_ch
            .combine( reads2_ch, by: 0 )
            .set { reads_ch }

        CAT_TWO_LANES(
            reads_ch
        )

        FASTQC(
            CAT_TWO_LANES.out.read_pairs
        )

        cat_reads_ch = CAT_TWO_LANES.out.read_pairs
        pre_fastqc_ch = FASTQC.out.qc
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
            .set { reads_ch }

        CAT_FOUR_LANES(
            reads_ch
        )

        FASTQC(
            CAT_FOUR_LANES.out.read_pairs
        )

        cat_reads_ch = CAT_FOUR_LANES.out.read_pairs
        pre_fastqc_ch = FASTQC.out.qc
    }

    emit:

    reads      = cat_reads_ch  // Channel: [ val(sample_id), [ reads1, reads2 ] ]
    pre_fastqc = pre_fastqc_ch // channel: [ val(sample_id), logs ]
}