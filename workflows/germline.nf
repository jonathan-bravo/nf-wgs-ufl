#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


include { CAT_READS } from '../subworkflows/cat_reads'
include { ALIGN     } from '../subworkflows/align'
include { VCF       } from '../subworkflows/vcf'
include { MULTIQC   } from '../subworkflows/multiqc'
include { REPORTING } from '../subworkflows/reporting'


workflow GERMLINE {
    if ( params.germline ) {

        CAT_READS()

        ALIGN(
            CAT_READS.out.reads
        )

        VCF(
            ALIGN.out.bam
        )

        MULTIQC(
            CAT_READS.out.pre_fastqc,
            ALIGN.out.trim_log,
            ALIGN.out.fastqc,
            ALIGN.out.picard_md_qc,
            ALIGN.out.picard_qc,
            VCF.out.vcf_stats
        )

        REPORTING(
            VCF.out.vcf
        )
    }
    else if ( params.cat_reads ) {

        CAT_READS()
    }
    else if ( params.align ) {

        ALIGN(
            null
        )
    }
    else if ( params.vcf ) {

        VCF(
            null
        )
    }
    else if ( params.multiqc ) {

        MULTIQC(
            null,
            null,
            null,
            null,
            null,
            null
        )
    }
    else if ( params.report ) {

        REPORTING(
            null
        )
    }
}