#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MANTA_WGS      } from '../modules/manta/manta_wgs'                  addParams([*:params, "run_dir" : params.run_dir])
include { CALL_SNV_WGS   } from '../modules/strelka2/call_snv_wgs'            addParams([*:params, "run_dir" : params.run_dir])
include { CALL_SNV_WES   } from '../modules/strelka2/call_snv_wes'            addParams([*:params, "run_dir" : params.run_dir])
include { CALL_CNV       } from '../modules/cn_mops/call_cnv'                 addParams([*:params, "run_dir" : params.run_dir])
include { ANNOTATE_CNV   } from '../modules/ubuntu_python3/annotate_cnv'      addParams([*:params, "run_dir" : params.run_dir])
include { CALL_EH        } from '../modules/expansion_hunter/call_eh'         addParams([*:params, "run_dir" : params.run_dir])
include { ANNOTATE_EH    } from '../modules/ubuntu_python3/annotate_eh'       addParams([*:params, "run_dir" : params.run_dir])
include { MERGE_VCF      } from '../modules/bcftools_tabix/merge_vcf'         addParams([*:params, "run_dir" : params.run_dir])
include { MERGE_GVCF     } from '../modules/bcftools_tabix/merge_gvcf'        addParams([*:params, "run_dir" : params.run_dir])
include { ANNOTATE_VCF   } from '../modules/snpeff_tabix/annotate_vcf'        addParams([*:params, "run_dir" : params.run_dir])
include { GAUCHIAN       } from '../modules/gauchian/gauchian'                addParams([*:params, "run_dir" : params.run_dir])
include { CYRIUS         } from '../modules/cyrius/cyrius'                    addParams([*:params, "run_dir" : params.run_dir])
include { CALL_EHD       } from '../modules/expansion_hunter_denovo/call_ehd' addParams([*:params, "run_dir" : params.run_dir])
include { VCF_TO_PARQUET } from '../modules/ubuntu_python3/vcf_to_parquet'


// Setting up the sample bam files
params.sample_bams = "${params.run_dir}/*/alignment/*_md.bam"
params.sample_bais = "${params.run_dir}/*/alignment/*_md.bam.bai"


workflow VCF {
    take:

    bam //Channel: [ val(sample_id), bam, bai ]

    main:

    if ( bam == null ) {

        bam_ch = Channel.fromFilePairs([params.sample_bams, params.sample_bais]), flat: true
    }
    else {

        bam_ch = bam
    }

    bam_ch.view()

    if ( params.exome ) {

        CALL_SNV_WES(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            bam_ch
        )

        snv_vcf_ch  = CALL_SNV_WES.out.snv_vcf
        snv_gvcf_ch = CALL_SNV_WES.out.snv_gvcf
    }
    else {

        CALL_SNV_WGS(
            params.reference,
            params.ref_fai,
            params.ref_gzi,
            bam_ch
        )

        snv_vcf_ch  = CALL_SNV_WGS.out.snv_vcf
        snv_gvcf_ch = CALL_SNV_WGS.out.snv_gvcf
    }

    ANNOTATE_VCF(
        snv_vcf_ch
    )

    MANTA_WGS(
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )

    CALL_CNV(
        params.cnv_control,
        params.cnv_vcf_header,
        bam_ch
    )

    ANNOTATE_CNV(
        params.hs37d5_genes,
        CALL_CNV.out.cnv_vcf
    )

    CALL_EH(
        params.variant_catalog,
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )

    ANNOTATE_EH(
        params.variant_catalog,
        CALL_EH.out.eh_vcf
    )

    CALL_EHD(
        params.ehd_controls,
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )

    GAUCHIAN(
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )

    CYRIUS(
        params.reference,
        params.ref_fai,
        params.ref_gzi,
        bam_ch
    )

    ANNOTATE_VCF.out.snpeff_vcf
        .mix(ANNOTATE_CNV.out.cnv_ann, ANNOTATE_EH.out.eh_ann)
        .groupTuple()
        .set { vcf_concat_ch }

    MERGE_VCF(
        vcf_concat_ch
    )

    snv_gvcf_ch
        .mix(CALL_CNV.out.cnv_gvcf, CALL_EH.out.eh_gvcf)
        .groupTuple()
        .set { gvcf_concat_ch }

    MERGE_GVCF(
        gvcf_concat_ch
    )

    VCF_TO_PARQUET(
        MERGE_VCF.out.vcf
    )

    emit:

    vcf       = MERGE_VCF.out.vcf             // Channel: [ val(sample_id), concat_vcf, vcf_tbi ]
    vcf_stats = ANNOTATE_VCF.out.snpeff_stats // Channel: [ val(sample_id), snpeff_stats ]
}