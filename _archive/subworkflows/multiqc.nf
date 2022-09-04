#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MULTIQC_SAMPLE } from '../modules/multiqc/multiqc_sample' addParams([*:params, "run_dir" : params.run_dir])


workflow MULTIQC {
    take:

    pre_fastqc
    trim_log
    fastqc
    picard_md_qc
    picard_qc
    vcf_stats

    main:

    if ( trim_log == null && pre_fastqc == null && fastqc == null && picard_md_qc == null && picard_qc == null && vcf_stats == null ) {

        if ( params.exome ) {

            picard_qc_ch = Channel.fromPath("${params.run_dir}/*/hs_metrics/*.txt")
                .map{ path -> tuple(path.baseName[0..-25], path) }
        }
        else {

            picard_qc_ch = Channel.fromPath("${params.run_dir}/*/wgs_metrics/*.txt")
                .map{ path -> tuple(path.baseName[0..-26], path) }
        }

        if ( params.match == '_{1,2}.fq.gz' ) {
            pre_fastqc_ch = Channel.fromPath("${params.run_dir}/*/fastqc_*[0-9]_logs/*")
                .map{ path -> tuple(path.baseName[0..-10], path) }
        }
        else {
            pre_fastqc_ch = Channel.fromPath("${params.run_dir}/*/fastqc_*[0-9]_logs/*")
                .map{ path -> tuple(path.baseName[0..-11], path) }
        }

        trim_log_ch = Channel.fromPath("${params.run_dir}/*/Trimmomatic/*.log")
            .map{ path -> tuple(path.baseName[0..-10], path) }

        fastqc_ch = Channel.fromPath("${params.run_dir}/*/fastqc_*_trimmed_logs/*")
            .map{ path -> tuple(path.baseName[0..-21], path) }

        picard_md_qc_ch = Channel.fromPath("${params.run_dir}/*/alignment/*.txt")
            .map{ path -> tuple(path.baseName[0..-12], path) }

        snpeff_qc_ch = Channel.fromPath("${params.run_dir}/*/snpEff/*.csv")
            .map{ path -> tuple(path.baseName[0..-14], path) }

        picard_qc_ch
            .mix(trim_log_ch, pre_fastqc_ch, fastqc_ch, picard_md_qc_ch, snpeff_qc_ch)
            .groupTuple()
            .set{ multiqc_ch }
    }
    else {

        pre_fastqc
            .mix( trim_log, fastqc ,picard_md_qc ,picard_qc ,vcf_stats )
            .groupTuple()
            .set{ multiqc_ch }
    }

    MULTIQC_SAMPLE(
        multiqc_ch
    )
}