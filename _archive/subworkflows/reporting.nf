#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { GENERAL_REPORTING } from '../modules/ubuntu_python3/general_reporting' addParams([*:params, "run_dir" : params.run_dir])
include { PANEL_REPORTING   } from '../modules/ubuntu_python3/panel_reporting'   addParams([*:params, "run_dir" : params.run_dir])

workflow REPORTING {
    take:

    vcf // Channel: [ val(sample_id), vcf, tbi ]

    main:

    if ( vcf == null ) {

        def sample_panel_pairs = params.panel_reporting_list.split(',')

        panel_ch = Channel
            .from( sample_panel_pairs )
            .collate( 2 )
            .branch{
                general: it[1] == 'general'
                panel  : true
            }
        
        panel_ch
            .panel
            .map{ sample_id, panel -> tuple( sample_id, "${params.run_dir}/${sample_id}/variants/${sample_id}_concat.vcf.gz", "${params.panel_dir}/${panel}", "${params.run_dir}/${sample_id}/variants/${sample_id}_concat.vcf.gz.tbi" ) }
            .set{ vcf_panel_ch }

        panel_ch
            .general
            .map{ sample_id, panel -> tuple( sample_id, "${params.run_dir}/${sample_id}/variants/${sample_id}_concat.vcf.gz", "${params.run_dir}/${sample_id}/variants/${sample_id}_concat.vcf.gz.tbi" ) }
            .set{ vcf_ch }

        GENERAL_REPORTING(
            params.reporting_dir,
            vcf_ch
        )

        PANEL_REPORTING(
            params.reporting_dir,
            vcf_panel_ch
        )
    }
    else {

        GENERAL_REPORTING(
            params.reporting_dir,
            vcf
        )
    }

    
}