#!/usr/bin/env nextflow

nextflow.enable.dsl   = 2

params.apply_panel_params = [:]

include { APPLY_PANEL } from './modules/ubuntu_python3/apply_panel' addParams([*:params, "run_dir" : params.run_dir])
include { CREATE_TSV  } from './modules/ubuntu_python3/create_tsv'  addParams([*:params, "glue_dir" : params.glue_dir])
include { PARSE_INFO  } from './modules/ubuntu_python3/parse_info'


pairs_ch = Channel.from()

workflow APPLY_PANELS {
    APPLY_PANEL(
        pairs_ch,
        params.panels_dir,
        params.run_dir
    )

    PARSE_INFO(
        APPLE_PANEL.out.panel
    )
    
    CREATE_TSV(
        APPLE_PANEL.out.panel,
        PARSE_INFO.out.data
    )
}
