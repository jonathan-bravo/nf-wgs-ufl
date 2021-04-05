#!/usr/bin/env nextflow

nextflow.enable.dsl        = 2

params.single_lane         = ""
params.exome               = ""
params.match               = ""
params.bucket              = ""
params.run_id              = ""
params.run_dir             = ""
params.panels_dir          = ""
params.ref_dir             = "${params.bucket}/Pipeline/Reference"
params.reference           = "${params.ref_dir}/hg19/hg19.fa"
params.bwa_amb             = "${params.reference}.amb"
params.bwa_ann             = "${params.reference}.ann"
params.bwa_bwt             = "${params.reference}.bwt"
params.bwa_pac             = "${params.reference}.pac"
params.bwa_sa              = "${params.reference}.sa"
params.ref_fai             = "${params.reference}.fai"
params.male_cnv_control    = "${params.ref_dir}/cnv/male_sim_control.RData"
params.female_cnv_control  = "${params.ref_dir}/cnv/female_sim_control.RData"
params.cnv_vcf_header      = "${params.ref_dir}/cnv/cnv_vcf_header.tsv"
params.bait                = "${params.ref_dir}/exome_targets/bait.interval_list"
params.target              = "${params.ref_dir}/exome_targets/target.interval_list"
params.dbnsfp              = "${params.ref_dir}/dbnsfp/dbNSFP4.1a.txt.gz"
params.dbnsfp_tbi          = "${params.dbnsfp}.tbi"
params.dbnsfp_dt           = "${params.dbnsfp}.data_types"
params.outdir              = "${params.bucket}/Pipeline_Output"
params.glue_dir            = "${params.bucket}/Pipeline_Output/_SampleTSV"

germline_params = [
    *:params,
    "single_lane"         : params.single_lane,
    "exome"               : params.exome,
    "match"               : params.match,
    "bucket"              : params.bucket,
    "run_id"              : params.run_id,
    "ref_dir"             : params.ref_dir,
    "reference"           : params.reference,
    "bwa_amb"             : params.bwa_amb,
    "bwa_ann"             : params.bwa_ann,
    "bwa_bwt"             : params.bwa_bwt,
    "bwa_pac"             : params.bwa_pac,
    "bwa_sa"              : params.bwa_sa,
    "ref_fai"             : params.ref_fai,
    "male_cnv_control"    : params.male_cnv_control,
    "female_cnv_control"  : params.female_cnv_control,
    "cnv_vcf_header"      : params.cnv_vcf_header,
    "bait"                : params.bait,
    "target"              : params.target,
    "dbnsfp"              : params.dbnsfp,
    "dbnsfp_tbi"          : params.dbnsfp_tbi,
    "dbnsfp_dt"           : params.dbnsfp_dt,
    "outdir"              : params.outdir
]

multiqc_params = [
    *:params,
    "run_id"  : params.run_id,
    "run_dir" : params.run_dir
]

// apply_panel_params = [
//     *:params,
//     "run_dir"    : params.run_dir,
//     "panels_dir" : params.panels_dir,
//     "bucket"     : params.bucket,
//     "glue_dir"   : params.glue_dir
// ]


workflow {
    include { GERMLINE } from './ufl-germline' addParams( germline_params: germline_params )
    include { MULTIQC  } from './ufl-multiqc'  addParams( multiqc_params:  multiqc_params  )
    GERMLINE ()
    MULTIQC ()
}
workflow.onComplete {}