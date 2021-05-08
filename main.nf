#!/usr/bin/env nextflow

nextflow.enable.dsl        = 2

params.pipeline            = ""
params.single_lane         = ""
params.exome               = ""
params.match               = ""
params.bucket              = ""
params.run_id              = ""
params.run_dir             = "${params.bucket}/Pipeline_Output/${params.run_id}"
params.ref_dir             = "${params.bucket}/Pipeline/Reference"
params.reference           = "${params.ref_dir}/hg19/hg19.fa"
params.bwa_amb             = "${params.reference}.amb"
params.bwa_ann             = "${params.reference}.ann"
params.bwa_bwt             = "${params.reference}.bwt"
params.bwa_pac             = "${params.reference}.pac"
params.bwa_sa              = "${params.reference}.sa"
params.ref_fai             = "${params.reference}.fai"
params.hg19_genes          = "${params.ref_dir}/hg19/hg19_genes.bed"
params.cnv_control         = "${params.ref_dir}/cnv/wgs_cnv_controls.RData"
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
    "hg19_genes"          : params.hg19_genes,
    "cnv_control"         : params.cnv_control,
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

workflow {
    if (params.pipeline == "GERMLINE") {
        include { GERMLINE } from './ufl-germline' addParams( germline_params: germline_params )
        GERMLINE ()
    }
    if (params.pipeline == "MULTIQC") {
        include { MULTIQC  } from './ufl-multiqc'  addParams( multiqc_params:  multiqc_params  )
        MULTIQC ()
    }
}