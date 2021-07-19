#!/usr/bin/env nextflow

nextflow.enable.dsl             = 2

params.pipeline                 = ""
params.lanes                    = ""
params.exome                    = ""
params.match                    = ""
params.bucket                   = ""
params.run_id                   = ""
params.run_dir                  = "${params.bucket}/Pipeline_Output/${params.run_id}"
params.ref_dir                  = "${params.bucket}/Pipeline/Reference"
params.trim_adapters            = "${params.ref_dir}/trim/NEBNext.fa"
params.reference                = "${params.ref_dir}/hg19/hg19.fa"
params.bwa_amb                  = "${params.reference}.amb"
params.bwa_ann                  = "${params.reference}.ann"
params.bwa_bwt                  = "${params.reference}.bwt"
params.bwa_pac                  = "${params.reference}.pac"
params.bwa_sa                   = "${params.reference}.sa"
params.ref_fai                  = "${params.reference}.fai"
params.hg19_genes               = "${params.ref_dir}/hg19/hg19_genes.bed"
params.cnv_control              = "${params.ref_dir}/cnv/wgs_cnv_controls.RData"
params.cnv_vcf_header           = "${params.ref_dir}/cnv/cnv_vcf_header.tsv"
params.bait                     = "${params.ref_dir}/exome_targets/bait.interval_list"
params.target                   = "${params.ref_dir}/exome_targets/target.interval_list"
params.variant_catalog          = "${params.ref_dir}/expansion_hunter/variant_catalog.json"
params.research_variant_catalog = "${params.ref_dir}/expansion_hunter/variant_catalog_research.json"
params.outdir                   = "${params.bucket}/Pipeline_Output"
params.glue_dir                 = "${params.bucket}/Pipeline_Output/_SampleTSV"

germline_params = [
    *:params,
    "lanes"                    : params.lanes,
    "exome"                    : params.exome,
    "match"                    : params.match,
    "bucket"                   : params.bucket,
    "run_id"                   : params.run_id,
    "ref_dir"                  : params.ref_dir,
    "trim_adapters"            : params.trim_adapters,
    "reference"                : params.reference,
    "bwa_amb"                  : params.bwa_amb,
    "bwa_ann"                  : params.bwa_ann,
    "bwa_bwt"                  : params.bwa_bwt,
    "bwa_pac"                  : params.bwa_pac,
    "bwa_sa"                   : params.bwa_sa,
    "ref_fai"                  : params.ref_fai,
    "hg19_genes"               : params.hg19_genes,
    "cnv_control"              : params.cnv_control,
    "cnv_vcf_header"           : params.cnv_vcf_header,
    "bait"                     : params.bait,
    "target"                   : params.target,
    "variant_catalog"          : params.variant_catalog,
    "research_variant_catalog" : params.research_variant_catalog
    "outdir"                   : params.outdir
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

// Add Fastqc step
// Add option for single sample runs