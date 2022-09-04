# Usage
## Work in progress
This pipeline is a work in progress. Currently the focus is on running the workflow in AWS.

&nbsp;

## AWS config
Additional configuration is included in the sub_config directory. Files other
than the local config in this directory are not included in the repo as they
would/do contain more sensative account information for running the pipelines.
The only additional file present here currently is `aws_id.config`. If this
pipeline is to be run in an AWS environment then a config file should be created
for doing so. Below is a skeleton of what our file looks like with specific
information stripped.

We have different labels for each docker container as well as different labels
for each compute environment/queue. This allowed us to optomise cost without
using SPOT instancing. You could create fewer queues, up the `errorStrategy`
number and potentially use SPOT instancing if desired.

```
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HaKmonkey/nf-wgs-ufl additional Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sensative AWS config options for all compute environments
-------------------------------------------------------------------------------
*/

plugins {
    id 'nf-amazon'
}

process {
    executor = 'awsbatch'
    errorStrategy = { task.attempt <= 2 ? 'retry' : 'finish' }
    maxRetries = 2

    withLabel: small_process {
        cpus = 2
        memory = 2.GB
        queue = ''
    }

    withLabel: medium_process {
        cpus = 8
        memory = 14.GB
        queue = ''
    }

    withLabel: high_mem {
        cpus = 8
        memory = 62.GB
        queue = ''
    }

    withLabel: alignment {
        cpus = 32
        memory = 62.GB
        queue = ''
    }

    withLabel: reporting {
        cpus = 16
        memory = 120.GB
        queue = ''
    }

    withLabel: ubuntu_python3 {
        container = ""
    }

    withLabel: bcftools_tabix {
        container = ''
    }

    withLabel: bwa {
        container = ''
    }

    withLabel: expansion_hunter {
        container = ''
    }

    withLabel: fastqc {
        container = ''
    }

    withLabel: cn_mops {
        container = ''
    }

    withLabel: picard {
        container = ''
    }

    withLabel: samtools {
        container = ''
    }

    withLabel: snpeff_tabix {
        container = ''
    }

    withLabel: strelka2 {
        container = ''
    }

    withLabel: trimmomatic {
        container = ''
    }

    withLabel: multiqc {
        container = ''
    }

    withLabel: manta {
        container = ''
    }

    withLabel: gauchian {
        container = ''
    }

    withLabel: cyrius {
        container = ''
    }
}

aws {
    region = ''

    client {
        maxErrorRetry = 4
    }

    batch {
        maxTransferAttempts = 8
        delayBetweenAttempts = 300
        maxParallelTransfers = 5
    }
}
```

You can also set the AWS access key, secret key, and default region here, but
we currently have these values set up in AWS secrets manager.

&nbsp;

## AWS bucket setup
Our AWS bucket is currently set up with the following directory tree:

```
pipeline-bucket/
  ├ CNV_Controls/
  ⎮   ├ Fastqs/
  ⎮   ⎮   ├ _Processed/
  ⎮   ⎮   ⎮   ├ 1_1.fastq.gz
  ⎮   ⎮   ⎮   ├ 1_2.fastq.gz
  ⎮   ⎮   ⎮   └ ...
  ⎮   ⎮   ├ 2_1.fastq.gz
  ⎮   ⎮   ├ 2_2.fastq.gz
  ⎮   ⎮   └ ...
  ⎮   └ RData/
  ⎮       ├ 1_cnv_control.RData
  ⎮       ├ 2_cnv_control.RData
  ⎮       └ ...
  ├ Exome_Fastqs/
  ⎮   ├ _Processed/
  ⎮   ⎮   ├ 1_1.fastq.gz
  ⎮   ⎮   ├ 1_2.fastq.gz
  ⎮   ⎮   └ ...
  ⎮   ├ 2_1.fastq.gz
  ⎮   ├ 2_2.fastq.gz
  ⎮   └ ...
  ├ Fastqs/
  ⎮   ├ _Processed/
  ⎮   ⎮   ├ 1_1.fastq.gz
  ⎮   ⎮   ├ 1_2.fastq.gz
  ⎮   ⎮   └ ...
  ⎮   ├ 2_1.fastq.gz
  ⎮   ├ 2_2.fastq.gz
  ⎮   └ ...
  ├ Pipeline_Output/
  ⎮   ├ run-1/
  ⎮   ⎮   ├ MultiQC/
  ⎮   ⎮   ⎮   └ ...
  ⎮   ⎮   ├ sample1/
  ⎮   ⎮   ⎮   └ ...
  ⎮   ⎮   ├ sample2/
  ⎮   ⎮   ⎮   └ ...
  ⎮   ⎮   └ ...
  ⎮   ├ run-2/
  ⎮   ⎮   └ ...
  ⎮   └ ...
  └ Pipeline/
      ├ Reference/
      ⎮   ├ cnv/
      ⎮   ⎮   ├ cnv_vcf_header.tsv
      ⎮   ⎮   └ wgs_cnv_controls.RData
      ⎮   ├ exome_targets/
      ⎮   ⎮   ├ bait.interval_list
      ⎮   ⎮   ├ target.interval_list
      ⎮   ⎮   ├ xgen-exome-research-panel-v2-probes-hg19.bed
      ⎮   ⎮   └ xgen-exome-research-panel-v2-targets-hg19.bed
      ⎮   ├ expansion_hunter_denovo/
      ⎮   ⎮   ├ manifest.tsv
      ⎮   ⎮   ├ control1.str_profile.json
      ⎮   ⎮   ├ control2.str_profile.json
      ⎮   ⎮   └ ...
      ⎮   ├ expansion_hunter/
      ⎮   ⎮   └ variant_catalog.json
      ⎮   ├ hs37d5/
      ⎮   ⎮   ├ hs37d5_genes.bed
      ⎮   ⎮   ├ hs37d5.fa.gz
      ⎮   ⎮   ├ hs37d5.fa.gz.amb
      ⎮   ⎮   ├ hs37d5.fa.gz.ann
      ⎮   ⎮   ├ hs37d5.fa.gz.bwt
      ⎮   ⎮   ├ hs37d5.fa.gz.fai
      ⎮   ⎮   ├ hs37d5.fa.gz.gzi
      ⎮   ⎮   ├ hs37d5.fa.gz.pac
      ⎮   ⎮   └ hs37d5.fa.gz.sa
      ⎮   ├ panels/
      ⎮   ⎮   ├ gene_panel_1
      ⎮   ⎮   ├ gene_panel_2
      ⎮   ⎮   └ ...
      ⎮   └ trim/
      ⎮      └ NEBNext.fa
      └ Reporting/
          ├ cadd/
          ⎮   ├ cadd.chr1.0.pkl
          ⎮   ├ cadd.chr1.1.pkl
          ⎮   ├ ...
          ⎮   ├ cadd.chr2.0.pkl
          ⎮   ├ cadd.chr2.1.pkl
          ⎮   ├ ...
          ⎮   └ cadd.chr22.19.pkl
          ├ clinvar/
          ⎮   ├ clinvar.chr1.pkl
          ⎮   ├ clinvar.chr2.pkl
          ⎮   ├ ...
          ⎮   └ clinvar.chrY.pkl
          ├ gnomad/
          ⎮   ├ gnomad.chr1.pkl
          ⎮   ├ gnomad.chr2.pkl
          ⎮   ├ ...
          ⎮   └ gnomad.chrY.pkl
          ├ omim/
          ⎮   └ omim_2_gene.tsv
          ├ pmc/
          ⎮   └ PMC-ids.csv
          └ revel/
              ├ revel.chr1.pkl
              ├ revel.chr2.pkl
              ├ ...
              └ revel.chrY.pkl
```

The `Pipeline_Ouput/` directory will be explaned in [Output](output.md).

&nbsp;

### Pipeline
***
The pipeline directory contains various controls or input files for the
processes in the pipeline. These are organized into different directories so,
later, I can identify which tools they belong to and add revisions, or clean
up old/ unused data.

### Reporting
***
The reporting directory contains sub-directories for each of the annotation
databases we are using during report generation.

The pickle files in the various directories were created from the vcf files
provided on the data home pages. The pickle files are the data split by
chromosome and turned into python dictionaries. The CADD data was split further
into 20 chunks each, for each chromosome because the data was so large.

I am still looking at potentially more efficient ways of applying the
annotations potentially using pyspark and the same data sources, but in parquet
format, or by creating a variant database.

### Fastqs and Exome Fastqs
***
These directories contain the fastq files split into lanes. Usually you will
see your reads split into 2 files, but because we are getting our data from
[Illumina's BaseSpace](https://www.illumina.com/products/by-type/informatics-products/basespace-sequence-hub.html)
and we use a larger flowcell, we get our files in pairs across 4 lanes.

Once the reads are processed through the pipeline they are archived in the
`_Processed` directory. Here the files are changed to the Glacier S3 instant
retrevial type.

### CNV Controls
***
This directory contains fastqs specifically selected to server as controls for
`cn.mops`. We generate the controls separately because this allows us to
pass a much smaller RData file into the module during the pipeline run.

Also here is a `_Processed` directory that functions similarly to the one in
the fastq data directories. In the event that new cnv controls are added the
control RData can be generated and then merged into the
`wgs_cnv_controls.RData`.

## Using the gatorgenome tkinter app
