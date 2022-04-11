# nf-wgs-ufl

Additional configuration is included in the sub_config directory. Files in this
directory are not included in the repo as they would/do contain more sensative
account information for running the pipelines. The only file present here
currently is `aws_id.config`. If this pipeline is to be run by someone else
then they should create this file and set up something that looks like this.

We have different labels for each docker container as well as different labels
for each compute environment/queue. This allowed us to optomise cost without
using SPOT instancing. You could create fewer queues, up the `errorStrategy`
number and potentially use SPOT instancing if desired.

```
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    HaKmonkey/nf-wgs-ufl additional Nextflow config file
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Sensative AWS config options for all compute environments
----------------------------------------------------------------------------------------
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

## Tools

This pipeline is a compilation of many open source tools:

- [TrimmomaticPE 0.39](https://doi.org/10.1093/bioinformatics/btu170)
- [BWA MEM 0.7.17-r1188](https://doi.org/10.1093/bioinformatics/btp324)
- [Samtools 1.10](https://doi.org/10.1093/gigascience/giab008)
- [Strelka2 2.9.10](https://doi.org/10.1038/s41592-018-0051-x)
- [cn.MOPS 1.36.0](https://doi.org/10.1093/nar/gks003)
- [ExpansionHunter 4.0.2](https://doi.org/10.1186/s13059-020-02017-z)
- [Manta 1.6.0](https://doi.org/10.1093/bioinformatics/btv710)
- [snpEff 5.0c](https://doi.org/10.4161/fly.19695)
- [Bcftools 1.10.2](https://doi.org/10.1093/gigascience/giab008)
- [Picard 2.23.8](http://broadinstitute.github.io/picard/)
- [FastQC 0.11.9](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
- [MultiQC 1.9](https://doi.org/10.1093/bioinformatics/btw354)
- [Nextflow 20.10](https://doi.org/10.1038/nbt.3820)
- [Genome Build hs37d5](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence)
- [ClassifyCNV 1.1.1](https://doi.org/10.1038/s41598-020-76425-3)
- [CADD v1.6](https://cadd.gs.washington.edu/)
- [REVEL 1.0](http://dx.doi.org/10.1016/j.ajhg.2016.08.016)
- [PMC](https://www.ncbi.nlm.nih.gov/pmc/)
- [OMIM](https://www.omim.org/)
- [gnomAD 2.1.1](https://doi.org/10.1038/s41586-020-2308-7)
- [ClinVar GRcH37_2021-04-18](https://doi.org/10.1093/nar/gkx1153)

## Workflow

This pipeline is used through the `pipeline_webpage` that is hosted in AWS
through CloudFront. All data is in AWS S3 and all computation occurs in AWS.
This is an ongoing project that will continue to evolve as the needs of our
lab evolves.

Version 1.0.0 of the pipeline runs as follows:

![](https://drive.google.com/uc?id=13KB4RFRjMlIpkyO3jEtgJTPS9hawK_qo)
