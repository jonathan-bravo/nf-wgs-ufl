# TO DD LIST

## DOCUMENT

- TODO Fill out SNP scoring table
- TODO Need data on CNV detection
- 5000 bp windows (5x depth @ 150bp reads is ~ 166.67 reads per window)
- Windows have 500 bp overlap so that all windows after the first and last
  - Will only have 4000 bp of unique information
- I believe we will have a high level of detection for even smaller CNVs because we are using a simulated data set as a result we must filter
- How are we filtering the CNVs down to pathogenic hits and how are we sure?
> NOTE: For all pathogenic hits may want to include the copy number so we can double check if that is within a normal range.
- TODO Figure out what the max SV length detected by Strelka2 is
- What about DeepVariant? DeepVariant is bound to be better and is part of the plan for v2
- TODO Document the CNV results generated


## CNV COMPARE

- TODO Download agregate giab bams
- TODO generate CNV data for agregate giab
- TODO Compare 30x agregate samples to 30x giab samples
- TODO RUN 30X NIST SAMPLE


## PIPELINE

- TODO Run MultiQC after germline pipeline automatically
- FIXME Why did the GE sample stop before alignment?
- TODO Remake the Nextflow EC2 instance
- TODO Put SNPID into pipeline
- TODO Make sure that `classify_vcf.py` is run in AWS
  - we want all data to remain in AWS, the less files that have to come down the better
- TODO Build a simple python gui for running pipeline
- TODO Translate gui to web app that will be hosted in AWS S3
- TODO Write a script for transfering files from Illumina BaseSpace to our AWS S3 bucket
- TODO Database the correct output from the pipeline
- TODO Request an increase of the vcpu limit to 2500 so that we don't encounter an upper limit when we expand


### CLASSIFY VCF

- TODO Possible build a docker image to run this script
- How do we want to run this in AWS?


### SNPID

- TODO Pull specific SNPS from genomic VCF (using PGX bed files)
- TODO Make sure that GENE information is included in the pulled SNPs
- TODO Build nextflow module for pulling the SNPID


### DATABASE

- TODO Make a database of vairants from the samples we have run
- What do we want to do with the data? This will change which data we need. We might be able to get away with just a datalake, might need multiple databases.
- Will want to filter all results to passing and then simplify some of the fields. Once again this depends on what we want to do with the data.
- There is some worry about determination of homo/heterozygosity at lower coverages. This is one of the fields we will probably want to simplify.


## OTHER

TODO Show Matt the CNV plot
