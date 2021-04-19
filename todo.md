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

## PIPELINE

- TODO Get files to automatically upload from NovaSeq to AWS
- TODO Put SNPID into pipeline
  - TODO Pull specific SNPS from GVCF (using PGX bed files)
  - TODO Make sure that GENE information is included in the pulled SNPs
  - TODO Build nextflow module for pulling the SNPID
- TODO Make sure that `classify_vcf.py` is run in AWS (runnabel image)
  - TODO Also need to hard code the paths to ClassifyCNV.py and the PMC_id.csv and remove them as options
  - We want all data to remain in AWS, the less files that have to come down the better
  - Image will plot the CNVs
  - Output file will be 'sample_id_panel_report.json'
  - TODO add excel output
- TODO add varClass to website
  - TODO allow users to select multiple panels
  - TODO send files to user
- TODO Create clinical bucket/ environment
- TODO Database the correct output from the pipeline
  - Frequency of variant startified by ethnicity
  - Position of variant
  - Zygosity of variant
  - Counts of the variants (tied into frequency?)
  - Affected vs Unaffected
- TODO Create cloud launch/ public version of pipeline
- TODO Request an increase of the vcpu limit to 2500 so that we don't encounter an upper limit when we expand
