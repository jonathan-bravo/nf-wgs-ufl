# To Do List

## Documentation

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

## General Pipeline

- TODO Get files to automatically upload from NovaSeq to AWS
- TODO Put SNPID into pipeline
    - TODO Pull specific SNPS from GVCF (using PGX bed files)
    - TODO Make sure that GENE information is included in the pulled SNPs
    - TODO Build nextflow module for pulling the SNPID
- TODO Create clinical bucket/ environment
- TODO Database the correct output from the pipeline
    - Frequency of variant startified by ethnicity
    - Position of variant
    - Zygosity of variant
    - Counts of the variants (tied into frequency?)
    - Affected vs Unaffected
- TODO Create cloud launch/ public version of pipeline

## Website

- TODO Check reporting module multiselect
    - When samples above others aren't checked then the parent selector doesn't work
- TODO Animate radio and check boxes
    - Checkbox check isn't showing up
- TODO Fix names in 'Reporting' module
    - The names of things in the 'Request' module have the last letter cut off
- TODO Fix back buttons for 'Request' module
- Potentialy add a VCF browser?
- Potentially add a GVCF database browser?
- TODO Finish the 'Request' module
  - [ ] make 'variants' file read as 'VCF'
- TODO Make login session expire after ~ 15 min of inactivity
- TODO Fix multi-lane back button
- TODO Grab Cgnito role from JWT Token
- TODO Look at pagify (JS) for all s3 queried lists

# varClass.py

- TODO Change the 'CNV Interactions' tab name to 'Compund Variants'
- TODO Make the 'Compound Variants' tab always have the SNP, SV, EXP, and CNV columns
- TODO Fix 'Supporting Lit' tab name spelling error
- TODO Include additional items in the 'metadata.tools' part of the report
    - [ ] Include clinvar version
    - [ ] Include genome build
    - [ ] Include Nextflow version
    - [ ] Include gnomAD version
    - [ ] All database versions included in current version of dbNSFP?
- TODO Update Expansions tab
    - [ ] Include ref allele from expansion hunter
    - [ ] Alt:Ref -> Ref Length
    - [ ] Include "Depths"
- TODO Include gnomAD SV frequencies in variant determination and report
- TODO Include clinvar data in variant determination and report
- TODO Make all titles start with a cap
