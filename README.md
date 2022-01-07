# nf-wgs-ufl
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

This pipeline is used through the `pipeline_webpage` that is hosted in AWS through CloudFront. All data is in AWS S3 and all computation occurs in AWS. This is an ongoing project that will continue to evolve as the needs of our lab evolves.

Version 1.0.0 of the pipeline runs as follows:

![](https://drive.google.com/uc?id=13KB4RFRjMlIpkyO3jEtgJTPS9hawK_qo)
