# Pipeline Tasks

### ToDo
 
- Split the WES and WGS pipelines
- Remove the MultiQC pipeline
- Add GRCh38.p13 reference genome
- MAYBE add ExpansionHunterDenovo
- MAYBE add DeepVariant

### In Progress

- Add a log to the `bs-to-aws.sh` script so that if any fastqs are skipped there is a record of it
- Add Gauchian to pipeline for GBA
- Add Cyrius to pipeline for CYP2D6
- Adding Nirvana to pipeline (in tandem with snpeff for now)
- Update README to include a 'methods like' section

### Done ✓

- Add Manta to the germline pipeline
  - For more robust structural variant calling and CNV validation
- Fix report and sample selection issue in web page
