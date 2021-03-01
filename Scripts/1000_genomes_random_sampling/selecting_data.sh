#!/bin/bash

###################
# Filtering 1000G #
###################

# Filtering the 1000g sample sheet by samples that are low-coverage WGS
grep -Fwf \
1000_genomes_phase3_samples.txt \
1000G_SampleSheetFull.csv > \
1000_genomes_filtered_phase3.csv

####################
# Population Files #
####################

# Split into separate files for each population
awk -F ',' '{print > ("1000g_"$3".csv")}' \
1000_genomes_filtered_phase3.csv

mkdir 1000_genomes_phase3_populations

mv 1000g* 1000_genomes_phase3_populations/

mkdir 1000_genomes_phase3_populations_subsamples

#############
# Subsample #
#############

# Get a subsample from each population
subsample -n 10 \
./1000_genomes_phase3_populations/1000g_AFR.csv > \
1000g_AFR_SubSample.csv

subsample -n 10 \
./1000_genomes_phase3_populations/1000g_AMR.csv > \
1000g_AMR_SubSample.csv

subsample -n 10 \
./1000_genomes_phase3_populations/1000g_EAS.csv > \
1000g_EAS_SubSample.csv

subsample -n 10 \
./1000_genomes_phase3_populations/1000g_EUR.csv > \
1000g_EUR_SubSample.csv

subsample -n 10 \
./1000_genomes_phase3_populations/1000g_SAS.csv > \
1000g_SAS_SubSample.csv

mv 1000g* 1000_genomes_phase3_populations_subsamples/

# Grabbing the sample names
awk -F ',' '{print $1}' \
1000_genomes_phase3_populations_subsamples/1000g_AFR_SubSample.csv > \
AFR_samples.txt

awk -F ',' '{print $1}' \
1000_genomes_phase3_populations_subsamples/1000g_AMR_SubSample.csv > \
AMR_samples.txt

awk -F ',' '{print $1}' \
1000_genomes_phase3_populations_subsamples/1000g_EAS_SubSample.csv > \
EAS_samples.txt

awk -F ',' '{print $1}' \
1000_genomes_phase3_populations_subsamples/1000g_EUR_SubSample.csv > \
EUR_samples.txt

awk -F ',' '{print $1}' \
1000_genomes_phase3_populations_subsamples/1000g_SAS_SubSample.csv > \
SAS_samples.txt

mkdir sample_lists

mv *_samples.txt sample_lists/