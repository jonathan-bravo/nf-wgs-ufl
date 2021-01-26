#!/bin/bash

set -e

PREFIX="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/chromosomes"
CHROM=("1" "2" "3" "4" "5" "6" "7" "8" "9" "10" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "X" "Y" "M")

if [[ -e hg19.fa ]]; then
  echo "hg19.fa already exists, delete to redownload"
  exit 0
fi

for x in ${CHROM[@]}; do
  echo "wget ${PREFIX}/chr${x}.fa.gz -O chr${x}.fa.gz"
  wget ${PREFIX}//chr${x}.fa.gz -O chr${x}.fa.gz
  gunzip -c "chr${x}.fa.gz" >> hg19.fa
  rm "chr${x}.fa.gz" 
done

sed -i "s/>chr/>/g" hg19.fa

sed -i "s/>M/>MT/" hg19.fa

samtools faidx hg19.fa

https://www.ncbi.nlm.nih.gov/nuccore/NC_000001.10