#!/usr/bin/env sh

sample_id="$1"

grep -v '^#' ${sample_id}_filtered_cnv.vcf | \
awk '{print $3}' | \
awk -F':' '{print $2":"$3}' > ${sample_id}_cnv_calls.csv

grep -v '^#' ${sample_id}_filtered_cnv.vcf | \
awk '{print $8}' | \
awk -F';' '{print $4}' | \
awk -F'=' '{print $2}' | \
cut -b 3 > ${sample_id}_cnv_class.csv

sed -i '${/^$/d;}' ${sample_id}_cnv_calls.csv
sed -i '${/^$/d;}' ${sample_id}_cnv_class.csv