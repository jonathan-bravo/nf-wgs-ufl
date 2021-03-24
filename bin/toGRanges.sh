#!/usr/bin/env sh

sample_id="$1"

grep '^chr' ${sample_id}_cnv.vcf | \
awk '{print $3}' | \
awk -F':' '{print $2":"$3}' > ${sample_id}_cnv_calls.csv

grep '^chr' ${sample_id}_cnv.vcf | \
awk '{print $8}' | \
awk -F';' '{print $3}' | \
awk -F'=' '{print $2}' | \
cut -b 3 > ${sample_id}_cnv_class.csv

sed -i '$ d' ${sample_id}_cnv_calls.csv
sed -i '$ d' ${sample_id}_cnv_class.csv