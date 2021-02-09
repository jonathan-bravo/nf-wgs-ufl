#!/usr/bin/env sh

sample_id="$1"
panel="$2"
panel_dir="$3"
sample_path="$4"

GENES=$(awk -F'\t' 'NR>1 {print $1}' ${panel_dir}/${panel} | tr -d '\r' | tr '\n' '|')
GENES=${GENES%?}

zgrep '#' ${sample_path}/${sample_id}/variants/${sample_id}_concat_snpsift.vcf.gz > ${sample_id}_${panel}.vcf

zgrep '^chr' ${sample_path}/${sample_id}/variants/${sample_id}_concat_snpsift.vcf | awk -v g="${GENES}" '$0 ~ g' - >> ${sample_id}_${panel}.vcf

cat ${sample_id}_${panel}.vcf | /snpEff/scripts/vcfEffOnePerLine.pl > ${sample_id}_${panel}_OPL.vcf