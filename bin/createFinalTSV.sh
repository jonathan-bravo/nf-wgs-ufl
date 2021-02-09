#!/usr/bin/env sh

sample_id="$1"
panel="$2"

ann=$(zgrep '^##INFO=<ID=ANN' ${sample_id}_${panel}_OPL.vcf | cut -c75-316 | sed -e 's|\s||g' | sed -e 's/|/\t/g')

awk -F'\t' '{ print $1 }' ${sample_id}_${panel}.info.tsv  | sed 's/|/\t/g' > ${sample_id}_${panel}.ann.tsv
sed -i "s|ANN$|$ann|" ${sample_id}_${panel}.ann.tsv
awk 'BEGIN { FS = OFS = "\t" } { for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "." }; 1' ${sample_id}_${panel}.ann.tsv | awk  'BEGIN { FS = OFS = "\t" } {if(NF==15){$16="."}; print $0}' > ${sample_id}_${panel}.ann.final.tsv

cut -f2- ${sample_id}_${panel}.info.tsv > ${sample_id}_${panel}.info.final.tsv

paste ${sample_id}_${panel}.data ${sample_id}_${panel}.ann.final.tsv ${sample_id}_${panel}.info.final.tsv > ${sample_id}_${panel}.final.tsv

sed -r -i 's/([0-9a-zA-Z.-_]+=)//g' ${sample_id}_${panel}.final.tsv