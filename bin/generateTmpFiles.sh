#!/usr/bin/env sh

sample_id="$1"
panel="$2"

zgrep '^#CHROM' ${sample_id}_${panel}_OPL.vcf > ${sample_id}_${panel}.data
sed -i 's/\tINFO//' ${sample_id}_${panel}.data
sed -i 's/#//' ${sample_id}_${panel}.data

zgrep '^chr' ${sample_id}_${panel}_OPL.vcf | awk -F'\t' '{
    print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$9"\t"$10;
}' >> ${sample_id}_${panel}.data

sed -rn 's/^##INFO=<ID=([0-9a-zA-Z.-_\;]+),.*/\1/gp' ${sample_id}_${panel}_OPL.vcf | tr '\n' '\t' | sed -e 's/\t$/\n/' > ${sample_id}_${panel}.info.txt

zgrep '^chr' ${sample_id}_${panel}_OPL.vcf | awk -F'\t' '{print $8}' | sed -e 's/\;/\t/g' >> ${sample_id}_${panel}.info.txt