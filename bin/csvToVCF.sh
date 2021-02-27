#!/usr/bin/env sh

sample_id="$1"
header="$2"

touch ${sample_id}_cnv.vcf

while read meta; do
    echo ${meta} >> ${sample_id}_cnv.vcf
done < ${header}

TODAY=$(date '+%D')

sed -i "s|##fileDate=|##fileDate=$TODAY|" ${sample_id}_cnv.vcf
sed -i "s|SAMPLE1|${sample_id}-sort|g" ${sample_id}_cnv.vcf

cut -d " " -f 3,5,6,7,8,9,10,11,12,13 ${sample_id}_cnv_table.csv | \
awk 'NR>1 {
    printf $1"\t"$3"\tpanelcn.MOPS:"$1":"$3"-"$4"\tN\t";
    if ($10 ~ "CN0" || $10 ~ "CN1") printf "<DEL>\t";
    else if ($10 ~ "CN3" || $10 ~ "CN4") printf "<DUP>\t";
    else printf ".\t";
    if ($9 == "lowQual") printf "lowQual\t";
    else printf "PASS\t";
    printf "SVTYPE=CNV,END="$4",CNCLASS="$10",CYTO="$2"\tRC:MRC:RCN:MRCN\t"$5":"$6":"$7":"$8"\n";
}' >> ${sample_id}_cnv.vcf

sed -i 's/"//g' ${sample_id}_cnv.vcf
sed -i 's/;/,GENE=/g' ${sample_id}_cnv.vcf

# Generate a version of the VCF file that only containes DEL/DUP
grep '^#' ${sample_id}_cnv.vcf > ${sample_id}_filtered_cnv.vcf
grep '^chr' ${sample_id}_cnv.vcf | \
awk '{if ($5 !=".")print}' - >> ${sample_id}_filtered_cnv.vcf
