#!/usr/bin/env sh

sample_id="$1"
header="$2"

touch ${sample_id}.head
touch ${sample_id}.tmp

while read meta; do
    echo ${meta} >> ${sample_id}.head
done < ${header}

TODAY=$(date '+%D')

sed -i "s|##fileDate=|##fileDate=$TODAY|" ${sample_id}.head
sed -i '$ s/\s/\t/g' ${sample_id}.head

cut -d "," -f 2,3,4,5,8,9,10 ${sample_id}_cnvs.csv | \
awk -F',' 'NR>1 {
    printf $1"\t"$2"\tcn.MOPS:"$1":"$2"-"$3"\tN\t";
    if ($7 ~ "CN0" || $7 ~ "CN1") printf "DEL\t";
    else if ($7 ~ "CN3" || $7 ~ "CN4" || $7 ~ "CN5" || $7 ~ "CN6" || $7 ~ "CN7" || $7 ~ "CN8") printf "DUP\t";
    else printf ".\t";
    printf ".\tPASS\tSVTYPE=CNV;END="$3";LENGTH="$4";CNCLASS="$7";GENES=.\tMED:MEAN\t"$5":"$6"\n";
}' >> ${sample_id}.tmp

sed -i 's/",//g' ${sample_id}.tmp
sed -i 's/"//g' ${sample_id}.tmp

cat ${sample_id}.head ${sample_id}.tmp > ${sample_id}_cnv.vcf

rm ${sample_id}.head
rm ${sample_id}.tmp

# Generate a version of the VCF file that only containes DEL/DUP
grep '^#' ${sample_id}_cnv.vcf > ${sample_id}_filtered_cnv.vcf
grep -v '^#' ${sample_id}_cnv.vcf | \
awk '{if ($5 !=".")print}' - >> ${sample_id}_filtered_cnv.vcf
