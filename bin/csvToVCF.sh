#!/usr/bin/env sh

#sample_id="$1"
header="$1"
tmp_head="$2"
tmp_data="$3"
csv_file="$4"
vcf_file="$5"

touch ${tmp_head}
touch ${tmp_data}

while read meta; do
    echo ${meta} >> ${tmp_head}
done < ${header}

TODAY=$(date '+%D')

sed -i "s|##fileDate=|##fileDate=$TODAY|" ${tmp_head}
sed -i '$ s/\s/\t/g' ${tmp_head}

cut -d "," -f 2,3,4,5,8,9,10 ${csv_file} | \
awk -F',' 'NR>1 {
    printf $1"\t"$2"\tcn.MOPS:"$1":"$2"-"$3"\tN\t";
    if ($7 ~ "CN0" || $7 ~ "CN1") printf "DEL\t";
    else if ($7 ~ "CN3" || $7 ~ "CN4" || $7 ~ "CN5" || $7 ~ "CN6" || $7 ~ "CN7" || $7 ~ "CN8") printf "DUP\t";
    else printf ".\t";
    printf ".\tPASS\tSVTYPE=CNV;END="$3";LENGTH="$4";CNCLASS="$7";GENES=.\tMED:MEAN\t"$5":"$6"\n";
}' >> ${tmp_data}

sed -i 's/",//g' ${tmp_data}
sed -i 's/"//g' ${tmp_data}

cat ${tmp_head} ${tmp_data} > ${vcf_file}