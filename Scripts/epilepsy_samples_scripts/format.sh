#!/bin/bash

awk '{if ($0 ~ /^#/ || $0 ~ /^I/ || $0 ~ /^$/ || $0 ~ /^-/) printf $0""; else print "\n" $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7;print $8}' gene_check_out.txt | sed 's/.*;ANN\=//' | tr ';' '\n' | tr ',' '\n' | sed '/^[ATCG]/s/^/ANN=/' > annotations.txt

mkdir annotations
cd annotations/

csplit -k ../annotations.txt  '/^[#*]/' '{*}'
rm ../annotations.txt
rm xx00
for FILE in ./*; do
  NAME=$("grep" '#' ${FILE} | "cut" -c3-)
  mv ${FILE} ./${NAME}.txt
done