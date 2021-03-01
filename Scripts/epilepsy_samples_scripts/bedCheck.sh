#!/bin/bash

printf_new() {
 str=$1
 num=$2
 v=$(printf "%-${num}s" "$str")
 echo "${v// /#}"
}

touch gene_check_out.txt

ISEC="./isec_files"
BED="./bed_files/hg19_intersect.bed"
OUT="./gene_check_out.txt"
GENES="CACNA1E|GABRB2|COL4A2|KIAA2022|SETBP1|KDM5C|ATP1A2|KCNT1|CPT1A|CDKL5|SCN2A|KCNH2|KCNT1|PIK3CA|TUBB3|SLC2A1|CPT1A|KCNB1"
A="NC37"
B="BS37"

for DIR in ${ISEC}/*
do
  NAME=$(echo ${DIR} | awk -F'[./_]' '{print $6}')
  printf_new "#" $((${#NAME} + 4)) >> ${OUT}
  echo "# ${NAME} #" >> ${OUT}
  printf_new "#" $((${#NAME} + 4)) >> ${OUT}
  echo >> ${OUT}
	for FILE in ${DIR}/*; do
		if [[ "${FILE}" == *"0000.vcf" || "${FILE}" == *"0001.vcf" || "${FILE}" == *"0002.vcf" ]]; then
			if [[ "${FILE}" == *"0000.vcf" ]]; then
				echo "In ${A} only" >> ${OUT}
        echo "--------------------------------------------------" >> ${OUT}
			elif [[ "${FILE}" == *"0001.vcf" ]]; then
				echo "In ${B} only" >> ${OUT}
        echo "--------------------------------------------------" >> ${OUT}
      elif [[ "${FILE}" == *"0002.vcf" ]]; then
        echo "In BOTH" >> ${OUT}
        echo "--------------------------------------------------" >> ${OUT}
			fi
			bedtools intersect -a ${BED} -b  ${FILE} > ${FILE}_genes
			awk -v g="${GENES}" '$0 ~ g' ${FILE}_genes >> ${OUT}
      echo >> ${OUT}
		fi
	done
done