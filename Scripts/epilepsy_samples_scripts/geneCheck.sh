#!/bin/bash

touch gene_check_out.txt

ISEC="./isec_files"
OUT="./gene_check_out.txt"
GENE_FILE="./test_genes.txt"
GENES=$(tr -d '\r' <${GENE_FILE} | tr '\n' '|')
A="NC37"
B="BS37"

for DIR in ${ISEC}/*; do
  NAME=$(echo ${DIR} | awk -F'[./_]' '{print $6}')
  echo "# ${NAME}" >> ${OUT}
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
      java -Xmx32g -jar ~/snpEff/snpEff.jar -noStats GRCh37.75 ${FILE} > ${FILE}_genes
      awk -v g="${GENES}" '$0 ~ g' ${FILE}_genes >> ${OUT}
      rm ${FILE}_genes
      echo >> ${OUT}
    fi
  done
done