#!/bin/bash

# Making temp folder to hold isec output
mkdir isec_files

# Making file to hold final analysis
touch NC_v_BS.csv

# Making header of file
printf "SAMPLE_A_VS_SAMPLE_B, " >> NC_v_BS.csv
printf "%%_CONCORDANCE, " >> NC_v_BS.csv
printf "IN_A_ONLY, " >> NC_v_BS.csv
printf "IN_B_ONLY, " >> NC_v_BS.csv
echo "IN_BOTH" >> NC_v_BS.csv

# Paths to files containing files of interest
HOME="/mnt/d/Epilepsy_Samples/Analysis"
NC37="./neurocode_vcfs_20_subset"
BS37="./basespace_neurocode_hgc37_hard_filtered_vcfs"
##BS38="./basespace_neurocode_hg38_hard_filtered_vcfs"
ISEC="./isec_files"

# IDT hg19 and Agilent hg19
BEDS37="./bed_files/hg19_intersect.bed"

# IDT hg38 and Agilent hg19
##BEDS38="./bed_files/hg38_intersect.bed"

for FILE1 in ${NC37}/*
do
    # File variables
    BUFFER=`echo ${FILE1} | cut -c 33-1000`
    FILE2="BS37_${BUFFER}"
    ##FILE3="BS38_${BUFFER}"

    # Indexing files
    bcftools index --threads 12 ${FILE1}
    bcftools index --threads 12 ${BS37}/${FILE2}
    ##bcftools index --threads 10 ${BS38}/${FILE3}

    # Filtering files with intersected bed files
    bcftools view -R ${BEDS37} -i'QUAL>=20 && INFO/DP>=10' --threads 12 ${FILE1} -o ${FILE1}_hg19_filtered.vcf
    ##bcftools view -R ${BEDS38} -i'QUAL>=20 && DP>=10' --threads 10 ${FILE1} -o ${FILE1}_hg38_filtered.vcf
    bcftools view -R ${BEDS37} -i'QUAL>=20 && INFO/DP>=10' --threads 12 ${BS37}/${FILE2} -o ${FILE2}_hg19_filtered.vcf
    ##bcftools view -R ${BEDS38} -i'QUAL>=20 && DP>=10' --threads 10 ${BS38}/${FILE3} -o ${FILE3}_hg38_filtered.vcf

    # Zipping files again
    bgzip -@ 12 ${FILE1}_hg19_filtered.vcf
    ##bgzip -@ 10 ${FILE1}_hg38_filtered.vcf
    bgzip -@ 12 ${FILE2}_hg19_filtered.vcf
    ##bgzip -@ 10 ${FILE3}_hg38_filtered.vcf

    # Indexing the new filtered files
    bcftools index --threads 12 ${FILE1}_hg19_filtered.vcf.gz
    ##bcftools index --threads 10 ${FILE1}_hg38_filtered.vcf.gz
    bcftools index --threads 12 ${FILE2}_hg19_filtered.vcf.gz
    ##bcftools index --threads 10 ${FILE3}_hg38_filtered.vcf.gz

    # Comparing the files
    bcftools isec -p ${FILE1}_v_${FILE2}_hg19_isec --threads 12 ${FILE1}_hg19_filtered.vcf.gz ${FILE2}_hg19_filtered.vcf.gz
    ##bcftools isec -p ${FILE1}_v_${FILE3}_hg38_isec --threads 10 ${FILE1}_hg38_filtered.vcf.gz ${FILE3}_hg38_filtered.vcf.gz
done

# Cleaning up index files
rm ${NC37}/*.csi
rm ${BS37}/*.csi
##rm ${BS38}/*.csi
rm ./*.csi

# Cleaning up intermediate filtered .vcf files
rm ${NC37}/*_filtered.vcf.gz
rm *_filtered.vcf.gz

# Moving `bcftools isec` output into folder
mv ${NC37}/*_isec ./isec_files/

for FILE in ${ISEC}/*
do
    a=`grep -v "^#" ${FILE}/0000.vcf | wc -l`
    b=`grep -v "^#" ${FILE}/0001.vcf | wc -l`
    c=`grep -v "^#" ${FILE}/0002.vcf | wc -l`
    d=`expr $a + $b + $c`
    e=`printf %.3f "$((10**3 * ${c}/${d} * 100))e-3"`
    printf "${FILE}, " >> NC_v_BS.csv
    printf "${e}%%, " >> NC_v_BS.csv
    printf "${a}, " >> NC_v_BS.csv
    printf "${b}, " >> NC_v_BS.csv
    echo "${d}" >> NC_v_BS.csv
done

# rm -r ${ISEC}

printf "Average of %% Concordance, " >> NC_v_BS.csv
echo "=AVERAGE(B2:B21)" >> NC_v_BS.csv
printf "Meadian of %% Concordance, " >> NC_v_BS.csv
echo "=MEDIAN(B2:B21)" >> NC_v_BS.csv