#!/usr/bin/env bash

RUN=$1

FILES=( $(aws s3 ls s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/ | colrm 1 31) )

mkdir ${RUN}_multiqc_sample_data

for SAMPLE in ${FILES[@]};
do
    SAMPLEID=${SAMPLE%/};
    aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/alignment/${SAMPLEID}_md_metrics.txt ${RUN}_multiqc_sample_data/;
    aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/Trimmomatic/${SAMPLEID}_trim_out.log ${RUN}_multiqc_sample_data/;
    aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/variants/${SAMPLEID}_snpeff_stats.csv ${RUN}_multiqc_sample_data/;
    aws s3 sync s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/fastqc_${SAMPLEID}_logs/ ${RUN}_multiqc_sample_data/;
    aws s3 sync s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/fastqc_${SAMPLEID}_trimmed_logs/ ${RUN}_multiqc_sample_data/;
    aws s3 sync s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/${SAMPLEID}/wgs_metrics/ ${RUN}_multiqc_sample_data/;
done

multiqc -n ${RUN}_multiqc_report ./${RUN}_multiqc_sample_data/

aws s3 cp ${RUN}_multiqc_report.html s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/MultiQC/
aws s3 cp --recursive ${RUN}_multiqc_report_data/ s3://hakmonkey-genetics-lab/Pipeline_Output/${RUN}/MultiQC/${RUN}_multiqc_report_data/

rm -rf ${RUN}_multiqc_sample_data
rm -rf ${RUN}_multiqc_report_data
rm ${RUN}_multiqc_report.html
