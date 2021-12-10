#!/usr/bin/env bash

RUN=$1

FILES=( $(aws s3 ls s3://hakmonkey-genetics-lab/Fastqs/ | colrm 1 31) )

for FASTQ in ${FILES[@]};
do
    if [[ ${FASTQ} == ${RUN}* ]]
    then
        aws s3 mv s3://hakmonkey-genetics-lab/Fastqs/${FASTQ} s3://hakmonkey-genetics-lab/Fastqs/_Processed/${RUN}/
    fi
done
