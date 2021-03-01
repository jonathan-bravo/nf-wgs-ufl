#!/bin/bash

############
# Get Data #
############

for file in ${PWD}/1000_genomes/sample_lists/*; do
    while read sample_id; do
        sudo aws s3 cp \
        s3://1000genomes/phase3/data/${sample_id}/alignment/ \
        ./ \
        --no-sign-request \
        --recursive
    done < ${file}
done