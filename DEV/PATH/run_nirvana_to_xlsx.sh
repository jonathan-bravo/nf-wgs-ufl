#!/usr/bin/env bash

RUN=$1

SAVEIFS=$IFS
IFS=$(echo -en "\n\b")

for FOLDER in /Volumes/PATH/DRL/Molecular/NGS/NGS\ Development/Dragen/GOAL-235856625/${RUN}*;
do
    FOLDERNAME=$(basename -- ${FOLDER});
    FILENAME=$(basename -- ${FOLDER}/${RUN}*);
    SAMPLE_ID=${FILENAME%-ds.*};
    if [ -f /Volumes/PATH/DRL/Molecular/NGS/NGS\ Development/Dragen/GOAL-235856625/${FOLDERNAME}/${FILENAME}/${SAMPLE_ID}.hard-filtered.annotations.json.gz ];
    then
        echo ${FILENAME}
        /Volumes/PATH/DRL/Molecular/NGS/NGS\ Development/Dragen/GOAL-235856625/nirvana_to_xlsx.py -s ${SAMPLE_ID} -p /Volumes/PATH/DRL/Molecular/NGS/NGS\ Development/Dragen/GOAL-235856625/${FOLDERNAME}/${FILENAME}
        mv ${SAMPLE_ID}.xlsx /Volumes/PATH/DRL/Molecular/NGS/NGS\ Development/Dragen/GOAL-235856625/${FOLDERNAME}/${SAMPLE_ID}.xlsx
    fi
done

IFS=$SAVEIFS
