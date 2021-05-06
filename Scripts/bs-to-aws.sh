#!/usr/bin/env bash

TYPE="WGS";
SPACED="";

usage () {
    cat << EOF

    How to use:
    
    -h -- to call this help function
    -r -- [REQUIRED] is the run id 'NQ-##-##'
    -b -- [REQUIRED] the destination s3 bucket
    -n -- [REQUIRED] the number of the file
    -s -- indicates that there are spaces in the desired file names
    -e -- indicates that the desired files are exome sequences

    EXAMPLE: ./bs-to-aws.sh -r NQ-20-10 -b awsbucket
    EXAMPLE: ./bs-to-aws.sh -r NQ-21-10 -b awsbucket -se

EOF
    exit 1;
}

while getopts "hser:b:" opt; do
    case $opt in
        h)
            usage
            ;;
        r)
            RUN=$OPTARG
            ;;
        b)
            BUCKET="s3://"$OPTARG
            ;;
        n)
            NUM=$OPTARG
            ;;
        s)
            SPACED="SPACED"
            ;;
        e)
            TYPE="WES"
            ;;
        \?)
            usage
            ;;
    esac
done

if [[ "${RUN}" == '' || "${BUCKET}" == '' ]];
then
    printf "    Missing required argument\n";
    usage;
    exit 1;
elif [[ ${SPACED} == "SPACED" && ${NUM} == '' ]];
then
    printf "    Missing required argument\n";
    usage;
    exit 1;
fi

case ${TYPE}" "${SPACED} in
    "WGS SPACED")
        for f in ~/BaseSpace/Projects/WGS/Samples/${RUN}*\ \(${NUM}\)/Files/;
        do
            aws s3 sync "${f}" ${BUCKET}/Fastqs/;
        done
        ;;
    "WGS ")
        for f in ~/BaseSpace/Projects/WGS/Samples/${RUN}*/Files/;
        do
            aws s3 sync "${f}" ${BUCKET}/Fastqs/;
        done
        ;;
    "WES SPACED")
        for f in ~/BaseSpace/Projects/WES/Samples/${RUN}*\ \(${NUM}\)/Files/;
        do
            aws s3 sync "${f}" ${BUCKET}/Exome_Fastqs/;
        done
        ;;
    "WES ")
        for f in ~/BaseSpace/Projects/WES/Samples/${RUN}*/Files/;
        do
            aws s3 sync "${f}" ${BUCKET}/Exome_Fastqs/;
        done
        ;;
esac