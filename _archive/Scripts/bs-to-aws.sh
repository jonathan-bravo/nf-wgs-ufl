#!/usr/bin/env bash

PROJECT="WGS";

usage () {
    cat << EOF

    How to use:
    
    -h -- to call this help function
    -r -- [REQUIRED] is the run id 'NQ-##-##'
    -b -- [REQUIRED] the destination s3 bucket
    -e -- indicates that the desired files are exome sequences

    EXAMPLE: ./bs-to-aws.sh -r NQ-20-10 -b awsbucket
    EXAMPLE: ./bs-to-aws.sh -r NQ-21-10 -b awsbucket -e

EOF
    exit 1;
}

while getopts "her:b:" opt; do
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
        e)
            PROJECT="WES"
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
fi

PROJECTID=$(./bs list projects | \
grep ${PROJECT} | \
awk -F'|' '{print $3}');

mkdir ${RUN};

./bs -v download project \
-i ${PROJECTID} \
-o ${RUN} \
--exclude '*' \
--include ${RUN}'*fastq.gz';

touch ${RUN}_transfer_log.txt

case ${PROJECT} in
    "WGS")
        MINSIZE=4200000000
        for f in ${RUN}/*;
        do
            if [[ ${f} == *.json ]];
            then
                mv ${f} ${RUN}/WGS_${RUN}.json
                aws s3 cp ${RUN}/WGS_${RUN}.json ${BUCKET}/Fastqs/_json_logs/
            else
                BIG_ENOUGH=true
                for FASTQ in ${f}/*;
                do
                    FILESIZE=$(stat -c%s "${FASTQ}")
                    if (( FILESIZE < MINSIZE ));
                    then
                        BIG_ENOUGH=false
                    fi
                done
                if [ "${BIG_ENOUGH}" = true ];
                then
                    aws s3 sync ${f} ${BUCKET}/Fastqs/
                else
                    echo ${f} >> ${RUN}_transfer_log.txt
                fi
            fi
        done
        ;;
    "WES")
        MINSIZE=1200000000
        for f in ${RUN}/*;
        do
            if [[ ${f} == *.json ]];
            then
                mv ${f} ${RUN}/WES_${RUN}.json
                aws s3 cp ${RUN}/WES_${RUN}.json ${BUCKET}/Fastqs/_json_logs/
            else
                BIG_ENOUGH=true
                for FASTQ in ${f}/*;
                do
                    FILESIZE=$(stat -c%s "${FASTQ}")
                    if (( FILESIZE < MINSIZE ));
                    then
                        BIG_ENOUGH=false
                    fi
                done
                if [ "${BIG_ENOUGH}" = true ];
                then
                    aws s3 sync ${f} ${BUCKET}/Exome_Fastqs/
                else
                    echo ${f} >> ${RUN}_transfer_log.txt
                fi
            fi
        done
        ;;
esac

aws s3 cp ${RUN}_transfer_log.txt ${BUCKET}/Fastqs/_json_logs/
