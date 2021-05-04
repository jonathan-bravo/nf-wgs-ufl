#!/usr/bin/env bash

DNS=$(
    dns_list=$(
        aws ec2 describe-instances --instance-ids i-0a21f1098dcb1575a --query 'Reservations[].Instances[].PublicDnsName'
    )
    echo ${dns_list} | 
    cut -c4-$(expr ${#dns_list} - 7)
);
STATUS=$(
    status_list=$(
        aws ec2 describe-instance-status --instance-ids i-0a21f1098dcb1575a | 
        jq '.InstanceStatuses[0].InstanceState.Name'
    )
    echo ${status_list} | 
    cut -c2-$(expr ${#status_list} - 1)
);

usage () {
    cat << EOM

    How to use:
    
    -h -- to call this help function
    -r -- [REQUIRED] is the run id 'NQ-##-##'
    -b -- [REQUIRED] the destination s3 bucket
    -s -- indicates that there are spaces in the desired file names
    -e -- indicates that the desired files are exome sequences

    EXAMPLE: ./bs-to-aws.sh -r NQ-20-10 -b awsbucket
    EXAMPLE: ./bs-to-aws.sh -r NQ-21-10 -b awsbucket -se

EOM
    exit 1;
}

to_aws () {
    TYPE="WGS";
    SPACED="";

    case ${TYPE}" "${SPACED} in
        "WGS SPACED")
            for f in ~/BaseSpace/Projects/WGS/Samples/${RUN}*\ \(2\)/Files/;
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
            for f in ~/BaseSpace/Projects/WES/Samples/${RUN}*\ \(2\)/Files/;
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
}

mount_bs () {
    basemount BaseSpace/
    echo "Mounting BaseSpace folder..."
    sleep 2 &
    PID=$!
    i=1
    sp="/-\|"
    echo -n ' '
    while [ -d /proc/$PID ]
    do
        printf "\b${sp:i++%${#sp}:1}"
    done
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
fi

if [[ ${STATUS} != "running" ]];
then 
    # aws ec2 start-instances --instance-ids i-0a21f1098dcb1575a;
    echo "Starting...";
    sleep 5 &
    PID=$!
    i=1
    sp="/-\|"
    echo -n ' '
    while [ -d /proc/$PID ]
    do
        printf "\b${sp:i++%${#sp}:1}"
    done
    ssh -i "~/Documents/j.bravo.pem" ubuntu@${DNS} "$(typeset -f); mount_bs; to_aws; exit"
    echo "stopping"
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a
else
    echo "already on"
    ssh -i "~/Documents/j.bravo.pem" ubuntu@${DNS} "$(typeset -f); to_aws; exit"
    echo "stopping"
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a
fi