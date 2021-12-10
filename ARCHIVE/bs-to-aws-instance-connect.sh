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
        aws ec2 describe-instances --instance-ids i-0a21f1098dcb1575a --query  'Reservations[].Instances[].State.Name'
    )
    echo ${status_list} | 
    cut -c4-$(expr ${#status_list} - 7)
);

wait_running () {
    while [[ ${STATUS} != "running" ]];
    do
        echo "Checking if status is running...";
        STATUS=$(
            status_list=$(
                aws ec2 describe-instances --instance-ids i-0a21f1098dcb1575a --query  'Reservations[].Instances[].State.Name'
            )
            echo ${status_list} | 
            cut -c4-$(expr ${#status_list} - 7)
        );
        echo "Status is "${STATUS}"...";
        sleep 10s;
    done
    sleep 10;
    DNS=$(
        dns_list=$(
            aws ec2 describe-instances --instance-ids i-0a21f1098dcb1575a --query 'Reservations[].Instances[].PublicDnsName'
        )
        echo ${dns_list} | 
        cut -c4-$(expr ${#dns_list} - 7)
    );
}

wait_stopped () {
    while [[ ${STATUS} != "stopped" ]];
    do
        echo "Checking if status is stopped...";
        STATUS=$(
            status_list=$(
                aws ec2 describe-instances --instance-ids i-0a21f1098dcb1575a --query  'Reservations[].Instances[].State.Name'
            )
            echo ${status_list} | 
            cut -c4-$(expr ${#status_list} - 7)
        );
        echo "Status is "${STATUS}"...";
        sleep 10s;
    done
}

if [[ ${STATUS} == "stopped" ]];
then 
    aws ec2 start-instances --instance-ids i-0a21f1098dcb1575a;
    echo "Starting...";
    wait_running;
    ssh -i "~/Documents/j.bravo.pem" ubuntu@${DNS};
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a;
elif [[ ${STATUS} == "stopping" ]];
then
    echo "Waiting for instance to be ready to start...";
    wait_stopped;
    aws ec2 start-instances --instance-ids i-0a21f1098dcb1575a;
    echo "Starting...";
    wait_running;
    ssh -i "~/Documents/j.bravo.pem" ubuntu@${DNS};
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a;
else
    echo "Already on";
    ssh -i "~/Documents/j.bravo.pem" ubuntu@${DNS};
    aws ec2 stop-instances --instance-ids i-0a21f1098dcb1575a;
fi