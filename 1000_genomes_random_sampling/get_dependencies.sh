#!/bin/bash

############################
# Downloading Dependencies #
############################

## pip3
if [[ $(which pip3) ]]
then
    pip3 -V
else
    sudo apt-get install -y python3-pip
fi

## subsample
if [[ $(pip3 list | grep -F subsample) ]]
then
    pip3 list subsample
else
    pip3 install subsample
fi

## curl
if [[ $(which curl) ]]
then
    curl -V
else
    sudo apt-get install -y curl
fi

## awscli
if [[ $(which aws) ]]
then
    aws --version
else
    curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" \
    -o "awscliv2.zip"
    unzip awscliv2.zip
    sudo ./aws/install
    rm awscliv2.zip
fi