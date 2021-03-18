FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python python3-pip curl wget unzip

RUN pip3 install \
awscli \
boto3 \
tqdm \
numpy \
scipy \
pandas \
scikit-learn \
matplotlib \
subsample \
pysam

ENV AWS_DEFAULT_REGION=us-east-1
