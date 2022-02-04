FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python \
python3-pip \
curl \
wget \
unzip \
libbz2-dev \
liblzma-dev

RUN pip3 install awscli \
pysam \
pyarrow \
pandas

ENV AWS_DEFAULT_REGION=us-east-1
