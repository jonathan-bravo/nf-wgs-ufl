FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y \
curl \
wget \
unzip \
default-jre \
python \
python3-pip \
libbz2-dev \
liblzma-dev

# RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
# RUN unzip awscli-exe-linux-x86_64.zip
# RUN rm awscli-exe-linux-x86_64.zip
# RUN /aws/install

RUN pip3 install awscli

ENV AWS_DEFAULT_REGION=us-east-1

RUN mkdir /data/
COPY *.nf /data/
COPY workflows/ /data/workflows/
COPY nextflow.config /data/
COPY bin/ /data/bin/
COPY modules/ /data/modules/
COPY dist/docker /usr/local/bin/docker

RUN wget https://github.com/nextflow-io/nextflow/releases/download/v21.10.4/nextflow -P /data/
RUN chmod +x /data/nextflow