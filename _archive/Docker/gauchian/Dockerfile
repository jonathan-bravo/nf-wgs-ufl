### Gauchian Dockerfile
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python \
curl \
wget \
unzip \
python3-pip

RUN pip install --upgrade pip
RUN pip install numpy \
scipy \
pysam \
statsmodels

RUN wget https://github.com/Illumina/Gauchian/archive/refs/tags/v1.0.2.zip; \
unzip v1.0.2.zip; \
rm v1.0.2.zip

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip; \
unzip awscli-exe-linux-x86_64.zip; \
rm awscli-exe-linux-x86_64.zip; \
/aws/install

WORKDIR /Gauchian-1.0.2

RUN python3 setup.py install

WORKDIR /

ENV AWS_DEFAULT_REGION=us-east-1