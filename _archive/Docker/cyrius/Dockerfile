### Cyrius Dockerfile
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

RUN wget https://github.com/Illumina/Cyrius/archive/refs/tags/v1.1.1.zip; \
unzip v1.1.1.zip; \
rm v1.1.1.zip

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip; \
unzip awscli-exe-linux-x86_64.zip; \
rm awscli-exe-linux-x86_64.zip; \
/aws/install

ENV AWS_DEFAULT_REGION=us-east-1
