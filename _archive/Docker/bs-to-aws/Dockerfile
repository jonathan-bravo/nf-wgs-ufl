FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python \
curl \
wget \
unzip

RUN wget "https://launch.basespace.illumina.com/CLI/latest/amd64-linux/bs"
RUN chmod +x ./bs

COPY .basespace/ /root/.basespace/
COPY bs-to-aws.sh /

RUN wget "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip"; \
unzip awscli-exe-linux-x86_64.zip; \
rm awscli-exe-linux-x86_64.zip; \
./aws/install; \
rm -rf ./aws/

ENV AWS_DEFAULT_REGION=us-east-1
