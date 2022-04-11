FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y curl \
wget \
unzip \
python3-pip \
bcftools

RUN pip3 install numpy \
scipy

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1

RUN wget https://github.com/HaKmonkey/ExpansionHunterDenovo/releases/download/v0.9.1/ExpansionHunterDenovo.zip
RUN unzip ExpansionHunterDenovo.zip && rm ExpansionHunterDenovo.zip
