FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y default-jre default-jdk curl wget unzip

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1

RUN wget https://github.com/RealTimeGenomics/rtg-tools/releases/download/3.12/rtg-tools-3.12-linux-x64.zip
RUN unzip rtg-tools-3.12-linux-x64.zip
RUN rm rtg-tools-3.12-linux-x64.zip
