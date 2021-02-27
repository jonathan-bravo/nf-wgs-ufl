FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y samtools curl wget unzip

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1
