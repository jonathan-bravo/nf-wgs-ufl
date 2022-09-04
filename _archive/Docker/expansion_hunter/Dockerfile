FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y curl wget unzip

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1

RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/v4.0.2/ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
RUN tar -xvzf ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
RUN rm ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
