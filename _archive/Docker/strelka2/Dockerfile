FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y python default-jre default-jdk curl wget unzip

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
RUN bzip2 -d strelka-2.9.10.centos6_x86_64.tar.bz2
RUN tar -xvf strelka-2.9.10.centos6_x86_64.tar
RUN rm strelka-2.9.10.centos6_x86_64.tar

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1
