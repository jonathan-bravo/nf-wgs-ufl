FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get install -y python
RUN apt-get install -y default-jre
RUN apt-get install -y default-jdk
RUN apt-get install -y curl wget unzip

RUN wget https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2
RUN bzip2 -d manta-1.6.0.centos6_x86_64.tar.bz2
RUN tar -xvf manta-1.6.0.centos6_x86_64.tar
RUN rm manta-1.6.0.centos6_x86_64.tar

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1
