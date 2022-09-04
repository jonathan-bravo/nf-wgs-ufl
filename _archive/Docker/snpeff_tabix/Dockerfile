FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y default-jre default-jdk curl wget unzip tabix

RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip
RUN unzip snpEff_latest_core.zip
RUN rm snpEff_latest_core.zip

RUN java -jar /snpEff/snpEff.jar download -v hg19
RUN java -jar /snpEff/snpEff.jar download -v hg38
RUN java -jar /snpEff/snpEff.jar download -v GRCh37.75
RUN java -jar /snpEff/snpEff.jar download -v GRCh38.99

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1
