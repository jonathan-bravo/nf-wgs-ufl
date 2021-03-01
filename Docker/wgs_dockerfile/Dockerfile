FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y

RUN apt-get install -y python python3-pip default-jre default-jdk fastqc trimmomatic bwa samtools curl libcurl4-openssl-dev r-base tabix bcftools

RUN apt-get install -y wget
RUN wget https://github.com/broadinstitute/picard/releases/download/2.23.8/picard.jar

RUN wget https://github.com/Illumina/strelka/releases/download/v2.9.10/strelka-2.9.10.centos6_x86_64.tar.bz2
RUN bzip2 -d strelka-2.9.10.centos6_x86_64.tar.bz2
RUN tar -xvf strelka-2.9.10.centos6_x86_64.tar
RUN rm strelka-2.9.10.centos6_x86_64.tar

RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('panelcn.mops')"

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

## I don't think these variables are necessary for running the pipeline
#ENV AWS_ACCESS_KEY_ID=XXX
#ENV AWS_SECRET_ACCESS_KEY=XXX
ENV AWS_DEFAULT_REGION=us-east-1

RUN pip3 install multiqc

RUN wget https://github.com/Illumina/ExpansionHunter/releases/download/v4.0.2/ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
RUN tar -xvzf ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
RUN rm ExpansionHunter-v4.0.2-linux_x86_64.tar.gz
