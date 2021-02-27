FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && apt-get install -y curl libcurl4-openssl-dev r-base wget libxml2-dev libxmlsec1-dev unzip

RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('panelcn.mops')"
RUN R -e "BiocManager::install('KaryoploteR')"
RUN R -e "BiocManager::install('CopyNumberPlots')"

RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip
RUN unzip awscli-exe-linux-x86_64.zip
RUN rm awscli-exe-linux-x86_64.zip
RUN /aws/install

ENV AWS_DEFAULT_REGION=us-east-1
