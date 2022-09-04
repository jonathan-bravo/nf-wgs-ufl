### how to setup nirvana
FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python \
curl \
wget \
unzip \
git

RUN wget https://packages.microsoft.com/config/ubuntu/20.04/packages-microsoft-prod.deb -O packages-microsoft-prod.deb && \
dpkg -i packages-microsoft-prod.deb && \
rm packages-microsoft-prod.deb

RUN apt-get update -y && \
apt-get install -y apt-transport-https && \
apt-get update -y && \
apt-get install -y dotnet-sdk-3.1

RUN dotnet tool install -g Amazon.Lambda.Tools

RUN mkdir -p Nirvana/Data
WORKDIR /Nirvana
RUN wget https://github.com/Illumina/Nirvana/releases/download/v3.17.0/Nirvana-3.17.0-dotnet-3.1.0.zip && \
unzip Nirvana-3.17.0-dotnet-3.1.0.zip && \
rm Nirvana-3.17.0-dotnet-3.1.0.zip

## awscli
WORKDIR /
RUN wget https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip && \
unzip awscli-exe-linux-x86_64.zip && \
rm awscli-exe-linux-x86_64.zip && \
/aws/install

ENV AWS_DEFAULT_REGION=us-east-1