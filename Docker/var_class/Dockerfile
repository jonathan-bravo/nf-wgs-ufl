FROM ubuntu:20.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y && \
apt-get install -y python \
python3-pip \
curl \
wget \
unzip \
libcurl4-openssl-dev \
r-base \
wget \
libxml2-dev \
libxmlsec1-dev \
bedtools \
bcftools \
tabix

RUN pip3 install \
awscli \
boto3 \
pandas \
pysam \
argparse \
openpyxl

RUN R -e "install.packages('RCurl')"
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('panelcn.mops')"
RUN R -e "BiocManager::install('KaryoploteR')"
RUN R -e "BiocManager::install('CopyNumberPlots')"

RUN mkdir Resources/
COPY varClass.py /
COPY ClassifyCNV.py /
COPY parse_clingen_tsv.py /
COPY resources.py /
COPY Resources /Resources/
COPY PMC-ids.csv.gz /
COPY CNV_json_plot.R /
COPY g_ranges.py /
COPY json_to_csv.py /
COPY clinvar.vcf.gz /
COPY clinvar.vcf.gz.tbi /
COPY hg19_regions.bed /

RUN gunzip PMC-ids.csv.gz