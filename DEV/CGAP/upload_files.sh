#!/usr/bin/env bash

BUCKET="ufl-to-cgap-msa"
XLSX_FILENAME=$1

# Install s3fs-fuse for mounting S3 buckets
aws s3 cp s3://hakmonkey-genetics-lab/home/j.bravo/cgap-msa/.cgap-keys.json .
aws s3 cp s3://hakmonkey-genetics-lab/home/j.bravo/cgap-msa/${XLSX_FILENAME} .
sudo yum update -y
sudo amazon-linux-extras install epel -y
sudo yum install s3fs-fuse -y
pip install --upgrade pip

# Mount buckets to ~/upload_files directory
mkdir upload_files
s3fs ${BUCKET} ~/upload_files/ -o iam_role

# Create virtual env for package installation
python3 -m venv ~/cgap_submission
source ~/cgap_submission/bin/activate

# Run SubmitCGAP with mounted files
pip install submit_cgap
submit-metadata-bundle ${XLSX_FILENAME} -s https://cgap-msa.hms.harvard.edu -u ~/upload_files/ -nq -sf