#!/usr/bin/env python3

import boto3
import argparse
from os import mkdir

def parse_args():
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-r',
        metavar = '--RUN_ID',
        type = str,
        help = 'The run reports you want',
        required = True
    )
    args = parser.parse_args()
    return args

def get_s3_resource():
    session = boto3.Session()
    s3 = session.resource('s3')
    return s3

def get_reports(bucket, s3, run_id):
    mkdir(f'{run_id}_reports/')
    my_bucket = s3.Bucket(bucket)
    data_source = my_bucket.objects.filter(Prefix = f'Pipeline_Output/{run_id}')
    for objects in data_source:
        report_file = (
            '_report.xlsx' in objects.key
            and 'low-qc_report.xlsx' not in objects.key
        )
        if report_file:
            my_bucket.download_file(objects.key, f"{run_id}_reports/{objects.key.split('/')[-1]}")

def main():
    args = parse_args()
    s3 = get_s3_resource()
    get_samples('hakmonkey-genetics-lab', s3, args.r)

if __name__ == '__main__':
	main()