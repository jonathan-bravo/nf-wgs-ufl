#!/usr/bin/env python3

import boto3
import argparse
from os import makedirs

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
    makedirs(f'{run_id}/Reports/', exist_ok = True)
    makedirs(f'{run_id}/MultiQC/multiqc_report_data/', exist_ok = True)

    my_bucket = s3.Bucket(bucket)

    data_source = my_bucket.objects.filter(Prefix = f'Pipeline_Output/{run_id}')
    
    for objects in data_source:
        copy_source = {
            'Bucket': my_bucket.name,
            'Key': objects.key
        }
        if objects.key.split('/')[2] == 'MultiQC':
            filename = objects.key.split('/')[-1]
            if '.html' in filename:
                my_bucket.download_file(objects.key, f"{run_id}/MultiQC/{filename}")
            else:
                my_bucket.download_file(objects.key, f"{run_id}/MultiQC/multiqc_report_data/{filename}")
        report_file = (
            '_report.xlsx' in objects.key
            and 'low-qc_report.xlsx' not in objects.key
        )
        if report_file:
            my_bucket.download_file(objects.key, f"{run_id}/Reports/{objects.key.split('/')[-1]}")
        s3.meta.client.copy(
            copy_source,
            my_bucket.name,
            f'Pipeline_Output/_GLACIER_ARCHIVE/{run_id}/{"/".join(objects.key.split("/")[2:])}',
            ExtraArgs = {
                'StorageClass': 'DEEP_ARCHIVE',
                'MetadataDirective': 'COPY'
            }
        )
        objects.delete()

def main():
    args = parse_args()
    s3 = get_s3_resource()
    get_reports('hakmonkey-genetics-lab', s3, args.r)

if __name__ == '__main__':
	main()