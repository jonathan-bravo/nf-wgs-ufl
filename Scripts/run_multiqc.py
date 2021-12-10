#!/usr/bin/env python3

import argparse
import boto3
import multiqc
from os import system, listdir


def parse_args():
    """
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-r',
        metavar = '--RUN_ID',
        type = str,
        help = '',
        required = True
    )
    args = parser.parse_args()
    return args


def get_s3_resource():
    """
    """
    session = boto3.Session()
    s3 = session.resource('s3')
    return s3


def get_qc_files(run_id, s3):
    """
    """
    system(f'mkdir {run_id}_multiqc_sample_data')
    my_bucket = s3.Bucket('hakmonkey-genetics-lab')
    for objects in my_bucket.objects.filter(Prefix = f'Pipeline_Output/{run_id}/'):
        qc_file = (
            '_md_metrics.txt' in objects.key
            or '_trim_out.log' in objects.key
            or '_snpeff_stats.csv' in objects.key
            or '_fastqc.html' in objects.key
            or '_fastqc.zip' in objects.key
            or '_wgs_metrics.txt' in objects.key
        )
        if qc_file:
            file_name = objects.key.split('/')[-1]
            my_bucket.download_file(objects.key, f'{run_id}_multiqc_sample_data/{file_name}')
            print(f'Downloaded: {file_name}')


def run_multiqc(run_id):
    """
    """
    multiqc.run(
        analysis_dir = f'{run_id}_multiqc_sample_data',
        filename = f'{run_id}_multiqc_report'
    )


def upload_multiqc(run_id, s3):
    """
    """
    s3.meta.client.upload_file(
        f'{run_id}_multiqc_report.html',
        'hakmonkey-genetics-lab',
        f'Pipeline_Output/{run_id}/MultiQC/{run_id}_multiqc_report.html'
    )
    for file in listdir(f'{run_id}_multiqc_report_data/'):
        s3.meta.client.upload_file(
            f'{run_id}_multiqc_report_data/{file}',
            'hakmonkey-genetics-lab',
            f'Pipeline_Output/{run_id}/MultiQC/{run_id}_multiqc_report_data/{file}'
        )
    system(f'rm -rf {run_id}_multiqc_sample_data {run_id}_multiqc_report_data {run_id}_multiqc_report.html')
    print(f'MultiQC data for {run_id} has been uploaded.')


def main():
    """
    """
    args = parse_args()
    s3 = get_s3_resource()
    get_qc_files(run_id = args.r, s3 = s3)
    run_multiqc(run_id = args.r)
    upload_multiqc(run_id = args.r, s3 = s3)


if __name__ == '__main__':
    main()