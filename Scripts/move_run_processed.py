#!/usr/bin/env python3

import argparse
import boto3


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


def move_fastqs(run_id, s3):
    """
    """
    my_bucket = s3.Bucket('hakmonkey-genetics-lab')
    for objects in my_bucket.objects.filter(Prefix = f'Fastqs/{run_id}'):
        copy_source = {
            'Bucket': my_bucket.name,
            'Key': objects.key
        }
        fastq = objects.key[7:]
        print(f'Moving: {fastq}')
        s3.meta.client.copy(
            copy_source,
            my_bucket.name,
            f'Fastqs/_Processed/{run_id}/{fastq}',
            ExtraArgs = {
                'StorageClass': 'DEEP_ARCHIVE',
                'MetadataDirective': 'COPY'
            }
        )
        objects.delete()
        print(f'Finished moving: {fastq}')


def main():
    """
    """
    args = parse_args()
    s3 = get_s3_resource()
    move_fastqs(run_id = args.r, s3 = s3)


if __name__ == '__main__':
    main()