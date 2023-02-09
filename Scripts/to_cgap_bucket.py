#!/usr/bin/env python3

import argparse
import boto3
import json
from itertools import chain
from concurrent.futures import ProcessPoolExecutor as ppe

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', required = True)
    parser.add_argument('--threads', type = int, default = 1)
    parser.add_argument('--source', default = 'hakmonkey-genetics-lab')
    parser.add_argument('--dest', default = 'ufl-to-cgap-msa')
    return parser.parse_args()

def get_s3_resource():
    session = boto3.Session()
    return session.resource('s3')

def grab_file_keys(data, source):
    s3 = get_s3_resource()
    source_bucket = s3.Bucket(source)
    objects = source_bucket.objects.filter(Prefix = f'Pipeline_Output/{data[0]}/')
    sam = lambda f: any(sample in f.key for sample in data[1])
    trm = lambda f: f.key.endswith('-p_trimmed.fastq.gz')
    key = lambda f: f.key
    return list(map(key, filter(trm, filter(sam, objects))))

def copy_file(f, source, dest):
    s3 = get_s3_resource()
    dest_bucket = s3.Bucket(dest)
    new_name = f.split('/')[-1].replace('-p_trimmed', '')
    copy_source = {'Bucket': source, 'Key': f}
    print(f'Copying {new_name} over now...')
    dest_bucket.copy(copy_source, new_name)

def main():
    args = parse_args()

    with open(args.json) as f:
        data = json.load(f)

    source_list = [args.source]*len(data)

    with ppe(args.threads) as p:
        files = chain(*list(p.map(grab_file_keys, data.items(), source_list)))

    source_list = [args.source]*sum(len(data[i]) for i in data)
    dest_list = [args.dest]*sum(len(data[i]) for i in data)

    with ppe(args.threads) as p:
        p.map(copy_file, files, source_list, dest_list)

if __name__ == '__main__':
    main()