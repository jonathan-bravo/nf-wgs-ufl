#!/usr/bin/env python3

import argparse
import boto3
import json
import pandas as pd
from itertools import chain
from datetime import date
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor as ppe

@dataclass
class SampleEntry:
    analysis_id:str = ''
    unique_analysis_id:str = ''
    individual_id:str = ''
    sex:str = 'M/F'
    birth_year:int = 1900
    relation_to_proband:str = 'Proband'
    specimen_type:str = 'Whole Blood'
    specimen_id:str = ''
    report_required:str = 'Yes'
    test_requested:str = 'WGS'
    sequencing_date:date = date.today().year
    tags:str = 'MSA/ MSA_control'
    files:str = ''

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

def gen_row(sample, files):
    f_read = files[0].split('/')[-1].replace('-p_trimmed', '')
    r_read = files[1].split('/')[-1].replace('-p_trimmed', '')
    return SampleEntry(
        analysis_id = sample,
        unique_analysis_id = f_read.replace('_R1.fastq.gz', ''),
        individual_id = sample,
        specimen_id = sample,
        files = f'{f_read}, {r_read}'
    )

def gen_xlsx(rows, outfile):
    with pd.ExcelWriter(outfile, mode='w') as writer:
        rows.to_excel(writer, index=False)

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

    samples = list(chain(*[data[run] for run in data]))
    outfile = f'cgap-msa_accessioning_{len(samples)}-samples_{date.today()}.xlsx'

    source_list = [args.source]*len(data)

    with ppe(args.threads) as p:
        files = list(chain(*list(p.map(grab_file_keys, data.items(), source_list))))

    file_pairs = [(files[i], files[i+1]) for i in range(0, len(files), 2)]

    rows = pd.DataFrame(list(map(gen_row, samples, file_pairs)))
    rows.columns = [
        'Analysis ID*',
        'Unique Analysis ID*',
        'Individual ID*',
        'Sex*',
        'Birth Year',
        'Relation to Proband*',
        'Specimen Type*',
        'Specimen ID*',
        'Report Required*',
        'Test Requested*',
        'Sequencing Date',
        'Tags',
        'Files'
    ]

    gen_xlsx(rows, outfile)

    source_list = [args.source]*len(files)
    dest_list = [args.dest]*len(files)

    with ppe(args.threads) as p:
        p.map(copy_file, files, source_list, dest_list)

if __name__ == '__main__':
    main()