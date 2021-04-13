#!/usr/bin/env python3

import boto3
import re
import argparse
from time import sleep

def parse_args():
    """
    """
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-b',
        metavar  = '--BUCKET',
        type     = str,
        help     = 'which S3 bucket to use when launching the workflow',
        required = True
    )
    parser.add_argument(
        '--exome',
        action = argparse.BooleanOptionalAction,
        default = False
    )
    parser.add_argument(
        '--single_lane',
        action = argparse.BooleanOptionalAction,
        default = False
    )
    parser.add_argument(
        '--multiqc',
        action = argparse.BooleanOptionalAction,
        default = False
    )
    args = parser.parse_args()
    return args


def get_data(bucket, prefix):
    """Grabs the correct data objects from the correct s3 bucket.
    
    Keyword arguments:

    bucket -- s3 bucket housing pipeline components and samples
    prefix -- dir that contains the data we are looking for

    Return
    
    data -- list of all data we are looking for [runs, panels]
    """
    client = boto3.client('s3')
    result = client.list_objects_v2(
        Bucket = bucket,
        Prefix = prefix,
        Delimiter = '/'
    )
    data = []
    if 'panels' in prefix or 'Fastqs' in prefix:
        for obj in result.get('Contents'):
            data.append(obj.get('Key').replace(str(prefix), ""))
    else:
        for obj in result.get('CommonPrefixes'):
            x = re.search('^Pipeline_Output/_', obj.get('Prefix'))
            if (x == None):
                data.append(obj.get('Prefix').replace(str(prefix), ""))
    return data


def get_choice(choices):
    """
    """
    selection = int(input("\n Select the desired index: "))
    try:
        print("\n Selected Object: " + choices[selection][1] + "\n")
        return choices[selection][1]
    except IndexError:
        print("\n Please select an option from the list.")
        get_choice(choices)


def get_runs(data):
    """
    """
    data = [sample[:8] for sample in data]
    data = list(set(data))
    if '' in data: data.remove('')
    data.sort()
    return data


def get_run_id(run_ids):
    """
    """
    avail_runs = list(enumerate(run_ids))
    print("\n" + "#" * 18)
    print("# Available Runs #")
    print("#" * 18 + "\n")
    print("(Index, Run ID)")
    for run in avail_runs: print(run)
    run = get_choice(avail_runs)
    return run


def yes_no(prompt):
    """
    """
    while True:
        answer = input(prompt).upper()
        if not answer in ["N","Y","YES","NO"]:
            print("\nPlease select either (Y)es or (N)o.\n")
            continue
        else:
            if answer == "N": answer = "NO"
            elif answer == "Y": answer = "YES"
            return answer


def get_match():
    """
    """
    match_choices = ["_{R1,R2}_001.fastq.gz", "_{1,2}.fq.gz"]
    while True:
        match_index = int(input("[0]_{R1,R2}_001.fastq.gz or [1]_{1,2}.fq.gz: "))
        if match_index > 1:
            print("\nPlease select either [0] or [1].\n")
            continue
        else:
            return match_choices[match_index]


def run_commands(bucket, out_dir, exome, single_lane, multiqc):
    """
    """
    if exome:
        samples_dir = "Exome_Fastqs/"
        exome = "YES"
    else:
        samples_dir = "Fastqs/"
        exome = "NO"
    if single_lane:
        match = get_match()
        single_lane = "YES"
    else:
        match = ""
        single_lane = "NO"
    if not multiqc: pipeline = "GERMLINE"
    else: pipeline = "MULTIQC"
    all_samples = get_data(bucket, samples_dir)
    all_run_ids = get_runs(all_samples)
    run = get_run_id(all_run_ids)
    client = boto3.client('batch')
    response = client.submit_job(
        jobName='ufl-germline',
        jobQueue='hakmonkey-nextflow',
        jobDefinition='nextflow-ufl-germline:10',
        containerOverrides={
            'command': [
                "bash",
                "-c",
                f"nextflow run /data/main.nf -work-dir s3://{bucket}/{out_dir}_work/ --bucket s3://{bucket} --run_id {run} --single_lane {single_lane} --match {match} --exome {exome} --pipeline {pipeline}"
            ],
        }
    )
    return response


def main():
    args = parse_args()
    bucket = args.b
    out_dir = 'Pipeline_Output/'
    response = run_commands(bucket, out_dir, args.exome, args.single_lane, args.multiqc)
    #print(response)


if __name__ == '__main__':
    main()
