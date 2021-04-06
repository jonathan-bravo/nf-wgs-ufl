#!/usr/bin/env python3

from os import system
import boto3
import re
import argparse

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
        '-p',
        metavar  = '--PEM',
        type     = str,
        help     = '',
        required = True
    )
    args = parser.parse_args()
    return args


def get_data(bucket, prefix):
    """
    This function grabs the correct data objects from the correct s3 bucket
    :param bucket: string, s3 bucket housing pipeline components
    :param prefix: string, dir that contains the data we are looking for
    :return data: list, list of all data we are looking for [runs, panels]
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


def select_run(runs):
    """
    """
    print("\n" + "#" * 18)
    print("# Available Runs #")
    print("#" * 18 + "\n")
    print("(Index, Run ID)")
    avail_runs = list(enumerate(runs))
    for run in avail_runs:
        print(run)
    run_id = get_choice(avail_runs)
    return run_id


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
    for run in avail_runs:
        print(run)
    run = get_choice(avail_runs)
    return run


def germline_nextflow(pem, bucket, out_dir, run, exome):
    """
    """
    while True:
        single_lane = input("Is this a single lane run? [Y/n]: ").upper()
        if not single_lane in ["N","Y","YES","NO"]:
            print("\nPlease select either (Y)es or (N)o.\n")
            continue
        else:
            break
    if single_lane == "N":
        single_lane = "NO"
    elif single_lane == "Y":
        single_lane = "YES"
    if single_lane == "NO":
        match = ""
    elif single_lane == "YES":
        match_index = int(input("[0]_{R1,R2}_001.fastq.gz or [1]_{1,2}.fq.gz: "))
        match_choices = ["_{R1,R2}_001.fastq.gz", "_{1,2}.fq.gz"]
        match = match_choices[match_index]

    ec2 = boto3.client('ec2')
    ec2.start_instances(InstanceIds=['i-0671758033a9db6fd'])
    ssm_client = boto3.client('ssm') # Need your credentials here
    commands = [
        'git -C /nf-wgs-ufl/ pull',
        f"nextflow run /nf-wgs-ufl/main.nf -work-dir s3://{bucket}/{out_dir}/_work/ --bucket 's3://{bucket}' --run_id '{run}' --single_lane '{single_lane}' --match '{match}' --exome '{exome}' --run_dir 's3://{bucket}/{out_dir}/{run}' -resume"
    ]
    execute_commands_on_linux_instances(ssm_client, commands)


def execute_commands_on_linux_instances(client, commands):
    """Runs commands on remote linux instances
    :param client: a boto/boto3 ssm client
    :param commands: a list of strings, each one a command to execute on the instances
    :param instance_ids: a list of instance_id strings, of the instances on which to execute the command
    :return: the response from the send_command function (check the boto3 docs for ssm client.send_command() )
    """

    resp = client.send_command(
        DocumentName="AWS-RunShellScript", # One of AWS' preconfigured documents
        Parameters={'commands': commands},
        InstanceIds=['i-0671758033a9db6fd']
    )
    return resp






def germline(pem, bucket, out_dir):
    """
    """
    while True:
        exome_data = input("Is this a WES run? [Y/n]: ")
        if not exome_data.upper() in ["N","Y","YES","NO"]:
            print("\nPlease select either (Y)es or (N)o.\n")
            continue
        else:
            break
    if exome_data.upper() == "N": exome_data = "NO"
    elif exome_data.upper() =="Y": exome_data = "YES"
    if exome_data.upper() == "YES":
        samples_dir = "Exome_Fastqs/"
        exome = "YES"
    elif exome_data.upper() == "NO":
        samples_dir = "Fastqs/"
        exome = "NO"
    all_samples = get_data(bucket, samples_dir)
    all_run_ids = get_runs(all_samples)
    run = get_run_id(all_run_ids)
    germline_nextflow(pem, bucket, out_dir.strip('/'), run, exome)


def main():
    """
    """
    args = parse_args()
    pem = args.p
    bucket = args.b # bucket = "hakmonkey-genetics-lab"
    out_dir = 'Pipeline_Output/'
    germline(pem, bucket, out_dir)


if __name__ == '__main__':
    main()
