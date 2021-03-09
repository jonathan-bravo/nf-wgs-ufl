#!/usr/bin/env python3

import os
import boto3
import re
import argparse
from tqdm import tqdm


parser = argparse.ArgumentParser(description='Append genes onto cytoband bed file')
parser.add_argument(
    '-w',
    metavar  = '--WORKFLOW',
    type     = str,
    help     = 'which workflow to launch',
    required = True
)
parser.add_argument(
    '-b',
    metavar  = '--BUCKET',
    type     = str,
    help     = 'which S3 bucket to use when launching the workflow',
    required = True
)
parser.add_argument(
    '-a',
    metavar  = '--AWS_ACCESS_KEY_ID',
    type     = str,
    help     = 'aws access key',
    required = True
)
parser.add_argument(
    '-s',
    metavar  = '--AWS_SECRET_ACCESS_KEY',
    type     = str,
    help     = 'aws secret access key',
    required = True
)


ARGS = parser.parse_args()
BUCKET = ARGS.b # bucket = "hakmonkey-genetics-lab"
WORKFLOW = ARGS.w
ACCESS_KEY = ARGS.a
SECREY_KEY = ARGS.s
WORKFLOW_LIST = ["GERMLINE", "MULTIQC", "APPLY_PANELS"]
OUT_DIR = 'Pipeline_Output/'
PROCESSED_DIR = "_Processed/"
PANELS_DIR = 'Pipeline/Reference/panels/'
TRASH = "sudo nextflow clean -f"

# From data_ops
def usage():
    """
    """


def check_workflow(workflow):
    if workflow.upper() not in WORKFLOW_LIST:
        print("please select one of these options: ", WORKFLOW_LIST)
        exit()


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
            data.append(obj.get('Key')
                           .replace(str(prefix), ""))
    else:
        for obj in result.get('CommonPrefixes'):
            x = re.search('^Pipeline_Output/_', obj.get('Prefix'))
            if (x == None):
                data.append(obj.get('Prefix')
                               .replace(str(prefix), ""))

    return data


def enumerate_data(data):
    """
    """

    selections = (list(enumerate(data)))

    return selections


def list_data(data):
    """
    """
    
    for datum in data:
        print(datum)


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

    avail_runs = enumerate_data(data = runs)
    list_data(data = avail_runs)
    run_id = get_choice(choices = avail_runs)

    return run_id


def get_runs(data):
    """
    """

    data = [sample[:8] for sample in data]
    data = list(set(data))
    if '' in data:
        data.remove('')
    data.sort()

    return data


def get_run_id(run_ids):
    """
    """

    avail_runs = enumerate_data(data = run_ids)

    print("\n" + "#" * 18)
    print("# Available Runs #")
    print("#" * 18 + "\n")
    print("(Index, Run ID)")

    list_data(data = avail_runs)

    run = get_choice(choices = avail_runs)

    return run


def archive_fastqs(bucket, processed_dir, samples_dir, run):
    """
    """

    s3 = boto3.resource('s3')
    dest_bucket = s3.Bucket(bucket)

    all_samples = get_data(bucket = bucket, prefix = samples_dir)
    for sample in tqdm(all_samples):
        if run in sample:
            copy_source = {
                'Bucket': bucket,
                'Key': samples_dir + sample
            }
            obj = dest_bucket.Object(samples_dir + processed_dir + sample)
            obj.copy(copy_source)
            src = dest_bucket.Object(samples_dir + sample)
            src.delete()
            print(sample + " archived")


def germline_nextflow(bucket, out_dir, run, exome):
    """
    """

    while True:
        single_lane = input("Is this a single lane run? [Y/n]: ")
        if not single_lane.upper() in ["N","Y","YES","NO"]:
            print("\nPlease select either (Y)es or (N)o.\n")
            continue
        else:
            break

    if single_lane.upper() == "N":
        single_lane = "NO"
    elif single_lane.upper() == "Y":
        single_lane = "YES"

    if single_lane.upper() == "NO":
        match = ""
    elif single_lane.upper() == "YES":
        match_index = int(input("[0]_{R1,R2}_001.fastq.gz or [1]_{1,2}.fq.gz: "))
        match_choices = ["_{R1,R2}_001.fastq.gz", "_{1,2}.fq.gz"]
        match = match_choices[match_index]


    launch = "nextflow run main.nf -work-dir s3://{bucket}/{out_dir}/_work/ --pipeline 'GERMLINE' --bucket 's3://{bucket}' --run_id '{run}' --single_lane '{laneage}' --match '{match_lane}' --exome '{exome}' -resume".format(
        bucket = bucket,
        out_dir = out_dir,
        run = run,
        laneage = single_lane.upper(),
        match_lane = match,
        exome = exome
    )

    #os.system(launch)
    submit_job(launch)

    #os.system(TRASH)


def germline():
    """
    """

    while True:
        exome_data = input("Is this a WES run? [Y/n]: ")
        if not exome_data.upper() in ["N","Y","YES","NO"]:
            print("\nPlease select either (Y)es or (N)o.\n")
            continue
        else:
            break

    if exome_data.upper() == "N":
        exome_data = "NO"
    elif exome_data.upper() =="Y":
        exome_data = "YES"

    if exome_data.upper() == "YES":
        samples_dir = "Exome_Fastqs/"
        exome = "YES"
    elif exome_data.upper() == "NO":
        samples_dir = "Fastqs/"
        exome = "NO"

    all_samples = get_data(bucket = BUCKET, prefix = samples_dir)

    all_run_ids = get_runs(data = all_samples)

    run = get_run_id(run_ids = all_run_ids)

    germline_nextflow(
        bucket = BUCKET,
        out_dir = OUT_DIR.strip('/'),
        run = run,
        exome = exome
    )

    # archive_fastqs(
    #     bucket = BUCKET,
    #     processed_dir = PROCESSED_DIR,
    #     samples_dir = samples_dir,
    #     run = run
    # )


def multiqc_nextflow(bucket, run_id, output_dir):
    """
    """

    launch = "nextflow run main.nf -work-dir s3://{bucket}/{output_dir}/_work/ --pipeline 'MULTIQC' --run_id '{run_id}' --run_dir 's3://{bucket}/{output_dir}/{run_id}'".format(
        run_id = run_id,
        bucket = bucket,
        output_dir = output_dir)

    #os.system(launch)
    submit_job(launch)

    #os.system(TRASH)


def multiqc():
    """
    """
    
    run = get_data(bucket = BUCKET, prefix = OUT_DIR)
    run_id = select_run(runs = run)

    multiqc_nextflow(
        bucket = BUCKET,
        run_id = run_id.strip('/'),
        output_dir = OUT_DIR.strip('/'))


def select_panels(sample_list, panel_list):
    """
    """

    test_pairs = []

    avail_samples = enumerate_data(data = sample_list)
    avail_panels = enumerate_data(data = panel_list)

    for sample in avail_samples:
        if not "MultiQC/" in sample:
            print("\n Sample:")
            print(sample[1] + "\n")

            test_num = int(input("How many panels would you like to run?: "))

            if test_num != 0:
                print("\n" + "#" * 20)
                print("# Available Panels #")
                print("#" * 20 + "\n")
                print("(Index, Panel)")
                list_data(data = avail_panels)

            for i in range(test_num):
                print("panel number " + str(i+1))

                test = get_choice(choices = avail_panels)
                test_pairs.append([sample[1].strip('/'), test])

    return test_pairs


def apply_panels_nextflow(test_pairs, bucket, runs_dir, run_id, panels_dir):
    """
    """

    pairs = "sed -i 's/pairs_ch = Channel.from()/pairs_ch = Channel.from({test_lists})/' ufl-apply_panels.nf".format(test_lists = test_pairs)

    os.system(pairs)

    quotes = "sed -ri \"s/\\[([a-zA-Z0-9_.-]+),\\s([a-zA-Z0-9_.-]+)\\]/\\['\\1','\\2'\\]/g\" ufl-apply_panels.nf"

    os.system(quotes)

    launch = "sudo nextflow run main.nf -work-dir s3://{bucket}/{runs_dir}/_work/ --pipeline 'APPLY_PANELS' --run_dir 's3://{bucket}/{runs_dir}/{run_id}' --panels_dir 's3://{bucket}/{panels_dir}' --bucket 's3://{bucket}' -resume".format(
        bucket = bucket,
        runs_dir = runs_dir,
        run_id = run_id,
        panels_dir = panels_dir)

    os.system(launch)

    reset = "sed -i 's/^pairs_ch = Channel.from(.*$/pairs_ch = Channel.from()/' reporting-ufl.nf"

    os.system(reset)

    os.system(TRASH)


def apply_panels():
    """
    """

    ################################
    # Listing & Getting Run Choice #
    ################################

    runs = get_data(bucket = BUCKET, prefix = OUT_DIR)

    run_id = select_run(runs = runs)

    ################################
    # Creating Sample/ Panel Pairs #
    ################################

    samples = get_data(bucket = BUCKET, prefix = OUT_DIR + run_id)
    
    panels = get_data(bucket = BUCKET, prefix = PANELS_DIR)
    while('' in panels):
        panels.remove('')

    test_pairs = select_panels(sample_list = samples, panel_list = panels)

    #################################
    # Running The Nextflow Pipeline #
    #################################

    apply_panels_nextflow(
        test_pairs = test_pairs,
        bucket = BUCKET,
        runs_dir = OUT_DIR.strip('/'),
        run_id = run_id.strip('/'),
        panels_dir = PANELS_DIR.strip('/'))


def submit_job(command, workflow=WORKFLOW, access=ACCESS_KEY, secret=SECREY_KEY):
    """
    """

    work = "bash -c 'cd /data/;{command}'".format(command = command)

    client = boto3.client('batch')

    response = client.submit_job(
        jobName=workflow,
        jobQueue='hakmonkey-nextflow',
        jobDefinition='nextflow-ufl-germline:3',
        containerOverrides={
            'command': [
                work
            ],
            'environment': [
            {
                'name': 'AWS_ACCESS_KEY_ID',
                'value': access
            },
            {
                'name': 'AWS_SECRET_ACCESS_KEY',
                'value': secret
            }
        ],
        }
    )

    print(response)
    print(command)



def main():
    """
    """
    # check workflow to make sure that it is a valid input, else exit
    check_workflow(WORKFLOW)
    # depending on the workflow selected, prompt user for correct inpout as before
    if WORKFLOW.upper() == "GERMLINE":
        germline()
    elif WORKFLOW.upper() == "MULTIQC":
        multiqc()
    elif WORKFLOW.upper() == "APPLY_PANELS":
        apply_panels()


if __name__ == '__main__':
    main()
