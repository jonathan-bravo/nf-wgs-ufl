#!/usr/bin/env python3

import os
import csv
import boto3

def usage():
    """
    """


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

    for obj in result.get('CommonPrefixes'):
        data.append(obj.get('Prefix')
                       .replace(str(prefix), ""))

    return data


def list_data(data):
    """
    """

    selections = (list(enumerate(data)))

    for datum in selections:
        print(datum)

    return selections


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


def select_panels(sample_ids):
    """
    """

    for sample in sample_ids:
        print(sample[1])


def main():
    """
    """

    bucket = 'hakmonkey-genetics-lab'
    runs_dir = 'Pipeline_Output/'
    panels_dir = 'Pipeline/Reference/panels/'

    runs = get_data(bucket = bucket, prefix = runs_dir)

    print("\n" + "#" * 18)
    print("# Available Runs #")
    print("#" * 18 + "\n")
    print("(Index, Run ID)")

    avail_runs = list_data(data = runs)

    run_id = get_choice(choices = avail_runs)

    samples = get_data(bucket = bucket, prefix = runs_dir + run_id)

    avail_samples = list_data(data = samples)

    select_panels(sample_ids = avail_samples)


if __name__ == '__main__':
    """
    """

    main()






## temp comment block

# pairs = []

# with open('barg.csv') as csv_file:
#     csv_reader = csv.reader(csv_file, delimiter='\t')
#     for row in csv_reader:
#         panels = row[1].split(",")
#         for panel in panels:
#             pairs.append([row[0], panel])

# # Will need to reverse this after running the nextflow pipeline

# os.system("sed -i 's/pairs_ch = Channel.from()/pairs_ch = Channel.from(%s)/' test-tuple.nf"%pairs)

# cmd = "sed -ri \"s/\\[([a-zA-Z0-9_.-]+),\\s([a-zA-Z0-9_.-]+)\\]/\\['\\1','\\2'\\]/g\" test-tuple.nf"

# os.system(cmd)

