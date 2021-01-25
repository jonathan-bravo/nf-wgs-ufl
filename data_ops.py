#!/usr/bin/env python3

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

    if 'panels' in prefix or 'Fastqs' in prefix:
        for obj in result.get('Contents'):
            data.append(obj.get('Key')
                           .replace(str(prefix), ""))
    else:
        for obj in result.get('CommonPrefixes'):
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