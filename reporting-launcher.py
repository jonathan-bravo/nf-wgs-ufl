#!/usr/bin/env python3

import os
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

    if 'panels' in prefix:
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


def select_panels(sample_list, panel_list):
    """
    """

    test_pairs = []

    avail_samples = enumerate_data(data = sample_list)
    avail_panels = enumerate_data(data = panel_list)

    for sample in avail_samples:
        print("\n" + "#" * 20)
        print("# Available Panels #")
        print("#" * 20 + "\n")
        print("(Index, Panel)")
        list_data(data = avail_panels)

        print("\n Sample:")
        print(sample[1] + "\n")

        test_num = int(input("How many panels would you like to run?: "))

        for i in range(test_num):
            print("panel number " + str(i+1))

            test = get_choice(choices = avail_panels)
            test_pairs.append([sample[1].strip('/'), test])

    return test_pairs


def launch_nextflow(test_pairs, bucket, runs_dir, run_id, panels_dir):
    """
    """

    pairs = "sed -i 's/pairs_ch = Channel.from()/pairs_ch = Channel.from({test_lists})/' reporting-ufl.nf".format(test_lists = test_pairs)

    os.system(pairs)

    quotes = "sed -ri \"s/\\[([a-zA-Z0-9_.-]+),\\s([a-zA-Z0-9_.-]+)\\]/\\['\\1','\\2'\\]/g\" reporting-ufl.nf"

    os.system(quotes)

    launch = "sudo nextflow run reporting-ufl.nf -c ~/Documents/nextflow.config.bk -work-dir s3://{bucket}/{runs_dir}/work/ --run_dir 's3://{bucket}/{runs_dir}/{run_id}' --panels_dir 's3://{bucket}/{panels_dir}'".format(
        bucket = bucket,
        runs_dir = runs_dir.strip('/'),
        run_id = run_id.strip('/'),
        panels_dir = panels_dir.strip('/'))

    os.system(launch)

    trash = "sudo nextflow clean -f"

    os.system(trash)

    reset = "sed -i 's/^pairs_ch = Channel.from(.*$/pairs_ch = Channel.from()/' reporting-ufl.nf"

    os.system(reset)


def main():
    """
    """

    bucket = 'hakmonkey-genetics-lab'
    runs_dir = 'Pipeline_Output/'
    panels_dir = 'Pipeline/Reference/panels/'

    ################################
    # Listing & Getting Run Choice #
    ################################

    runs = get_data(bucket = bucket, prefix = runs_dir)

    run_id = select_run(runs = runs)

    ################################
    # Creating Sample/ Panel Pairs #
    ################################

    samples = get_data(bucket = bucket, prefix = runs_dir + run_id)
    
    panels = get_data(bucket = bucket, prefix = panels_dir)
    while('' in panels):
        panels.remove('')

    test_pairs = select_panels(sample_list = samples, panel_list = panels)

    #################################
    # Running The Nextflow Pipeline #
    #################################

    launch_nextflow(
        test_pairs = test_pairs,
        bucket = bucket,
        runs_dir = runs_dir,
        run_id = run_id,
        panels_dir = panels_dir)


if __name__ == '__main__':
    main()