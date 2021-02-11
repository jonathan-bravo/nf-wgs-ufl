#!/usr/bin/env python3

import os
import boto3
from data_ops import enumerate_data
from data_ops import list_data
from data_ops import get_choice
from data_ops import get_data
from data_ops import select_run


def usage():
    """
    """


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


def launch_nextflow(test_pairs, bucket, runs_dir, run_id, panels_dir):
    """
    """

    pairs = "sed -i 's/pairs_ch = Channel.from()/pairs_ch = Channel.from({test_lists})/' reporting-ufl.nf".format(test_lists = test_pairs)

    os.system(pairs)

    quotes = "sed -ri \"s/\\[([a-zA-Z0-9_.-]+),\\s([a-zA-Z0-9_.-]+)\\]/\\['\\1','\\2'\\]/g\" reporting-ufl.nf"

    os.system(quotes)

    launch = "sudo nextflow run reporting-ufl.nf -work-dir s3://{bucket}/{runs_dir}/_work/ --run_dir 's3://{bucket}/{runs_dir}/{run_id}' --panels_dir 's3://{bucket}/{panels_dir}' --bucket 's3://{bucket}' -resume".format(
        bucket = bucket,
        runs_dir = runs_dir,
        run_id = run_id,
        panels_dir = panels_dir)

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
        runs_dir = runs_dir.strip('/'),
        run_id = run_id.strip('/'),
        panels_dir = panels_dir.strip('/'))


if __name__ == '__main__':
    main()