#!/usr/bin/env python3

import os
import boto3
from data_ops import get_data
from data_ops import select_run


def usage():
    """
    """


def launch_nextflow(bucket, run_id, output_dir):
    """
    """

    cmd = "sudo nextflow run multiqc-ufl.nf -work-dir s3://hakmonkey-genetics-lab/Pipeline_Output/work/ --run_id '{run_id}' --run_dir 's3://{bucket}/{runs_dir}/{run_id}'".format(
        run_id = run_id,
        bucket = bucket,
        runs_dir = output_dir)

    os.system(cmd)


def main():
    """
    """
    
    bucket = 'hakmonkey-genetics-lab'
    output_dir = 'Pipeline_Output/'

    run = get_data(bucket = bucket, prefix = output_dir)
    run_id = select_run(runs = run)

    launch_nextflow(
        bucket = bucket,
        run_id = run_id,
        output_dir = output_dir)


if __name__ == '__main__':
    main()