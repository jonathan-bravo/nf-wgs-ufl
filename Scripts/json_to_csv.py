#!/usr/bin/env python3

import argparse
import json
import csv
import pandas as pd
from pprint import pprint

def parse_args():
    """Parse input arguments.

    Keyword arguments:
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputs for checking TP, FP, TN, FN values of cnv data'
    )
    parser.add_argument(
        '-j',
        metavar = '--JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    file_name = args.j.split('.')[0]
    with open(args.j, 'r') as jason_file:
        data = json.load(jason_file)
        p_snp = data['snp']['pathogenic_snp']
        l_snp = data['snp']['likely_pathogenic_snp']
        p_sv = data['sv']['pathogenic_sv']
        l_sv = data['sv']['likely_pathogenic_sv']
        p_cnv = data['cnv']['pathogenic_cnv']
        l_cnv = data['cnv']['likely_pathogenic_cnv']
        exp = data['exp']['all_expansions']


if __name__ == '__main__':
    main()
