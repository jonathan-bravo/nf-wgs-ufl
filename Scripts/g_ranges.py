#!/usr/bin/env python3

import argparse
import json

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
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = '',
        required = True
    )
    args = parser.parse_args()
    return args


def main():
    calls = []
    classes = []
    args = parse_args()
    with open(args.j, 'r') as jason_file:
        data = json.load(jason_file)
        for classification in data['cnv']:
            for cnv in data['cnv'][classification]:
                calls.append(f"{cnv['chrom']}:{cnv['start']}-{cnv['stop']}")
                if 'DEL' in cnv['alt']: classes.append('0')
                elif 'DUP' in cnv['alt']: classes.append('4')
                else: classes.append('2')
    with open(f'{args.s}_cnv_class.csv', 'w') as cnv_class:
        for c in classes:
            cnv_class.write(f'{c}\n')
    with open(f'{args.s}_cnv_calls.csv', 'w') as cnv_calls:
        for c in calls:
            cnv_calls.write(f'{c}\n')


if __name__ == '__main__':
    main()
