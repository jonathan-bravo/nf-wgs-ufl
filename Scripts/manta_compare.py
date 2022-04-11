#!/usr/bin/env python3

import argparse
import pandas as pd
from pysam import VariantFile
from pprint import pprint
from json_to_csv import read_json

def parse_args():
    """Parse input arguments.

    Keyword arguments:
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-j',
        metavar = '--JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    parser.add_argument(
        '-v',
        metavar = '--MANTA_VCF',
        type = str,
        help = 'the input Manta vcf for the sample',
        required = True
    )
    args = parser.parse_args()
    return args


def parse_manta(vcf):
    """
    """
    contigs = (
        '1','2','3','4','5','6','7','8',
        '9','10','11','12','13','14','15','16',
        '17','18','19','20','21','22','X','Y'
    )
    cnvs = []
    other_data = []
    for variant in vcf.fetch():
        svtype = variant.info['SVTYPE']
        if variant.contig in contigs and ('DEL' in svtype or 'DUP' in svtype):
            cnvs.append({
                'Chrom': f'chr{variant.contig}',
                'Start': variant.start,
                'Stop': variant.stop,
                'Alt Allele': svtype,
                'Length': variant.info['SVLEN'][0]
            })
        else:
            other_data.append({
                'Chrom': f'chr{variant.contig}',
                'Start': variant.start,
                'Stop': variant.stop,
                'Ref Allele': variant.ref,
                'Alt Allele': variant.alts[0],
                'Type': svtype
            })
    return (pd.DataFrame(cnvs), pd.DataFrame(other_data))


def compare_manta_to_cnmops(contig, start, alt, length, manta_cnvs):
    """
    """
    cnv_range = range(length - 7000, length + 7000)
    start_range = range(start - 7000, start + 7000)
    for i, manta_cnv in manta_cnvs.iterrows():
        match = (
            manta_cnv.Chrom == contig
            and manta_cnv['Alt Allele'] == alt
            and manta_cnv.Start in start_range
            and abs(manta_cnv.Length) in cnv_range
        )
        if match:
            return (
                manta_cnv.Chrom,
                manta_cnv.Start,
                manta_cnv.Stop,
                manta_cnv['Alt Allele'],
                manta_cnv.Length
            )
    return 'NA'


def main():
    args = parse_args()
    file_name = args.j.split('.')[0]
    json_data = read_json(args.j)
    manta_cnvs, manta_data = parse_manta(VariantFile(args.v))
    cnmops_data = json_data[2]
    results = [compare_manta_to_cnmops(a,b,c,d,manta_cnvs) for a,b,c,d in zip(cnmops_data['Chrom'], cnmops_data['Start'], cnmops_data['Alt Allele'], cnmops_data['Length'])]
    cnmops_data['Manta Data'] = results
    with pd.ExcelWriter(f'{file_name}_CNVS.xlsx') as writer:
        cnmops_data.to_excel(writer, sheet_name = 'CNVs', index = False)
        manta_data.to_excel(writer, sheet_name = 'Manta', index = False)


if __name__ == '__main__':
    main()