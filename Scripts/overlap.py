#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile



parser = argparse.ArgumentParser(description='Compare overlaps with a CNV vcf file and a standard SNV vcf file')
parser.add_argument(
    '-c',
    metavar  = '--CNV',
    type     = str,
    help     = 'name of the CNV vcf',
    required = True
)
parser.add_argument(
    '-s',
    metavar  = '--SNV',
    type     = str,
    help     = 'name of the SNV file',
    required = True
)
parser.add_argument(
    '-o',
    metavar  = '--OUT',
    type     = str,
    help     = 'name of the outputfile',
    required = True
)

args = parser.parse_args()
cnv = VariantFile(args.c) # vcf = pysam.VariantFile('NQ-21-01-BC708501-385_S25_filtered_cnv.vcf.gz')
snv = VariantFile(args.s) # vcf2 = pysam.VariantFile('strelka2_variants.vcf.gz')


for entry in vcf.fetch():
    if entry.info["GENE"] != '.' and entry.info["CNCLASS"] != "CN2":
        test.append((entry.chrom, entry.start, entry.stop, entry.info["CNCLASS"], entry.info["CYTO"], entry.info["GENE"]))

def main():
    """
    """


if __name__ == '__main__':
    main()