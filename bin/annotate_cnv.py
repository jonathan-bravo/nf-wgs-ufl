#!/usr/bin/env python3

import argparse
import pandas as pd
from pysam import VariantFile
from pysam import bcftools
from pysam import tabix_compress


def parse_args():
    """
    """
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-v',
        metavar = '--VCF',
        type = str,
        help = 'The input sample vcf',
        required = True
    )
    parser.add_argument(
        '-b',
        metavar = '--BED',
        type = str,
        help = 'hg19 genes bed file',
        required = True
    ),
    parser.add_argument(
        '-o',
        metavar = '--OUT_FILE',
        type = str,
        help = 'annotated vcf output',
        required = True
    )
    args = parser.parse_args()
    return args


def index_input_vcf(vcf):
    """
    """
    tabix_compress(vcf, f'{vcf}.gz')
    bcftools.index('--tbi', f'{vcf}.gz')


def load_bed(bed):
    """
    """
    headers = ['Contig', 'Start', 'Stop', 'Gene']
    df = pd.read_csv(bed, sep='\t', names=headers)
    return df


def get_genes(variant, bed):
    """
    """
    genes = bed.loc[
        (bed.Contig == variant.contig) & (
            (
                (bed.Start >= variant.start) & 
                (bed.Start <= variant.stop)
            ) | (
                (bed.Stop >= variant.start) & 
                (bed.Stop <= variant.stop)
            )
        )
    ].Gene.tolist()
    return ",".join(genes)


def write_updated_vcf(out_file, vcf, bed):
    """
    """
    with open(out_file, 'w') as out_vcf:
        out_vcf.write(str(vcf.header))
        for variant in vcf.fetch():
            genes = get_genes(variant, bed)
            variant.info.update({'GENES': genes})
            out_vcf.write(str(variant))            


def main():
    """
    """
    args = parse_args()
    #sample_id = args.v.split('.')[0]
    index_input_vcf(args.v)
    vcf = VariantFile(f'{args.v}.gz')
    bed = load_bed(args.b)
    write_updated_vcf(args.o, vcf, bed)    


if __name__ == '__main__':
    main()