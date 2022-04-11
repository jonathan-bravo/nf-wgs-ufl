#!/usr/bin/env python3

import argparse
from datetime import date
from os import system

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
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = '',
        required = True
    )
    parser.add_argument(
        '-t',
        metavar = '--TSV',
        type = str,
        help = '',
        required = True
    )
    args = parser.parse_args()
    return args


def get_locus_data(tsv_file, sample_id):
    """
    """
    with open(tsv_file, 'r') as dataset:
        locus_data = []
        for locus in dataset:
            loci = locus.split('\t')
            counts = loci[6].split(',')
            include = False
            for count in counts:
                if sample_id in count.split(':')[0]:
                    include = True
                    loci[6] = count
            if include: locus_data.append(loci)
    return locus_data


def to_vcf(sample_id, header, contigs, locus_data):
    """
    """
    with open(f'{sample_id}_ehd.vcf', 'w') as vcf:
        for entry in header:
            vcf.write(f'{entry}\n')
        for locus in locus_data:
            contig = locus[0]
            start = locus[1]
            end = locus[2]
            motif = locus[3]
            pvalue = locus[4]
            bonf_pvalue = locus[5]
            counts = locus[6].split(':')[1].strip()
            if contig in contigs:
                vcf.write(f'{contig}\t{start}\t.\tN\t{motif}\t{pvalue}\tPASS\tEND={end}\tBP:CT\t{bonf_pvalue}:{counts}\n')


def sort_vcf(sample_id):
    """
    """
    bcftools_command = f'bcftools sort -o {sample_id}_ehd_sorted.vcf {sample_id}_ehd.vcf'
    clean_command = f'rm {sample_id}_ehd.vcf'
    adjust_name_command = f'mv {sample_id}_ehd_sorted.vcf {sample_id}_ehd.vcf'
    system(bcftools_command)
    system(clean_command)
    system(adjust_name_command)


def main():
    """
    """
    HEADER = (
        '##fileformat=VCFv4.1',
        f'##fileDate={date.today()}',
        '##source=ExpansionHunterDenovo',
        '##source_version=0.90',
        '##reference=hs37d5.fa.gz',
        '##contig=<ID=1,length=249250621>',
        '##contig=<ID=2,length=243199373>',
        '##contig=<ID=3,length=198022430>',
        '##contig=<ID=4,length=191154276>',
        '##contig=<ID=5,length=180915260>',
        '##contig=<ID=6,length=171115067>',
        '##contig=<ID=7,length=159138663>',
        '##contig=<ID=8,length=146364022>',
        '##contig=<ID=9,length=141213431>',
        '##contig=<ID=10,length=135534747>',
        '##contig=<ID=11,length=135006516>',
        '##contig=<ID=12,length=133851895>',
        '##contig=<ID=13,length=115169878>',
        '##contig=<ID=14,length=107349540>',
        '##contig=<ID=15,length=102531392>',
        '##contig=<ID=16,length=90354753>',
        '##contig=<ID=17,length=81195210>',
        '##contig=<ID=18,length=78077248>',
        '##contig=<ID=19,length=59128983>',
        '##contig=<ID=20,length=63025520>',
        '##contig=<ID=21,length=48129895>',
        '##contig=<ID=22,length=51304566>',
        '##contig=<ID=X,length=155270560>',
        '##contig=<ID=Y,length=59373566>',
        '##content=ExpansionHunterDenovo locus calls',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the variant described in this record">',
        '##FORMAT=<ID=CT,Number=1,Type=Float,Description="depth-normalized count of reads originating inside each STR">',
        '##FORMAT=<ID=BP,Number=1,Type=Float,Description="MED: median individual call for the copy number segment">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1'
    )
    CONTIGS = (
        '1','2','3','4','5','6','7','8','9','10',
        '11','12','13','14','15','16','17','18','19','20',
        '21','22','X','Y'
    )
    args = parse_args()
    sample_id = args.s
    tsv_file = args.t
    locus_data = get_locus_data(tsv_file, sample_id)
    to_vcf(sample_id, HEADER, CONTIGS, locus_data)
    sort_vcf(sample_id)


if __name__ == '__main__':
    main()