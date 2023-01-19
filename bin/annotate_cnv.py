#!/usr/bin/env python3

import argparse
import pandas as pd
from datetime import datetime
from pysam import VariantFile
from pysam import bcftools
from pysam import tabix_compress

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'The input sample id',
        required = True
    )
    parser.add_argument(
        '-c',
        metavar = '--CSV',
        type = str,
        help = 'The input sample csv from cn.mops',
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
    return parser.parse_args()

def csv_to_vcf(csv, sample_id):
    header = '##fileformat=VCFv4.1\n' \
        f'##fileDate={datetime.today()}\n' \
        '##source=cn.MOPS\n' \
        '##source_version=1.36.0\n' \
        '##reference=hs37d5.fa.gz\n' \
        '##contig=<ID=1,length=249250621>\n' \
        '##contig=<ID=2,length=243199373>\n' \
        '##contig=<ID=3,length=198022430>\n' \
        '##contig=<ID=4,length=191154276>\n' \
        '##contig=<ID=5,length=180915260>\n' \
        '##contig=<ID=6,length=171115067>\n' \
        '##contig=<ID=7,length=159138663>\n' \
        '##contig=<ID=8,length=146364022>\n' \
        '##contig=<ID=9,length=141213431>\n' \
        '##contig=<ID=10,length=135534747>\n' \
        '##contig=<ID=11,length=135006516>\n' \
        '##contig=<ID=12,length=133851895>\n' \
        '##contig=<ID=13,length=115169878>\n' \
        '##contig=<ID=14,length=107349540>\n' \
        '##contig=<ID=15,length=102531392>\n' \
        '##contig=<ID=16,length=90354753>\n' \
        '##contig=<ID=17,length=81195210>\n' \
        '##contig=<ID=18,length=78077248>\n' \
        '##contig=<ID=19,length=59128983>\n' \
        '##contig=<ID=20,length=63025520>\n' \
        '##contig=<ID=21,length=48129895>\n' \
        '##contig=<ID=22,length=51304566>\n' \
        '##contig=<ID=X,length=155270560>\n' \
        '##contig=<ID=Y,length=59373566>\n' \
        '##content=cn.MOPS cnv calls\n' \
        '##ALT=<ID=DEL,' \
        'Description="Deletion relative to the reference">\n' \
        '##ALT=<ID=DUP,' \
        'Description="Region of elevated copy number relative to the reference">\n' \
        '##INFO=<ID=SVTYPE,Number=1,Type=String,' \
        'Description="Type of structural variant">\n' \
        '##INFO=<ID=END,Number=1,Type=Integer,' \
        'Description="End position of the variant described in this record">\n' \
        '##INFO=<ID=LENGTH,Number=1,Type=Integer,' \
        'Description="The length of the called copy number variant">\n' \
        '##INFO=<ID=CNCLASS,Number=1,Type=String,' \
        'Description="Class given to the I vector value">\n' \
        '##INFO=<ID=GENES,Number=1,Type=String,' \
        'Description="The target gene(s)">\n' \
        '##FORMAT=<ID=MED,Number=1,Type=Float,' \
        'Description="MED: median individual call for the copy number segment">\n' \
        '##FORMAT=<ID=MEAN,Number=1,Type=Float,' \
        'Description="MEAN: individual call for the copy number segment">\n' \
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE1\t\n'
    cnv_strings = []
    with open(csv) as f:
        next(f)
        for line in f:
            cn_data = line.split(',')[1:5] + line.split(',')[7:]
            if '0' in cn_data[6] or '1' in cn_data[6]:
                cn = 'DEL'
            elif '2' in cn_data[6]:
                cn = '.'
            else:
                cn = 'DUP'
            if '.' not in cn:
                cnv_strings.append(
                    f'{cn_data[0]}' \
                    f'\t{cn_data[1]}' \
                    f'\tcn.MOPS:{cn_data[0]}:{cn_data[1]}-{cn_data[2]}' \
                    '\tN' \
                    f'\t{cn}' \
                    '\t.' \
                    '\tPASS' \
                    '\tSVTYPE=CNV;' \
                    f'END={cn_data[2]};' \
                    f'LENGTH={cn_data[3]};' \
                    f'CNCLASS={cn_data[6]};' \
                    'GENES=.' \
                    '\tMED:MEAN' \
                    f'\t{cn_data[4]}:{cn_data[5]}\n'
                )
        with open(f'{sample_id}.cnv.vcf' ,'w') as v:
            for cnv in cnv_strings:
                v.write(cnv)
        tabix_compress(f'{sample_id}.cnv.vcf', f'{sample_id}.cnv.vcf.gz')
        bcftools.index('--tbi', f'{sample_id}.cnv.vcf.gz')

def load_bed(bed):
    headers = ['Contig', 'Start', 'Stop', 'Gene']
    df = pd.read_csv(bed, sep='\t', names=headers)
    return df

def get_genes(variant, bed):
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
    with open(out_file, 'w') as out_vcf:
        out_vcf.write(str(vcf.header))
        for variant in vcf.fetch():
            genes = get_genes(variant, bed)
            variant.info.update({'GENES': genes})
            out_vcf.write(str(variant))            

def main():
    args = parse_args()
    csv_to_vcf(args.c, args.s)
    vcf = VariantFile(f'{args.s}.cnv.vcf.gz')
    bed = load_bed(args.b)
    write_updated_vcf(args.o, vcf, bed)    

if __name__ == '__main__':
    main()