#!/usr/bin/env python3

import argparse
from pysam import VariantFile


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
    )
    args = parser.parse_args()
    return args


def main():
    """
    """
    args = parse_args()
    vcf = VariantFile(args.v)
    outfile = args.v.split('.')[0]
    with open(outfile + "_ann.vcf", "w") as out:
        out.write(str(vcf.header))
        for variant in vcf.fetch():
            if 'CNCLASS' in variant.info.keys():
                gene_list = []
                with open(args.b) as f:
                    for line in f:
                        gene = line.split('\t')
                        check = (
                            gene[0] == variant.contig
                            and (
                                (
                                    int(gene[1]) >= variant.start
                                    and int(gene[1]) <= variant.stop
                                ) or
                                (
                                    int(gene[2]) >= variant.start
                                    and int(gene[2]) <= variant.stop
                                )
                            )
                        )
                        if check:
                            gene_list.append(gene[3].strip())
                genes = ","
                variant.info['GENES'] = genes.join(gene_list)
            out.write(str(variant))
    out.close()


if __name__ == '__main__':
    main()