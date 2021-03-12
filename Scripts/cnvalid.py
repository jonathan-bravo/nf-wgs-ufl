#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile

from pprint import pprint


def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -b, --BENCHMARK         -- name of the benchmark VCF file
    -c, --CNV_CONTIG_REPORT -- input report generated from `cnv_contigs.py`
    -s, --SAMPLE_ID         -- sample id to be included in report name

    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputds for checking TP, FP, TN, FN values of cnv data'
    )
    parser.add_argument(
        '-b',
        metavar = '--BENCHMARK',
        type = str,
        help = 'name of the benchmark VCF file',
        required = True
    )
    parser.add_argument(
        '-c',
        metavar = '--CNV_CONTIG_REPORT',
        type = str,
        help = 'the input report generated from `cnv_contigs.py`',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'sample id to be included in report name',
        required = True
    )
    args = parser.parse_args()
    return args


def parse_cnv_contigs(report):
    """This function will parses the cnv_contigs_report.

    First a list is created that will contain a tuple for each CNV contig in
    the report file. Each contig is pulled out and only the chromosome,
    start, stop, and CNV type is kept in the list so that the data matches
    the data generated in `parse_vcf()`.

    Keyword arguments:

    report -- the input report generated from `cnv_contigs.py`

    Return:

    contigs -- the chrom, start, stop, and CNV type from our contig report
    """
    contigs = []
    with open(report) as r:
        for entry in r:
            entry = entry.strip()
            entry = tuple(entry.split('\t'))
            contigs.append((
                entry[0].replace('CHROM: ', ''),
                int(entry[3].replace('START: ', '')),
                int(entry[4].replace('STOP: ', '')),
                entry[2].replace('CNV: ', '')
            ))
    return contigs


def parse_vcf(vcf):
    """This function parses the benchmark VCF.

    First a list is created that will contain a tuple for each CNV in the
    benchmark file. Then each entry in the benchmark VCF is read and pulled
    into the list if it is a true CNV (deletion or duplication). The returned
    list of tuples has the chromosome, start position, stop position, and CNV
    type for all CNVs in the benchmark file.

    Keyword arguments:

    vcf -- the input benchmark VCF

    Return:

    bench_cnvs -- the chrom, start, stop, and CNV type from our benchmark VCF
    """
    bench_cnvs = []
    for entry in vcf.fetch():
        if entry.info["SVTYPE"] == "DEL" or entry.info["SVTYPE"] == "DUP":
            bench_cnvs.append((
                'chr' + entry.chrom,
                entry.start,
                entry.stop,
                entry.info["SVTYPE"]
            ))
    return bench_cnvs



def main():
    """
    """
    args = parse_args()
    benchmark = parse_vcf(VariantFile(args.b))
    contigs = parse_cnv_contigs(args.c)
    sample_id = args.s
    pprint(benchmark)
    pprint(contigs)


if __name__ == '__main__':
    main()
