#!/usr/bin/env python3

import argparse
from pysam import VariantFile
from concurrent.futures import ProcessPoolExecutor
from numpy import setdiff1d
from pprint import pprint

def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -v, --BENCHMARK_VCF     -- name of the benchmark VCF file
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
        metavar = '--BENCHMARK_VCF',
        type = str,
        help = 'the benchmark VCF file',
        required = True
    )
    parser.add_argument(
        '-v',
        metavar = '--SAMPLE_VCF',
        type = str,
        help = '',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'sample id to be included in report name',
        required = False
    )
    parser.add_argument(
        '-t',
        metavar = '--THREADS',
        type = str,
        help = '',
        required = False
    )
    args = parser.parse_args()
    return args


def parse_vcf(vcf, chrom = None):
    """
    """
    cnvs = []
    for cnv in vcf.fetch():
        cnvs.append((cnv.contig, cnv.start, cnv.stop, cnv.alts))
    return cnvs


def compare_cnvs(chrom, bench_cnvs, sample_cnvs):
    """
    """
    bench_cnv_for_chr = [x for x in bench_cnv if x[0] == chrom]
    sample_cnv_for_chr = [x for x in sample_cnv if x[0] == chrom]
    tp = 0
    fp = 0
    fp_list = []
    for cnv in sample_cnv_for_chr:
        if cnv in bench_cnv_for_chr: tp += 1
        else:
            fp += 1
            fp_list.append(cnv)
    return (tp, fp, fp_list)


def chunk_compare(chrom_tup, bench_list, sample_list, cpus):
    """
    """
    if cpus == None: cpus = 1
    with ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(compare_cnvs, chrom_tup, bench_list, sample_list)
    return results


def main():
    """
    """
    chrom_tup = (
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr22', 'chr21', 'chrX', 'chrY', 'chrM'
    )
    args = parse_args()
    bench_vcf = parse_vcf(VariantFile(args.b))
    sample_vcf = parse_vcf(VariantFile(args.s))
    bench_list = [bench_vcf] * len(chrom_tup)
    sample_list = [sample_vcf] * len(chrom_tup)
    cpus = int(args.t)
    results = chunk_compare(chrom_tup, bench_list, sample_list, cpus)

    tp = 0
    fp = 0
    fp_list = []
    for result in results:
        tp += result[0]
        fp += result[1]
        for fp_cnv in result[2]:
            fp_list.append(fp_cnv)

    bench_cnvs = parse_vcf(VariantFile(args.b))
    sample_cnvs = parse_vcf(VariantFile(args.v))
    fn_list = setdiff1d(bench_cnvs, sample_cnvs)
    fn = len(fn_list)

    print(f'tp: {tp}, fp: {fp}, fn: {fn}')


if __name__ == '__main__':
    main()
