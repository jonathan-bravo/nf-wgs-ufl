#!/usr/bin/env python3

import argparse
from pysam import VariantFile
from concurrent.futures import ProcessPoolExecutor
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


def parse_vcf(vcf, chrom):
    """
    """
    cnvs = []
    for cnv in vcf.fetch():
        if cnv.contig == chrom:
            cnvs.append((cnv.contig, cnv.start, cnv.stop, cnv.alts[0]))
    return cnvs


def compare_cnvs(chrom, bench, sample):
    """
    """
    bench_vcf = parse_vcf(VariantFile(bench), chrom)
    sample_vcf = parse_vcf(VariantFile(sample), chrom)
    tp = 0
    fp = 0
    fp_list = []
    fn_list = [x for x in bench_vcf if x not in sample_vcf]
    for cnv in sample_vcf:
        if cnv in bench_vcf: tp += 1
        else:
            fp += 1
            fp_list.append(cnv)
    return (tp, fp, fp_list, fn_list)


def chunk_compare(chrom_tup, bench, sample, cpus):
    """
    """
    if cpus == None: cpus = 1
    with ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(compare_cnvs, chrom_tup, bench, sample)
    return results


def parse_results(results, bench_vcf):
    """
    """
    tp = 0
    fp = 0
    tp_base = 0
    fp_list = []
    fn_list = []
    for _ in VariantFile(bench_vcf).fetch(): tp_base += 1
    
    for result in results:
        tp += result[0]
        fp += result[1]
        for fp_cnv in result[2]: fp_list.append(fp_cnv)
        for fn_cnv in result[3]: fn_list.append(fn_cnv)
    fn = len(fn_list)
    ppv = tp / (tp + fp)
    tpr = tp_base / (tp_base + fn)
    return(tp, fp, fn, ppv, tpr, fp_list, fn_list)


def make_outfile(parsed_results, bench, sample):
    """
    """
    bench_name = bench.split('_', 0)
    sample_name = sample.split('_', 0)
    f = open(f'{sample_name}_vs_{bench_name}.txt', "w")
    f.write('VALUES\n\n')
    f.write(f'True Positive: {parsed_results[0]}\n')
    f.write(f'False Positive: {parsed_results[1]}\n')
    f.write(f'False Negatives: {parsed_results[2]}\n')
    f.write(f'Precision: {parsed_results[3]}\n')
    f.write(f'Sensitivity: {parsed_results[4]}\n\n')
    f.write(f'FALSE POSITIVES\n\n')
    for fp in parsed_results[5]: f.write(f'{fp}\n')
    f.write('\nFALSE NEGATIVES\n\n')
    for fn in parsed_results[6]: f.write(f'{fn}\n')
    f.close()


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
    bench = [args.b] * len(chrom_tup)
    sample = [args.v] * len(chrom_tup)
    cpus = int(args.t)
    results = chunk_compare(chrom_tup, bench, sample, cpus)
    parsed_results = parse_results(results, args.b)
    make_outfile(parsed_results, args.b, args.v)


if __name__ == '__main__':
    main()
