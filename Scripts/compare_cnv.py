#!/usr/bin/env python3

import argparse
from pysam import VariantFile
from concurrent.futures import ProcessPoolExecutor

def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -b, --BENCHMARK_VCF -- name of the benchmark VCF file
    -v, --SAMPLE_VCF    -- name of the sample VCF file
    -t, --THREADS       -- number of cpus to use, Default[1]
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputs for checking TP, FP, TN, FN values of cnv data'
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
        help = 'the sample VCF file',
        required = True
    )
    parser.add_argument(
        '-t',
        metavar = '--THREADS',
        type = str,
        help = 'number of cpus',
        required = False
    )
    args = parser.parse_args()
    return args


def parse_vcf(vcf, chrom):
    """Parse input VCF files.

    This function grabs the contig, start, stop, and alternate information
    for each cnv entry in the input VCF file, and returns a list of tuples
    that contain this data.

    Keyword arguments:

    vcf   -- input vcf file to parse
    chrom -- the current chromosome being checked in `compare_cnvs()`

    Return:

    cnvs -- a list of tuples of the cnv data
    """
    cnvs = []
    for cnv in vcf.fetch():
        if cnv.contig == chrom:
            if cnv.alts == None: alt = '.'
            else: alt = cnv.alts[0]
            cnvs.append((cnv.contig, cnv.start, cnv.stop, alt))
    return cnvs


def compare_cnvs(chrom, bench, sample):
    """Comparing the benchmark and sample VCF cnv files.

    This function takes the two input VCF files and the specific chromosome
    from the multiprocessing and getting the initial accuracy metrics.

    Keyword arguments:

    chrom  -- the current chromosome being parsed
    bench  -- the input benchmark VCF file
    sample -- the input sample VCF file

    Return:

    tp      -- the number of true positives
    fp      -- the number of false positives
    fp_list -- a list of the false positive entries
    fn_list -- a list of the false negatice entries
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
    """Using multiprocessing to leverage multi-core cpus.

    This function takes the given input lists and a number of cpus and chunks
    the data by chromosome, and returns a list of lists, each containing
    the returned result from `compare_cnvs()`.

    Keyword arguments:

    chrom_tup -- tuple of chromosomes defined in `main()`
    bench     -- list of length of `chrom_tup` of args.b
    sample    -- list of length of `chrom_tup` of args.v
    cpus      -- the number of cpus to use; Default[1]

    Return:

    results -- list of results from multiprocessing
    """
    if cpus == None: cpus = 1
    with ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(compare_cnvs, chrom_tup, bench, sample)
    return results


def parse_results(results, bench_vcf):
    """Create precision and sensativity values.

    This function takes the returned results from `chunk_compare()` and
    merges them into single values. It also generate the ppv, or precision
    and tpr, or sensativity for generating the report.

    Keyword arguments:

    results   -- the returned results from `chunk_compare()`
    bench_vcf -- the input benchmark vcf

    Return:

    tp      -- the number of true positives
    fp      -- the number of false positives
    fn      -- the number of false negatives
    ppv     -- the percision
    tpr     -- the sensativity
    fp_list -- a list of the false positive entries
    fn_list -- a list of the false negatice entries
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
    ppv = round(tp / (tp + fp), 3)
    tpr = round(tp_base / (tp_base + fn), 3)
    return(tp, fp, fn, ppv, tpr, fp_list, fn_list)


def make_outfile(parsed_results, bench, sample):
    """Generate the outfile.

    Keyword arguments:

    parsed_results -- the returned value from `parse_results()`
    bench          -- the input benchmark vcf
    sample         -- the input sample vcf
    """
    bench_name = bench.split('_')[0]
    sample_name = sample.split('_')[0]
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
