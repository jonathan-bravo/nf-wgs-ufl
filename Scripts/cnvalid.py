#!/usr/bin/env python3

import argparse
from os import path
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
        '-v',
        metavar = '--BENCHMARK_VCF',
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


def parse_bench_vcf(vcf):
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


def create_file(sample_id):
    """Check presence of file and create file name.

    This function checks if a report exists for a given sample and then
    quits if `True`. If the report is not present then the file name
    is returned so a new report can be generated.

    Keyword arguments:

    sample_id -- the specified sample id from `args.s`

    Return:

    file_name -- the name of the report to be written
    """
    file_name = f'{sample_id}_cnvalid_report.vcf'
    if path.exists(file_name):
        print("The report file `" + str(file_name) + "` already exists.")
        quit()
    return file_name


def check_overlap(contig_chr, contig_range, bench):
    """This function checks if two CNVs overlap.

    First we check if the to CNVs are in the same chromosome. Then we check
    that the CNV from a benchmark file is contained within the CNV contig
    from the cnv_contig_report. The hit length is set to 80 percent because
    if 80% of a hit from the file overlaps our contig we want to count it as
    a hit.

    Keyword Arguments:

    contig_chr   -- the chromosome of the CNV contig
    contig_range -- the range of the CNV contig
    bench        -- the current entry from the benchmark VCF file

    Return:

    success -- a boolean value of either True or False
    """
    if bench[0] == contig_chr:
        hit_len = (bench[2] - bench[1]) * 0.8
        check = (
            bench[1] in contig_range
            and bench[1] + hit_len in contig_range
            or bench[2] - hit_len in contig_range
            and bench[2] in contig_range
            or bench[1] in contig_range
            and bench[2] in contig_range
        )
        if check:
            return True
    return False


def get_overlaps(contig, benchmark):
    """
    First we check if the called CNVs are both in the same chromosome. If this
    is true we then check to see if the CNV from the second is contained within
    the first. If all the checks are passed then we return a boolean value of
    `True` or `False`.
    """
    overlaps = []
    contig_chr = contig[0]
    contig_range = range(int(contig[1]), int(contig[2])+1)
    for i in range(len(benchmark)):
        if check_overlap(contig_chr, contig_range, benchmark[i]):
            overlaps.append(benchmark[i])
    return overlaps


def get_positives(contig, overlaps):
    """
    """
    # yes, we need to be able to report what we got wrong (false positives)
    # and what we missed (false negatives). The other benchmarking tools
    # produced a vcf of those things. The same would be great.
    cnv_type = contig[3]
    tp = 0
    fp = 0
    fp_list = []
    for i in range(len(overlaps)):
        if overlaps[i][3] == cnv_type:
            tp += 1
        else:
            fp += 1
            fp_list.append(overlaps[i])
    return (
        tp,
        fp,
        fp_list
    )


def get_false_negative(benchmark, all_overlaps):
    """
    """
    fn_list = list(set(benchmark) - set(all_overlaps))
    fn = len(fn_list)
    return (
        fn,
        fn_list
    )


def get_validation_values(contigs, benchmark):
    """
    The main loop
    """
    true_pos = 0
    false_pos = 0
    misses = []

    all_overlaps = []
    for contig in contigs:
        overlaps = get_overlaps(contig, benchmark)
        tp, fp, fp_list = get_positives(contig, overlaps)
        true_pos += tp
        false_pos += fp
        for miss in fp_list:
            misses.append(miss)
        for overlap in overlaps:
            all_overlaps.append(overlap)
    
    false_neg, fn_list = get_false_negative(benchmark, all_overlaps)
    for miss in fn_list:
        misses.append(miss)

    true_pos_base = len(benchmark)
    return (
        true_pos,
        false_pos,
        true_pos_base,
        false_neg,
        misses
    )


def data_tuple(tp, fp, fn, tp_base):
    """
    """
    ppv = tp / (tp + fp)
    tpr = tp_base / (tp_base + fn)
    data = (tp, fp, fn, ppv, tpr)
    return data


def generate_report(file_name, results, misses):
    """
    """
    pprint(misses)
    pprint(file_name)
    pprint(results)


def main():
    """
    """
    args = parse_args()
    file_name = create_file(args.s)
    benchmark = parse_bench_vcf(VariantFile(args.v))
    contigs = parse_cnv_contigs(args.c)
    tp, fp, tp_base, fn, misses = get_validation_values(contigs, benchmark)
    results = data_tuple(tp, fp, fn, tp_base)
    generate_report(file_name, results, misses)


if __name__ == '__main__':
    main()


#(1259, 1516, 34647, 0.45369369369369367, 0.5191856673004066)