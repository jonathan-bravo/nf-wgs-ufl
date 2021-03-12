#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile


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


def parse_bench_bed():
    """
    """
    return tmp


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


def generate_report(file_name, X):
    """
    """


def check_overlap(contig_chr, contig_range, bench):
    """This function checks if two CNVs overlap.

    First we check if the to CNVs are in the same chromosome. Then we check
    that the CNV from a benchmark file is contained within the CNV contig
    from the cnv_contig_report.

    Keyword Arguments:

    contig_chr   -- the chromosome of the CNV contig
    contig_range -- the range of the CNV contig
    bench        -- the current entry from the benchmark VCF file

    Return:

    success -- a boolean value of either True or False
    """
    if bench[0] == contig_chr:
        if bench[1] in contig_range and bench[2] in contig_range:
            return True
    return False


def get_overlaps(contig, benchmark, start_index):
    """
    First we check if the called CNVs are both in the same chromosome. If this
    is true we then check to see if the CNV from the second is contained within
    the first. If all the checks are passed then we return a boolean value of
    `True` or `False`.
    """
    overlaps = []
    contig_chr = contig[0]
    contig_range = range(int(contig[1]), int(contig[2])+1)
    for i in range(start_index, len(benchmark)+1):
        end_index = i
        if check_overlap(contig_chr, contig_range, benchmark[i]):
            overlaps.append(benchmark[i])
        else:
            end_index += 1
            break
    return overlaps, end_index


def generate_percision(contig, overlaps):
    """
    """
    # Precision (PPV): TP-call/(TP-call+FP)
    # Sensitivity (TPR): true positive rate = TP-base/(TP-base+FN)
    ## as we did for the SNVs
    # yes, we need to be able to report what we got wrong (false positives)
    # and what we missed (false negatives). The other benchmarking tools
    # produced a vcf of those things. The same would be great.
    cnv_type = contig[3]
    true_pos = 0
    false_pos = 0
    false_pos_list = []
    for i in range(len(overlaps)):
        if overlaps[i][3] == cnv_type:
            true_pos += 1
        else:
            false_pos += 1
            false_pos_list.append(overlaps[i])
    percision = true_pos/(true_pos+false_pos)
    return (true_pos, false_pos, false_pos_list, percision)


def generate_sensitivity(contig, overlaps):
    """
    """
    cnv_type = contig[3]
    true_neg = 0
    false_neg = 0
    false_neg_list = []
    for i in range(len(overlaps)):
        if overlaps[i][3] == cnv_type:
            true_pos += 1
        else:
            false_pos += 1
            false_pos_list.append(overlaps[i])
    sensitivity = true_neg/(true_neg+false_neg)
    return (true_neg, false_neg, false_neg_list, sensitivity)


def validation(contigs, benchmark):
    """
    The main loop
    """
    start_index = 0
    for contig in contigs:
        overlaps, start_index = get_overlaps(contig, benchmark, start_index)
        true_pos, false_pos, false_pos_list, percision = generate_percision(contig, overlaps)
        true_neg, false_neg, false_neg_list, sensativity = generate_sensitivity(contig, overlaps)
    return temp_value


def main():
    """
    """
    args = parse_args()
    benchmark = parse_bench_vcf(VariantFile(args.b))
    confidence = parse_bench_bed()
    contigs = parse_cnv_contigs(args.c)
    file_name = create_file(args.s)
    temp_value = validation(contigs, benchmark)
    generate_report(file_name, temp_value)


if __name__ == '__main__':
    main()
