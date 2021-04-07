#!/usr/bin/env python3

import argparse
from pysam import VariantFile
from concurrent.futures import ProcessPoolExecutor

def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -c, --CHILD_VCF  -- VCF file of the child
    -m, --MOTHER_VCF -- VCF file of the mother
    -f, --FATHER_VCF -- VCF file of the father
    -t, --THREADS    -- number of cpus to use, Default[1]
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputs for doing cnv trio analysis'
    )
    parser.add_argument(
        '-c',
        metavar = '--CHILD_VCF',
        type = str,
        help = 'the VCF file of the child',
        required = True
    )
    parser.add_argument(
        '-m',
        metavar = '--MOTHER_VCF',
        type = str,
        help = 'the VCF file of the mother',
        required = True
    )
    parser.add_argument(
        '-f',
        metavar = '--FATHER_VCF',
        type = str,
        help = 'the VCF file of the father',
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


def compare_cnvs(chrom, child, mother, father):
    """Comparing the benchmark and sample VCF cnv files.

    This function takes the two input VCF files and the specific chromosome
    from the multiprocessing and getting the initial accuracy metrics.

    Keyword arguments:

    chrom  -- the current chromosome being parsed
    child  -- the input VCF file of the child
    mother -- the input VCF file of the mother
    father -- the input VCF file of the father

    Return:

    both_parents -- number of cnvs present in both parents
    only_mother  -- number of cnvs present in only the mother
    only_father  -- number of cnvs present in only the father
    neither      -- number of cnvs present in only the child
    both_list    -- list of cnvs in both parents
    mother_list  -- list of cnvs in only mother
    father_list  -- list of cnvs in only father
    neither_list -- list of cnvs in only child
    """
    child_vcf = parse_vcf(VariantFile(child), chrom)
    mother_vcf = parse_vcf(VariantFile(mother), chrom)
    father_vcf = parse_vcf(VariantFile(father), chrom)
    both_parents = 0
    only_mother = 0
    only_father = 0
    neither = 0
    both_list = []
    mother_list = []
    father_list = []
    neither_list = []
    for cnv in child_vcf:
        if cnv in mother_vcf and cnv in father_vcf: 
            both_parents += 1
            both_list.append(cnv)
        elif cnv in mother_vcf:
            only_mother += 1
            mother_list.append(cnv)
        elif cnv in father_vcf:
            only_father += 1
            father_list.append(cnv)
        else:
            neither += 1
            neither_list.append(cnv)
    return (
        both_parents, only_mother, only_father, neither,
        both_list, mother_list, father_list ,neither_list
    )


def chunk_compare(chrom_tup, child, mother, father, cpus):
    """Using multiprocessing to leverage multi-core cpus.

    This function takes the given input lists and a number of cpus and chunks
    the data by chromosome, and returns a list of lists, each containing
    the returned result from `compare_cnvs()`.

    Keyword arguments:

    chrom_tup -- tuple of chromosomes defined in `main()`
    child     -- list of length of `chrom_tup` of args.c
    mother    -- list of length of `chrom_tup` of args.m
    father    -- list of length of `chrom_tup` of args.f
    cpus      -- the number of cpus to use; Default[1]

    Return:

    results -- list of results from multiprocessing
    """
    if cpus == None: cpus = 1
    with ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(compare_cnvs, chrom_tup, child, mother, father)
    return results


def parse_results(results):
    """Create merged result values for reporting.

    This function takes the returned results from `chunk_compare()` and
    merges them into single values.

    Keyword arguments:

    results   -- the returned results from `chunk_compare()`

    Return:

    both_parents -- number of cnvs present in both parents
    only_mother  -- number of cnvs present in only the mother
    only_father  -- number of cnvs present in only the father
    neither      -- number of cnvs present in only the child
    both_list    -- list of cnvs in both parents
    mother_list  -- list of cnvs in only mother
    father_list  -- list of cnvs in only father
    neither_list -- list of cnvs in only child
    """
    both_parents = 0
    only_mother = 0
    only_father = 0
    neither = 0
    both_list = []
    mother_list = []
    father_list = []
    neither_list = []
    for result in results:
        both_parents += result[0]
        only_mother += result[1]
        only_father += result[2]
        neither += result[3]
        for cnv in result[4]: both_list.append(cnv)
        for cnv in result[5]: mother_list.append(cnv)
        for cnv in result[6]: father_list.append(cnv)
        for cnv in result[7]: neither_list.append(cnv)
    return (
        both_parents, only_mother, only_father, neither,
        both_list, mother_list, father_list ,neither_list
    )


def make_outfile(parsed_results, child, mother, father):
    """Generate the outfile.

    Keyword arguments:

    parsed_results -- the returned value from `parse_results()`
    child          -- the input vcf of the child
    mother         -- the input vcf of the mother
    father         -- the input vcf of the father
    """
    child_name = child.split('_')[0]
    mother_name = mother.split('_')[0]
    father_name = father.split('_')[0]
    b = open(f'both_{mother_name}_and_{father_name}.txt', "w")
    b.write(f'CNVs in both {mother_name} and {father_name}\n\n')
    b.write(f'Number of CNVs: {parsed_results[0]}\n\n')
    for cnv in parsed_results[4]: b.write(f'{cnv}\n')
    b.close()

    m = open(f'only_{mother_name}.txt', "w")
    m.write(f'CNVs in only {mother_name}\n\n')
    m.write(f'Number of CNVs: {parsed_results[1]}\n\n')
    for cnv in parsed_results[5]: m.write(f'{cnv}\n')
    m.close()

    f = open(f'only_{father_name}.txt', "w")
    f.write(f'CNVs in only {father_name}\n\n')
    f.write(f'Number of CNVs: {parsed_results[2]}\n\n')
    for cnv in parsed_results[6]: f.write(f'{cnv}\n')
    f.close()

    n = open(f'only_{child_name}.txt', "w")
    n.write(f'CNVs in only {child_name}\n\n')
    n.write(f'Number of CNVs: {parsed_results[3]}\n\n')
    for cnv in parsed_results[7]: n.write(f'{cnv}\n')
    n.close()


def main():
    chrom_tup = (
        'chr1', 'chr2', 'chr3', 'chr4', 'chr5',
        'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
        'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
        'chr22', 'chr21', 'chrX', 'chrY', 'chrM'
    )
    args = parse_args()
    child = [args.c] * len(chrom_tup)
    mother = [args.m] * len(chrom_tup)
    father = [args.f] * len(chrom_tup)
    cpus = int(args.t)
    results = chunk_compare(chrom_tup, child, mother, father, cpus)
    parsed_results = parse_results(results)
    make_outfile(parsed_results, args.c, args.m, args.f)


if __name__ == '__main__':
    main()
