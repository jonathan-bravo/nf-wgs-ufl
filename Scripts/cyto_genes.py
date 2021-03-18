#!/usr/bin/env python3

import argparse
import concurrent.futures
from os import path

CHROM_LIST = [
    'chr1',
    'chr2',
    'chr3',
    'chr4',
    'chr5',
    'chr6',
    'chr7',
    'chr8',
    'chr9',
    'chr10',
    'chr11',
    'chr12',
    'chr13',
    'chr14',
    'chr15',
    'chr16',
    'chr17',
    'chr18',
    'chr19',
    'chr20',
    'chr21',
    'chr22',
    'chrX',
    'chrY',
    'chrM'
]


def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -w, --WINDOWS -- input file containing the count windows
    -g, --GENE    -- input bed file containing the ucsc hg19 genes
    -t, --THREADS -- the number of cpus to use
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Append genes onto cytoband bed file'
    )
    parser.add_argument(
        '-w',
        metavar  = '--WINDOWS',
        type     = str,
        help     = 'bed file of the count windows',
        required = True
    )
    parser.add_argument(
        '-g',
        metavar  = '--GENE',
        type     = str,
        help     = 'bed file of the genes from ucsc table browser',
        required = True
    )
    parser.add_argument(
        '-t',
        metavar  = '--THREADS',
        type     = int,
        help     = 'number of the cpus to use for multi-processing',
        required = False
    )
    args = parser.parse_args()
    return args


def create_file():
    """Check presence of file and create file name.

    This function checks if a report exists for a given sample and then
    quits if `True`. If the report is not present then the file name
    is returned so a new report can be generated.

    Return:

    file_name -- the name of the report to be written
    """
    file_name = 'cyto_genes.bed'
    if path.exists(file_name):
        print("The report file `" + str(file_name) + "` already exists.")
        quit()
    return file_name


def get_list(item_file):
    """
    """
    item_list = []
    with open(item_file) as f:
        for item in f:
            item = item.strip()
            item = tuple(item.split('\t'))
            item_list.append(item)
    return item_list


def get_chrom_data(chrom, windows, gene_list):
    """
    """
    windows_for_chr = []
    genes_for_chr = []
    for w in windows:
        if w[0] == chrom:
            windows_for_chr.append(w)
    for g in gene_list:
        if g[0] == chrom:
            genes_for_chr.append(g)
    gene_matches = gene_map(windows_for_chr, genes_for_chr)
    results = []
    for match in gene_matches:
        results.append(match)
    return gene_matches


def check_gene(window, gene):
    """
    """
    window_range = range(int(window[1]), int(window[2]))
    gene_range = range(int(gene[1]), int(gene[2]))
    if len(window_range) > len(gene_range):
        if int(gene[1]) in window_range or int(gene[2]) in window_range:
            return True
        return False
    elif len(gene_range) > len(window_range):
        if int(window[1]) in gene_range or int(window[2]) in gene_range:
            return True
        return False


def gene_map(windows, gene_list):
    """
    """
    gene_matches = [None] * len(windows)
    for index, window in enumerate(windows):
        gene_list = clean_gene_list(window, gene_list)
        gene_overlaps = []
        for gene in gene_list:
            if check_gene(window, gene):
                gene_overlaps.append(gene[3])
        if not gene_overlaps:
            gene_matches[index] = '.'
        else: 
            gene_matches[index] = gene_overlaps[:]
    return gene_matches


def clean_gene_list(window, gene_list):
    """
    """
    for gene in gene_list:
        if int(window[1]) > int(gene[2]):
            gene_list.remove(gene)
        else:
            return gene_list
    return gene_list


def chunk_data(chrom_list, windows_list, gene_list_list, cpus):
    """
    """
    if cpus == None: cpus = 1
    with concurrent.futures.ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(
            get_chrom_data,
            chrom_list,
            windows_list,
            gene_list_list
        )
    return results


def create_report(file_name, windows, results):
    """
    """
    genes = []
    for result in results:
        for g in result:
            genes.append(g)
    f = open(file_name, "w")
    for i, window in enumerate(windows):
        chrom = window[0]
        start = window[1]
        stop = window[2]
        cyto = window[3]
        f.write(f'{chrom}\t{start}\t{stop}\t\"{cyto};')
        for gene in genes[i]:
            if gene == genes[i][-1]:
                f.write(gene)
            else:
                f.write(f'{gene},')
        f.write('\"\n')
    f.close()


def main():
    """
    """
    args = parse_args()
    file_name = create_file()
    windows = [get_list(args.w)] * len(CHROM_LIST)
    gene_list = [get_list(args.g)] * len(CHROM_LIST)
    results = chunk_data(CHROM_LIST, windows, gene_list, args.t)
    create_report(file_name, windows[0], results)


if __name__ == '__main__':
    main()
