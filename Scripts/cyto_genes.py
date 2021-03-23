#!/usr/bin/env python3

import argparse
from concurrent.futures import ProcessPoolExecutor
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
    """Parses input files.

    This function takes an input file and creates a list of tuples that
    contain the chromosome, start, stop, and cyto/ gene information of a
    given window or gene.

    Keyword arguments:

    item_file -- the input file to use.

    Return:

    item_list -- the list of tuples from the input file
    """
    item_list = []
    with open(item_file) as f:
        for item in f:
            item = item.strip()
            item = tuple(item.split('\t'))
            item_list.append(item)
    return item_list


def get_chrom_data(chrom, windows, gene_list):
    """Grab data for a specific chromosome.

    This function takes all the windows and all the genes, then grabs
    data for only the chromosome specified. After, `gene_map()` is called.
    `gene_map()` will give use the list of gene matches for the selected
    windows.

    Keyword arguments:

    chrom     -- the selected chromosome for the chunk
    windows   -- the list of all windows
    gene_list -- the list of all genes

    Return:

    gene_matches -- a list of genes that maps onto the windows
    """
    print(f'get_chrom_data: {chrom}')
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
    """Determine window/ gene overlap.

    This function takes the range of the selected window and selected gene,
    compares them to eachother to determine which is larger, then does the
    appropriate check to see if there is any overlap present at all.

    Ketword arguments:

    window -- the selected window from the list of windows
    gene   -- the selected gene from the list of genes

    Return:

    success -- boolean value that indicates overlap
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
    """Map all gene overlaps onto count windows.

    This function is looping through all the windows and genes. For each
    window we check the list of genes for an overlap. At the start of each
    new window we check to see if any of the stop positions of the genes are
    less than the starting position of the current window. If this is true
    we remove the gene from our gene list so as we loop through are windows
    we reduce the size of the gene list.

    Keyword arguments:

    windows   -- list of all the windows
    gene_list -- list of all the genes

    Return:

    gene_matches -- a list of genes, mapped to the windows
    """
    gene_matches = ['.'] * len(windows)
    for index, window in enumerate(windows):
        gene_list = clean_gene_list(window, gene_list)
        gene_overlaps = []
        for gene in gene_list:
            if check_gene(window, gene):
                gene_overlaps.append(gene)
        if not gene_overlaps:
            gene_matches[index] = '.'
        else: 
            gene_matches[index] = gene_overlaps[:]
    return gene_matches


def clean_gene_list(window, gene_list):
    """Clean up the gene list.

    This function takes the current window and the whole gene list to loop
    through and checks if the start of the current window is greater than the
    stop of any of the genes in the gene list. If this is true these genes
    will no longer map to any of the remaining windows we are checking, so
    the genes are removed from the list so our loop is smaller.

    Keyword arguments:

    window    -- the current window being checked
    gene_list -- the list of genes we are looping through

    Return:

    gene_list -- the shortened/ cleaned up gene_list
    """
    for gene in gene_list:
        if int(window[1]) > int(gene[2]):
            gene_list.remove(gene)
        else:
            return gene_list
    return gene_list


def chunk_data(chrom_list, windows_list, gene_list, cpus):
    """Split data data into chunks.

    This function uses multiprocessing to execute each chunk of data in a
    different cpu, for as many cpus as selected (Default: 1).

    Keyword arguments:

    chrom_list   -- the list of chromosomes (default for human genome)
    windows_list -- a list the length of the CHROM_LIST that contains the
                    same repeated list of count windows
    gene_list    -- a list the length of the CHROM_LIST that contains the
                    same repeated list of genes
    cpus         -- the number of cpus to use (Default: 1)

    Return:

    results -- the generator object from ProcessPoolExecutor
    """
    if cpus == None: cpus = 1
    with ProcessPoolExecutor(max_workers = cpus) as executor:
        results = executor.map(
            get_chrom_data,
            chrom_list,
            windows_list,
            gene_list
        )
    return results


def create_report(file_name, windows, results):
    """Creating the report.

    Taking all the windows and the mapped gene matches and formatting them
    so they can be used as count windows with the R library `panelcn.mops`.

    Keyword arguments:

    file_name -- the name of the output file
    windows   -- the list of all the windows
    results   -- the list of gene_matches that maps to the windows
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
                f.write(str(gene))
            else:
                f.write(f'{gene},')
        f.write('\"\n')
    f.close()


def main():
    """The main function.

    First we parse the arguments, check that the output file doesn't already
    exist, then parse both our windows and genes lists. This data is then
    passed in the `chunk_data()` function and then the 'report' is
    created.
    """
    args = parse_args()
    file_name = create_file()
    windows = [get_list(args.w)] * len(CHROM_LIST)
    gene_list = [get_list(args.g)] * len(CHROM_LIST)
    results = chunk_data(CHROM_LIST, windows, gene_list, args.t)
    create_report(file_name, windows[0], results)


if __name__ == '__main__':
    main()
