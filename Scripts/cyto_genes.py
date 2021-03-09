#!/usr/bin/env python3

import argparse

parser = argparse.ArgumentParser(description='Append genes onto cytoband bed file')
parser.add_argument(
    '-b',
    metavar  = '--BAND',
    type     = str,
    help     = 'bed file of the cytobands from ucsc table browser',
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
    '-o',
    metavar  = '--OUT',
    type     = str,
    help     = 'name of the outputfile'
)

def get_list(item_file):
    item_list = []
    with open(item_file) as f:
        for item in f:
            item = item.strip()
            item = tuple(item.split('\t'))
            item_list.append(item)
    return item_list


def check_gene(band_chr, band_range, gene_chr, gene_start, gene_end):
    if gene_chr == band_chr:
        if gene_start in band_range and gene_end in band_range:
            return True
        else:
            return False
    else:
        return False

def make_output(out_file, band_list, gene_matches):
    f = open(out_file, "w")
    for i in range(len(band_list)):
        f.write(band_list[i][0])
        f.write('\t')
        f.write(band_list[i][1])
        f.write('\t')
        f.write(band_list[i][2])
        f.write('\t')
        f.write('\"')
        f.write(band_list[i][3])
        f.write(';')
        for gene in gene_matches[i]:
            if gene == gene_matches[i][-1]:
                f.write(gene)
            else:
                f.write(gene)
                f.write(',')
        f.write('\"')
        f.write('\n')
    f.close()


def main():
    args = parser.parse_args()
    band_list = get_list(args.b)
    gene_list = get_list(args.g)
    if args.o == None:
        out_file = "cyto_genes.bed"
    else:
        out_file = args.o
    gene_matches = [None] * len(band_list)

    for i, band in enumerate(band_list):
        band_range = range(int(band[1]), int(band[2])+1)
        band_chr = band[0]
        genes = []
        for j, gene in enumerate(gene_list):
            gene_chr = gene[0]
            gene_start = gene[1]
            gene_end = gene[2]
            if check_gene(band_chr, band_range, gene_chr, int(gene_start), int(gene_end)):
                genes.append(gene[3])
        if not genes:
            gene_matches[i] = '.'
        else: 
            gene_matches[i] = genes[:]
                
    make_output(out_file, band_list, gene_matches)

if __name__ == '__main__':
    main()
