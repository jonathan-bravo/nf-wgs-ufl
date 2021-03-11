#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile


def parse_args():
    """This function parses the input arguments.

    Keyword Arguments:

    -c, --CNV:
        the input CNV VCF file filtered to just deletions or duplications
    
    -s, --SAMPLE_NAME:
        the desired sample name to be included in the report name

    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputs for generating multiple reports for CNV contigs'
    )
    parser.add_argument(
        '-c',
        metavar  = '--CNV',
        type     = str,
        help     = 'name of the CNV VCF file',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_NAME',
        type = str,
        help = 'the sample id/name'
    )
    args = parser.parse_args()
    return args


def get_entries(vcf):
    """This function parses the input CNV VCF file.

    This functions takes a CNV VCF file that has been filtered to only
    deletions and duplications, and stores each hit as an entry.

    The structure of an entry is as follows for any given entry at index `i`:

        * entries[i][0] = chrom
        * entries[i][1] = start
        * entries[i][2] = stop
        * entries[i][3] = cyto
        * entries[i][4] = alt
        * entries[i][5] = genes

    Keyword Arguments:

    vcf -- the input VCF file from `args.c`

    Return:

    entries -- a list of all the hits from the VCF file
    """
    entries = []
    for entry in vcf.fetch():
        entries.append((
            entry.chrom,
            entry.start,
            entry.stop,
            entry.info["CYTO"],
            entry.alts,
            entry.info["GENE"]
        ))
    return entries


def check_overlap(first, second):
    """This function checks if two CNVs overlap.

    First we check if the called CNVs are both in the same cytoband. If this
    is true we then check to see if the CNVs are of the same type, either
    a deletion or duplication. Then we check to see if the starting position
    of the second CNV is contained within the first CNV. If all the checks
    are passed then we return a boolean value of `True` or `False` followed
    by the genes from overlapping CNV.

    Keyword Arguments:

    first  -- the CNV at the end of the current contig
    second -- the CNV that we are looking for an overlap in

    Return:

    success   -- a boolean value of either True or False
    gene_list -- a list of genes from `second` or `None` if `success` is False
    """
    if second[3] == first[3]:
        if second[4] == first[4]:
            if second[1]+1 < first[2]:
                return (True, second[5])
    return (False, None)


def contig(start_index, entries):
    """Generating CNV contigs that are returned in a final report.

    We create a list that will contain all the information for the CNV
    contig and a list that will contain all the genes from all the CNVs
    in the contig. First we append the chromosome, cytoband, CNV type, and
    start of the first CNV into the contig list. Then we loop through the
    entries until we either hit the end of the list or a CNV that does not
    overlap with our 'current contig'. At the end we return the CNV
    contig, which will either be a single CNV or multiple CNVs merged into
    a single contig and the `end_index`, which is the starting index for
    the next CNV contig.

    Keyword arguments:

    start_index -- the index of the CNV that starts the contig
    entries     -- the entire list of VCF entries

    Return:

    contig    -- a list that is the CNV contig info
    end_index -- this is the index of the start of the next CNV contig
    """
    for i in range(len(entries)):
        contig = []
        genes = []
        genes.append(str(entries[start_index][5]))
        chrom = "CHROM: " + str(entries[start_index][0])
        cyto = "CYTO: " + str(entries[start_index][3])
        cnv = "CNV: " + str(entries[start_index][4])
        start = "START: " + str(entries[start_index][1]+1)
        contig.append(chrom)
        contig.append(cyto)
        contig.append(cnv)
        contig.append(start)
        for i in range(start_index, len(entries)):
            end_index = i
            if entries[i] == entries[-1]:
                stop = "STOP: " + str(entries[i][2])
                contig.append(stop)
                end_index += 1
                break
            else:
                success, gene_list = check_overlap(entries[i], entries[i+1])
                if success:
                    genes.append(gene_list)
                    continue
                else:
                    stop = "STOP: " + str(entries[i][2])
                    contig.append(stop)
                    end_index += 1
                    break
        genes[:] = [x for x in genes if x != '.']
        gene = "GENES: " + str(genes)
        contig.append(gene)
        return contig, end_index


def strip_cnv(cnv):
    """This function cleans up the CNV contig CNType.

    This function takes the CNV type, either deletion or duplication from
    each CNV contig and strips it of the unnecessary characters.

    Keyword arguments:

    cnv -- is the CNV type of the contig

    Return:

    strip_cnv -- the returned string that is the CNV type of the contig
              -- with extranious characters removed
    """
    strip_cnv = cnv.replace('\'', '')
    strip_cnv = strip_cnv.replace('(', '')
    strip_cnv = strip_cnv.replace(')', '')
    strip_cnv = strip_cnv.replace(',', '')
    return strip_cnv


def strip_gene(gene):
    """This function cleans up the CNV contig gene list.

    This function takes the gene list from each CNV contig and strips it of
    the unnecessary characters.

    Keyword arguments:

    gene -- is the list of genes for the CNV contig

    Return:

    strip_gene -- the returned string that is the list of genes of the contig
               -- with extranious characters removed
    """
    strip_gene = gene.replace('\"', '')
    strip_gene = strip_gene.replace('\'', '')
    strip_gene = strip_gene.replace('(', '')
    strip_gene = strip_gene.replace(')', '')
    strip_gene = strip_gene.replace('[', '')
    strip_gene = strip_gene.replace(']', '')
    return strip_gene


def generate_report(file_name, cnvs):
    """Generates the final reports.

    This function takes the CNV contigs and prints them out into four different
    reports. There is a report for all duplications, all deletions, only
    duplications that contain genes, and only deletions that contain genes.

    Keyword arguments:

    cnvs -- the cnv contigs generated from `contig()`

    Output:

    {sample_id}_cnv_contigs_report.tsv
    """
    f = open(file_name, "w")
    for i in range(len(cnvs)):
        chrom = cnvs[i][0]
        cyto  = cnvs[i][1]
        cnv   = strip_cnv(cnvs[i][2])
        start = cnvs[i][3]
        stop  = cnvs[i][4]
        genes = strip_gene(cnvs[i][5])
        entry = f'{chrom}\t{cyto}\t{cnv}\t{start}\t{stop}\t{genes}\n'
        f.write(entry)
    f.close()


def main():
    """ Create CNV contigs from sliding windows CNV calls.

    This script looks at 25kb sliding window CNV calls from a VCF file and
    merges the CNV windows that are on the same cytoband, are the same
    CNV type (DEL, DUP), and are contiguous. The resulting report will
    allow better comparisons with CGH Array data and tools already offered
    from companies such as NxClinical and XXX.

    Output:

    file_name -- the name of the resulting report 
    """
    args = parse_args()
    vcf = VariantFile(args.c)
    sample_id = args.s

    file_name = f'{sample_id}_cnv_contigs_report.tsv'

    if path.exists(file_name):
        print("The report file `" + str(file_name) + "` already exists.")
        quit()

    cnvs = []
    entries = get_entries(vcf)
    index = 0

    while index < len(entries):
        cnv, index = contig(index, entries)
        cnvs.append(cnv)

    generate_report(file_name, cnvs)


if __name__ == '__main__':
    main()
