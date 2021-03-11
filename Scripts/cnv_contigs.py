#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile

def parse_args():
    """
    This function parses the input arguments.

    Keyword Arguments:

    -c, --CNV:
        the input CNV VCF filtered to just deletions or duplications
    
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
        help     = 'name of the CNV vcf',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_NAME',
        type = str,
        help = 'the sample name'
    )
    args = parser.parse_args()
    return args

def get_entries(vcf):
    """
    This function parses the input CNV file.

    This functions takes a CNV variant call file that has been filtered to only
    deletions and duplications, and stores each hit as an entry.

    The structure of an entry is as follows for any given entry at index `i`:

        * entries[i][0] = chrom
        * entries[i][1] = start
        * entries[i][2] = stop
        * entries[i][3] = cyto
        * entries[i][4] = alt
        * entries[i][5] = genes

    Keyword Arguments:

    vcf -- the input VCF from `args.c`

    Return:

    entries -- a list of all the hits from the VCF

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

# This function will be for checking if the second entry is contiguous with
# the first entry
# ex: first = entries[0], second = entries[1]
def check_overlap(first, second):
    """
    This function checks if two CNVs overlap.

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
    gene_list -- a list of genes from `second` or None if `success` is False

    """
    if second[3] == first[3]:
        if second[4] == first[4]:
            if second[1]+1 < first[2]:
                return (True, second[5])
    return (False, None)
            


def contig(start_index, entries):
    for i in range(len(entries)):
        # This is the list that will hold the info for the CNV contig
        contig = []
        # This is a list that will contain all genes for the CNV contig
        genes = []
        genes.append(str(entries[start_index][5]))
        # FIRST we need to put the information from the first entry into a list
        chrom = "CHROM: " + str(entries[start_index][0])
        cyto = "CYTO: " + str(entries[start_index][3])
        cnv = "CNV: " + str(entries[start_index][4])
        start = "START: " + str(entries[start_index][1]+1)
        contig.append(chrom)
        contig.append(cyto)
        contig.append(cnv)
        contig.append(start)
        # we loop, starting with a specified index, so we can iterate through
        # the whole list of entries without duplication.
        for i in range(start_index, len(entries)):
            end_index = i
            # If we are at the last entry we are done matching 
            # (there is nothing left in the vcf file)
            if entries[i] == entries[-1]:
                stop = "STOP: " + str(entries[i][2])
                contig.append(stop)
                end_index += 1
                break
            else:
                # We get a boolean value of True or False and a tuple of
                # genes from the entry that belongs in our contig
                success, gene_list = check_overlap(entries[i], entries[i+1])
                if success:
                    genes.append(gene_list)
                    # Since we were successful we want to continue our for loop
                    continue
                else:
                    # Stop is determined from the last contig that overlaps
                    stop = "STOP: " + str(entries[i][2])
                    contig.append(stop)
                    end_index += 1
                    # Then we break from the for loop because the next entry
                    # does not overlap or is not contiguous
                    break
        # We then set the value of the genes equal to all the genes from the
        # contigs, but first we will want to strip any '.' (NA) values from our
        # list
        genes[:] = [x for x in genes if x != '.']
        gene = "GENES: " + str(genes)
        contig.append(gene)
        # we are returning the entire contig with:
        # CHROM, CYTO, CNV (DEL/DUP), START, STOP, GENES
        # We also need to return the index of the last entry so we can start
        # at the next entry
        return contig, end_index



def generate_report(sample_id, cnvs):
    """
    Generates the final reports.

    This function takes the CNV contigs and prints them out into four different
    reports. There is a report for all duplications, all deletions, only
    duplications that contain genes, and only deletions that contain genes.

    Keyword arguments:

    cnvs -- the cnv contigs generated from `contig()`

    Output:

    {sample_id}_cnv_contigs_report.tsv
    
    """
    file_name = f'{sample_id}_cnv_contigs_report.tsv'
    f = open(file_name, "w")

    for i in range(len(cnvs)):
        chrom = cnvs[i][0]
        cyto  = cnvs[i][1]
        cnv   = cnvs[i][2].replace('\'', '').replace('(', '').replace(')', '').replace(',', '')
        start = cnvs[i][3]
        stop  = cnvs[i][4]
        genes = cnvs[i][5].replace('\"', '').replace('\'', '').replace('(', '').replace(')', '').replace('[', '').replace(']', '')
        
        entry = f'{chrom}\t{cyto}\t{cnv}\t{start}\t{stop}\t{genes}\n'
        
        f.write(entry)
    
    f.close()



def main():
    args = parse_args()
    vcf = VariantFile(args.c)
    sample_id = VariantFile(args.s)

    cnvs = []
    entries = get_entries(vcf)
    index = 0

    while index < len(entries):
        cnv, index = contig(index, entries)
        cnvs.append(cnv)

    generate_report(sample_id, cnvs)

        

if __name__ == '__main__':
    main()
