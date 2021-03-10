#!/usr/bin/env python3

import argparse
from os import path
from pysam import VariantFile

##  ~/Documents/nf-wgs-ufl/Scripts/contiguous.py -c NQ-21-01-BC708501-385_S25_filtered_cnv.vcf.gz

parser = argparse.ArgumentParser(description='Compare overlaps with a CNV vcf file and a standard SNV vcf file')
parser.add_argument(
    '-c',
    metavar  = '--CNV',
    type     = str,
    help     = 'name of the CNV vcf',
    required = True
)
parser.add_argument(
    '-o',
    metavar  = '--OUT',
    type     = str,
    help     = 'name of the outputfile'
)



# This function will pull all entries from the vcf file into a list and
# return the list
def get_entries(vcf):
    # a list that will contain a tuple for every contig from the filtered
    # cnv vcf file
    entries = []
    # a for loop that will append all the vcf entries from the filtered cnv vcf
    # as a tuple into the contig list
    for entry in vcf.fetch():
        entries.append((
            entry.chrom,
            entry.start,
            entry.stop,
            entry.info["CYTO"],
            entry.alts,
            entry.info["GENE"]
        ))
    # Structure of entry:
    # entries[i][0] = chrom
    # entries[i][1] = start
    # entries[i][2] = stop
    # entries[i][3] = cyto
    # entries[i][4] = alt
    # entries[i][5] = genes
    # Returning the list of entries from the vcf file
    return entries



# This function will be for checking if the second entry is contiguous with
# the first entry
# ex: first = entries[0], second = entries[1]
def check_overlap(first, second):
    # check that the cytoband matched
    if second[3] == first[3]:
        # check that the alt matchs (DEL or DUP)
        if second[4] == first[4]:
            # check that the second entry overlaps with the first
            # second START is before the first END
            # If everything passes we teturn a tuple that containes our answer
            # and the genes from the second entry to append to our genes list
            if second[1]+1 < first[2]:
                return (True, second[5])
    # If the second entry is not contiguous return false. The return is a
    # tuple of False, and None to maintain the correct formatting throughout
    # the function 
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



def generate_report(cnvs):
    """Generates the final reports

    This function takes the CNV contigs and prints them out into four different
    reports. There is a report for all duplications, all deletions, only
    duplications that contain genes, and only deletions that contain genes.
    """
        
    f = open("report.txt", "w")
    for i in range(len(cnvs)):
        entry = '{chrom}\t{cyto}\t{cnv}\t{start}\t{stop}\t{genes}\n'.format(
            chrom = cnvs[i][0],
            cyto  = cnvs[i][1],
            cnv   = cnvs[i][2]
                .replace('\'', '')
                .replace('(', '')
                .replace(')', '')
                .replace(',', ''),
            start = cnvs[i][3],
            stop  = cnvs[i][4],
            genes = cnvs[i][5]
                .replace('\"', '')
                .replace('\'', '')
                .replace('(', '')
                .replace(')', '')
                .replace('[', '')
                .replace(']', '')
        )
        f.write(entry)



def main():
    args = parser.parse_args()
    vcf = VariantFile(args.c)

    cnvs = []
    entries = get_entries(vcf)
    index = 0

    while index < len(entries):
        cnv, index = contig(index, entries)
        cnvs.append(cnv)

    generate_report(cnvs)

        

if __name__ == '__main__':
    main()
