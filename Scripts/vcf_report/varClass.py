#!/usr/bin/env python3

import argparse
from os import system
from pysam import VariantFile
from pprint import pprint

def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -b, --VCF       -- input sorted merged vcf file
    -s, --SAMPLE_ID -- sample id

    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'generation of json for variant classification for clinical reporting'
    )
    parser.add_argument(
        '-v',
        metavar = '--VCF',
        type = str,
        help = 'The input sample vcf',
        required = True
    )
    parser.add_argument(
        '-p',
        metavar = '--GENE_PANEL',
        type = str,
        help = '',
        required = True
    )
    parser.add_argument(
        '-t',
        metavar = '--THREADS',
        type = int,
        help = '',
        required = False
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'sample id to be included in report name',
        required = False
    )
    args = parser.parse_args()
    return args


def get_panel(panel):
    """
    """
    panel_genes = []
    with open(panel) as f:
        for line in f:
            panel_genes.append(line.split('\t')[0])
    panel_genes.pop(0)
    return panel_genes


def check_overlap(first, second):
    """
    """
    if second.info['CYTO'] == first.info['CYTO']:
        if second.alts[0] == first.alts[0]:
            if second.start < first.stop:
                return (True, second.info['GENE'])
    return (False, None)


def contig(start_index, cnv_list):
    """
    """
    contig = []
    genes = []
    genes.append(cnv_list[start_index].info['GENE'])
    contig.append(cnv_list[start_index].contig)
    contig.append(cnv_list[start_index].info['CYTO'])
    contig.append(cnv_list[start_index].alts[0])
    contig.append(cnv_list[start_index].start)
    for i in range(start_index, len(cnv_list)):
        end_index = i
        if cnv_list[i] == cnv_list[-1]:
            contig.append(cnv_list[i].stop)
            end_index += 1
            break
        else:
            success, gene_list = check_overlap(cnv_list[i], cnv_list[i+1])
            if success:
                genes.append(gene_list)
                continue
            else:
                contig.append(cnv_list[i].stop)
                end_index += 1
                break
    genes_tmp = []
    for gene in genes:
        if isinstance(gene, tuple):
            for x in gene:
                genes_tmp.append(x)
        else:
            genes_tmp.append(gene)
    genes_final = []
    [genes_final.append(x) for x in genes_tmp if x not in genes_final]
    contig.append(genes_final)
    return contig, end_index


def check_entries(cnv_list):
    """
    """
    cnvs = []
    index = 0
    while index < len(cnv_list):
        cnv, index = contig(index, cnv_list)
        cnvs.append(cnv)
    return cnvs


def cnv_bed():
    """
    """
    
def filter_vcf(vcf, panel):
    """This function parses the input VCF file.

    This functions takes a CNV VCF file that has been filtered to only
    deletions and duplications, and stores each hit as an entry.

    The structure of an entry is as follows for any given entry at index `i`:

    Keyword Arguments:

    vcf -- the input VCF file from `args.v`

    Return:

    entries -- a list of all the hits from the VCF file
    """
    snp_list = []
    sv_list = []
    cnv_list = []
    exp_list = []
    for variant in vcf:
        # CADD_phred >= 30 is %0.1 of population and highly causative
        # phastCons >= 0.999 is extremely likely to be deleterious (%0.1 of population and highly causative)
        # MutationTaster == A is disease causing automatic - i.e. known to be deleterious
        # LRT_pred == D is a predicition of deleteriousness
        snp = (
            'SNVHPOL' in variant.info.keys()
            and len(variant.info.keys()) > 3
            and str(variant.info['ANN']).split('|')[3] in panel
            and 'dbNSFP_CADD_phred' in variant.info.keys()
            and variant.info['dbNSFP_CADD_phred'][0] > 30
            and 'LOF' in variant.info.keys()
            and 'dbNSFP_phastCons100way_vertebrate' in variant.info.keys()
            and variant.info['dbNSFP_phastCons100way_vertebrate'][0] >= 0.999
            and 'dbNSFP_MutationTaster_pred' in variant.info.keys()
            and 'A' in variant.info['dbNSFP_MutationTaster_pred'][0]
            and 'dbNSFP_LRT_pred' in variant.info.keys()
            and 'D' in variant.info['dbNSFP_LRT_pred'][0]
        )
        sv = (
            'CIGAR' in variant.info.keys()
            and str(variant.info['ANN']).split('|')[3] in panel
            and 'LOF' in variant.info.keys()
        )
        cnv = (
            'SVTYPE' in variant.info.keys()
            and variant.info['GENE'] in panel
        )
        cnv_multi_gene = (
            'SVTYPE' in variant.info.keys()
            and isinstance(variant.info['GENE'], tuple)
            and any(item in variant.info['GENE'] for item in panel)
        )
        exp = (
            'VARID' in variant.info.keys()
            and variant.info['VARID'] in panel
        )
        if 'PASS' in variant.filter.keys():
            if snp:
                snp_list.append(variant)
            elif sv:
                sv_list.append(variant)
            elif cnv or cnv_multi_gene:
                cnv_list.append(variant)
            elif exp:
                exp_list.append(variant)

    cnv_contigs = check_entries(cnv_list)
    print("number of CNV Contigs: ", len(cnv_contigs))
    print("number of CNVs: ", len(cnv_list))
    print("number of SNPs: ", len(snp_list))
    print("number of SVs: ", len(sv_list))
    print("number of EXPs: ", len(exp_list))


def main():
    """
    """
    args = parse_args()
    panel = get_panel(args.p)
    vcf = VariantFile(args.v)
    filter_vcf(vcf, panel)


if __name__ == '__main__':
    main()
