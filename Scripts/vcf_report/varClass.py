#!/usr/bin/env python3

import argparse
from json import dump
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
        required = True
    )
    parser.add_argument(
        '-c',
        metavar = '--CLASSIFY_CNV',
        type = str,
        help = '',
        required = True
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
            if second.start < first.stop + 4500:
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
    contig.append(cnv_list[start_index].start+1)
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


def cnv_bed(contigs, sample_id):
    """
    """
    file_name = f'{sample_id}_cnv.bed'
    f = open(file_name, "w")
    for contig in contigs:
        f.write(f'{contig[0]}\t{contig[3]}\t{contig[4]}\t{contig[2]}\n')
    f.close()


def call_classify_cnv(cpus, sample_id, script):
    """
    """
    infile = f'{sample_id}_cnv.bed'
    outfile = f'{sample_id}_ClassifyCNV_out'
    launch = f'{script} --infile {infile} --GenomeBuild hg19 --cores {cpus} --outdir {outfile} --precise'
    system(launch)
    clean = f'rm {infile}; rm -rf {outfile}/Intermediate_files'
    system(clean)


def get_cnv_determination(sample_id):
    """
    """
    cnvs = []
    infile = f'{sample_id}_ClassifyCNV_out/Scoresheet.txt'
    with open(infile) as f:
        for line in f:
            entry = line.split('\t')
            if entry[5] == 'Likely pathogenic' or entry[5] == 'Pathogenic':
                cnvs.append((entry[1], entry[2], entry[3], entry[5]))
    clean = f'rm -rf {sample_id}_ClassifyCNV_out'
    system(clean)
    return cnvs


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
    return (snp_list, sv_list, cnv_list, exp_list)


def filter_cnv(cnv_contigs, cnv_contig_determinations):
    """
    """
    path_cnvs = []
    for contig in cnv_contigs:
        for determination in cnv_contig_determinations:
            match = (
                determination[0] in contig[0]
                and contig[3] == int(determination[1])
                and contig[4] == int(determination[2])
            )
            if match:
                path_cnvs.append((
                    contig[0],
                    contig[1],
                    contig[2],
                    contig[3],
                    contig[4],
                    contig[5],
                    determination[3]
                ))
    return path_cnvs


def make_json(snp_list, sv_list, exp_list, path_cnvs):
    """
    """
    data = {}
    data['snp'] = []
    data['sv'] = []
    data['exp'] = []
    data['cnv'] = []
    for snp in snp_list:
        data['snp'].append({
            'chrom': snp.contig,
            'start': snp.start,
            'stop': snp.stop,
            'alt': snp.alts,
            'lof': snp.info['LOF'],
            'ann': snp.info['ANN']
        })
    for sv in sv_list:
        data['sv'].append({
            'chrom': sv.contig,
            'start': sv.start,
            'stop': sv.stop,
            'alt': sv.alts,
            'lof': sv.info['LOF']
        })
    for exp in exp_list:
        data['exp'].append({
            'chrom': exp.contig,
            'start': exp.start,
            'stop' : exp.stop,
            'gene' : exp.info['VARID'],
            'alt' : exp.alts
        })
    for cnv in path_cnvs:
        data['cnv'].append({
            'chrom': cnv[0],
            'cyto': cnv[1],
            'start': cnv[2],
            'stop': cnv[3],
            'alt': cnv[4],
            'genes': cnv[5],
            'variant class': cnv[6]
        })
    with open('data.json', 'w') as outfile:
        dump(data, outfile, indent = 2)


def main():
    """
    """
    args = parse_args()
    panel = get_panel(args.p)
    vcf = VariantFile(args.v)
    sample_id = args.s
    cpus = args.t
    script = args.c
    if cpus == None: cpus = 1
    snp_list, sv_list, cnv_list, exp_list = filter_vcf(vcf, panel)
    cnv_contigs = check_entries(cnv_list)
    cnv_bed(cnv_contigs, sample_id)
    call_classify_cnv(cpus, sample_id, script)
    cnv_contig_determinations = get_cnv_determination(sample_id)
    path_cnvs = filter_cnv(cnv_contigs, cnv_contig_determinations)
    make_json(snp_list, sv_list, exp_list, path_cnvs)


if __name__ == '__main__':
    main()
