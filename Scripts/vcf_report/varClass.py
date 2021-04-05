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
        required = False
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


def map_score(x, low, high, a, b):
    """
    """
    return((((x - low) / (high - low)) * (a - b)) + b)


def decision_tree(info, result_list):
    """
    """
    score = 0.0
    test = [x for x in info if x != '.']
    value = 1 / len(test)
    if len(result_list) == 4:
        for result in test:
            if result_list[0] in result: score += value
            elif result_list[1] in result: score += value * 0.5
            elif result_list[2] in result: score -= value * 0.5
            elif result_list[3] in result: score -= value
    elif len(result_list) == 3:
        for result in test:
            if result_list[0] in result: score += value
            elif result_list[1] in result: score += value * 0.75
            elif result_list[2] in result: score -= value * 0.5
    elif len(result_list) == 2:
        for result in test:
            if result_list[0] in result: score += value
            elif result_list[1] in result: score -= value
    return score


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
            if entry[5] == 'Pathogenic':
                cnvs.append((entry[1], entry[2], entry[3], entry[6]))
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
        if panel != None:
            snp = (
                'SNVHPOL' in variant.info.keys()
                and len(variant.info.keys()) > 3
                and str(variant.info['ANN']).split('|')[3] in panel
            )
            sv = (
                'CIGAR' in variant.info.keys()
                and str(variant.info['ANN']).split('|')[3] in panel
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
        else:
            snp = ('SNVHPOL' in variant.info.keys() and len(variant.info.keys()) > 3)
            sv = ('CIGAR' in variant.info.keys())
            cnv = ('SVTYPE' in variant.info.keys())
            cnv_multi_gene = ('SVTYPE' in variant.info.keys())
            exp = ('VARID' in variant.info.keys())
        if 'PASS' in variant.filter.keys():
            if snp:
                snp_score = 0.0
                ann = str(variant.info['ANN']).split('|')[2]
                if 'HIGH' in ann: snp_score += 1.0
                elif 'MODERATE' in ann: snp_score += 0.5
                if 'LOF' in variant.info.keys(): snp_score += float(
                    variant.info['LOF'][0].split('|')[3].strip(')')
                )
                if 'dbNSFP_CADD_phred' in variant.info.keys():
                    cadd = variant.info['dbNSFP_CADD_phred'][0]
                    if cadd >= 15: snp_score += map_score(
                        10**(cadd / -10),
                        10**(15 / -10),
                        10**(99 / -10),
                        1.0,
                        0.01
                    )
                if 'dbNSFP_phastCons100way_vertebrate' in variant.info.keys():
                    info = variant.info['dbNSFP_phastCons100way_vertebrate'][0]
                    snp_score += map_score(info, 0.0, 1.0, 1.0, -1.0)
                if 'dbNSFP_MutationTaster_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_MutationTaster_pred']
                    snp_score += decision_tree(info, ['A', 'D', 'N', 'P'])
                if 'dbNSFP_MutationAssessor_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_MutationAssessor_pred']
                    snp_score += decision_tree(info, ['H', 'M', 'L', 'N'])
                if 'dbNSFP_LRT_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_LRT_pred']
                    snp_score += decision_tree(info, ['D', 'N'])
                if 'dbNSFP_FATHMM_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_FATHMM_pred']
                    snp_score += decision_tree(info, ['D', 'T'])
                if 'dbNSFP_MetaSVM_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_MetaSVM_pred']
                    snp_score += decision_tree(info, ['D', 'T'])
                if 'dbNSFP_PROVEAN_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_PROVEAN_pred']
                    snp_score += decision_tree(info, ['D', 'N'])
                if 'dbNSFP_Polyphen2_HVAR_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_Polyphen2_HVAR_pred']
                    snp_score += decision_tree(info, ['D', 'P', 'B'])
                if 'dbNSFP_SIFT_pred' in variant.info.keys():
                    info = variant.info['dbNSFP_SIFT_pred']
                    snp_score += decision_tree(info, ['D', 'T'])
                if 'dbNSFP_GERP___RS' in variant.info.keys():
                    gerp = variant.info['dbNSFP_GERP___RS'][0]
                    if gerp >= 0: snp_score += map_score(gerp, 0.0, 6.18, 1.0, 0.01)
                    elif 0 > gerp: snp_score += map_score(gerp, -10.7, 0.0, 0.0, -1.0)
                n_score = map_score(snp_score, -9.5, 13.0, 1.0, 0.0)
                if n_score >= 0.87: snp_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.alts,
                    variant.info['ANN'],
                    n_score
                ))
            elif sv:
                sv_score = 0.0
                ann = str(variant.info['ANN']).split('|')[2]
                if 'HIGH' in ann: sv_score += 1.0
                elif 'MODERATE' in ann: sv_score += 0.5
                if 'LOF' in variant.info.keys():
                    sv_score += float(variant.info['LOF'][0].split('|')[3].strip(')'))
                if sv_score >= 1.0: sv_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.alts,
                    variant.info['ANN'],
                    sv_score
                ))
            elif cnv or cnv_multi_gene:
                cnv_list.append(variant)
            elif exp:
                alt_length = []
                for alt in variant.alts:
                    length = int(alt.strip('<STR').strip('>'))
                    alt_length.append(length)
                if any(x >= 70 for x in alt_length): exp_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.alts,
                    variant.info['VARID'],
                    alt_length
                ))
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


def get_literature(panel, genes):
    """
    """
    lit = []
    with open(panel) as f:
        for line in f:
            entry = line.split('\t')
            if entry[0] in genes: 
                for pub in entry[12].split(';'):
                    if pub != '': lit.append(pub)
    lit = set(lit)
    lit_list = []
    with open('PMC-ids.csv') as f:
        for line in f:
            pub = line.split(',') # 7 = DOI, 9 = PMID
            if pub[9] in lit: lit_list.append(f'{pub[8]}: {pub[7]}')
    return lit_list


def make_json(panel, gene_panel, snp_list, sv_list, exp_list, path_cnvs, sample_id):
    """
    """
    genes = []
    data = {}
    data['snp'] = {
        'pathogenic_snp':[],
        'likely_pathogenic_snp':[]
    }
    data['sv'] = {
        'pathogenic_sv':[],
        'likely_pathogenic_sv':[]
    }
    data['cnv'] = {
        'pathogenic_cnv':[],
        'likely_pathogenic_cnv':[]
    }
    data['exp'] = {
        'pathogenic_exp':[],
        'likely_pathogenic_exp':[]
    }
    data['metadata'] = {
        'genes_in_panel': gene_panel,
        'supporting_literature':[],
        'pipeline':[
            {
                'name': 'TrimmomaticPE',
                'version': '0.39',
                'purpose': 'trim reads',
                'citation': '10.1093/bioinformatics/btu170'
            },
            {
                'name': 'BWA MEM',
                'version': '0.7.17-r1188',
                'purpose': 'paired read alignment',
                'citation': ' 10.1093/bioinformatics/btp324'
            },
            {
                'name': 'Samtools',
                'version': '1.10',
                'purpose': 'SAM to BAM, sorting BAM, and indexing BAM',
                'citation': '10.1093/gigascience/giab008'
            },
            {
                'name': 'Strelka2',
                'version': '2.9.10',
                'purpose': 'detection of SNPs and SVs',
                'citation': '10.1038/s41592-018-0051-x'
            },
            {
                'name': 'Panelcn.MOPS',
                'version': '1.8.0',
                'purpose': 'detection of CNVs within specified windows',
                'citation': '10.1002/humu/23237'
            },
            {
                'name': 'ExpansionHunter',
                'version': '4.0.2',
                'purpose': 'detection of repeat expansions',
                'citation': '10.1186/s13059-020-02017-z'
            },
            {
                'name': 'snpEff',
                'version': '5.0c',
                'purpose': 'snp annotation',
                'citation': '10.4161/fly.19695'
            },
            {
                'name': 'Bcftools',
                'version': '1.10.2',
                'purpose': 'merging and indexing vcf files',
                'citation': '10.1093/gigascience/giab008'
            },
            {
                'name': 'Picard',
                'version': '2.23.8',
                'purpose': 'generating bam metrics',
                'citation': 'http://broadinstitute.github.io/picard/',
            },
            {
                'name': 'FastQC',
                'version': '0.11.9',
                'purpose': 'generating fastq metrics',
                'citation': 'https://www.bioinformatics.babraham.ac.uk/projects/fastqc/'
            },
            {
                'name': 'MultiQC',
                'version': '1.9',
                'purpose': 'creating merged qc page',
                'citation': '10.1093/bioinformatics/btw354'
            }
        ]
    }
    for snp in snp_list:
        genes.append(str(snp[4]).split('|')[3])
        snp_dict = {
            'chrom': snp[0],
            'start': snp[1],
            'stop': snp[2],
            'alt': snp[3],
            'ann': snp[4]
        }
        if float(snp[5]) >= 0.9: data['snp']['pathogenic_snp'].append(snp_dict)
        else: data['snp']['likely_pathogenic_snp'].append(snp_dict)
    for sv in sv_list:
        genes.append(str(sv[4]).split('|')[3])
        sv_dict = {
            'chrom': sv[0],
            'start': sv[1],
            'stop': sv[2],
            'alt': sv[3],
            'ann': sv[4]
        }
        if float(sv[5]) >= 1.5: data['sv']['pathogenic_sv'].append(sv_dict)
        else: data['sv']['likely_pathogenic_sv'].append(sv_dict)
    for cnv in path_cnvs:
        genes.append(cnv[5])
        cnv_dict = {
            'chrom': cnv[0],
            'cyto': cnv[1],
            'start': cnv[3],
            'stop': cnv[4],
            'alt': cnv[2],
            'genes': cnv[5]
        }
        if float(cnv[6]) > 1.0: data['cnv']['pathogenic_cnv'].append(cnv_dict)
        else: data['cnv']['likely_pathogenic_cnv'].append(cnv_dict)
    for exp in exp_list:
        genes.append(exp[4])
        exp_dict = {
            'chrom': exp[0],
            'start': exp[1],
            'stop' : exp[2],
            'alt' : exp[3],
            'gene': exp[4]
        }
        if any(x >= 115 for x in exp[5]): data['exp']['pathogenic_exp'].append(exp_dict)
        else: data['exp']['likely_pathogenic_exp'].append(exp_dict)
    if panel != None: data['metadata']['supporting_literature'] = get_literature(panel, genes)
    else: data['metadata']['supporting_literature'] =  None
    with open(f'{sample_id}_report.json', 'w') as outfile:
        dump(data, outfile, indent = 4)


def main():
    """
    """
    args = parse_args()
    if args.p != None: panel = get_panel(args.p)
    else: panel = None
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
    make_json(args.p, panel, snp_list, sv_list, exp_list, path_cnvs, sample_id)


if __name__ == '__main__':
    main()
