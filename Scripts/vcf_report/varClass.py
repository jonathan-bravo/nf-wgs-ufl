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
    parser.add_argument(
        '-l',
        metavar = '--PMC_IDS',
        type = str,
        help = '',
        required = True
    )
    parser.add_argument(
        '--low_coverage',
        action = argparse.BooleanOptionalAction,
        default = False
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
    rcn = 0
    mrcn = 0
    genes.append(cnv_list[start_index].info['GENE'])
    contig.append(cnv_list[start_index].contig)
    contig.append(cnv_list[start_index].info['CYTO'])
    contig.append(cnv_list[start_index].alts[0])
    contig.append(cnv_list[start_index].start+1)
    rcn += cnv_list[start_index].samples['SAMPLE1'].get('RCN')
    mrcn += cnv_list[start_index].samples['SAMPLE1'].get('MRCN')
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
                rcn += cnv_list[start_index].samples['SAMPLE1'].get('RCN')
                mrcn += cnv_list[start_index].samples['SAMPLE1'].get('MRCN')
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
    contig.append(rcn)
    contig.append(mrcn)
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


def filter_vcf(vcf, panel, low_coverage):
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
                'CNCLASS' in variant.info.keys()
                and variant.info['GENE'] in panel
            )
            cnv_multi_gene = (
                'CNCLASS' in variant.info.keys()
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
            cnv = ('CNCLASS' in variant.info.keys())
            cnv_multi_gene = ('SVTYPE' in variant.info.keys())
            exp = ('VARID' in variant.info.keys())
        if 'PASS' in variant.filter.keys():
            if snp and not low_coverage:
                snp_score = 0.0
                ann = str(variant.info['ANN']).split('|')[2]
                if 'HIGH' in ann: snp_score += 1.0
                elif 'MODERATE' in ann: snp_score += 0.5
                if 'LOF' in variant.info.keys(): 
                    snp_score += float(
                        variant.info['LOF'][0].split('|')[3].strip(')')
                    )
                    lof = variant.info['LOF'][0]
                else: lof = '.'
                if 'NMD' in variant.info.keys():
                    nmd = variant.info['NMD'][0]
                else: nmd = '.'
                if 'dbNSFP_CADD_phred' in variant.info.keys():
                    cadd = variant.info['dbNSFP_CADD_phred'][0]
                    if cadd >= 15: snp_score += map_score(
                        10**(cadd / -10),
                        10**(15 / -10),
                        10**(99 / -10),
                        1.0,
                        0.01
                    )
                else: cadd = '.'
                if 'dbNSFP_phastCons100way_vertebrate' in variant.info.keys():
                    phastcons = variant.info['dbNSFP_phastCons100way_vertebrate'][0]
                    snp_score += map_score(phastcons, 0.0, 1.0, 1.0, -1.0)
                if 'dbNSFP_MutationTaster_pred' in variant.info.keys():
                    mutaste = variant.info['dbNSFP_MutationTaster_pred']
                    snp_score += decision_tree(mutaste, ['A', 'D', 'N', 'P'])
                if 'dbNSFP_MutationAssessor_pred' in variant.info.keys():
                    muassess = variant.info['dbNSFP_MutationAssessor_pred']
                    snp_score += decision_tree(muassess, ['H', 'M', 'L', 'N'])
                if 'dbNSFP_LRT_pred' in variant.info.keys():
                    lrt = variant.info['dbNSFP_LRT_pred']
                    snp_score += decision_tree(lrt, ['D', 'N'])
                if 'dbNSFP_FATHMM_pred' in variant.info.keys():
                    fathmm = variant.info['dbNSFP_FATHMM_pred']
                    snp_score += decision_tree(fathmm, ['D', 'T'])
                if 'dbNSFP_MetaSVM_pred' in variant.info.keys():
                    metasvm = variant.info['dbNSFP_MetaSVM_pred']
                    snp_score += decision_tree(metasvm, ['D', 'T'])
                if 'dbNSFP_PROVEAN_pred' in variant.info.keys():
                    provean = variant.info['dbNSFP_PROVEAN_pred']
                    snp_score += decision_tree(provean, ['D', 'N'])
                if 'dbNSFP_Polyphen2_HVAR_pred' in variant.info.keys():
                    polyphen = variant.info['dbNSFP_Polyphen2_HVAR_pred']
                    snp_score += decision_tree(polyphen, ['D', 'P', 'B'])
                if 'dbNSFP_SIFT_pred' in variant.info.keys():
                    sift = variant.info['dbNSFP_SIFT_pred']
                    snp_score += decision_tree(sift, ['D', 'T'])
                if 'dbNSFP_GERP___RS' in variant.info.keys():
                    gerp = variant.info['dbNSFP_GERP___RS'][0]
                    if gerp >= 0:
                        snp_score += map_score(gerp, 0.0, 6.18, 1.0, 0.01)
                    elif 0 > gerp:
                        snp_score += map_score(gerp, -10.7, 0.0, 0.0, -1.0)
                if 'dbNSFP_ExAC_AF' in variant.info.keys():
                    gnomad = variant.info['dbNSFP_ExAC_AF'][0]
                    snp_score += map_score(
                        10**(gnomad / -10),
                        10**(1.0 / -10),
                        10**(0.0 / -10),
                        1.0,
                        0.0
                    )
                else: 
                    snp_score += 1.0
                    gnomad = '.'
                if 'dbNSFP_1000Gp3_AF' in variant.info.keys():
                    gp3 = variant.info['dbNSFP_1000Gp3_AF'][0]
                    snp_score += map_score(
                        10**(gp3 / -10),
                        10**(1.0 / -10),
                        10**(0.0 / -10),
                        1.0,
                        0.0
                    )
                else:
                    snp_score += 1.0
                    gp3 = '.'
                gt = variant.samples['SAMPLE1'].get('GT')
                if gt != None and len(gt) > 1:
                    if gt[0] == gt[1]: genotype = 'heterozygous'
                    elif gt[0] != gt[1]: genotype = 'homozygous'
                else: genotype == 'None'
                adf = variant.samples['SAMPLE1'].get('ADF')
                adr = variant.samples['SAMPLE1'].get('ADR')
                fr = f'{adf[0]}:{adr[0]}, {adf[1]}:{adr[1]}'
                if isinstance(gnomad, float): gnomad = f'{gnomad:.3e}'
                if isinstance(gp3, float): gp3 = f'{gp3:.3e}'
                n_score = map_score(snp_score, -9.5, 15.0, 1.0, 0.0)
                if n_score >= 0.87: snp_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.ref,
                    variant.alts,
                    variant.info['ANN'],
                    lof,
                    nmd,
                    genotype,
                    fr,
                    round(cadd, 3),
                    gnomad,
                    gp3,
                    n_score
                ))
            elif sv and not low_coverage:
                sv_score = 0.0
                ann = str(variant.info['ANN']).split('|')[2]
                if 'HIGH' in ann: sv_score += 1.0
                elif 'MODERATE' in ann: sv_score += 0.5
                if 'LOF' in variant.info.keys():
                    sv_score += float(variant.info['LOF'][0].split('|')[3].strip(')'))
                    lof = variant.info['LOF'][0]
                else: lof = '.'
                if 'NMD' in variant.info.keys():
                    nmd = variant.info['NMD'][0]
                else: nmd = '.'
                gt = variant.samples['SAMPLE1'].get('GT')
                if gt != None and len(gt) > 1:
                    if gt[0] == gt[1]: genotype = 'heterozygous'
                    elif gt[0] != gt[1]: genotype = 'homozygous'
                else: genotype == 'None'
                adf = variant.samples['SAMPLE1'].get('ADF')
                adr = variant.samples['SAMPLE1'].get('ADR')
                fr = f'{adf[0]}:{adr[0]}, {adf[1]}:{adr[1]}'
                if sv_score >= 1.0: sv_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.ref,
                    variant.alts,
                    variant.info['ANN'],
                    lof,
                    nmd,
                    genotype,
                    fr,
                    sv_score
                ))
            elif cnv or cnv_multi_gene:
                cnv_list.append(variant)
            elif exp:
                # alt_length = []
                # for alt in variant.alts:
                #     length = int(alt.strip('<STR').strip('>'))
                #     alt_length.append(length)
                # if any(x >= 70 for x in alt_length):
                gt = variant.samples['SAMPLE1'].get('GT')
                if gt != None and len(gt) > 1:
                    if gt[0] == gt[1]: genotype = 'heterozygous'
                    elif gt[0] != gt[1]: genotype = 'homozygous'
                else: genotype == 'None'
                repcn = variant.samples['SAMPLE1'].get('REPCN').split('/')
                exp_list.append((
                    variant.contig,
                    variant.start,
                    variant.stop,
                    variant.alts,
                    variant.info['VARID'],
                    genotype,
                    repcn
                    #alt_length
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
                    contig[6],
                    contig[7],
                    determination[3]
                ))
    return path_cnvs


def get_literature(panel, genes, pmc):
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
    with open(pmc) as f:
        for line in f:
            pub = line.split(',') # 7 = DOI, 9 = PMID
            if pub[9] in lit: lit_list.append(f'{pub[8]}: {pub[7]}')
    return lit_list


def check_interactions(path_cnvs, snp_list, sv_list, exp_list):
    """
    """
    overlaps = []
    for cnv in path_cnvs:
        cnv_dict = {
            'chrom': cnv[0],
            'cyto': cnv[1],
            'start': cnv[3],
            'stop': cnv[4],
            'alt': cnv[2],
            'genes': cnv[5],
            'RCN:MRCN': f'{cnv[6]}:{cnv[7]}'
        }
        overlap = {
            'cnv': cnv_dict,
            'snps': [],
            'svs': [],
            'exps': []
        }
        cnv_range = range(cnv[3], cnv[4]+1)
        for snp in snp_list:
            snp_overlap = (snp[0] == cnv[0] and snp[1] in cnv_range)
            if snp_overlap:
                snp_dict = snp_dict = {
                    'chrom': snp[0],
                    'start': snp[1],
                    'stop': snp[2],
                    'ref': snp[3],
                    'alt': snp[4],
                    'annotation': snp[5],
                    'lof': snp[6],
                    'nmd': snp[7],
                    'genotype': snp[8],
                    'F:R_ref_F:R_alt': snp[9],
                    'cadd': snp[10],
                    'gnomad_freq': snp[11],
                    '1000_genomes_freq': snp[12]
                }
                overlap['snps'].append(snp_dict)
        for sv in sv_list:
            sv_overlap = (sv[0] == cnv[0] and sv[1] in cnv_range)
            if sv_overlap:
                sv_dict = {
                    'chrom': sv[0],
                    'start': sv[1],
                    'stop': sv[2],
                    'ref': sv[3],
                    'alt': sv[4],
                    'ann': sv[5],
                    'lof': sv[6],
                    'nmd': sv[7],
                    'genotype': sv[8],
                    'F:R_ref_F:R_alt': sv[9]
                }
                overlap['svs'].append(sv_dict)
        for exp in exp_list:
            exp_overlap = (exp[0] == cnv[0] and exp[1] in cnv_range)
            if exp_overlap:
                exp_dict = {
                    'chrom': exp[0],
                    'start': exp[1],
                    'stop' : exp[2],
                    'alt' : exp[3],
                    'gene': exp[4],
                    'genotype': exp[5],
                    'Alt:Ref': f'{exp[6][1]}:{exp[6][0]}'
                }
                overlap['exps'].append(exp_dict)
        empty = (
            not overlap['snps']
            and not overlap['svs']
            and not overlap['exps']
        )
        if not empty: overlaps.append(overlap)
    return overlaps


def make_json(panel, gene_panel, snp_list, sv_list, exp_list, path_cnvs, sample_id, pmc, interactions):
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
        'all_expansions':[],
    }
    data['cnv_interactions'] = {
        'all_interactions': interactions
    }
    data['metadata'] = {
        'genes_in_panel': gene_panel,
        'supporting_literature':[],
        'snp_qc_metrics': [
            {
                'name': 'snpeff_impact',
                'max_score': '1',
                'min_score': '0',
                'HIGH': '1',
                'MEDIUM': '0.5',
                'all_else': '0'
            },
            {
                'name': 'snpeff_lof',
                'max_score': '1',
                'min_score': '0',
                'percentage_of_lof': '0...1'
            },
            {
                'name': 'CADD_phred',
                'max_score': '1',
                'min_score': '0',
                'scale': 'log10',
                'CADD': '0...99',
                'CADD < 15': '0',
                'CADD >= 15': '0.01...1'
            },
            {
                'name': 'phastCons100way_vertebrate',
                'max_score': '1',
                'min_score': '-1',
                'scale': '0...1 to -1...1'
            },
            {
                'name': 'MutationTaster_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'A': '+  value',
                'D': '+  value * 0.5',
                'N': '-  value * 0.5',
                'P': '-  value'
            },
            {
                'name': 'MutationAssessor_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'H': '+  value',
                'M': '+  value * 0.5',
                'L': '-  value * 0.5',
                'N': '-  value'
            },
            {
                'name': 'LRT_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'N': '-  value'
            },
            {
                'name': 'FATHMM_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'T': '-  value'
            },
            {
                'name': 'MetaSVM_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'T': '-  value'
            },
            {
                'name': 'PROVEAN_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'N': '-  value'
            },
            {
                'name': 'Polyphen2_HVAR_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'P': '+  value * 0.75',
                'B': '-  value'
            },
            {
                'name': 'SIFT_pred',
                'max_score': '1',
                'min_score': '-1',
                'value': '1/ length of result list',
                'D': '+  value',
                'T': '-  value'
            },
            {
                'name': 'GERP___RS',
                'max_score': '1',
                'min_score': '-1',
                'GERP >= 0': '0...1',
                'GERP < 0': '0...-1'
            },
            {
                'name': 'dbNSFP_ExAC_AF',
                'max_score': '1',
                'min_score': '0',
                'scale': 'log10',
                'gnomad': '0...1',
                'no_value': '1'
            },
            {
                'name': '1000Gp3_AF',
                'max_score': '1',
                'min_score': '0',
                'scale': 'log10',
                'gp3': '0...1',
                'no_value': '1'
            },
            {
                'name': 'final_impact',
                'max_score': '1',
                'min_score': '0',
                'scale': '-9.5...15 to 0...1',
                'accepted_values': '>= 0.87'
            }
        ],
        'sv_qc_metrics': [
                        {
                'name': 'snpeff_impact',
                'max_score': '1',
                'min_score': '0',
                'HIGH': '1',
                'MEDIUM': '0.5',
                'all_else': '0'
            },
            {
                'name': 'snpeff_lof',
                'max_score': '1',
                'min_score': '0',
                'percentage_of_lof': '0...1'
            },
            {
                'name': 'final_impact',
                'max_score': '2',
                'min_score': '0',
                'accepted_values': '>= 1'
            }
        ],
        'cnv_qc_metrics': [
            {
                'name': 'ClassifyCNV',
                'max_score': '1.8',
                'min_score': '-2.6',
                'any_promoters_or_enhancers': [
                    {
                        'answer': 'yes',
                        'score': '0'
                    },
                    {
                        'answer': 'no',
                        'score': '-.6'
                    }
                ],
                'number_of_genes': [
                    {
                        'number_of_genes': '< 2',
                        'score': '0'
                    },
                    {
                        'number_of_genes': '< 3',
                        'score': '0.45'
                    },
                    {
                        'number_of_genes': '> 3',
                        'score': '0.9'
                    }
                ],
                'region_overlap_del': [
                    {
                        'overlap_type': 'no overlap',
                        'inside_benign': [
                            {
                                'answer': 'yes',
                                'score': '-1'
                            },
                            {
                                'answer': 'no',
                                'multiple_predictors_genes': [
                                    {
                                        'answer': 'yes',
                                        'score': '0.15'
                                    },
                                    {
                                        'answer': 'no',
                                        'pop_freq_high': [
                                            {
                                                'answer': 'yes',
                                                'score': '-1'
                                            },
                                            {
                                                'answer': 'no',
                                                'score': '0'
                                            }
                                        ]
                                    }
                                ]
                            }
                        ]
                    },
                    {
                        'overlap_type': 'complete overlap',
                        'score': '1'
                    },
                    {
                        'overlap_type': 'partial overlap with the 5\' end, coding sequence involved',
                        'score': '0.9'
                    },
                    {
                        'overlap_type': 'only the last exon involved',
                        'score': '0.3'
                    },
                    {
                        'overlap_type': 'last exon and other exons involved, nonsense-mediated decay expected',
                        'score': '0.9'
                    },
                    {
                        'overlap_type': 'both CNV breakpoints are withint the same genes + PVS1',
                        'score': '0.9'
                    }
                ],
                'region_overlap_dup': [
                    {
                        'overlap_type': 'CNV fully contained inside triplosensative region',
                        'score': '1'
                    },
                    {
                        'overlap_type': 'no overlap',
                        'score': '0',
                        'popultaion_fre_high': '-1'
                    },
                    {
                        'overlap_type': 'CNV identical to a benign region',
                        'score': '-1',
                        'popultaion_fre_high': '-1'
                    },
                    {
                        'overlap_type': 'larger than a benign region, no other protein-coding genes included',
                        'score': '-1',
                        'popultaion_fre_high': '-1'
                    },
                    {
                        'overlap_type': 'smaller than a benign region, breakpoints don\'t interrupt protein-coding genes',
                        'score': '-1',
                        'popultaion_fre_high': '-1'
                    }
                ]
            },
            {
                'name': 'final_impact',
                'max_score': '1.8',
                'min_score': '-2.6',
                'accepted_values': '>= 1'
            }
        ],
        'all_variants_qc': [
            {
                'filter': 'PASS'
            }
        ],
        'pipeline': [
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
        genes.append(str(snp[5]).split('|')[3])
        snp_dict = {
            'chrom': snp[0],
            'start': snp[1],
            'stop': snp[2],
            'ref': snp[3],
            'alt': snp[4],
            'annotation': snp[5],
            'lof': snp[6],
            'nmd': snp[7],
            'genotype': snp[8],
            'F:R_ref_F:R_alt': snp[9],
            'cadd': snp[10],
            'gnomad_freq': snp[11],
            '1000_genomes_freq': snp[12]
        }
        if float(snp[13]) >= 0.9: data['snp']['pathogenic_snp'].append(snp_dict)
        else: data['snp']['likely_pathogenic_snp'].append(snp_dict)
    for sv in sv_list:
        genes.append(str(sv[5]).split('|')[3])
        sv_dict = {
            'chrom': sv[0],
            'start': sv[1],
            'stop': sv[2],
            'ref': sv[3],
            'alt': sv[4],
            'ann': sv[5],
            'lof': sv[6],
            'nmd': sv[7],
            'genotype': sv[8],
            'F:R_ref_F:R_alt': sv[9]
        }
        if float(sv[10]) >= 1.5: data['sv']['pathogenic_sv'].append(sv_dict)
        else: data['sv']['likely_pathogenic_sv'].append(sv_dict)
    for cnv in path_cnvs:
        genes.append(cnv[5])
        cnv_dict = {
            'chrom': cnv[0],
            'cyto': cnv[1],
            'start': cnv[3],
            'stop': cnv[4],
            'alt': cnv[2],
            'genes': cnv[5],
            'RCN:MRCN': f'{cnv[6]}:{cnv[7]}'
        }
        if float(cnv[8]) > 1.0: data['cnv']['pathogenic_cnv'].append(cnv_dict)
        else: data['cnv']['likely_pathogenic_cnv'].append(cnv_dict)
    for exp in exp_list:
        genes.append(exp[4])
        exp_dict = {
            'chrom': exp[0],
            'start': exp[1],
            'stop' : exp[2],
            'alt' : exp[3],
            'gene': exp[4],
            'genotype': exp[5],
            'Alt:Ref': f'{exp[6][1]}:{exp[6][0]}'
        }
        data['exp']['all_expansions'].append(exp_dict)
        # if any(x >= 115 for x in exp[5]): data['exp']['pathogenic_exp'].append(exp_dict)
        # else: data['exp']['likely_pathogenic_exp'].append(exp_dict)
    if panel != None: data['metadata']['supporting_literature'] = get_literature(panel, genes, pmc)
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
    pmc = args.l
    if cpus == None: cpus = 1
    snp_list, sv_list, cnv_list, exp_list = filter_vcf(vcf, panel, args.low_coverage)
    cnv_contigs = check_entries(cnv_list)
    cnv_bed(cnv_contigs, sample_id)
    call_classify_cnv(cpus, sample_id, script)
    cnv_contig_determinations = get_cnv_determination(sample_id)
    path_cnvs = filter_cnv(cnv_contigs, cnv_contig_determinations)
    interactions = check_interactions(path_cnvs, snp_list, sv_list, exp_list)
    make_json(args.p, panel, snp_list, sv_list, exp_list, path_cnvs, sample_id, pmc, interactions)


if __name__ == '__main__':
    main()
