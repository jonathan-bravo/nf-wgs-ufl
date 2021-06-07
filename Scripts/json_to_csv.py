#!/usr/bin/env python3

import argparse
import json
import pandas as pd
from pprint import pprint

def parse_args():
    """Parse input arguments.

    Keyword arguments:
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'inputs for checking TP, FP, TN, FN values of cnv data'
    )
    parser.add_argument(
        '-j',
        metavar = '--JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    args = parser.parse_args()
    return args


def parse_variant(json, var_type):
    """
    """
    path = pd.DataFrame(json[var_type][f'pathogenic_{var_type}'])
    likely = pd.DataFrame(json[var_type][f'likely_pathogenic_{var_type}'])
    path_status = ['Pathogenic'] * len(path)
    likely_status = ['Likely Pathogenic'] * len(likely)
    if var_type == 'snp' or var_type == 'sv':
        vus = pd.DataFrame(json[var_type][f'vus_{var_type}'])
        vus_status = ['VUS'] * len(vus) # test
        status = path_status + likely_status + vus_status
        result = pd.concat([path, likely, vus])
    else:
        status = path_status + likely_status
        result = pd.concat([path, likely])
    result.insert(0, 'Status', status)
    return result


def parse_qc(json):
    """
    """
    var_types = ('snp', 'sv', 'cnv')
    qc = []
    for var_type in var_types:
        name = []
        max_score = []
        min_score = []
        scale = []
        for metric in json['metadata'][f'{var_type}_qc_metrics']:
            name.append(metric['name'])
            max_score.append(metric['max_score'])
            min_score.append(metric['min_score'])
            if 'scale' in metric.keys(): scale.append(metric['scale'])
            elif 'value' in metric.keys(): scale.append(metric['value'])
            else: scale.append('NaN')
        tag = [var_type] * len(name)
        qc.append(pd.DataFrame({'QC Type': tag, 'Name': name, 'Max Score': max_score, 'Min Score': min_score, 'Scale': scale}))
    return pd.concat(qc)


def parse_cnv_interactions(json):
    """
    """
    interactions = []
    index = 0
    if json['cnv_interactions']['all_interactions'] != []:
        for interaction in json['cnv_interactions']['all_interactions']:
            cnv = str((
                interaction['CNV']['Chrom'],
                interaction['CNV']['Start'],
                interaction['CNV']['Stop'],
                interaction['CNV']['Alt']
            ))
            snps = interaction['SNPs']
            for snp in snps:
                current_snp = str((
                    snp['Chrom'],
                    snp['Start'],
                    snp['Stop'],
                    snp['Gene']
                ))
                interactions.append(pd.DataFrame({'CNV': cnv, 'SNP': current_snp}, index = [index]))
                index += 1
            svs = interaction['SVs']
            for sv in svs:
                current_sv = str((
                    sv['Chrom'],
                    sv['Start'],
                    sv['Stop'],
                    sv['Gene']
                ))
                interactions.append(pd.DataFrame({'CNV': cnv, 'SV': current_sv}, index = [index]))
                index += 1
            exps = interaction['EXPs']
            for exp in exps:
                current_exp = str((
                    exp['Chrom'],
                    exp['Start'],
                    exp['Stop'],
                    exp['Gene']
                ))
                interactions.append(pd.DataFrame({'CNV': cnv, 'EXP': current_exp}, index = [index]))
                index += 1
        return pd.concat(interactions)
    else: return pd.DataFrame()

def read_json(json_file):
    """
    """
    with open(json_file, 'r') as jason_file:
        data = json.load(jason_file)
        snp = parse_variant(data, 'snp')
        sv = parse_variant(data, 'sv')
        cnv = parse_variant(data, 'cnv')
        exp = pd.DataFrame(data['exp']['all_expansions'])
        interactions = parse_cnv_interactions(data)
        genes = pd.DataFrame(data['metadata']['genes_in_panel'], columns = ['genes'])
        lit = pd.DataFrame(data['metadata']['supporting_literature'], columns = ['literature'])
        tools = pd.DataFrame(data['metadata']['pipeline'])
        qc = parse_qc(data)
    return (snp, sv, cnv, exp, interactions, genes, lit, tools, qc)


def write_xlsx(data, file_name):
    """
    """
    with pd.ExcelWriter(f'{file_name}.xlsx') as writer:
        data[0].to_excel(writer, sheet_name = 'SNPs', index = False)
        data[1].to_excel(writer, sheet_name = 'SVs', index = False)
        data[2].to_excel(writer, sheet_name = 'CNVs', index = False)
        data[3].to_excel(writer, sheet_name = 'Expansions', index = False)
        data[4].to_excel(writer, sheet_name = 'Compound Variants', index = False)
        data[5].to_excel(writer, sheet_name = 'Genes in Panel', index = False)
        data[6].to_excel(writer, sheet_name = 'Supporting Literature', index = False)
        data[7].to_excel(writer, sheet_name = 'Tools in Pipeline', index = False)
        data[8].to_excel(writer, sheet_name = 'QC', index = False)


def main():
    args = parse_args()
    file_name = args.j.split('.')[0]
    data = read_json(args.j)
    write_xlsx(data, file_name)
            

if __name__ == '__main__':
    main()
