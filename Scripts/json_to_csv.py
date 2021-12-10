#!/usr/bin/env python3

import argparse
import json
import pandas as pd

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
                interaction['CNV']['Alt'],
                interaction['CNV']['Determination']
            ))
            snps = interaction['SNPs']
            for snp in snps:
                current_snp = str((
                    snp['Chrom'],
                    snp['Start'],
                    snp['Stop'],
                    snp['Gene'],
                    snp['REVEL'],
                    snp['CADD']
                ))
                interactions.append(pd.DataFrame({'CNV': cnv, 'SNP': current_snp}, index = [index]))
                index += 1
            svs = interaction['SVs']
            for sv in svs:
                current_sv = str((
                    sv['Chrom'],
                    sv['Start'],
                    sv['Stop'],
                    sv['Gene'],
                    sv['CADD']
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
    with open(json_file, 'r') as j_file:
        data = json.load(j_file)
        snp = pd.DataFrame(data['snp']['all_snps'])
        sv = pd.DataFrame(data['sv']['all_svs'])
        cnv = pd.DataFrame(data['cnv']['all_cnvs'])
        exp = pd.DataFrame(data['exp']['all_expansions'])
        interactions = parse_cnv_interactions(data)
        genes = pd.DataFrame(data['metadata']['genes_in_panel'], columns = ['genes'])
        lit = pd.DataFrame(data['metadata']['supporting_literature'], columns = ['literature'])
        tools = pd.DataFrame(data['metadata']['pipeline'])
    return (snp, sv, cnv, exp, interactions, genes, lit, tools)


def write_xlsx(data, file_name):
    """
    """
    if data[6].empty:
        with pd.ExcelWriter(f'{file_name}.xlsx') as writer:
            data[0].to_excel(writer, sheet_name = 'SNVs', index = False)
            data[1].to_excel(writer, sheet_name = 'SVs', index = False)
            data[2].to_excel(writer, sheet_name = 'CNVs', index = False)
            data[3].to_excel(writer, sheet_name = 'Expansions', index = False)
            data[4].to_excel(writer, sheet_name = 'Compound Variants', index = False)
            data[7].to_excel(writer, sheet_name = 'Tools in Pipeline', index = False)
    else:
        with pd.ExcelWriter(f'{file_name}.xlsx') as writer:
            data[0].to_excel(writer, sheet_name = 'SNPs', index = False)
            data[1].to_excel(writer, sheet_name = 'SVs', index = False)
            data[2].to_excel(writer, sheet_name = 'CNVs', index = False)
            data[3].to_excel(writer, sheet_name = 'Expansions', index = False)
            data[4].to_excel(writer, sheet_name = 'Compound Variants', index = False)
            data[5].to_excel(writer, sheet_name = 'Genes in Panel', index = False)
            data[6].to_excel(writer, sheet_name = 'Supporting Literature', index = False)
            data[7].to_excel(writer, sheet_name = 'Tools in Pipeline', index = False)


def main():
    args = parse_args()
    file_name = args.j.split('.')[0]
    data = read_json(args.j)
    write_xlsx(data, file_name)
            

if __name__ == '__main__':
    main()
