#!/usr/bin/env python3

import argparse
from json import load, dump


def parse_args():
    """Parse input arguments.

    Keyword arguments:
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-p',
        metavar = '--PROBAND_JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    parser.add_argument(
        '-m',
        metavar = '--MATERNAL_JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    parser.add_argument(
        '-f',
        metavar = '--PATERNAL_JSON',
        type = str,
        help = 'the input JSON report',
        required = True
    )
    args = parser.parse_args()
    return args


def read_proband_json(proband):
    """
    """
    variants = {}
    with open(proband, 'r') as p:
        data = load(p)
        snvs = data['snp']['all_snps']
        svs = data['sv']['all_svs']
        cnvs = data['cnv']['all_cnvs']
        exps = data['exp']['all_expansions']
    for snv in snvs:
        variants[f"{snv['Chrom']};{snv['Start']};{snv['Stop']};{snv['Alt Allele']};{snv['Ref Allele']}"] = ['proband']
    for sv in svs:
        variants[f"{sv['Chrom']};{sv['Start']};{sv['Stop']};{sv['Alt Allele']};{sv['Ref Allele']}"] = ['proband']
    for cnv in cnvs:
        variants[f"{cnv['Chrom']};{cnv['Start']};{cnv['Stop']};{snv['Alt Allele']}"] = ['proband']
    for exp in exps:
        variants[f"{exp['Chrom']};{exp['Start']};{exp['Stop']};{exp['Alt Allele']};{exp['Ref Allele']}"] = ['proband']
    return variants
        

def check_inheritance(proband_variants, parental_json, trio):
    """
    """
    with open(parental_json, 'r') as j:
        data = load(j)
        snvs = data['snp']['all_snps']
        svs = data['sv']['all_svs']
        cnvs = data['cnv']['all_cnvs']
        exps = data['exp']['all_expansions']
    for snv in snvs:
        try: proband_variants[f"{snv['Chrom']};{snv['Start']};{snv['Stop']};{snv['Alt Allele']};{snv['Ref Allele']}"].append(trio)
        except KeyError: continue
    for sv in svs:
        try: proband_variants[f"{sv['Chrom']};{sv['Start']};{sv['Stop']};{sv['Alt Allele']};{sv['Ref Allele']}"].append(trio)
        except KeyError: continue
    for cnv in cnvs:
        try: proband_variants[f"{cnv['Chrom']};{cnv['Start']};{cnv['Stop']};{snv['Alt Allele']}"].append(trio)
        except KeyError: continue
    for exp in exps:
        try: proband_variants[f"{exp['Chrom']};{exp['Start']};{exp['Stop']};{exp['Alt Allele']};{exp['Ref Allele']}"].append(trio)
        except KeyError: continue
    return proband_variants


def convert_trio(trio):
    """
    """
    if trio == ['proband']: inheritance = 'de novo'
    elif trio == ['proband', 'maternal']: inheritance = 'maternal'
    elif trio == ['proband', 'paternal']: inheritance = 'paternal'
    elif trio == ['proband', 'maternal', 'paternal']: inheritance = 'both'
    return inheritance


def make_new_proband_json(proband_json, proband_variants):
    """
    """
    with open(proband_json, 'r') as p:
        data = load(p)
        snvs = data['snp']['all_snps']
        svs = data['sv']['all_svs']
        cnvs = data['cnv']['all_cnvs']
        exps = data['exp']['all_expansions']
    for snv in snvs:
        trio = proband_variants[f"{snv['Chrom']};{snv['Start']};{snv['Stop']};{snv['Alt Allele']};{snv['Ref Allele']}"]
        inheritance = convert_trio(trio)
        snv['Inheritance'] = inheritance
    for sv in svs:
        trio = proband_variants[f"{sv['Chrom']};{sv['Start']};{sv['Stop']};{sv['Alt Allele']};{sv['Ref Allele']}"]
        inheritance = convert_trio(trio)
        sv['Inheritance'] = inheritance
    for cnv in cnvs:
        trio = proband_variants[f"{cnv['Chrom']};{cnv['Start']};{cnv['Stop']};{snv['Alt Allele']}"]
        inheritance = convert_trio(trio)
        cnv['Inheritance'] = inheritance
    for exp in exps:
        trio = proband_variants[f"{exp['Chrom']};{exp['Start']};{exp['Stop']};{exp['Alt Allele']};{exp['Ref Allele']}"]
        inheritance = convert_trio(trio)
        exp['Inheritance'] = inheritance
    data['snp']['all_snps'] = snvs
    data['sv']['all_svs'] = svs
    data['cnv']['all_cnvs'] = cnvs
    data['exp']['all_expansions'] = exps
    data['small_var']['small_variants'] = []
    for snv in snvs:
        data['small_var']['small_variants'].append(snv)
    for sv in svs:
        data['small_var']['small_variants'].append(sv)
    with open(f'new_proband_report.json', 'w') as outfile:
            dump(data, outfile, indent = 4)


def main():
    """
    """
    args = parse_args()
    proband = read_proband_json(args.p)
    proband = check_inheritance(proband, args.m, 'maternal')
    proband = check_inheritance(proband, args.f, 'paternal')
    make_new_proband_json(args.p, proband)


if __name__ == '__main__':
    main()