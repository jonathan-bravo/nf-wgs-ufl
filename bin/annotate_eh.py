#!/usr/bin/env python3

import argparse
import json
from pysam import VariantFile
from pysam import bcftools
from pysam import tabix_compress

def parse_args():
    """Parse input arguments.
    """
    parser = argparse.ArgumentParser(
        description = ''
    )
    parser.add_argument(
        '-v',
        metavar = '--VCF',
        type = str,
        help = 'The input sample vcf',
        required = True
    )
    parser.add_argument(
        '-c',
        metavar = '--VARIANT_CATALOG',
        type = str,
        help = '',
        required = False
    )
    args = parser.parse_args()
    return args


def index_input_vcf(vcf):
    """
    """
    tabix_compress(vcf, f'{vcf}.gz')
    bcftools.index('--tbi', f'{vcf}.gz')


def load_catalog(json_file):
    """
    """
    catalog = {}
    with open(json_file, 'r') as f:
        data = json.load(f)
    for variant_entry in data:
        locus_id = variant_entry['LocusId']
        disease = variant_entry['Disease']
        inheritance_mode = variant_entry['InheritanceMode']
        norm_max = int(variant_entry['NormalMax'])
        path_min = int(variant_entry['PathologicMin'])
        source = variant_entry['SourceDisplay']
        source_id = variant_entry['SourceId']
        ref_region = variant_entry['ReferenceRegion']
        try:
            var_id = variant_entry['VariantId']
            path_region = variant_entry['PathologicRegion']
        except:
            var_id = None
            path_region = ref_region
        if var_id:
            for i, vid in enumerate(var_id):
                catalog_entry_dict = {
                    'disease': disease,
                    'inheritance_mode': inheritance_mode,
                    'norm_max': norm_max,
                    'path_min': path_min,
                    'source': source,
                    'source_id': source_id,
                    'path_region': path_region,
                }
                catalog_entry_dict['ref_region'] = ref_region[i]
                catalog[vid] = catalog_entry_dict
        else:
            catalog_entry_dict = {
                'disease': disease,
                'inheritance_mode': inheritance_mode,
                'norm_max': norm_max,
                'path_min': path_min,
                'source': source,
                'source_id': source_id,
                'path_region': path_region,
            }
            catalog_entry_dict['ref_region'] = ref_region
            catalog[locus_id] = catalog_entry_dict
    return catalog


def update_vcf_header(vcf):
    """
    """
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"DISEASE"),
            ('Number',1),
            ('Type','String'),
            ('Description','Disease associated with repeat expansion')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"IM"),
            ('Number',1),
            ('Type','String'),
            ('Description','Inheritance mode')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"NM"),
            ('Number',1),
            ('Type','Integer'),
            ('Description','Normal expansion range max value')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"PM"),
            ('Number',1),
            ('Type','Integer'),
            ('Description','Pathologic expansion range min value')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"SOURCE"),
            ('Number',1),
            ('Type','String'),
            ('Description','Reference literature resource type, eg GeneReviews or PubMed')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"SOURCEID"),
            ('Number',1),
            ('Type','String'),
            ('Description','PMID or GeneReviews book ID for references')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"REFREG"),
            ('Number',1),
            ('Type','String'),
            ('Description','Reference region for the repeat expansion')
        ]
    )
    vcf.header.add_meta(
        'INFO', items=[
            ('ID',"PATHREG"),
            ('Number',1),
            ('Type','String'),
            ('Description','Pathogenic region for the repeat expansion')
        ]
    )


def write_updated_vcf(sample_id, vcf, catalog):
    """
    """
    with open(f'{sample_id}_ann.vcf', 'w') as out_vcf:
        out_vcf.write(str(vcf.header))
        for variant in vcf.fetch():
            var_id = variant.info['VARID']
            variant.info.update({'DISEASE': catalog[var_id]['disease']})
            variant.info.update({'IM': catalog[var_id]['inheritance_mode']})
            variant.info.update({'NM': catalog[var_id]['norm_max']})
            variant.info.update({'PM': catalog[var_id]['path_min']})
            variant.info.update({'SOURCE': catalog[var_id]['source']})
            variant.info.update({'SOURCEID': catalog[var_id]['source_id']})
            variant.info.update({'REFREG': catalog[var_id]['ref_region']})
            variant.info.update({'PATHREG': catalog[var_id]['path_region']})
            out_vcf.write(str(variant))


def main():
    """
    """
    args = parse_args()
    sample_id = args.v.split('.')[0]
    index_input_vcf(args.v)
    vcf = VariantFile(f'{args.v}.gz')
    update_vcf_header(vcf)
    catalog = load_catalog(args.c)
    write_updated_vcf(sample_id, vcf, catalog)


if __name__ == '__main__':
    main()