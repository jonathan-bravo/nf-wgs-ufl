#!/usr/bin/env python3

from argparse import ArgumentParser
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from pysam import VariantFile


def parse_args():
    parser = ArgumentParser()
    parser.add_argument(
        '-v',
        metavar = '--VCF',
        nargs = '+',
        help = 'The input sample vcfs',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'sample id to for the output parquet file',
        required = True
    )
    return parser.parse_args()


def get_header_length(vcf):
    pysam_vcf = VariantFile(vcf)
    header_count = len(list(pysam_vcf.header.records))
    del(pysam_vcf)
    return header_count


def read_vcf(vcf, header_count):
    col_names = [
        'Contig',
        'Position',
        'ID',
        'Ref Allele',
        'Alt Allele',
        'Quality',
        'Filter',
        'Info',
        'Format',
        'Sample'
    ]
    data_types = {
        'Contig': str,
        'Position': 'int64',
        'ID': str,
        'Ref Allele': str,
        'Alt Allele': str,
        'Quality': str,
        'Filter': str,
        'Info': str,
        'Format': str,
        'Sample': str
    }
    vcf_stream = pd.read_csv(
        vcf,
        compression='gzip',
        skiprows=header_count+1,
        sep='\t',
        chunksize=100_000,
        low_memory=False,
        names=col_names,
        dtype=data_types
    )
    return vcf_stream


def write_parquet(sample_id, vcf_type, vcf_stream):
    fields = [
        pa.field('Sample ID', pa.string()),
        pa.field('Contig', pa.string()),
        pa.field('Position', pa.int64()),
        pa.field('ID', pa.string()),
        pa.field('Ref Allele', pa.string()),
        pa.field('Alt Allele', pa.string()),
        pa.field('Quality', pa.string()),
        pa.field('Filter', pa.string()),
        pa.field('Info', pa.string()),
        pa.field('Format', pa.string()),
        pa.field('Sample', pa.string()),
    ]
    for i, chunk in enumerate(vcf_stream):
        print("Chunk", i)
        if i == 0:
            parquet_schema = pa.schema(fields)
            parquet_writer = pq.ParquetWriter(
                f'{sample_id}.{vcf_type}.parquet',
                parquet_schema,
                compression='snappy'
            )
        chunk['Sample ID'] = f'{sample_id}'
        table = pa.Table.from_pandas(chunk, schema=parquet_schema)
        parquet_writer.write_table(table)
    parquet_writer.close()


def main():
    args = parse_args()

    for vcf in args.v:
        if 'variants.' in vcf:
            vcf_type = 'SNV'
        elif '.cnv.' in vcf:
            vcf_type = 'CNV'
        elif '.en.' in vcf:
            vcf_type = 'EXP'
        elif 'diploidSV.' in vcf:
            vcf_type = 'INDEL'
        header_count = get_header_length(vcf)
        vcf_stream = read_vcf(vcf, header_count)
        write_parquet(args.s, vcf_type, vcf_stream)


if __name__ == '__main__':
    main()