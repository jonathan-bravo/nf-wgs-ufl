#!/usr/bin/env python3

import argparse
import json
import pandas as pd
import gzip


def parse_sample_id_args():
    """Parse input arguments.

    The input sample id will be used to grab all required qc csv files and all
    the required json data files.

    Keyword arguments:

    -s, --SAMPLE_ID  -- sample id
    -p, --PATH -- path to files
    
    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'sample id for json data files and qc files'
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'the input sample id',
        required = True
    )
    parser.add_argument(
        '-p',
        metavar = '--PATH',
        type = str,
        help = 'the path to all files',
        required = True
    )
    args = parser.parse_args()
    return args


def parse_transcripts(transcripts):
    """Parses the first transcript for the variant

    Each variant in a position has a list of transcripts. We are looping through
    the dictionary to see if any of the transcripts have our desired source,
    refseq, and then we parse the additional info we want.

    Keyword arguments:

    transcripts -- the transcript dictionary for the variant
    
    Return:

    transcript_name -- the name of the transcript
    source          -- the source of the data
    bio_type        -- the type of transcript
    hgnc            -- the hgnc value of the transcript
    hgvsc           -- the hgvsc value of the transcript
    hgvsp           -- the hgvsp value of the transcript
    """
    if transcripts != 'NA':
        for transcript in transcripts:
            try: source = transcript['source']
            except KeyError: source = 'NA'
            if source == 'RefSeq':
                try: transcript_name = transcript['transcript']
                except KeyError: transcript_name = 'NA'
                try: bio_type = transcript['bioType']
                except KeyError: bio_type = 'NA'
                try: hgnc = transcript['hgnc']
                except KeyError: hgnc = "NA"
                try: hgvsc = transcript['hgvsc']
                except KeyError: hgvsc = "NA"
                try: hgvsp = transcript['hgvsp']
                except KeyError: hgvsp = "NA"
                return (
                    transcript_name,
                    source,
                    bio_type,
                    hgnc,
                    hgvsc,
                    hgvsp
                )
            else:
                transcript_name = 'NA'
                source = 'NA'
                bio_type = 'NA'
                hgnc = 'NA'
                hgvsc = 'NA'
                hgvsp = 'NA'
    else:
        transcript_name = transcripts
        source = transcripts
        bio_type = transcripts
        hgnc = transcripts
        hgvsc = transcripts
        hgvsp = transcripts
    return (
        transcript_name,
        source,
        bio_type,
        hgnc,
        hgvsc,
        hgvsp
    )


def parse_clinvar(clinvar):
    """Parses the clinvar data for a variant

    Each variant has some data from ClinVar that we flatten here. Four values
    are returned as a tuple. If the values are not present then a value of NA
    is returned.

    Keyword arguments:

    clinvar -- the first clinvar dictionary for the variant

    Return:

    clinvar_id            -- the id value for the clinvar entry
    clinvar_review_status -- the review status for the clinvar entry
    clinvar_phenotypes    -- the listed phenotypes for the clinvar entry
    clinvar_significance  -- the significance in the clinvar entry
    """
    if clinvar != 'NA':
        try: clinvar_id = clinvar['id']
        except KeyError: clinvar_id = 'NA'
        try: clinvar_review_status = clinvar['reviewStatus']
        except KeyError: clinvar_review_status = 'NA'
        try: clinvar_phenotypes = clinvar['phenotypes']
        except KeyError: clinvar_phenotypes = 'NA'
        try: clinvar_significance = clinvar['significance']
        except KeyError: clinvar_significance = 'NA'
    else:
        clinvar_id = clinvar
        clinvar_review_status = clinvar
        clinvar_phenotypes = clinvar
        clinvar_significance = clinvar
    return (
        clinvar_id,
        clinvar_review_status,
        clinvar_phenotypes,
        clinvar_significance
    )


def parse_clingen(clingen):
    """Parses the clingen data for a SV variant

        Each structural variant has some data from ClinGen that we flatten here.
        Three values are returned as a tuple. If the values are not present then
        a value of NA is returned.

        Keyword arguments:

        clingen -- the first clingen dictionary for the variant

        Return:

        clingen_id              -- the id value for the clingen entry
        clingen_interperitation -- the interperitation for the clingen entry
        clingen_phenotypes      -- the listed phenotypes for the clingen entry
    """
    if clingen != 'NA':
        try: clingen_id = clingen['id']
        except KeyError: clingen_id = 'NA'
        try: clingen_interperitation = clingen['clinicalInterpretation']
        except KeyError: clingen_interperitation = 'NA'
        try: clingen_phenotypes = clingen['phenotypes']
        except KeyError: clingen_phenotypes = 'NA'
    else:
        clingen_id = clingen
        clingen_interperitation = clingen
        clingen_phenotypes = clingen
    return (
        clingen_id,
        clingen_interperitation,
        clingen_phenotypes
    )


def parse_variant(position):
    """Parse the variants for a given position

    Each variant contains many data points we want to extract, the main ones
    being the contig, start position, stop position, reference allele,
    altternate allele, and variant type. On top of these data points are
    additional ones from public datasets like gnomAD and clinical datasets like
    ClinVar and ClinGen. This function takes a position and parses each variant
    listed. Here we only take a small amount of the potential data for each
    variant, so it would be easy to add an additional `try, except` block and
    dictionary key, entry pair for a new piece of information as needed.
    
    Keyword arguments:

    position -- the position dictionary

    Return:

    variants -- a list of dictionaries, one for each variant parsed
    """
    variants = []
    for variant in position['variants']:
        contig = variant['chromosome']
        start = variant['begin']
        stop = variant['end']
        ref = variant['refAllele']
        alt = variant['altAllele']
        var_type = variant['variantType']
        try: 
            hgvsg_vid = variant['hgvsg']
        except KeyError: 
            try:
                hgvsg_vid = variant['vid']
            except KeyError:
                hgvsg_vid = 'NA'
        try: clinvar = parse_clinvar(clinvar = variant['clinvar'][0])
        except KeyError: clinvar = parse_clinvar('NA')
        try: dbsnp = variant['dbsnp'][0]
        except KeyError: dbsnp = 'NA'
        try: global_minor_allele_freq = variant['globalAllele']['globalMinorAlleleFrequency']
        except KeyError: global_minor_allele_freq = 'NA'
        try: gnomad = variant['gnomad']['allAf']
        except KeyError: gnomad = 'NA'
        try: onekg = variant['oneKg']['allAf']
        except KeyError: onekg = 'NA'
        try: revel = variant['revel']['score']
        except KeyError: revel = 'NA'
        try: topmed = variant['topmed']['allAf']
        except KeyError: topmed = 'NA'
        try: transcript = parse_transcripts(transcripts = variant['transcripts'])
        except KeyError: transcript = parse_transcripts(transcripts = 'NA')
        variants.append({
            'Contig': contig,
            'Start': start,
            'Stop': stop,
            'Ref Allele': ref,
            'Alt Allele': alt,
            'Variant Type': var_type,
            'HGVSG/ VID': hgvsg_vid,
            'Clinvar ID': clinvar[0],
            'Clinvar Review Status': clinvar[1],
            'Clinvar Phenotypes': clinvar[2],
            'Clinvar Significance': clinvar[3],
            'dbSNP': dbsnp,
            'Global Minor Allele Freq': global_minor_allele_freq,
            'gnomAD': gnomad,
            'oneKG': onekg,
            'REVEL': revel,
            'topMED': topmed,
            'Transcript': transcript[0],
            'Transcript Source': transcript[1],
            'Transcript Bio Type': transcript[2],
            'Transcript HGNC': transcript[3],
            'Transcript HGVSC': transcript[4],
            'Transcript HGVSP': transcript[5]
        })
    return variants


def parse_position(data, var_type):
    """Parse the positions for a sample

    Each position contains some basic information like sample info, and position
    filtering, but they also contain lists for the variants, and some other
    quality metrics. This function parses each position for an input json file
    and returns a list of dictionaries, one dictionary for each position.
    
    Keyword arguments:

    data     -- the loaded json file/ data
    var_type -- 'SV' or 'SNV'

    Return:

    positions -- a list of dictionaries, one for each position parsed
    """
    positions = []
    for position in data['positions']:
        sample = position['samples'][0]
        try: filter = position['filters']
        except KeyError: filter = 'NA'
        if var_type == 'SNV':
            try: mapping_quality = position['mappingQuality']
            except KeyError: mapping_quality = 'NA'
            try: variant_freq = sample['variantFrequencies']
            except KeyError: variant_freq = 'NA'
            try: total_depth = sample['totalDepth']
            except KeyError: total_depth = 'NA'
            try: allele_depths = sample['alleleDepths']
            except KeyError: allele_depths = 'NA'
            try: somatic_quality = sample['somaticQuality']
            except KeyError: somatic_quality = 'NA'
            variants = parse_variant(position = position)
            positions.append({
                'Filter': filter,
                'Mapping Quality': mapping_quality,
                'Variant Frequency': variant_freq,
                'Total Depth': total_depth,
                'Allele Depths': allele_depths,
                'Somatic Quality': somatic_quality,
                'Variants': variants
            })
        elif var_type == 'SV':
            try: split_read_counts = sample['splitReadCounts']
            except KeyError: split_read_counts = 'NA'
            try: paired_end_read_counts = sample['pairedEndReadCounts']
            except KeyError: paired_end_read_counts = 'NA'
            try: clingen = parse_clingen(position['clingen'][0])
            except KeyError: clingen = parse_clingen('NA')
            variants = parse_variant(position = position)
            positions.append({
                'Filter': filter,
                'Split Read Counts': split_read_counts,
                'Paired End Read Counts': paired_end_read_counts,
                'Clingen ID': clingen[0],
                'Clingen Interperitation': clingen[1],
                'Clingen Phenotypes': clingen[2],
                'Variants': variants
            })
    return positions


def parse_hits(positions, var_type):
    """Take the parsed positions and create flattened 'hits'

    This function takes the parsed positions from `parse_positions()` and
    creates a new list called hits. These hits are just the 'flattened' data
    from the positions and variants within those positions. This function is
    basically just cleaning up the data so it can be made into a pandas
    DataFrame with little effort.

    Keyword arguments:

    positions -- the list of parsed positions from `parse_positions()`
    var_type  -- 'SV' or 'SNV'

    Return:

    df -- the pd DataFrame of the positions
    """
    hits = []
    for position in positions:
        for variant in position['Variants']:
            var_dict = {
                'Contig': variant['Contig'],
                'Start': variant['Start'],
                'Stop': variant['Stop'],
                'Ref Allele': variant['Ref Allele'],
                'Alt Allele': variant['Alt Allele'],
                'Variant Type': variant['Variant Type'],
                'Transcript': variant['Transcript'],
                'Transcript Source': variant['Transcript Source'],
                'Transcript Bio Type': variant['Transcript Bio Type'],
                'Transcript HGNC': variant['Transcript HGNC'],
                'Transcript HGVSC': variant['Transcript HGVSC'],
                'Transcript HGVSP': variant['Transcript HGVSP'],
                'HGVSG/ VID': variant['HGVSG/ VID'],
                'Filter': position['Filter'],
                'Clinvar ID': variant['Clinvar ID'],
                'Clinvar Review Status': variant['Clinvar Review Status'],
                'Clinvar Phenotypes': variant['Clinvar Phenotypes'],
                'Clinvar Significance': variant['Clinvar Significance'],
                'dbSNP': variant['dbSNP'],
                'Global Minor Allele Freq': variant['Global Minor Allele Freq'],
                'gnomAD': variant['gnomAD'],
                'oneKG': variant['oneKG']
            }
            if var_type == 'SNV':  
                var_dict['Mapping Quality'] = position['Mapping Quality']
                var_dict['Variant Frequency'] = position['Variant Frequency']
                var_dict['Total Depth'] = position['Total Depth']
                var_dict['Allele Depths'] = position['Allele Depths']
                var_dict['Somatic Quality'] = position['Somatic Quality']
                var_dict['REVEL'] = variant['REVEL']
                var_dict['topMED'] = variant['topMED']
            elif var_type == 'SV':
                var_dict['Split Read Counts'] = position['Split Read Counts']
                var_dict['Paired End Read Counts'] = position['Paired End Read Counts']
                var_dict['Clingen ID'] = position['Clingen ID']
                var_dict['Clingen Interperitation'] = position['Clingen Interperitation']
                var_dict['Clingen Phenotypes'] = position['Clingen Phenotypes']
    df = pd.DataFrame(hits)
    return df


def parse_json(json_file, var_type):
    """Load in the json file and parse the data

    This function is rather simple. It is loading in the gzipped json data,
    parsing the positions for the input file in `parse_positions()`, and then
    flattening the data in `parse_hits()`. At the end the returned value is the
    same DataFrame that is returned from `parse_hits()`.

    Keyword arguments:

    json_file -- the sv or snv json file to be parsed
    var_type  -- 'SV' or 'SNV'

    Return:

    df -- the DataFrame created in `parse_hits()`
    """
    with gzip.open(json_file, 'r') as j:
        data = json.load(j)
        positions = parse_position(data = data, var_type = var_type)
    df = parse_hits(positions = positions, var_type = var_type)
    return df


def parse_qc_csv(csv_file, qc_type):
    """Parse the qc metrics csv files

    This function handles reading in all the extra qc csv files. The three
    options are the tmb, summary, or coverage csv file.

    Keyword arguments:

    csv_file -- the input qc csv file to be parsed
    qc_type  -- 'TMB', 'SUMMARY', or 'COVERAGE'

    Return:

    df -- a DataFrame of the input csv file
    """
    if qc_type == 'TMB':
        df = pd.read_csv(
            csv_file,
            usecols=[2,3],
            names = ['TMB Summary', 'Value']
        )
    elif qc_type == 'SUMMARY':
        df = pd.read_csv(
            csv_file,
            skiprows=4,
            names = ['DRAGEN Enrichment Summary Report', 'Value']
        )
    elif qc_type == 'COVERAGE':
        df = pd.read_csv(
            csv_file,
            usecols=[2,3],
            names = ['Coverage Summary', 'Value']
        )
    return df


def parse_bed(bed_file):
    """Parse the input bed file

    This function parses the input bed file and returns a DataFrame containing
    only the entries from the bed file where the total coverage for the region
    is at or below 500.

    Keyword arguments:

    bed_file -- the input bed file

    Return:

    df -- a DataFrame of the low coverage regions from the bed file
    """
    low_cov_regions = []
    with open(bed_file, 'r') as bed:
        lines = bed.readlines()[1:]
        for line in lines:
            row = line.split('\t')
            if int(row[5]) <= 500:
                low_cov_regions.append({
                    'Contig': row[0],
                    'Start': row[1],
                    'End': row[2],
                    'Name': row[3],
                    'Gene ID': row[4],
                    'Total Coverage': row[5],
                    'Read1 Coverage': row[6],
                    'Read2 Coverage': row[7]
                })
    df = pd.DataFrame(low_cov_regions)
    return df


def parse_cnv_seg(cnv_seg):
    """
    """
    cnvs = []
    with open(cnv_seg, 'r') as seg:
        for line in seg:
            entry = line.split('\t').strip()
            sample = entry[0]
            chrom = entry[1]
            start = entry[2]
            end = entry[3]
            num_targets = entry[4]
            seg_mean = entry[5]
            seg_call = entry[6]
            qual = entry[7]
            qual_filter = entry[8]
            copy_number = entry[9]
            ploidy = entry[10]
            imp_pairs = entry[11]
            cnvs.append({
                'Sample': sample,
                'Chromosome': chrom,
                'Start': start,
                'End': end,
                'Num_Targets': num_targets,
                'Segment_Mean': seg_mean,
                'Segment_Call': seg_call,
                'Qual': qual,
                'Filter': qual_filter,
                'Copy_Number': copy_number,
                'Ploidy': ploidy,
                'Improper_Pairs': imp_pairs
            })
    return cnvs


def parse_cnv_json(cnv_json):
    """
    """
    cnvs = []
    with gzip.open(cnv_json, 'r') as j:
        data = json.load(j)
        for position in data['positions']:
            contig = position['chromosome']
            start = position['position']
            stop = position['svEnd']
            sv_len = position['svLength']
            cytoband = position['cytogeneticBand']
            variant = position['variants'][0]
            var_type = variant['variantType']
            for transcript in variant['transcripts']:
                if transcript['source'] == 'RefSeq':
                    source = transcript['source']
                    try: transcript_name = transcript['transcript']
                    except KeyError: transcript_name = 'NA'
                    try: bio_type = transcript['bioType']
                    except KeyError: bio_type = 'NA'
                    try: hgnc = transcript['hgnc']
                    except KeyError: hgnc = "NA"
                    break
            if source != 'RefSeq':
                source = variant['transcripts'][0]['source']
                try: transcript_name = variant['transcripts'][0]['transcript']
                except KeyError: transcript_name = 'NA'
                try: bio_type = variant['transcripts'][0]['bioType']
                except KeyError: bio_type = 'NA'
                try: hgnc = variant['transcripts'][0]['hgnc']
                except KeyError: hgnc = "NA"
        cnvs.append({
            'Chromosome': contig,
            'Start': start,
            'End': stop,
            'svLength': sv_len,
            'cytogeneticBand': cytoband,
            'variantType': var_type,
            'transcript': transcript_name,
            'bioType': bio_type,
            'hgnc': hgnc
        })
    return cnvs


def parse_cnvs(seg_cnvs, json_cnvs):
    """
    """
    cnvs = []
    for cnv in seg_cnvs:
        for j in json_cnvs:
            svLength = j['svLength'],
            cytogeneticBand =  j['cytogeneticBand'],
            variantType =  j['variantType'],
            transcript =  j['transcript'],
            bioType =  j['bioType'],
            hgnc =  j['hgnc']
        cnv['svLength'] = svLength
        cnv['cytogeneticBand'] = cytogeneticBand
        cnv['variantType'] = variantType
        cnv['transcript'] = transcript
        cnv['bioType'] = bioType
        cnv['hgnc'] = hgnc
    cnvs.append(cnv)
    df = pd.DataFrame(cnvs)
    return df


def write_xlsx(data, sample_id):
    """Write the output xlsx file

    This function takes all the data generated from parsing the sv and snv
    json files, the tmb, summary, and coverage csv files, and the bed file.
    Each DataFrame is written to its own sheet, or tab, within an excel file/
    workbook.

    Keyword arguments:

    data      -- a tuple containing all the parsed file's data
    sample_id -- the sample id
    """
    with pd.ExcelWriter(f'{sample_id}.xlsx') as writer:
        data[0].to_excel(writer, sheet_name = 'SNVs', index = False)
        data[1].to_excel(writer, sheet_name = 'SVs', index = False)
        data[2].to_excel(writer, sheet_name = 'CNVs', index = False)
        data[3].to_excel(writer, sheet_name = 'TMB', index = False)
        data[4].to_excel(writer, sheet_name = 'SUMMARY', index = False)
        data[5].to_excel(writer, sheet_name = 'COVERAGE', index = False),
        data[6].to_excel(writer, sheet_name = 'LOW COVERAGE', index = False)


def main():
    """Main function that runs

    1. Parse the arguments to get the sample id as `args.s`
    2. Parse the SNV json file; var_type = SNV
    3. Parse the SV json file; var_type = SV
    4. Parse the CNV json file; var_type = CNV
    5. Parse the TMB csv file; qc_type = TMB
    6. Parse the SUMMARY csv file; qc_type = SUMMARY
    7. Parse the COVERAGE csv file; qc_type = COVERAGE
    8. Parse the bed file
    9. Pass all data into the xlsx writer
    10. ???
    11. Profit
    """
    args = parse_sample_id_args()
    sample_id = args.s
    data_path = args.p
    snv_json = f'{data_path}/{sample_id}.hard-filtered.annotations.json.gz'
    sv_json = f'{data_path}/{sample_id}.sv.annotations.json.gz'
    cnv_json = f'{data_path}/{sample_id}.cnv.annotations.json.gz'
    cnv_seg = f'{data_path}/{sample_id}.seg.called.merged'
    tmb_csv = f'{data_path}/{sample_id}.tmb.metrics.csv'
    summary_csv = f'{data_path}/Additional Files/{sample_id}.summary.csv'
    coverage_csv = f'{data_path}/{sample_id}.qc-coverage-region-1_coverage_metrics.csv'
    bed_file = f'{data_path}/{sample_id}.qc-coverage-region-1_read_cov_report.bed'
    snv_hits = parse_json(
        json_file = snv_json,
        var_type = 'SNV'
    )
    sv_hits = parse_json(
        json_file = sv_json,
        var_type = 'SV'
    )
    seg_cnvs = parse_cnv_seg(
        cnv_seg = cnv_seg
    )
    json_cnvs = parse_cnv_json(
        cnv_json = cnv_json
    )
    cnv_hits = parse_cnvs(
        seg_cnvs = seg_cnvs,
        json_cnvs = json_cnvs
    )
    tmb = parse_qc_csv(
        csv_file = tmb_csv,
        qc_type = 'TMB'
    )
    summary = parse_qc_csv(
        csv_file = summary_csv,
        qc_type = 'SUMMARY'
    )
    coverage = parse_qc_csv(
        csv_file = coverage_csv,
        qc_type = 'COVERAGE'
    )
    bed = parse_bed(
        bed_file = bed_file
    )
    write_xlsx(
        data = (
            snv_hits,
            sv_hits,
            cnv_hits,
            tmb,
            summary,
            coverage,
            bed
        ),
        sample_id = sample_id
    )
                

if __name__ == '__main__':
    main()
