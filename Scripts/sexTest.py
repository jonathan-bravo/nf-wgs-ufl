#!/usr/bin/env python3

import argparse
from concurrent.futures import ProcessPoolExecutor
from pysam import coverage

CHROM_LIST = [
    "-rchr1",
    "-rchr2",
    "-rchr3",
    "-rchr4",
    "-rchr5",
    "-rchr6",
    "-rchr7",
    "-rchr8",
    "-rchr9",
    "-rchr10",
    "-rchr11",
    "-rchr12",
    "-rchr13",
    "-rchr14",
    "-rchr15",
    "-rchr16",
    "-rchr17",
    "-rchr18",
    "-rchr19",
    "-rchr20",
    "-rchr21",
    "-rchr22",
    "-rchrX",
    "-rchrY"
]

def parse_args():
    """Parse input arguments.

    Keyword arguments:

    -b, --BAM       -- input sorted bam file
    -s, --SAMPLE_ID -- sample id
    -t, --THREADS   -- number of cpus to use (Default: 1)

    Return:

    args -- the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description = 'Test for the sex of a sample from sorted bam file'
    )
    parser.add_argument(
        '-b',
        metavar = '--BAM',
        type = str,
        help = 'The input sample bam',
        required = True
    )
    parser.add_argument(
        '-s',
        metavar = '--SAMPLE_ID',
        type = str,
        help = 'sample id to be included in report name',
        required = True
    )
    parser.add_argument(
        '-t',
        metavar = '--THREADS',
        type = int,
        help = 'The number of cpus to use for multiprocessing',
        required = False
    )
    args = parser.parse_args()
    return args


def get_coverage(chrom, bam_file):
    """Get coverage for the specified chromosome.

    Keyword arguments:

    chrom    -- the specified chromoeoms
    bam_file -- the sample bam file to get the information from

    Return:

    chr_coverage -- a float of the percentage of coverage
    """
    chr_coverage = float(coverage("-H", chrom, bam_file).split('\t')[5])
    return chr_coverage


def chunk_coverage(bam_file, cpus, chrom_list = CHROM_LIST):
    """Split data data into chunks.

    This function uses multiprocessing to execute each chunk of data in a
    different cpu, for as many cpus as selected (Default: 1).

    Keyword arguments:

    bam_file   -- input sample sorted bam file
    cpus       -- the number of cpus to use for multiprocessing (Default: 1)
    chrom_list -- list of chromosomes from CHROM_LIST

    Return:

    results -- the generator object from ProcessPoolExecutor
    """
    if cpus == None: cpus = 1
    sample = [bam_file] * len(chrom_list)
    with ProcessPoolExecutor(max_workers=cpus) as executor:
        results = executor.map(
            get_coverage,
            chrom_list,
            sample
        )
    return results


def get_sex(results):
    """Determine sex of sample.

    This function takes the generator from `chunk_coverage()` and gets
    the mean of coverage percent for the autosomal chromosomes. This mean,
    in conjuction with the sex chromosome concentration is used to
    determine a sample's sex.

    Keyword arguments:

    results -- the generator that was returned from `chunk_coverage()`

    Return:

    sex -- a string value that is either 'Male' or 'Female' 
    """
    generator = []
    for result in results:
        generator.append(result)
    acov_list = [generator[i] for i in range(22)]
    # xcov = generator[22]
    ycov = generator[23]
    cov_mean = sum(acov_list) / len(acov_list)
    if ycov/cov_mean >= 0.25:
        return 'Male'
    return 'Female'


def write_out(sample_id, sex):
    """Create the output file.

    Keyword arguments:

    sample_id -- the sample id to be included in the out file name
    sex       -- the string result to be written from `get_sex()`
    """
    file_name = f'{sample_id}_m_or_f.txt'
    f = open(file_name, "w")
    f.write(sex)
    f.close()


def main():
    """The main function.

    First we parse our arguments, execute chunking/ getting the chromosomal
    coverage percent, and then determining sex.
    """
    args = parse_args()
    results = chunk_coverage(args.b, args.t)
    sex = get_sex(results)
    print(sex)
    write_out(args.s, sex)


if __name__ == '__main__':
    main()
