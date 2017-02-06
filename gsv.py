"""
OOP this whole thing? (Util, Variant, Bam, Vcf), defining some kind of interface would make plugging objects in easy
"""

import pdb # TODO: debug

import pysam, argparse, vcf

def get_parsed_args ():
    """
    Sets up the commmand line parser.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run test cases')
    parser.add_argument('-b', '--input_bam', help='Path to input BAM file')
    parser.add_argument('-v', '--input_vcf', help='Path to input VCF file')
    parser.add_argument('--fetch_flank', type=int, default=20, help='How far around confidence interval to pull reads from BAM')
    return parser.parse_args()

# TODO: check chromosome length to make sure query bounds don't go too far
def get_query_bounds (svtype, position, conf_int, fetch_flank):
    """
    Returns the boudns to query the BAM file around the given position/confidence interval for the specified svtype and fetch_flank.
    >>> get_query_bounds('DEL', 10, (-4, 10), 1)
    (5, 22)
    >>> get_query_bounds('DEL', 100, (-10, 0), 20)
    (70, 121)
    """
    if svtype == 'DEL':
        return (max(position + conf_int[0] - fetch_flank, 0), position + conf_int[1] + fetch_flank + 1)

def bam_reads_in_range (bam, chrom_id, lower_bound, upper_bound):
    """
    Fetches the reads from the bam from the specified location/range.
    """
    yield bam.fetch(chrom_id, lower_bound, upper_boud)

def variants_from_vcf (input_vcf):
    """
    Return all the variants from the input_vcf file with fields in their order (svtype, pos, conf_int).
    """
    yield [('chr1', 'DEL', 5, (-2, 5)), ('chr1', 'DEL', 100, (-43, 25))]


if __name__ == '__main__':
    pdb.set_trace()    
    # Get the command line parser
    args = get_parsed_args()
    if args.test:
        import doctest
        doctest.testmod()
    else:
        # TODO: Parse this with vcf?
        input_vcf = args.input_vcf
        # Loop over all variants to get svtype, pos, conf_int from VCF file
        for variant in variants_from_vcf(input_vcf):
            chrom_id = variant[0]
            svtype = variant[1]
            pos = variant[2]
            conf_int = variant[3]
        
            # Set up the input BAM file
            input_bam = pysam.AlignmentFile(args.input_bam, 'rb')

            # Find where we want to query the input BAM
            lower_bound, upper_bound = get_query_bounds(svtype, pos, conf_int)

            for read in bam_reads_in_range(bam, lower_bound, upper_bound):
                print(read)
