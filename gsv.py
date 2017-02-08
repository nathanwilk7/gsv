"""
TODO: OOP this whole thing? (Util, Variant, Bam, Vcf), defining some kind of interface would make plugging objects in easy
"""

import pdb # TODO: debug

import pysam, argparse
import utility as util

def get_parsed_args ():
    """
    Sets up the commmand line parser.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run test cases')
    parser.add_argument('-b', '--input_bam', help='Path to input BAM file')
    parser.add_argument('-v', '--input_vcf', help='Path to input VCF file')
    parser.add_argument('--fetch_flank', type=int, default=util.FETCH_FLANK, help='How far around confidence interval to pull reads from BAM')
    return parser.parse_args()

# TODO: check chromosome length to make sure query bounds don't go too far
def get_query_bounds (svtype, position, conf_int, fetch_flank):
    """
    Returns the bounds to query the BAM file around the given position/confidence interval for the specified svtype and fetch_flank.
    >>> get_query_bounds(util.SVTYPE_DEL, 10, (-4, 10), 1)
    (5, 22)
    >>> get_query_bounds(util.SVTYPE_DEL, 100, (-10, 0), 20)
    (70, 121)
    """
    if svtype == util.SVTYPE_DEL:
        return max(position + conf_int[0] - fetch_flank, 0), position + conf_int[1] + fetch_flank + 1

def read_depth_at_pos (bam, chrom_id, pos):
    return sum(1 for _ in bam.fetch(chrom_id, pos, pos+1))

# TODO: Use svtyper logic for splitters?
def is_splitter (read):
    return True

def is_discordant (read):
    return True

def is_reference (read):
    return True

if __name__ == '__main__':
    pdb.set_trace() # TODO: debug
    # Get the command line parser
    args = get_parsed_args()
    if args.test:
        import doctest
        doctest.testmod()
    else:
        # Get args
        input_vcf_str = args.input_vcf
        input_bam_str = args.input_bam
        fetch_flank = args.fetch_flank
        
        
        input_vcf = pysam.VariantFile(input_vcf_str)

        # Loop over all variants to get svtype, pos, conf_int from VCF file
        for variant in input_vcf.fetch():
            # Set up the input BAM file
            # TODO: multiple input bams
            input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

            # Get variant info
            chrom_id = variant.chrom
            svtype = variant.info[util.VCF_SVTYPE]
            pos = variant.pos
            conf_int = variant.info[util.VCF_CIPOS]

            if svtype == util.SVTYPE_DEL:
                # Check read depth at start, end, and in the middle of the deletion
                start_read_depth = read_depth_at_pos(input_bam, chrom_id, pos)
                end_pos = variant.info[util.VCF_END]
                end_read_depth = read_depth_at_pos(input_bam, chrom_id, end_pos)
                mid_read_depth = read_depth_at_pos(input_bam, chrom_id, (pos + end_pos) // 2)

            # Get range from which to get BAM reads
            lower_bound, upper_bound = get_query_bounds(svtype, pos, conf_int, fetch_flank)

            # Pull out all the reads from the BAM file within the specified range
            reference_reads = []
            splitters = []
            discordant_reads = []

            # Loop through reads around SV site and classify them
            for read in input_bam.fetch(chrom_id, lower_bound, upper_bound):
                if svtype == util.SVTYPE_DEL:
                    read_is_splitter = is_splitter(read)
                    read_is_discordant = is_discordant(read)
                    read_is_reference = is_reference(read)
                    
