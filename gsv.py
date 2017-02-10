"""
TODO: OOP this whole thing? (Util, Variant, Bam, Vcf), defining some kind of interface would make plugging objects in easy
TODO: One off errors...I am certain they're everywhere
TODO: How to test funcs that depend on 
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
    # if svtype == util.SVTYPE_DEL:
    return max(position + conf_int[0] - fetch_flank, 0), position + conf_int[1] + fetch_flank + 1

def read_depth_at_pos (bam, chrom_id, pos):
    return sum(1 for _ in bam.fetch(chrom_id, pos, pos+1))

def percent_bases_aligned_in_range (read, start_pos, end_pos):
    return read.get_overlap(start_pos, end_pos) / (end_pos - start_pos)

def genotype_read_depth (start_read_depth, mid_read_depth, end_read_depth):
    avg_end_read_depth = (start_read_depth + end_read_depth) / 2
    if mid_read_depth <= avg_end_read_depth * util.MAX_FRAC_HOM_ALT_READ_DEPTH:
        return util.HOM_ALT
    elif mid_read_depth <= avg_end_read_depth * util.MAX_FRAC_HET_READ_DEPTH:
        return util.HET
    else:
        return util.HOM_REF

def genotype_splitters (splitter_count, reference_read_count):
    # TODO: Think about fixing this division by zero for real...
    if splitter_count == 0 and reference_read_count == 0:
        return util.HOM_ALT
    frac_splitters = splitter_count / (splitter_count + reference_read_count)
    if frac_splitters <= util.MAX_FRAC_HOM_REF_SPLITTERS:
        return util.HOM_REF
    elif frac_splitters <= util.MAX_FRAC_HET_SPLITTERS:
        return util.HET
    else:
        return util.HOM_ALT

def genotype_variants (input_vcf_str, input_bam_str, fetch_flank):
    input_vcf = pysam.VariantFile(input_vcf_str)

    # Loop over all variants to get svtype, pos, conf_int from VCF file
    # natee = 0
    for variant in input_vcf.fetch():
        # if not natee < 15:
            # break
        # natee += 1
        # Set up the input BAM file
        # TODO: multiple input bams
        input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

        # Get variant info
        chrom_id = variant.chrom
        svtype = variant.info[util.VCF_SVTYPE]
        pos = variant.pos
        conf_int = variant.info[util.VCF_CIPOS]
        end_pos = variant.info[util.VCF_END]

        # We can use read depth if it's a deletion
        if svtype == util.SVTYPE_DEL:
            # Check read depth at start, end, and in the middle of the deletion
            start_read_depth = read_depth_at_pos(input_bam, chrom_id, pos)
            mid_read_depth = read_depth_at_pos(input_bam, chrom_id, (pos + end_pos) // 2)
            end_read_depth = read_depth_at_pos(input_bam, chrom_id, end_pos)

        # Get range from which to get BAM reads
        lower_bound, upper_bound = get_query_bounds(svtype, pos, conf_int, fetch_flank)

        # Count vars for splitters vs reference
        splitter_count = 0
        reference_count = 0

        # Loop through reads around SV site and classify them
        for read in input_bam.fetch(chrom_id, lower_bound, upper_bound):
            if svtype == util.SVTYPE_DEL:
                frac_aligned = percent_bases_aligned_in_range(read, pos, end_pos)
                if frac_aligned < util.MIN_FRAC_ALIGN_REF:
                    splitter_count += 1
                else:
                    reference_count += 1


        read_depth_genotype = genotype_read_depth(start_read_depth, mid_read_depth, end_read_depth)
        splitter_genotype = genotype_splitters(splitter_count, reference_count)

        print("SVPOS: {}".format(pos))
        #print("Start Read Depth: {}".format(start_read_depth))
        #print("Mid Read Depth: {}".format(mid_read_depth))
        #print("End Read Depth: {}".format(end_read_depth))
        #print("Splitter Count: {}".format(splitter_count))
        #print("Reference Count: {}".format(reference_count))
        #print("GENOTYPE")
        print("Read Depth Genotype: {}".format(read_depth_genotype))
        print("Splitter Genotype: {}\n".format(splitter_genotype))

if __name__ == '__main__':
    # pdb.set_trace() # TODO: debug
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

        genotype_variants(input_vcf_str, input_bam_str, fetch_flank)
        exit()
        
        input_vcf = pysam.VariantFile(input_vcf_str)

        # Loop over all variants to get svtype, pos, conf_int from VCF file
        natee = 0
        for variant in input_vcf.fetch():
            if not natee < 15:
                break
            natee += 1
            # Set up the input BAM file
            # TODO: multiple input bams
            input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

            # Get variant info
            chrom_id = variant.chrom
            svtype = variant.info[util.VCF_SVTYPE]
            pos = variant.pos
            conf_int = variant.info[util.VCF_CIPOS]
            end_pos = variant.info[util.VCF_END]

            # We can use read depth if it's a deletion
            if svtype == util.SVTYPE_DEL:
                # Check read depth at start, end, and in the middle of the deletion
                start_read_depth = read_depth_at_pos(input_bam, chrom_id, pos)
                mid_read_depth = read_depth_at_pos(input_bam, chrom_id, (pos + end_pos) // 2)
                end_read_depth = read_depth_at_pos(input_bam, chrom_id, end_pos)

            # Get range from which to get BAM reads
            lower_bound, upper_bound = get_query_bounds(svtype, pos, conf_int, fetch_flank)

            # Count vars for splitters vs reference
            splitter_count = 0
            reference_count = 0

            # Loop through reads around SV site and classify them
            for read in input_bam.fetch(chrom_id, lower_bound, upper_bound):
                if svtype == util.SVTYPE_DEL:
                    frac_aligned = percent_bases_aligned_in_range(read, pos, end_pos)
                    if frac_aligned < util.MIN_FRAC_ALIGN_REF:
                        splitter_count += 1
                    else:
                        reference_count += 1


            read_depth_genotype = genotype_read_depth(start_read_depth, mid_read_depth, end_read_depth)
            splitter_genotype = genotype_splitters(splitter_count, reference_count)

            print("SVPOS: {}".format(pos))
            print("Start Read Depth: {}".format(start_read_depth))
            print("Mid Read Depth: {}".format(mid_read_depth))
            print("End Read Depth: {}".format(end_read_depth))
            print("Splitter Count: {}".format(splitter_count))
            print("Reference Count: {}".format(reference_count))
            print("GENOTYPE")
            print("Read Depth Genotype: {}".format(read_depth_genotype))
            print("Splitter Genotype: {}\n".format(splitter_genotype))
