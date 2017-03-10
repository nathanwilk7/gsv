"""
TODO: How to test funcs that depend on stuff, test directory with sample stuff? Or just create fake objects and add the needed fields
TODO: VCF output print for real and file match up and take care of fields with values like True (IMPRECISE), do VCF stuff for real
TODO: Deal with other svtypes besides DEL
TODO: Make this python2 and python3 compatible
"""

#from __future__ import print_function # TODO: Python2
import pdb # TODO: debug
import pysam, argparse, re

# Constants and adjustable parameters
# TODO: Make these easy to adjust via command line arguments

FETCH_FLANK = 200
MIN_ALIGNED = 200 # 200 good
MIN_PCT_ALIGNED = 0.80 # 0.85 worse, 0.80 good
SPLIT_SLOP = 300 # 300 good, 250 good, 150 gets worse genotypes correct
READ_DEPTH_SKIP = 50 # 50 good
READ_DEPTH_INTERVAL = 50 # 50 good
MISMATCH_SLOP = 100 # 100 good, 300 good
MISMATCH_PCT = 0.30
EXTRA_CLIP_SLOP = 1000

MIN_PCT_HET = 0.15

def get_parsed_args ():
    """
    Sets up the commmand line parser and returns the parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run test cases')
    parser.add_argument('-b', '--input_bam', help='Path to input BAM file')
    parser.add_argument('-v', '--input_vcf', help='Path to input VCF file')
    parser.add_argument('-o', '--output_vcf', help='Path to output VCF file')
    return parser.parse_args()

def read_header (input_vcf_str):
    """
    Returns the header of the input VCF file
    """
    header = ''
    with open(input_vcf_str, 'r') as f:
        for line in f:
            if line[:2] == '##':
                header += line
            else:
                break
    # TODO: Make the columns legit
    header += '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878\n'
    return header

def write_output_vcf_header (output_vcf_str, header):
    """
    Writes the output VCF header to the specified filepath or standard out.
    """
    if output_vcf_str is not None:
        with open(output_vcf_str, 'w') as f:
            f.write(header)
    else:
        print(header, end='')

def get_query_bounds (svtype, position, conf_int, fetch_flank, chrom_length):
    """
    Returns the bounds to query the BAM file around the given position/confidence 
     interval for the specified svtype and fetch_flank. Basically gets all the reads
     that are +/- the confidence interval and fetch flank around the position.
    >>> get_query_bounds('DEL', 10, (-4, 10), 1)
    (5, 22)
    >>> get_query_bounds('DEL', 100, (-10, 0), 20)
    (70, 121)
    """
    return max(position + conf_int[0] - fetch_flank, 0), min(position + conf_int[1] + fetch_flank + 1, chrom_length)

def fetch_one_loci_reads (input_bam, chrom, lower_bound, upper_bound):
    """
    Get the reads from the bam that are on the specified chromosome and 
     within the specified bounds.
    """
    reads = {}
    for read in input_bam.fetch(chrom, lower_bound, upper_bound):
        if read.is_unmapped or read.is_duplicate:
            continue
        if read.query_name in reads:
            reads[read.query_name].append(read)
        else:
            reads[read.query_name] = []
            reads[read.query_name].append(read)
    return reads

def add_reads (reads, reads_to_add):
    """
    Adds the read from reads_to_add to reads. Puts reads with the same name (splitters) 
     into the same read_list and avoid adding reads again if they are already in reads.
    """
    for read_list in reads_to_add.values():
        for read in read_list:
            if read.query_name in reads:
                do_add = True
                for read2 in reads[read.query_name]:
                    if read == read2:
                        do_add = False
                if do_add:
                    reads[read.query_name].append(read)
            else:
                reads[read.query_name] = []
                reads[read.query_name].append(read)

def combine_reads (reads, left_reads, right_reads):
    """
    Combines left_reads and right_reads into reads using the add_reads function.
    """
    add_reads(reads, left_reads)
    add_reads(reads, right_reads)
    return reads

def fetch_reads (input_bam, svtype, chrom, left_pos, left_conf_int, right_pos, right_conf_int, fetch_flank, chrom_length):
    """
    Fetches all the reads from the bam at the given location using the parameters specified and the get_query_bounds function.
    """
    left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, FETCH_FLANK, chrom_length)
    left_reads = fetch_one_loci_reads(input_bam, chrom, left_lower_bound, left_upper_bound)

    right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, FETCH_FLANK, chrom_length)
    right_reads = fetch_one_loci_reads(input_bam, chrom, right_lower_bound, right_upper_bound)

    reads = {}
    reads = combine_reads(reads, left_reads, right_reads)
    return reads

# TODO: Consider conf_int
def spans_breakpoint (read, pos, min_aligned, min_pct_aligned):
    """
    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = max(0, pos - min_aligned), pos + min_aligned
    return read.get_overlap(lower_bound, upper_bound) >= 2 * min_aligned * min_pct_aligned

# TODO: Consider conf_int
def split_by_breakpoint (read_list, left_pos, right_pos, split_slop):
    """
    Returns True if a split happens on the left or side within the split slop of the corresponding 
     position and if the other split aligns to the other side of the variant. Otherwise return False.
     Assumes that the read_list contains all split pieces that should be considered.
    """
    for read in read_list:
        split_left = read.reference_end >= left_pos - split_slop and read.reference_end <= left_pos + split_slop
        split_right = read.reference_start >= right_pos - split_slop and read.reference_start <= right_pos + split_slop
        if split_left or split_right:
            for read2 in read_list:
                if read == read2:
                    continue
                if split_left:
                    if read2.reference_start >= right_pos - split_slop and read2.reference_start<= right_pos + split_slop:
                        return True

                if split_right:
                    if read2.reference_end>= left_pos - split_slop and read2.reference_start<= left_pos + split_slop:
                        return True
    return False
"""
        if read.has_tag('SA'):
            split_left = read.reference_end >= left_pos - split_slop and read.reference_end <= left_pos + split_slop
            split_right = read.reference_start >= right_pos - split_slop and read.reference_start <= right_pos + split_slop
            if split_left or split_right:
                split_alignment_list = read.get_tag('SA').rstrip(';').split(';')
                mates = []
                for split_alignment in split_alignment_list:
                    sa = split_alignment.split(',')
                    mate = {}
                    mate['chrom'] = sa[0]
                    mate['pos'] = int(sa[1]) - 1
                    mate['is_reverse'] = sa[2] == '-'
                    #mate['cigar'] = cigarstring_to_tuple(sa[3])
                    mate['mapping_quality'] = int(sa[4])
                    mates.append(mate)

                if split_left:
                    for mate in mates:
                        if mate['pos'] >= right_pos - split_slop and mate['pos'] <= right_pos + split_slop:
                            return True

                if split_right:
                    for mate in mates:
                        if mate['pos'] >= left_pos - split_slop and mate['pos'] <= left_pos + split_slop:
                            return True

    return False # If none of the reads in the list were valid splitters
"""

def clipped_by_breakpoint (read, left_pos, right_pos, split_slop, extra_clip_slop):
    """
    Returns True if the read was clipped on it's right or left side by the variant. Otherwise False.
    """
    right_side_clipped = read.reference_end >= left_pos - split_slop and read.reference_end <= left_pos + split_slop
    left_side_clipped = read.reference_start >= right_pos - split_slop and read.reference_start <= right_pos + split_slop
    if right_side_clipped:
        cigar = read.cigartuples
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            right_clip_size = cigar[-1][1]
            if read.reference_end + right_clip_size >= left_pos + split_slop + extra_clip_slop:
                return True

    if left_side_clipped:
        cigar = read.cigartuples
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            left_clip_size = cigar[0][1]
            if read.reference_start - left_clip_size <= right_pos - split_slop - extra_clip_slop:
                return True

    return False

def get_coverage_diff_skip (read, left_pos, right_pos, chrom_length, read_depth_skip, read_depth_interval):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version skips the specified amount around the positions.
    """
    # TODO: conf_int might extend too far into SV, could be better to have read_depth_interval parameter around pos
    left_outside_outer = max(left_pos - read_depth_skip - read_depth_interval, 0)
    left_outside_inner = max(left_pos - read_depth_skip, 0)
    left_inside_inner = min(left_pos + read_depth_skip + 1, chrom_length)
    left_inside_outer = min(left_pos + read_depth_skip + read_depth_interval + 1, chrom_length)
    left_outside_overlap = read.get_overlap(left_outside_outer, left_outside_inner)
    left_inside_overlap = read.get_overlap(left_inside_inner, left_inside_outer)

    right_inside_inner = max(right_pos - read_depth_skip - read_depth_interval, 0)
    right_inside_outer = max(right_pos - read_depth_skip, 0)
    right_outside_inner = min(right_pos + read_depth_skip + 1, chrom_length)
    right_outside_outer = min(right_pos + read_depth_skip + read_depth_interval + 1, chrom_length)
    right_inside_overlap = read.get_overlap(right_inside_outer, right_inside_inner)
    right_outside_overlap = read.get_overlap(right_outside_inner, right_outside_outer)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        return 1.0
    return float(inside_overlap) / float(outside_overlap)

def get_coverage_diff (read, left_pos, right_pos, chrom_length, read_depth_interval):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version starts counting overlap from the positions.
    """
    left_outside = max(left_pos - read_depth_interval, 0)
    left_inside = min(left_pos + read_depth_interval, chrom_length)
    left_outside_overlap = read.get_overlap(left_outside, left_pos)
    left_inside_overlap = read.get_overlap(left_pos, left_inside)

    right_inside = max(right_pos - read_depth_interval, 0)
    right_outside = min(right_pos + read_depth_interval, chrom_length)
    right_inside_overlap = read.get_overlap(right_inside, right_pos)
    right_outside_overlap = read.get_overlap(right_pos, right_outside)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        return 1.0
    return float(inside_overlap) / float(outside_overlap) 

def get_coverage_diff_conf_int (read, left_pos, left_conf_int, right_pos, right_conf_int, chrom_length, fetch_flank):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version starts counting from the confidence intervals.
    """
    # TODO: conf_int might extend too far into SV, could be better to have read_depth_interval parameter around pos
    left_outside_outer = max(left_pos + left_conf_int[0] - fetch_flank, 0)
    left_outside_inner = max(left_pos + left_conf_int[0], 0)
    left_inside_inner = min(left_pos + left_conf_int[1] + 1, chrom_length)
    left_inside_outer = min(left_pos + left_conf_int[1] + fetch_flank + 1, chrom_length)
    left_outside_overlap = read.get_overlap(left_outside_outer, left_outside_inner)
    left_inside_overlap = read.get_overlap(left_inside_inner, left_inside_outer)

    right_inside_inner = max(right_pos + right_conf_int[0] - fetch_flank, 0)
    right_inside_outer = max(right_pos + right_conf_int[0], 0)
    right_outside_inner = min(right_pos + right_conf_int[1] + 1, chrom_length)
    right_outside_outer = min(right_pos + right_conf_int[1] + fetch_flank + 1, chrom_length)
    right_inside_overlap = read.get_overlap(right_inside_outer, right_inside_inner)
    right_outside_overlap = read.get_overlap(right_outside_inner, right_outside_outer)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        return 1.0
    return float(inside_overlap) / float(outside_overlap)

def is_mismatched_over_sv (read, left_pos, right_pos, mismatch_slop, mismatch_pct):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= left_pos + mismatch_slop
    enters_sv_from_right = read.reference_start >= left_pos + mismatch_slop
    if not enters_sv_from_left and not enters_sv_from_right:
        return False
    if enters_sv_from_left:
        length_over_sv = read.reference_end - (left_pos + mismatch_slop)
        length_of_sv = (right_pos - mismatch_slop) - (left_pos + mismatch_slop)
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(left_pos + mismatch_slop, left_pos + mismatch_slop + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
    
    if enters_sv_from_right:
        length_over_sv = (right_pos - mismatch_slop) - read.reference_start
        length_of_sv = (right_pos - mismatch_slop) - (left_pos + mismatch_slop)
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(right_pos - mismatch_slop - min_length, right_pos - mismatch_slop)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
            
    return False

def get_genotype (alt_support, ref_support, min_pct_het):
    """
    Returns the genotype based on the alternate and reference support.
    """
    total_support = alt_support + ref_support
    if total_support > 0:
        if ref_support > total_support * min_pct_het and alt_support > total_support * min_pct_het:
            genotype = '0/1'
        elif alt_support > ref_support:
            genotype = '1/1'
        else:
            genotype = '0/0'
    else:
        genotype = '1/1'
    return genotype

class Variant:
    """
    Used to store variant info for output VCF.
    """
    def __init__ (self):
        """
        Sets up fields that should be used to output this
        """
        self.chrom = None
        self.pos = None
        self.id = None
        self.ref = None
        self.alt = None
        self.qual = None
        self.filter = None
        self.info = {}
        self.format = {}

def write_variant (output_vcf_str, output_variant, genotype):
    """
    Write the variant to the output VCF or standard out.
    """
    if output_vcf_str is not None:
        with open(output_vcf_str, 'a') as f:
            f.write(output_variant.chrom + '\t')
            f.write(str(output_variant.pos) + '\t')
            f.write(output_variant.id + '\t')
            f.write(output_variant.ref + '\t')
            f.write(output_variant.alt + '\t')
            f.write(output_variant.qual + '\t')
            f.write(output_variant.filter + '\t')
            for k, v in output_variant.info.items():
                if isinstance(v, int) or isinstance(v, str):
                    f.write('{}={}'.format(k, v) + ';')
                else:
                    f.write(k + '=')
                    for i in range(len(v) - 1):
                        f.write(str(v[i]) + ',')
                    f.write(str(v[len(v)-1]) + ';')
            f.write('\t')
            temp_keys = output_variant.format.keys()
            num_keys = len(temp_keys)
            temp_count = 0
            for k in temp_keys:
                if temp_count >= num_keys - 1:
                    break
                temp_count += 1
                f.write(k + ':')
            f.write(k + '\t')
            temp_count = 0
            for k in temp_keys:
                if temp_count >= num_keys - 1:
                    break
                temp_count += 1
                f.write(str(output_variant.format[k]) + ':')
            f.write(str(output_variant.format[k]) + '\n')
    else:
        print(output_variant.chrom, end='\t')
        print(str(output_variant.pos), end='\t')
        print(output_variant.id, end='\t')
        print(output_variant.ref, end='\t')
        print(output_variant.alt, end='\t')
        print(output_variant.qual, end='\t')
        print(output_variant.filter, end='\t')
        for k, v in output_variant.info.items():
            if k == 'IMPRECISE':
                print(k, end=';')
            elif isinstance(v, int) or isinstance(v, str):
                print('{}={}'.format(k, v), end=';')
            else:
                print(k, end='=')
                for i in range(len(v) - 1):
                    print(str(v[i]), end=',')
                print(str(v[len(v)-1]), end=';')
        print('', end='\t')
        print('GT', end='\t')
        print(genotype)
        #temp_keys = output_variant.format.keys()
        #num_keys = len(temp_keys)
        #temp_count = 0
        #for k in temp_keys:
        #    if temp_count >= num_keys - 1:
        #        break
        #    temp_count += 1
        #    print(k, end=':')
        #print(k, end='\t')
        #temp_count = 0
        #for k in temp_keys:
        #    if temp_count >= num_keys - 1:
        #        break
        #    temp_count += 1
        #    print(output_variant.format[k], end=':')
        #print(output_variant.format[k], end='\n')    

def replace_none (field):
    """
    Replace fields that are None with periods.
    """
    if field == None:
        return '.'
    else:
        return field

def replace_empty (_list):
    """
    Replaces empty lists with a period.
    """
    if len(_list) == 0:
        return '.'
    else:
        return _list

def genotype_variants (input_vcf_str, input_bam_str, output_vcf_str):
    """
    Genotype the variants specified in the VCF file using data from the BAM file.
    """
    # Read header from input VCF and write it to the output VCF
    header = read_header(input_vcf_str)
    write_output_vcf_header(output_vcf_str, header)

    # Set up the input VCF and input BAM
    input_vcf = pysam.VariantFile(input_vcf_str)
    # TODO: multiple input bams/samples
    input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

    # Loop over all variants in input VCF file
    for variant in input_vcf.fetch():
        # Abstraction layer for variant info
        chrom = variant.chrom
        svtype = variant.info['SVTYPE']
        left_pos = variant.pos
        right_pos = variant.info['END']
        # TODO: Use the 95% conf int instead to limit the size of this thing?
        left_conf_int = variant.info['CIPOS']
        right_conf_int = variant.info['CIEND']
        # TODO: left and right on different chromosomes?
        chrom_length = input_bam.lengths[input_bam.gettid(chrom)]

        # Get the reads around the variant on both sides
        reads = fetch_reads (input_bam, svtype, chrom, left_pos, left_conf_int, 
                             right_pos, right_conf_int, FETCH_FLANK, chrom_length)
        
        # These numbers are used to genotype
        ref_support = 0
        alt_support = 0

        # TODO: Instead of ifs, get probabilities of support 
        # Extract support from individual reads
        for read_list in reads.values():
            for read in read_list:
                if spans_breakpoint(read, left_pos, MIN_ALIGNED, MIN_PCT_ALIGNED):
                    ref_support += 0.5

                if spans_breakpoint(read, right_pos, MIN_ALIGNED, MIN_PCT_ALIGNED):
                    ref_support += 0.5

                if clipped_by_breakpoint(read, left_pos, right_pos, SPLIT_SLOP, EXTRA_CLIP_SLOP):
                    alt_support += 0

                if is_mismatched_over_sv(read, left_pos, right_pos, MISMATCH_SLOP, MISMATCH_PCT):
                    alt_support += 0

                #coverage_diff = get_coverage_diff_skip(read, left_pos, right_pos, chrom_length, READ_DEPTH_SKIP, READ_DEPTH_INTERVAL)
                coverage_diff = get_coverage_diff(read, left_pos, right_pos, chrom_length, READ_DEPTH_INTERVAL)
                if coverage_diff >= 0.8:
                    ref_support += 0 # min(coverage_diff - 1, 2)
                elif coverage_diff < 0.1:
                    alt_support += 0 # (1 - coverage_diff) * 0.20

        # If there are splittlers that have pieces near both breakpoints, they're in the same read_list
        for read_list in reads.values():
            if split_by_breakpoint(read_list, left_pos, right_pos, SPLIT_SLOP):
                alt_support += 1

        # Genotype based on support for alternate and reference
        genotype = get_genotype(alt_support, ref_support, MIN_PCT_HET)

        # Set up output variant info
        output_variant = Variant()
        output_variant.chrom = variant.chrom
        output_variant.pos = variant.pos
        output_variant.id = variant.id
        output_variant.ref = variant.ref
        output_variant.alt = variant.alts[0]
        output_variant.qual = replace_none(variant.qual)
        output_variant.filter = replace_empty(variant.filter)
        for k in variant.info.keys():
            output_variant.info[k] = variant.info[k]

        # Write this variant and its genotype
        write_variant(output_vcf_str, output_variant, genotype)

if __name__ == '__main__':
    """
    This function is called when this program is called from the command line. It sets 
     up the arguments and gets everything started.
    """
    args = get_parsed_args()
    if args.test:
        import doctest
        doctest.testmod()
    else:
        input_vcf_str = args.input_vcf
        input_bam_str = args.input_bam
        output_vcf_str = args.output_vcf
        genotype_variants(input_vcf_str, input_bam_str, output_vcf_str)

### MISC ###

"""
# Added this per Brent's request but I'll need to do some more refactoring to make it work well, also as long as I'm prototyping it is a pain to maintain...
def get_section_bounds (svtype, pos_left, conf_int_left, chrom_length_left, pos_right, conf_int_right, chrom_length_right, fetch_flank, min_aligned):
    query_bounds_left = max(pos_left + conf_int_left[0] - fetch_flank, 0), min(pos_left + conf_int_left[1] + fetch_flank + 1, chrom_length_left)
    query_bounds_right = max(pos_right + conf_int_right[0] - fetch_flank, 0), min(pos_right + conf_int_right[1] + fetch_flank + 1, chrom_length_right)

    reference_bounds_left = max(0, pos_left - min_aligned), pos_left + min_aligned
    reference_bounds_right = max(0, pos_right - min_aligned), pos_right + min_aligned

    splitter_bounds_left = pos_left - split_slop, pos_left + split_slop
    splitter_bounds_right = pos_right - split_slop, pos_right + split_slop

    return query_bounds_left, query_bounds_right, reference_bounds_left, reference_bounds_right, splitter_bounds_left, splitter_bounds_right

left_query_bounds, right_query_bounds, left_reference_bounds, right_reference_bounds, left_splitter_bounds, right_splitter_bounds = get_section_bounds (svtype, left_pos, left_conf_int, chrom_length, right_pos, right_conf_int, chrom_length, FETCH_FLANK, MIN_ALIGNED)
left_lower_bound, left_upper_bound = left_query_bounds
right_lower_bound, right_upper_bound = right_query_bounds
left_lower_ref, left_upper_ref = left_reference_bounds
right_lower_ref, right_upper_ref = right_reference_bounds
left_lower_split, left_upper_split = left_splittler_bounds
right_lower_split, right_upper_split = right_splittler_bounds
"""
