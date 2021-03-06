import pdb # TODO: debug
import pysam, argparse, re
import cyvcf2

# Constants and adjustable parameters
# TODO: Make these easy to adjust via command line arguments, or avoid them altogether

FETCH_FLANK = 200
MIN_ALIGNED = 200 # 200 good, 250 same
MIN_PCT_ALIGNED = 0.80 # 0.85 worse, 0.80 good
SPLIT_SLOP = 300 # 300 good, 250 good, 150 gets worse genotypes correct
READ_DEPTH_SKIP = 50 # 50 good
READ_DEPTH_INTERVAL = 50 # 50 good
MISMATCH_SLOP = 100 # 100 good, 300 good
MISMATCH_PCT = 0.90
EXTRA_CLIP_SLOP = 1000
INNER_READ_FETCH_FLANK = 20

MIN_PCT_HET = 0.15 # .15 normally

def get_parsed_args ():
    """
    Sets up the commmand line parser and returns the parsed arguments.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test', action='store_true', help='Run test cases')
    parser.add_argument('-b', '--input_bam', help='Path to input BAM file')
    parser.add_argument('-v', '--input_vcf', help='Path to input VCF file')
    parser.add_argument('-o', '--output_vcf', help='Path to output VCF file')
    parser.add_argument('-f', '--fasta', help='Path to input FASTA file')
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
    >>> get_query_bounds('DEL', 10, (-4, 10), 1, 25)
    (5, 22)
    >>> get_query_bounds('DEL', 100, (-10, 0), 20, 150)
    (70, 121)
    """
    return max(position + conf_int[0] - fetch_flank, 0), min(position + conf_int[1] + fetch_flank + 1, chrom_length)

# TODO: Test with directory
def fetch_one_locus_reads (input_bam, chrom, lower_bound, upper_bound):
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
     into the same read_list and avoids adding reads again if they are already in reads.
     Depends on == being True if two reads are the same read and False otherwise
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
    left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, fetch_flank, chrom_length)
    left_reads = fetch_one_locus_reads(input_bam, chrom, left_lower_bound, left_upper_bound)

    right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, fetch_flank, chrom_length)
    right_reads = fetch_one_locus_reads(input_bam, chrom, right_lower_bound, right_upper_bound)

    reads = {}
    reads = combine_reads(reads, left_reads, right_reads)
    return reads

def get_inner_query_bounds (svtype, left_pos, left_conf_int, right_pos, right_conf_int, inner_read_fetch_flank):
    chrom_length = right_pos + 1
    left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, inner_read_fetch_flank, chrom_length)
    right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, inner_read_fetch_flank, chrom_length)
    
    # If they cross over, just return the middle
    if left_upper_bound > right_lower_bound:
        middle = int((left_upper_bound + right_lower_bound) / 2)
        return middle, middle + 1
    
    return left_upper_bound, right_lower_bound

def fetch_inner_reads (input_bam, svtype, chrom, left_pos, left_conf_int, right_pos, right_conf_int, inner_read_fetch_flank):
    """
    Fetches all the reads from the bam between the left and right sides
    """
    # Chromosome length shouldn't matter since we're looking in between two breakpoints, so just set it something irrelevant
    left_upper_bound, right_lower_bound = get_inner_query_bounds(svtype, left_pos, left_conf_int, right_pos, right_conf_int, inner_read_fetch_flank)
    inner_reads = fetch_one_locus_reads(input_bam, chrom, left_upper_bound, right_lower_bound)
    return inner_reads

# http://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
MATCH = 0  # M
INS = 1    # I
DEL = 2    # D
SKIP = 3   # N
SOFT = 4   # S
HARD = 5   # H
PAD = 6    # P
EQUAL = 7  # =
DIFF = 8   # X

def get_chrom(fasta_fh, chrom):
    """
    Return the chromosome sequence
    """
    return fasta_fh.fetch(chrom)

def expand_match(qry_seq, ref_seq):
    """
    Return an expanded list of CIGAR ops
    based upon nuceltotide matches (=, or 7)
    and mismatches (X, or 8)
    """
    prev_op = None
    curr_op = None
    length = 1
    for idx, q_nucl in enumerate(qry_seq):
        if q_nucl == ref_seq[idx]: curr_op = 7  # EQUAL (=)
        else: curr_op = 8  # DIFF (X)

        if curr_op == prev_op:
            length += 1
        elif prev_op is not None:
            yield (prev_op, length)
            length = 1
        prev_op = curr_op
    yield (curr_op, length)

def get_new_cigar (fa, bam, read, left_pos, right_pos):
    pdb.set_trace()
    curr_chrom_id = read.reference_id

    curr_chrom_seq = get_chrom(fa, bam.get_reference_name(curr_chrom_id))
    ref_pos = read.reference_start
    qry_pos = 0
    new_cigar = []
    for cigar_tuple in read.cigartuples:
        op = cigar_tuple[0]
        op_len = cigar_tuple[1]
        if op == EQUAL or op == DIFF or op == MATCH:
            if op == MATCH:
                qry_seq = read.query_sequence[qry_pos:qry_pos + op_len]
                ref_seq = curr_chrom_seq[ref_pos:ref_pos + op_len]
                if qry_seq == ref_seq:
                    new_cigar.append((7, len(qry_seq))) # EQUAL (=)
                else: # expand the M cigar op into X and = ops
                    for new_cigar_tuple in expand_match(qry_seq, ref_seq):
                        new_cigar.append(new_cigar_tuple)
            else:
                new_cigar.append(cigar_tuple)
            ref_pos += op_len
            qry_pos += op_len
        elif op == DEL or op == SKIP:
            ref_pos += op_len
            new_cigar.append(cigar_tuple)
        elif op == INS or op == SOFT:
            qry_pos += op_len
            new_cigar.append(cigar_tuple)
        elif op == HARD:
            new_cigar.append(cigar_tuple)
            
    return new_cigar

def get_num_matches (fa, bam, read, left_pos, right_pos):
    # to keep track of how many matches there are within the left_pos and right_pos
    num_matches = 0
    pdb.set_trace()
    curr_chrom_id = read.reference_id

    curr_chrom_seq = get_chrom(fa, bam.get_reference_name(curr_chrom_id))
    ref_pos = read.reference_start

    # how much further the left_pos is than the ref_pos, we need to skip this much of the cigar
    diff_ref_vs_left = 0

    if ref_pos > left_pos:
        # do nothing, the gap is mismatches
        a = 1
    elif ref_pos == left_pos:
        # do nothing
        a = 1
    elif ref_pos < left_pos:
        diff_ref_vs_left = left_pos - ref_pos

    #ref_pos += diff_ref_vs_left
    #qry_pos = diff_ref_vs_left
    qry_pos = 0

    # skip tuples until we reach the one before the left_pos
    i = 0
    #previous_tuple = None
    
    #next_tuple = read.cigartuples[1]
    while i < len(read.cigartuples):
        current_tuple = read.cigartuples[i]
        if current_tuple[0] == 5:
            i += 1
            continue
        if qry_pos + current_tuple[1] >= diff_ref_vs_left:
            break
        if (current_tuple[0] == 0 or current_tuple[0] == 1 or current_tuple[0] == 7 or current_tuple[0] == 8 or current_tuple[0] == 4):
            qry_pos += current_tuple[1]
        i += 1

    

    #for cigar_tuple in read.cigartuples:
    while i < len(read.cigartuples):
        current_tuple = read.cigartuples[i]
        op = cigar_tuple[0]
        op_len = cigar_tuple[1]
        if op == EQUAL or op == DIFF or op == MATCH:
            if op == MATCH:
                qry_seq = read.query_sequence[qry_pos:qry_pos + op_len]
                ref_seq = curr_chrom_seq[ref_pos:ref_pos + op_len]
                if qry_seq == ref_seq:
                    num_matches += len(qry_seq)
                else: # expand the M cigar op into X and = ops
                    for new_cigar_tuple in expand_match(qry_seq, ref_seq):
                        if new_cigar_tuple[0] == EQUAL:
                            num_matches += new_cigar_tuple[1]
            else:
                a = 1
            ref_pos += op_len
            qry_pos += op_len
        elif op == DEL or op == SKIP:
            ref_pos += op_len
        elif op == INS or op == SOFT:
            qry_pos += op_len
        elif op == HARD:
            a = 1
            
    return num_matches

def spans_breakpoint_mismatches (fa, bam, read, pos, conf_int, min_aligned, min_pct_aligned):
    """
    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = max(0, pos + conf_int[0] - min_aligned), pos + min_aligned + conf_int[1]
    conf_int_length = abs(conf_int[0] - conf_int[1])
    #new_cigar = get_new_cigar(fa, bam, read, lower_bound, upper_bound)
    num_matches = get_num_matches(fa, bam, read, lower_bound, upper_bound)
    return read.get_overlap(lower_bound, upper_bound) >= (conf_int_length + (2 * min_aligned)) * min_pct_aligned

def spans_breakpoint (read, pos, conf_int, min_aligned, min_pct_aligned):
    """
    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = max(0, pos + conf_int[0] - min_aligned), pos + min_aligned + conf_int[1]
    conf_int_length = abs(conf_int[0] - conf_int[1])
    return read.get_overlap(lower_bound, upper_bound) >= (conf_int_length + (2 * min_aligned)) * min_pct_aligned

def spans_breakpoint_old (read, pos, min_aligned, min_pct_aligned):
    """
    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = max(0, pos - min_aligned), pos + min_aligned
    return read.get_overlap(lower_bound, upper_bound) >= 2 * min_aligned * min_pct_aligned

def split_by_breakpoint (read_list, left_pos, left_conf_int, right_pos, right_conf_int, split_slop):
    """
    Returns True if a split happens on the left or side within the split slop of the corresponding 
     position and if the other split aligns to the other side of the variant. Otherwise return False.
     Assumes that the read_list contains all split pieces that should be considered.
    """
    for read in read_list:
        split_left = read.reference_end >= left_pos + left_conf_int[0] - split_slop and read.reference_end <= left_pos + left_conf_int[1] + split_slop
        split_right = read.reference_start >= right_pos + right_conf_int[0] - split_slop and read.reference_start <= right_pos + right_conf_int[1] + split_slop
        if split_left or split_right:
            for read2 in read_list:
                if read == read2:
                    continue
                if split_left:
                    if read2.reference_start >= right_pos + right_conf_int[0] - split_slop and read2.reference_start <= right_pos + right_conf_int[1] + split_slop:
                        return True

                if split_right:
                    if read2.reference_end >= left_pos + left_conf_int[0] - split_slop and read2.reference_start <= left_pos + left_conf_int[1] + split_slop:
                        return True
    return False

def split_by_breakpoint_old (read_list, left_pos, right_pos, split_slop):
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

def clipped_by_breakpoint (read, left_pos, left_conf_int, right_pos, right_conf_int, split_slop, extra_clip_slop):
    """
    Returns True if the read was clipped on it's right or left side by the variant. Otherwise False.
    """
    right_side_clipped = read.reference_end >= left_pos + left_conf_int[0] - split_slop and read.reference_end <= left_pos + left_conf_int[1] + split_slop
    left_side_clipped = read.reference_start >= right_pos + right_conf_int[0] - split_slop and read.reference_start <= right_pos + right_conf_int[1] + split_slop
    if right_side_clipped:
        cigar = read.cigartuples
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            right_clip_size = cigar[-1][1]
            if read.reference_end + right_clip_size >= left_pos + left_conf_int[1] + split_slop + extra_clip_slop:
                return True

    if left_side_clipped:
        cigar = read.cigartuples
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            left_clip_size = cigar[0][1]
            if read.reference_start - left_clip_size <= right_pos + right_conf_int[0] - split_slop - extra_clip_slop:
                return True

    return False

def clipped_by_breakpoint_old (read, left_pos, right_pos, split_slop, extra_clip_slop):
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

def get_coverage_diff_skip (read, left_pos, left_conf_int, right_pos, right_conf_int, chrom_length, read_depth_skip, read_depth_interval):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version skips the specified amount around the positions.
    """
    # TODO: conf_int might extend too far into SV, could be better to have read_depth_interval parameter around pos
    left_outside_outer = max(left_pos + left_conf_int[0] - read_depth_skip - read_depth_interval, 0)
    left_outside_inner = max(left_pos + left_conf_int[0] - read_depth_skip, 0)
    left_inside_inner = min(left_pos + left_conf_int[1] + read_depth_skip + 1, chrom_length)
    left_inside_outer = min(left_pos + left_conf_int[1] + read_depth_skip + read_depth_interval + 1, chrom_length)
    left_outside_overlap = read.get_overlap(left_outside_outer, left_outside_inner)
    left_inside_overlap = read.get_overlap(left_inside_inner, left_inside_outer)

    right_inside_inner = max(right_pos + right_conf_int[0] - read_depth_skip - read_depth_interval, 0)
    right_inside_outer = max(right_pos + right_conf_int[0] - read_depth_skip, 0)
    right_outside_inner = min(right_pos + right_conf_int[1] + read_depth_skip + 1, chrom_length)
    right_outside_outer = min(right_pos + right_conf_int[1] + read_depth_skip + read_depth_interval + 1, chrom_length)
    right_inside_overlap = read.get_overlap(right_inside_outer, right_inside_inner)
    right_outside_overlap = read.get_overlap(right_outside_inner, right_outside_outer)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        return 1.0
    return float(inside_overlap) / float(outside_overlap)

def get_coverage_diff_skip_old (read, left_pos, right_pos, chrom_length, read_depth_skip, read_depth_interval):
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

def get_coverage_diff (read, left_pos, left_conf_int, right_pos, right_conf_int, chrom_length, read_depth_interval):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version starts counting overlap from the positions.
    """
    left_outside = max(left_pos + left_conf_int[0] - read_depth_interval, 0)
    left_inside = min(left_pos + left_conf_int[1] + read_depth_interval, chrom_length)
    left_outside_overlap = read.get_overlap(left_outside, left_pos)
    left_inside_overlap = read.get_overlap(left_pos, left_inside)

    right_inside = max(right_pos + right_conf_int[0] - read_depth_interval, 0)
    right_outside = min(right_pos + right_conf_int[1] + read_depth_interval, chrom_length)
    right_inside_overlap = read.get_overlap(right_inside, right_pos)
    right_outside_overlap = read.get_overlap(right_pos, right_outside)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        return 1.0
    return float(inside_overlap) / float(outside_overlap) 

def get_coverage_diff_old (read, left_pos, right_pos, chrom_length, read_depth_interval):
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

def is_mismatched_over_sv (read, left_pos, left_conf_int, right_pos, right_conf_int, mismatch_slop, mismatch_pct):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= left_pos + left_conf_int[1] + mismatch_slop
    enters_sv_from_right = read.reference_start >= left_pos + left_conf_int[1] + mismatch_slop
    if not enters_sv_from_left and not enters_sv_from_right:
        return False
    if enters_sv_from_left:
        length_over_sv = read.reference_end - (left_pos + left_conf_int[1] + mismatch_slop)
        length_of_sv = (right_pos + right_conf_int[0] - mismatch_slop) - (left_pos + left_conf_int[1] + mismatch_slop)
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(left_pos + left_conf_int[1] + mismatch_slop, left_pos + left_conf_int[1] + mismatch_slop + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
    
    if enters_sv_from_right:
        length_over_sv = (right_pos + right_conf_int[0] - mismatch_slop) - read.reference_start
        length_of_sv = (right_pos + right_conf_int[0] - mismatch_slop) - (left_pos + left_conf_int[1] + mismatch_slop)
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(right_pos + right_conf_int[0] - mismatch_slop - min_length, right_pos + right_conf_int[0] - mismatch_slop)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
            
    return False

def is_mismatched_over_sv_old (read, left_pos, right_pos, mismatch_slop, mismatch_pct):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= left_pos  + mismatch_slop
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

def is_mismatched_over_sv_ver2 (read, left_pos, left_conf_int, right_pos, right_conf_int, mismatch_pct):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= left_pos + left_conf_int[1]
    enters_sv_from_right = read.reference_start > left_pos + left_conf_int[1]
    if not enters_sv_from_left and not enters_sv_from_right:
        return False
    if enters_sv_from_left:
        length_over_sv = read.reference_end - left_pos + left_conf_int[1]
        length_of_sv = (right_pos + right_conf_int[0]) - (left_pos + left_conf_int[1])
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(left_pos + left_conf_int[1], left_pos + left_conf_int[1] + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
    
    if enters_sv_from_right:
        length_over_sv = right_pos + right_conf_int[0] - read.reference_start
        length_of_sv = (right_pos + right_conf_int[0]) - (left_pos + left_conf_int[1])
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(right_pos + right_conf_int[0] - min_length, right_pos + right_conf_int[0])
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
            
    return False

def is_mismatched_over_sv_ver2_old (read, left_pos, right_pos, mismatch_pct):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= left_pos
    enters_sv_from_right = read.reference_start > left_pos
    if not enters_sv_from_left and not enters_sv_from_right:
        return False
    if enters_sv_from_left:
        length_over_sv = read.reference_end - left_pos
        length_of_sv = right_pos - left_pos
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(left_pos, left_pos + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
    
    if enters_sv_from_right:
        length_over_sv = right_pos - read.reference_start
        length_of_sv = right_pos - left_pos
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return False
        matches = read.get_overlap(right_pos - min_length, right_pos)
        mismatches = min_length - matches
        if mismatches / min_length > mismatch_pct:
            return True
            
    return False

def get_bases_possible_one_side (lower_bound, upper_bound):
    return upper_bound - lower_bound

def get_bases_possible (left_lower_bound, left_upper_bound, right_lower_bound, right_upper_bound):
    return get_bases_possible_one_side(left_lower_bound, left_upper_bound), get_bases_possible_one_side(right_lower_bound, right_upper_bound) 
    
def get_num_bases_covered (read, lower_bound, upper_bound):
    return read.get_overlap(lower_bound, upper_bound)

def get_cov_pct_avg_base (bases_covered_actual, bases_covered_possible, read_count):
    if bases_covered_possible > 0 and read_count > 0:
        coverage_pct = bases_covered_actual / (bases_covered_possible * read_count)
        avg_base_coverage = bases_covered_actual / bases_covered_possible
    else:
        coverage_pct = 0
        avg_base_coverage = 0
    return coverage_pct, avg_base_coverage

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

def cyvcf_genotype (genotype_str):
    if genotype_str == '0/0':
        return [[0, 0, False]]
    elif genotype_str == '0/1':
        return [[0, 1, False]]
    elif genotype_str == '1/1':
        return [[1, 1, False]]
    else:
        print('Genotype not implemented: ', genotype_str)
        exit()

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

def write_variant ():
    a = 1


"""
def write_variant (output_vcf_str, output_variant, genotype):
    "#""
    Write the variant to the output VCF or standard out.
    "#""
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
                if len(v) > 0:
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
    "#""
    Replace fields that are None with periods.
    "#""
    if field == None:
        return '.'
    else:
        return field

def replace_empty (_list):
    "#""
    Replaces empty lists with a period.
    "#""
    if len(_list) == 0:
        return '.'
    else:
        return _list
"""
def genotype_variants (input_vcf_str, input_bam_str, output_vcf_str, input_fasta_str):
    """
    Genotype the variants specified in the VCF file using data from the BAM file.
    """
    # Read header from input VCF and write it to the output VCF
    #header = read_header(input_vcf_str)
    #write_output_vcf_header(output_vcf_str, header)

    # Set up the input VCF and input BAM
    #input_vcf2 = pysam.VariantFile(input_vcf_str)
    input_vcf = cyvcf2.VCF(input_vcf_str)

    # TODO: multiple input bams/samples
    input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

    # Variants will be written to this
    output_vcf = cyvcf2.Writer(output_vcf_str, input_vcf)

    if input_fasta_str is not None:
        fa = pysam.FastaFile(input_fasta_str)
    else:
        fa = None

    # Loop over all variants in input VCF file
    #for variant2 in input_vcf2.fetch():
    for variant in input_vcf:
        # Get variant info into local vars so this can be swapped out with different libraries
        #chrom = variant2.chrom
        #svtype = variant2.info['SVTYPE']
        #left_pos = variant2.pos
        #right_pos = variant2.info['END']
        # TODO: Use the 95% conf int instead to limit the size of this thing?
        #left_conf_int = variant2.info['CIPOS']
        #right_conf_int = variant2.info['CIEND']
        # TODO: left and right on different chromosomes?
        #chrom_length = input_bam.lengths[input_bam.gettid(chrom)]

        chrom = variant.CHROM
        svtype = variant.INFO.get('SVTYPE')
        left_pos = variant.start
        right_pos = variant.end # variant.INFO.get('END')
        left_conf_int = variant.INFO.get('CIPOS')
        right_conf_int = variant.INFO.get('CIEND')
        chrom_length = input_bam.lengths[input_bam.gettid(chrom)]

        # Get the reads around the variant on both sides
        reads = fetch_reads(input_bam, svtype, chrom, left_pos, left_conf_int, 
                            right_pos, right_conf_int, FETCH_FLANK, chrom_length)

        # Get the reads in the middle of the variant
        inner_reads = fetch_inner_reads(input_bam, svtype, chrom, left_pos, left_conf_int, 
                                        right_pos, right_conf_int, INNER_READ_FETCH_FLANK)

        # Uses coverage and aligned bases to find coverage depth
        left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, FETCH_FLANK, chrom_length)
        right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, FETCH_FLANK, chrom_length)
        ref_bases_covered_left_possible, ref_bases_covered_right_possible = get_bases_possible(left_lower_bound, left_upper_bound, 
                                                                                                   right_lower_bound, right_upper_bound)
        ref_bases_covered_left_actual, ref_bases_covered_right_actual = 0, 0
        ref_left_read_count, ref_right_read_count = 0, 0

        # Save reads that span breakpoint
        reads_span_breakpoint = set()

        # These numbers are used to genotype
        # TODO: score struct, named tuple
        ref_support = 0
        alt_support = 0
        
        # TODO: Instead of ifs, get probabilities of support 
        # Extract support from individual reads
        for read_list in reads.values():
            for read in read_list:
                # TODO: Refactor into named tuples so funcs can take these in
                #options = {} # cmd options
                #params = {} # svnumbers
                #generic_feature_extractor(read, options, params)
                if spans_breakpoint(read, left_pos, left_conf_int, MIN_ALIGNED, MIN_PCT_ALIGNED):
                    ref_support += 0.5
                    reads_span_breakpoint.add(read)

                if spans_breakpoint(read, right_pos, right_conf_int, MIN_ALIGNED, MIN_PCT_ALIGNED):
                    ref_support += 0.5
                    reads_span_breakpoint.add(read)
                
                if fa is not None:
                    if spans_breakpoint_mismatches (fa, input_bam, read, left_pos, left_conf_int, MIN_ALIGNED, MIN_PCT_ALIGNED):
                        ref_support += 0
                    if spans_breakpoint_mismatches (fa, input_bam, read, right_pos, right_conf_int, MIN_ALIGNED, MIN_PCT_ALIGNED):
                        ref_support += 0

                if clipped_by_breakpoint(read, left_pos, left_conf_int, right_pos, right_conf_int, SPLIT_SLOP, EXTRA_CLIP_SLOP):
                    alt_support += 0

                #if is_mismatched_over_sv(read, left_pos, right_pos, MISMATCH_SLOP, MISMATCH_PCT):
                if is_mismatched_over_sv_ver2(read, left_pos, left_conf_int, right_pos, right_conf_int, MISMATCH_PCT):
                    alt_support += 0

                #coverage_diff = get_coverage_diff_skip(read, left_pos, right_pos, chrom_length, READ_DEPTH_SKIP, READ_DEPTH_INTERVAL)
                coverage_diff = get_coverage_diff(read, left_pos, left_conf_int, right_pos, right_conf_int, chrom_length, READ_DEPTH_INTERVAL)
                if coverage_diff >= 0.8:
                    ref_support += 0 # min(coverage_diff - 1, 2)
                elif coverage_diff < 0.1:
                    alt_support += 0 # (1 - coverage_diff) * 0.20

                # Find the number of bases covered by this read on the left breakpoint
                read_ref_bases_covered_left_actual = get_num_bases_covered(read, left_lower_bound, left_upper_bound)
                if read_ref_bases_covered_left_actual > 0:
                    ref_left_read_count += 1
                    ref_bases_covered_left_actual += read_ref_bases_covered_left_actual

                # Find the number of bases covered by this read on the right breakpoint
                read_ref_bases_covered_right_actual = get_num_bases_covered(read, right_lower_bound, right_upper_bound)
                if read_ref_bases_covered_right_actual > 0:
                    ref_right_read_count += 1
                    ref_bases_covered_right_actual += read_ref_bases_covered_right_actual

        # If there are splitters that have pieces near both breakpoints, they're in the same read_list
        for read_list in reads.values():
            if split_by_breakpoint(read_list, left_pos, left_conf_int, right_pos, right_conf_int, SPLIT_SLOP):
                alt_support += 1
                
                # If a read is classified as a splitter, then we don't want to give support to the reference for it
                for read in read_list:
                    if read in reads_span_breakpoint:
                        ref_support -= 0.0 # ignoring this does best on training set, 0.5 fives alt_het calls, but adds more het_alt calls

        # Find the avg coverage per base for the region between the breakpoints
        inner_lower_bound, inner_upper_bound = get_inner_query_bounds(svtype, left_pos, left_conf_int, 
                                                                      right_pos, right_conf_int, INNER_READ_FETCH_FLANK)
        alt_bases_covered_possible = get_bases_possible_one_side(inner_lower_bound,inner_upper_bound)
        alt_bases_covered_actual = 0
        alt_read_count = 0

        for read_list in inner_reads.values():
            for read in read_list:
                alt_bases_covered_actual += get_num_bases_covered(read, inner_lower_bound, inner_upper_bound)
                alt_read_count += 1

        ref_left_coverage_pct, ref_left_avg_base_coverage = get_cov_pct_avg_base(ref_bases_covered_left_actual, ref_bases_covered_left_possible, 
                                                                                      ref_left_read_count)
        ref_right_coverage_pct, ref_right_avg_base_coverage = get_cov_pct_avg_base(ref_bases_covered_right_actual, ref_bases_covered_right_possible, 
                                                                                      ref_right_read_count)
        ref_avg_base_coverage = (ref_left_avg_base_coverage + ref_right_avg_base_coverage) / 2

        alt_coverage_pct, alt_avg_base_coverage = get_cov_pct_avg_base(alt_bases_covered_actual, alt_bases_covered_possible, alt_read_count)

        # Genotype based on support for alternate and reference
        genotype = get_genotype(alt_support, ref_support, MIN_PCT_HET)

        # Save the genotype to this variant
        variant.genotypes = cyvcf_genotype(genotype)

        # Write this variant to the output VCF
        output_vcf.write_record(variant)
    output_vcf.close()
    input_vcf.close()

        #write_genotype(
        # Set up output variant info
        #output_variant = Variant()
        #output_variant.chrom = variant.CHROM
        #output_variant.pos = variant.start
        #output_variant.id = variant.ID
        #output_variant.ref = variant.REF
        #output_variant.alt = variant.ALT[0]
        #output_variant.qual = replace_none(variant.QUAL)
        #output_variant.filter = replace_empty([variant.FILTER])
        #output_variant.info = {}
        #for k in dict(variant.INFO).keys():
        #    output_variant.info[k] = variant.INFO[k]

        #output_variant.chrom = variant2.chrom
        #output_variant.pos = variant2.pos
        #output_variant.id = variant2.id
        #output_variant.ref = variant2.ref
        #output_variant.alt = variant2.alts[0]
        #output_variant.qual = replace_none(variant2.qual)
        #output_variant.filter = replace_empty(variant2.filter)
        #for k in variant2.info.keys():
        #    output_variant.info[k] = variant2.info[k]

        # Write this variant and its genotype
        #write_variant(output_vcf_str, output_variant, genotype)

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
        input_fasta_str = args.fasta
        genotype_variants(input_vcf_str, input_bam_str, output_vcf_str, input_fasta_str)
