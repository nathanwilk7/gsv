import pdb # TODO: debug
import pysam, argparse, re
import cyvcf2
from collections import namedtuple
import math

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
    parser.add_argument('-c', '--csv', help='Path to output CSV file of features')
    return parser.parse_args()

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

def fetch_reads (input_bam, sv_info, options): #svtype, chrom, left_pos, left_conf_int, right_pos, right_conf_int, fetch_flank, chrom_length):
    """
    Fetches all the reads from the bam at the given location using the parameters specified and the get_query_bounds function.
    """
    sv_info.left_lower_bound, sv_info.left_upper_bound = get_query_bounds(sv_info.svtype, sv_info.left_pos, sv_info.left_conf_int, 
                                                                          options.fetch_flank, sv_info.chrom_length)
    left_reads = fetch_one_locus_reads(input_bam, sv_info.chrom, sv_info.left_lower_bound, sv_info.left_upper_bound)

    sv_info.right_lower_bound, sv_info.right_upper_bound = get_query_bounds(sv_info.svtype, sv_info.right_pos, sv_info.right_conf_int, 
                                                                            options.fetch_flank, sv_info.chrom_length)
    right_reads = fetch_one_locus_reads(input_bam, sv_info.chrom, sv_info.right_lower_bound, sv_info.right_upper_bound)

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

def fetch_inner_reads (input_bam, sv_info, options): #svtype, chrom, left_pos, left_conf_int, right_pos, right_conf_int, inner_read_fetch_flank):
    """
    Fetches all the reads from the bam between the left and right sides
    """
    # Chromosome length shouldn't matter since we're looking in between two breakpoints, so just set it something irrelevant
    sv_info.inner_lower_bound, sv_info.inner_upper_bound = get_inner_query_bounds(sv_info.svtype, sv_info.left_pos, sv_info.left_conf_int, sv_info.right_pos, 
                                                                 sv_info.right_conf_int, options.inner_read_fetch_flank)
    inner_reads = fetch_one_locus_reads(input_bam, sv_info.chrom, sv_info.inner_lower_bound, sv_info.inner_upper_bound)
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

def spans_breakpoint_left (variant_features, read, sv_info, options):
    """
    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = (max(0, sv_info.left_pos + sv_info.left_conf_int[0] - options.min_aligned), 
                                sv_info.left_pos + options.min_aligned + sv_info.left_conf_int[1])
    conf_int_length = abs(sv_info.left_conf_int[0] - sv_info.left_conf_int[1])
    does_span_left_breakpoint = read.get_overlap(lower_bound, upper_bound) >= (conf_int_length + (2 * options.min_aligned)) * options.min_pct_aligned
    if does_span_left_breakpoint:
        variant_features.left_breakpoint_spanners += 1

def spans_breakpoint_right (variant_features, read, sv_info, options):
    """    Returns True if the read as at least the minimum percent aligned with the reference around the position 
     +/- the minimum number of aligned bases on each side. Otherwise, returns False.
    """
    lower_bound, upper_bound = (max(0, sv_info.right_pos + sv_info.right_conf_int[0] - options.min_aligned), 
                                sv_info.right_pos + options.min_aligned + sv_info.right_conf_int[1])
    conf_int_length = abs(sv_info.right_conf_int[0] - sv_info.right_conf_int[1])
    does_span_right_breakpoint = read.get_overlap(lower_bound, upper_bound) >= (conf_int_length + (2 * options.min_aligned)) * options.min_pct_aligned
    if does_span_right_breakpoint:
        variant_features.right_breakpoint_spanners += 1

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

def split_by_breakpoint (variant_features, read_list, sv_info, options):
    """
    Returns True if a split happens on the left or side within the split slop of the corresponding 
     position and if the other split aligns to the other side of the variant. Otherwise return False.
     Assumes that the read_list contains all split pieces that should be considered.
    """
    for read in read_list:
        split_left = (read.reference_end >= sv_info.left_pos + sv_info.left_conf_int[0] - options.split_slop and 
                      read.reference_end <= sv_info.left_pos + sv_info.left_conf_int[1] + options.split_slop)
        split_right = (read.reference_start >= sv_info.right_pos + sv_info.right_conf_int[0] - options.split_slop and 
                       read.reference_start <= sv_info.right_pos + sv_info.right_conf_int[1] + options.split_slop)
        if split_left or split_right:
            for read2 in read_list:
                if read == read2:
                    continue
                if split_left:
                    if (read2.reference_start >= sv_info.right_pos + sv_info.right_conf_int[0] - options.split_slop and 
                        read2.reference_start <= sv_info.right_pos + sv_info.right_conf_int[1] + options.split_slop):
                        variant_features.left_splitters += 1
                        return

                if split_right:
                    if (read2.reference_end >= sv_info.left_pos + sv_info.left_conf_int[0] - options.split_slop and 
                        read2.reference_start <= sv_info.left_pos + sv_info.left_conf_int[1] + options.split_slop):
                        variant_features.right_splitters += 1
                        return

def split_by_breakpoint_old2 (read_list, left_pos, left_conf_int, right_pos, right_conf_int, split_slop):
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

# Save this because i may need to go back to SA tags in the end...

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

def clipped_by_breakpoint (variant_features, read, sv_info, options):
    """
    Returns True if the read was clipped on it's right or left side by the variant. Otherwise False.
    """
    right_side_clipped = (read.reference_end >= sv_info.left_pos + sv_info.left_conf_int[0] - options.split_slop and 
                          read.reference_end <= sv_info.left_pos + sv_info.left_conf_int[1] + options.split_slop)
    left_side_clipped = (read.reference_start >= sv_info.right_pos + sv_info.right_conf_int[0] - options.split_slop and 
                         read.reference_start <= sv_info.right_pos + sv_info.right_conf_int[1] + options.split_slop)
    if right_side_clipped:
        cigar = read.cigartuples
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            right_clip_size = cigar[-1][1]
            if read.reference_end + right_clip_size >= sv_info.left_pos + sv_info.left_conf_int[1] + options.split_slop + options.extra_clip_slop:
                variant_features.right_clippers += 1
                return

    if left_side_clipped:
        cigar = read.cigartuples
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            left_clip_size = cigar[0][1]
            if read.reference_start - left_clip_size <= sv_info.right_pos + sv_info.right_conf_int[0] - options.split_slop - options.extra_clip_slop:
                variant_features.left_clippers += 1
                return

def clipped_by_breakpoint_old (read, left_pos, left_conf_int, right_pos, right_conf_int, split_slop, extra_clip_slop):
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

def get_coverage_diff (variant_features, read, sv_info, options):
    """
    Returns the ratio of the overlap with the reference on the inside of the variant versus the outside of the variant.
     This version starts counting overlap from the positions.
    """
    left_outside = max(sv_info.left_pos + sv_info.left_conf_int[0] - options.read_depth_interval, 0)
    left_inside = min(sv_info.left_pos + sv_info.left_conf_int[1] + options.read_depth_interval, sv_info.chrom_length)
    left_outside_overlap = read.get_overlap(left_outside, sv_info.left_pos)
    left_inside_overlap = read.get_overlap(sv_info.left_pos, left_inside)

    right_inside = max(sv_info.right_pos + sv_info.right_conf_int[0] - options.read_depth_interval, 0)
    right_outside = min(sv_info.right_pos + sv_info.right_conf_int[1] + options.read_depth_interval, sv_info.chrom_length)
    right_inside_overlap = read.get_overlap(right_inside, sv_info.right_pos)
    right_outside_overlap = read.get_overlap(sv_info.right_pos, right_outside)
    
    outside_overlap = left_outside_overlap + right_outside_overlap
    inside_overlap = left_inside_overlap + right_inside_overlap

    if outside_overlap <= 0:
        variant_features.coverage_diff.append(1.0)
    else:
        variant_features.coverage_diff.append(float(inside_overlap) / float(outside_overlap))

def get_coverage_diff_old (read, left_pos, left_conf_int, right_pos, right_conf_int, chrom_length, read_depth_interval):
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

def is_mismatched_over_sv_ver2 (variant_features, read, sv_info, options):
    """
    This function returns True if the read extends into the variant and has a mismatch 
     percentage higher than the specified one.
    """
    enters_sv_from_left = read.reference_start <= sv_info.left_pos + sv_info.left_conf_int[1]
    enters_sv_from_right = read.reference_start > sv_info.left_pos + sv_info.left_conf_int[1]

    if enters_sv_from_left:
        length_over_sv = read.reference_end - sv_info.left_pos + sv_info.left_conf_int[1]
        length_of_sv = (sv_info.right_pos + sv_info.right_conf_int[0]) - (sv_info.left_pos + sv_info.left_conf_int[1])
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return
        matches = read.get_overlap(sv_info.left_pos + sv_info.left_conf_int[1], sv_info.left_pos + sv_info.left_conf_int[1] + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > options.mismatch_pct:
            variant_features.mismatched_over_sv += 1
            return
    
    if enters_sv_from_right:
        length_over_sv = sv_info.right_pos + sv_info.right_conf_int[0] - read.reference_start
        length_of_sv = (sv_info.right_pos + sv_info.right_conf_int[0]) - (sv_info.left_pos + sv_info.left_conf_int[1])
        min_length = max(min(length_over_sv, length_of_sv), 0)
        if min_length == 0:
            return
        matches = read.get_overlap(sv_info.right_pos + sv_info.right_conf_int[0] - min_length, sv_info.right_pos + sv_info.right_conf_int[0])
        mismatches = min_length - matches
        if mismatches / min_length > options.mismatch_pct:
            variant_features.mismatched_over_sv += 1
            return

def is_mismatched_over_sv_ver2_old2 (read, left_pos, left_conf_int, right_pos, right_conf_int, mismatch_pct):
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

def get_left_bounds_for_coverage (sv_info, options):
    return (max(0, sv_info.left_pos + sv_info.left_conf_int[0] - options.outer_read_fetch_flank), 
            sv_info.left_pos + sv_info.left_conf_int[0])

def get_right_bounds_for_coverage (sv_info, options):
    return (max(0, sv_info.right_pos + sv_info.right_conf_int[1]), 
            sv_info.right_pos + sv_info.right_conf_int[1] + options.outer_read_fetch_flank)

def get_bases_possible_one_side (lower_bound, upper_bound):
    return upper_bound - lower_bound

def get_bases_possible (variant_features, sv_info, options):
    left_lower_bound_possible, left_upper_bound_possible = get_left_bounds_for_coverage(sv_info, options)
    right_lower_bound_possible, right_upper_bound_possible = get_right_bounds_for_coverage(sv_info, options)
    variant_features.ref_bases_covered_left_possible = get_bases_possible_one_side(left_lower_bound_possible, left_upper_bound_possible)
    variant_features.ref_bases_covered_right_possible = get_bases_possible_one_side(right_lower_bound_possible, right_upper_bound_possible)

def get_inner_bases_possible (variant_features, sv_info):
    variant_features.alt_bases_covered_possible = get_bases_possible_one_side(sv_info.inner_lower_bound, sv_info.inner_upper_bound)

def get_bases_possible_old (left_lower_bound, left_upper_bound, right_lower_bound, right_upper_bound):
    return get_bases_possible_one_side(left_lower_bound, left_upper_bound), get_bases_possible_one_side(right_lower_bound, right_upper_bound) 
    
def get_num_bases_covered_old (variant_features, read, sv_info, options):
    variant_features.ref_bases_covered_left_actual += read.get_overlap(sv_info.left_lower_bound, sv_info.left_upper_bound)
    variant_features.ref_bases_covered_right_actual += read.get_overlap(sv_info.right_lower_bound, sv_info.right_upper_bound)

def get_num_bases_covered (variant_features, read, sv_info, options):
    ref_bases_left_lower_bound, ref_bases_left_upper_bound = get_left_bounds_for_coverage(sv_info, options)
    ref_bases_right_lower_bound, ref_bases_right_upper_bound = get_right_bounds_for_coverage(sv_info, options)
    variant_features.ref_bases_covered_left_actual += read.get_overlap(ref_bases_left_lower_bound, ref_bases_left_upper_bound)
    variant_features.ref_bases_covered_right_actual += read.get_overlap(ref_bases_right_lower_bound, ref_bases_right_upper_bound)

def get_inner_num_bases_covered (variant_features, inner_read, sv_info, options):
    variant_features.alt_bases_covered_actual += inner_read.get_overlap(sv_info.inner_lower_bound, sv_info.inner_upper_bound)

def get_num_bases_covered_old (read, lower_bound, upper_bound):
    return read.get_overlap(lower_bound, upper_bound)

def get_cov_pct_avg_base (bases_covered_actual, bases_covered_possible, read_count):
    if bases_covered_possible > 0 and read_count > 0:
        coverage_pct = bases_covered_actual / (bases_covered_possible * read_count)
        avg_base_coverage = bases_covered_actual / bases_covered_possible
    else:
        coverage_pct = 0
        avg_base_coverage = 0
    return coverage_pct, avg_base_coverage

def csvify (sv_info, variant_features, pred_genotype):
    csv_str = str(sv_info.chrom) + ','
    csv_str += str(sv_info.svtype) + ','
    csv_str += str(sv_info.left_pos) + ','
    csv_str += str(sv_info.right_pos) + ','
    #csv_str += str(sv_info.left_conf_int) + ','
    #csv_str += str(sv_info.right_conf_int) + ','
    csv_str += str(variant_features.left_breakpoint_spanners) + ','
    csv_str += str(variant_features.right_breakpoint_spanners) + ','
    csv_str += str(variant_features.avg_coverage_diff) + ','
    csv_str += str(variant_features.right_clippers) + ','
    csv_str += str(variant_features.left_clippers) + ','
    csv_str += str(variant_features.mismatched_over_sv) + ','
    csv_str += str(variant_features.ref_bases_covered_left_possible) + ','
    csv_str += str(variant_features.ref_bases_covered_left_actual) + ','
    csv_str += str(variant_features.ref_bases_covered_right_possible) + ','
    csv_str += str(variant_features.ref_bases_covered_right_actual) + ','
    csv_str += str(variant_features.ref_avg_base_coverage) + ','
    csv_str += str(variant_features.ref_coverage_pct) + ','
    csv_str += str(variant_features.alt_bases_covered_possible) + ','
    csv_str += str(variant_features.alt_bases_covered_actual) + ','
    csv_str += str(variant_features.alt_avg_base_coverage) + ','
    csv_str += str(variant_features.alt_coverage_pct) + ','
    csv_str += str(variant_features.left_splitters) + ','
    csv_str += str(variant_features.right_splitters) + ','
    csv_str += str(variant_features.splitter_spanners) + ','
    csv_str += str(pred_genotype) + '\n'
    return csv_str

def get_ref_het_alt_support (variant_features):
    """
    Calculate alt and ref support from variant features
    """
    ref_support, het_support, alt_support = 0, 0, 0
    
    if variant_features.left_breakpoint_spanners > 18:
        ref_support += 1
    elif variant_features.left_breakpoint_spanners < 2:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.right_breakpoint_spanners > 16:
        ref_support += 1
    elif variant_features.right_breakpoint_spanners < 4:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.avg_coverage_diff > .7:
        ref_support += 1
    elif variant_features.avg_coverage_diff < .3:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.right_clippers < 2:
        ref_support += 1
    elif variant_features.right_clippers > 16:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.left_clippers < 5:
        ref_support += 1
    elif variant_features.left_clippers > 17:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.mismatched_over_sv < 2:
        ref_support += 1
    elif variant_features.mismatched_over_sv > 20:
        alt_support += 1
    else:
        het_support += 1
    
    if variant_features.alt_avg_base_coverage > 16:
        ref_support += 1
    elif variant_features.alt_avg_base_coverage < 6:
        alt_support += 1
    else:
        het_support += 1
    
    if variant_features.alt_coverage_pct > .75:
        ref_support += 1
    elif variant_features.alt_coverage_pct < .2:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.left_splitters < 2:
        ref_support += 1
    elif variant_features.left_splitters > 18:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.splitter_spanners > 20:
        ref_support += 1
    elif variant_features.splitter_spanners < 5:
        alt_support += 1
    else:
        het_support += 1
    
    if variant_features.spanners > 35:
        ref_support += 1
    elif variant_features.spanners < 5:
        alt_support += 1
    else:
        het_support += 1

    if variant_features.splitters < 2:
        ref_support += 1
    elif variant_features.splitters > 20:
        alt_support += 1
    else:
        het_support += 1

    return ref_support, het_support, alt_support

def get_genotype (ref_support, het_support, alt_support): #min_pct_het):
    """
    Returns the genotype based on the alternate and reference support.
    """
    max_support = max(ref_support, het_support, alt_support)
    if het_support == max_support:
        return '0/1'
    elif alt_support == max_support:
        return '1/1'
    else:
        return '0/0'
    #total_support = alt_support + ref_support
    #if total_support > 0:
    #    if ref_support > total_support * min_pct_het and alt_support > total_support * min_pct_het:
    #        genotype = '0/1'
    #    elif alt_support > ref_support:
    #        genotype = '1/1'
    #    else:
    #        genotype = '0/0'
    #else:
    #    genotype = '1/1'
    #return genotype

def genotype_eq (G, k, a, e):
    """
    G: 2 for HOM REF, 1 for ALT, 0 for HOM ALT
    k: Total number of reads
    a: Number of reads supporting the alternate
    e: Estimated error at loci
    
    http://biorxiv.org/content/biorxiv/early/2017/02/01/105080.full.pdf
    https://academic.oup.com/bioinformatics/article/27/21/2987/217423/A-statistical-framework-for-SNP-calling-mutation
    """
    return (-k * math.log(2) + (k - a) * math.log((2 - G) * e + G * (1 - e)) +
            a * math.log((2 - G) * (1 - e) + G * e))

def get_genotype_with_prob (num_ref_reads, num_alt_reads):
    """
    Returns the genotype as 0/0 (HOM ALT), 0/1 (HET), or 1/1 (HOM REF) based on most likely genotype.
    """
    k = num_ref_reads + num_alt_reads
    a = num_alt_reads
    e = 0.05
    log_prob_hom_ref = genotype_eq(2, k, a, e)
    log_prob_het = genotype_eq(1, k, a, e)
    log_prob_hom_alt = genotype_eq(0, k, a, e)

    # Call het if they're the same
    if log_prob_het >= log_prob_hom_ref and log_prob_het >= log_prob_hom_alt:
        return '0/1'
    elif log_prob_hom_alt >= log_prob_hom_ref:
        return '1/1'
    else:
        return '0/0'

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

class SVInfo:
    """
    Holds information about this structural variant
    """
    def __init__ (self):
        self.chrom = None
        self.svtype = None
        self.left_pos = None
        self.right_pos = None
        self.left_conf_int = None
        self.right_conf_int = None
        self.chrom_length = None
        self.left_lower_bound = None
        self.left_upper_bound = None
        self.right_lower_bound = None
        self.right_upper_bound = None
        self.inner_lower_bound = None
        self.inner_upper_bound = None

class VariantFeatures:
    """
    When features are extracted from read data, this object's fields keep track of the sum of features for the variant.
    """
    def __init__ (self):
        self.left_breakpoint_spanners = 0.0
        self.right_breakpoint_spanners = 0.0
        self.coverage_diff = []
        self.avg_coverage_diff = 0
        self.right_clippers = 0.0
        self.left_clippers = 0.0
        self.mismatched_over_sv = 0.0
        self.ref_bases_covered_left_possible = 0.0
        self.ref_bases_covered_left_actual = 0.0
        self.ref_bases_covered_right_possible = 0.0
        self.ref_bases_covered_right_actual = 0.0
        self.ref_avg_base_coverage = 0.0
        self.ref_coverage_pct = 0.0
        self.alt_bases_covered_possible = 0.0
        self.alt_bases_covered_actual = 0.0
        self.alt_avg_base_coverage = 0.0
        self.alt_coverage_pct = 0.0
        self.left_splitters = 0.0
        self.right_splitters = 0.0
        self.splitter_spanners = 0.0
        self.spanners = 0.0
        self.splitters = 0.0

def genotype_variants (args):
    """
    Genotype the variants specified in the VCF file using data from the BAM file. optional fasta file can be passed in.
    """
    # Set up the input VCF and input BAM
    input_vcf = cyvcf2.VCF(args.input_vcf)

    # TODO: multiple input bams/samples
    input_bam = pysam.AlignmentFile(args.input_bam, 'rb')

    # Variants will be written to this
    output_vcf = cyvcf2.Writer(args.output_vcf, input_vcf)

    # if a fasta file was provided, then read it in
    if args.fasta is not None:
        fasta = pysam.FastaFile(args.fasta)
    else:
        fasta = None

    if args.csv is not None:
        with open(args.csv, 'w') as f:
            #csv_header = ('chrom,svtype,left_pos,right_pos,left_conf_int,right_conf_int,left_breakpoint_spanners,right_breakpoint_spanners,avg_coverage_diff,'
            #              'right_clippers,left_clippers,mismatched_over_sv,ref_bases_covered_left_possible,ref_bases_covered_left_actual,'
            #              'ref_bases_covered_right_possible,ref_bases_covered_right_actual,ref_avg_base_coverage,ref_coverage_pct,'
            #              'alt_bases_covered_possible,alt_bases_covered_actual,alt_avg_base_coverage,alt_coverage_pct,'
            #              'left_splitters,right_splitters,splitter_spanners,pred_genotype')
            csv_header = ('chrom,svtype,left_pos,right_pos,left_breakpoint_spanners,right_breakpoint_spanners,avg_coverage_diff,'
                          'right_clippers,left_clippers,mismatched_over_sv,ref_bases_covered_left_possible,ref_bases_covered_left_actual,'
                          'ref_bases_covered_right_possible,ref_bases_covered_right_actual,ref_avg_base_coverage,ref_coverage_pct,'
                          'alt_bases_covered_possible,alt_bases_covered_actual,alt_avg_base_coverage,alt_coverage_pct,'
                          'left_splitters,right_splitters,splitter_spanners,pred_genotype')
            f.write(csv_header + '\n')

    Options = namedtuple('Options', ['fetch_flank', 'min_aligned', 'min_pct_aligned', 'split_slop', 'read_depth_skip', 'read_depth_interval', 
                                     'mismatch_slop', 'mismatch_pct', 'extra_clip_slop', 'inner_read_fetch_flank', 'outer_read_fetch_flank'])
    # TODO: set these via cmd args
    options = Options(fetch_flank=200, 
                      min_aligned=200, min_pct_aligned=0.80, 
                      split_slop=300, 
                      read_depth_skip=50, read_depth_interval=50, 
                      mismatch_slop=100, mismatch_pct=0.90, 
                      extra_clip_slop=1000, 
                      inner_read_fetch_flank=20, 
                      outer_read_fetch_flank=20)
        
    # Loop over all variants in input VCF file
    for variant in input_vcf:
        sv_info = SVInfo()
        sv_info.chrom = variant.CHROM
        sv_info.svtype = variant.INFO.get('SVTYPE')
        sv_info.left_pos = variant.start
        sv_info.right_pos = variant.end
        sv_info.left_conf_int = variant.INFO.get('CIPOS')
        sv_info.right_conf_int = variant.INFO.get('CIEND')
        sv_info.chrom_length = input_bam.lengths[input_bam.gettid(variant.CHROM)]

        variant_features = VariantFeatures()

        # Get the reads around the variant on both sides
        reads = fetch_reads(input_bam, sv_info, options)

        # Get the reads in the middle of the variant
        inner_reads = fetch_inner_reads(input_bam, sv_info, options)

        # Save reads that span breakpoint
        reads_span_breakpoint = set()

        # meta info when extracting features from reads
        prev_ref_bases_covered_right_actual = 0.0
        prev_ref_bases_covered_left_actual = 0.0
        ref_left_read_count = 0
        ref_right_read_count = 0
        prev_left_breakpoint_spanners = 0.0
        prev_right_breakpoint_spanners = 0.0
        alt_read_count = 0

        read_funcs = []
        read_funcs.append(spans_breakpoint_left)
        read_funcs.append(spans_breakpoint_right)
        read_funcs.append(get_coverage_diff)
        read_funcs.append(clipped_by_breakpoint)
        read_funcs.append(is_mismatched_over_sv_ver2)
        read_funcs.append(get_num_bases_covered)

        # functions that need the fasta file to work
        fasta_funcs = []

        # Extract support from individual reads
        for read_list in reads.values():
            for read in read_list:
                for read_func in read_funcs:
                    read_func(variant_features, read, sv_info, options)

                # increments read count if there was at least one base of the read was "near" the left/right breakpoint
                if variant_features.ref_bases_covered_left_actual > prev_ref_bases_covered_left_actual:
                    ref_left_read_count += 1
                    prev_ref_bases_covered_left_actual = variant_features.ref_bases_covered_left_actual

                if variant_features.ref_bases_covered_right_actual > prev_ref_bases_covered_right_actual:
                    ref_right_read_count += 1
                    prev_ref_bases_covered_right_actual = variant_features.ref_bases_covered_right_actual

                # keeps track of reads that are considered spanners on at least one side
                if variant_features.left_breakpoint_spanners > prev_left_breakpoint_spanners:
                    reads_span_breakpoint.add(read)
                    prev_left_breakpoint_spanners = variant_features.left_breakpoint_spanners

                if variant_features.right_breakpoint_spanners > prev_right_breakpoint_spanners:
                    reads_span_breakpoint.add(read)
                    prev_right_breakpoint_spanners = variant_features.right_breakpoint_spanners

                if fasta is not None:
                    for fasta_func in fasta_funcs:
                        fasta_func(fasta, variant_features, read, sv_info, options)

        read_list_funcs = []
        read_list_funcs.append(split_by_breakpoint)

        for read_list in reads.values():
            for read_list_func in read_list_funcs:
                read_list_func(variant_features, read_list, sv_info, options)
                for read in read_list:
                    if read in reads_span_breakpoint:
                        variant_features.splitter_spanners += 1
    
        inner_read_funcs = []
        inner_read_funcs.append(get_inner_num_bases_covered)
        for inner_read_list in inner_reads.values():
            for inner_read in inner_read_list:
                alt_read_count += 1
                for inner_read_func in inner_read_funcs:
                    inner_read_func(variant_features, inner_read, sv_info, options)

        # post processing
        if float(len(variant_features.coverage_diff)) != 0:
            variant_features.avg_coverage_diff = sum(variant_features.coverage_diff) / float(len(variant_features.coverage_diff))
        else:
            variant_features.avg_coverage_diff = 1.0
        get_bases_possible(variant_features, sv_info, options)
        get_inner_bases_possible(variant_features, sv_info)

        ref_left_coverage_pct, ref_left_avg_base_coverage = get_cov_pct_avg_base(variant_features.ref_bases_covered_left_actual, 
                                                                                 variant_features.ref_bases_covered_left_possible, ref_left_read_count)
        ref_right_coverage_pct, ref_right_avg_base_coverage = get_cov_pct_avg_base(variant_features.ref_bases_covered_right_actual, 
                                                                                   variant_features.ref_bases_covered_right_possible, ref_right_read_count)
        variant_features.ref_avg_base_coverage = (ref_left_avg_base_coverage + ref_right_avg_base_coverage) / 2.0
        variant_features.ref_coverage_pct = (ref_left_coverage_pct + ref_right_coverage_pct) / 2.0

        variant_features.alt_coverage_pct, variant_features.alt_avg_base_coverage = get_cov_pct_avg_base(variant_features.alt_bases_covered_actual, 
                                                                                                         variant_features.alt_bases_covered_possible, alt_read_count)

        variant_features.spanners = variant_features.left_breakpoint_spanners + variant_features.right_breakpoint_spanners
        variant_features.splitters = variant_features.left_splitters + variant_features.right_splitters

        # get alt/ref support from variant features
        #ref_support, het_support, alt_support = get_ref_het_alt_support(variant_features)

        # Genotype based on support for alternate and reference
        #genotype = get_genotype(ref_support, het_support, alt_support) #MIN_PCT_HET)

        genotype = get_genotype_with_prob(variant_features.spanners, variant_features.splitters)

        # Save the genotype to this variant
        variant.genotypes = cyvcf_genotype(genotype)

        # Write this variant to the output VCF
        output_vcf.write_record(variant)
        
        # write this record to csv
        if args.csv is not None:
            with open(args.csv, 'a') as f:
                f.write(csvify(sv_info, variant_features, genotype))

    output_vcf.close()
    input_vcf.close()

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
        genotype_variants(args)
