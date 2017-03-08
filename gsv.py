"""
TODO: How to test funcs that depend on stuff, test directory with sample stuff?
TODO: VCF output print for real and file match up and take care of fields with values like True (IMPRECISE)
"""

import pdb # TODO: debug

import pysam, argparse, re

class Variant:
    """

    """
    def __init__ (self):
        """

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

# adapted from Matt Shirley (http://coderscrowd.com/app/public/codes/view/171)
def cigarstring_to_tuple (cigarstring):
    """

    """
    cigar_dict = {'M':0, 'I':1,'D':2,'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
    pattern = re.compile('([MIDNSHPX=])')
    values = pattern.split(cigarstring)[:-1] ## turn cigar into tuple of values
    paired = (values[n:n+2] for n in range(0, len(values), 2)) ## pair values by twos
    return [(cigar_dict[pair[1]], int(pair[0])) for pair in paired]


# TODO: check chromosome length to make sure query bounds don't go too far if necessary, test it to find out...
def get_query_bounds (svtype, position, conf_int, fetch_flank, chrom_length):
    """
    Returns the bounds to query the BAM file around the given position/confidence 
    interval for the specified svtype and fetch_flank.
    >>> get_query_bounds('DEL', 10, (-4, 10), 1)
    (5, 22)
    >>> get_query_bounds('DEL', 100, (-10, 0), 20)
    (70, 121)
    """
    # if svtype == util.SVTYPE_DEL:
    return max(position + conf_int[0] - fetch_flank, 0), min(position + conf_int[1] + fetch_flank + 1, chrom_length)

def get_section_bounds (svtype, pos_left, conf_int_left, chrom_length_left, pos_right, conf_int_right, chrom_length_right, fetch_flank, min_aligned):
    query_bounds_left = max(pos_left + conf_int_left[0] - fetch_flank, 0), min(pos_left + conf_int_left[1] + fetch_flank + 1, chrom_length_left)
    query_bounds_right = max(pos_right + conf_int_right[0] - fetch_flank, 0), min(pos_right + conf_int_right[1] + fetch_flank + 1, chrom_length_right)

    reference_bounds_left = max(0, pos_left - min_aligned), pos_left + min_aligned
    reference_bounds_right = max(0, pos_right - min_aligned), pos_right + min_aligned

    splitter_bounds_left = pos_left - split_slop, pos_left + split_slop
    splitter_bounds_right = pos_right - split_slop, pos_right + split_slop

    return query_bounds_left, query_bounds_right, reference_bounds_left, reference_bounds_right, splitter_bounds_left, splitter_bounds_right

# TODO: Consider conf_int
def spans_breakpoint (read, pos, min_aligned, min_pct_aligned):
    """

    """
    lower_bound, upper_bound = max(0, pos - min_aligned), pos + min_aligned
    return read.get_overlap(lower_bound, upper_bound) >= 2 * min_aligned * min_pct_aligned

# TODO: Consider conf_int
def split_by_breakpoint (read_list, left_pos, right_pos, split_slop):
    """

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

def clipped_by_breakpoint (read, left_pos, right_pos, split_slop):
    """

    """
    right_side_clipped = read.reference_end >= left_pos - split_slop and read.reference_end <= left_pos + split_slop
    left_side_clipped = read.reference_start >= right_pos - split_slop and read.reference_start <= right_pos + split_slop
    if right_side_clipped:
        cigar = read.cigartuples
        if cigar[-1][0] == 4 or cigar[-1][0] == 5:
            right_clip_size = cigar[-1][1]
            if read.reference_end + right_clip_size >= left_pos + split_slop:
                return True

    if left_side_clipped:
        cigar = read.cigartuples
        if cigar[0][0] == 4 or cigar[0][0] == 5:
            left_clip_size = cigar[0][1]
            if read.reference_start - left_clip_size <= right_pos - split_slop:
                return True

    return False

def get_coverage_diff_skip (read, left_pos, right_pos, chrom_length, read_depth_skip, read_depth_interval):
    """

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

def is_mismatched_over_sv (read, left_pos, right_pos, split_slop, mismatched_pct):
    """

    """
    from_left = read.reference_end >= left_pos + split_slop 
    from_right = read.reference_start <= right_pos - split_slop
    if not from_left and not from_right:
        return False
    if from_left:
        length_over_sv = read.reference_end - (left_pos + split_slop)
        length_of_sv = (right_pos - split_slop) - (left_pos + split_slop)
        min_length = max(min(length_over_sv, length_of_sv), 0)
        #if min_length == 0:
        #    return False
        matches = read.get_overlap(left_pos + split_slop, left_pos + split_slop + min_length)
        mismatches = min_length - matches
        if mismatches / min_length > mismatched_pct:
            

def replace_none (field):
    """

    """
    if field == None:
        return '.'
    else:
        return field

def replace_empty (_list):
    """

    """
    if len(_list) == 0:
        return '.'
    else:
        return _list

def genotype_variants (input_vcf_str, input_bam_str, output_vcf_str):
    """
    Genotype the variants specified in the VCF file using data from the BAM file.
    """
    header = ''
    with open(input_vcf_str, 'r') as f:
        for line in f:
            if line[:2] == '##':
                header += line
            else:
                break
    header += '#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	NA12878\n'
    
    if output_vcf_str is not None:
        with open(output_vcf_str, 'w') as f:
            f.write(header)
    else:
        print(header, end='')

    input_vcf = pysam.VariantFile(input_vcf_str)
    # TODO: multiple input bams/samples
    input_bam = pysam.AlignmentFile(input_bam_str, 'rb')

    # Loop over all variants to get svtype, pos, conf_int from VCF file
    for variant in input_vcf.fetch():
        chrom = variant.chrom
        svtype = variant.info['SVTYPE']
        left_pos = variant.pos
        right_pos = variant.info['END']
        # TODO: Use the 95% conf int instead to limit the size of this thing?
        left_conf_int = variant.info['CIPOS']
        right_conf_int = variant.info['CIEND']
        # TODO: left and right on different chromosomes?
        chrom_length = input_bam.lengths[input_bam.gettid(chrom)]

        fetch_flank = 200
        min_aligned = 200
        min_pct_aligned = 0.80
        split_slop = 250 # 150 gets worse genotypes correct
        read_depth_skip = 50
        read_depth_interval = 50

        # Added this per Brent's request but I'll need to do some more refactoring to make it work well, also as long as I'm prototyping it is a pain to maintain...
        #left_query_bounds, right_query_bounds, left_reference_bounds, right_reference_bounds, left_splitter_bounds, right_splitter_bounds = get_section_bounds (svtype, left_pos, left_conf_int, chrom_length, right_pos, right_conf_int, chrom_length, fetch_flank, min_aligned)
        #left_lower_bound, left_upper_bound = left_query_bounds
        #right_lower_bound, right_upper_bound = right_query_bounds
        #left_lower_ref, left_upper_ref = left_reference_bounds
        #right_lower_ref, right_upper_ref = right_reference_bounds
        #left_lower_split, left_upper_split = left_splittler_bounds
        #right_lower_split, right_upper_split = right_splittler_bounds

        left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, fetch_flank, chrom_length)
        right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, fetch_flank, chrom_length)

        left_reads = {}
        for read in input_bam.fetch(chrom, left_lower_bound, left_upper_bound):
            if read.is_unmapped or read.is_duplicate:
                continue
            if read.query_name in left_reads:
                left_reads[read.query_name].append(read)
            else:
                left_reads[read.query_name] = []
                left_reads[read.query_name].append(read)

        right_reads = {}
        for read in input_bam.fetch(chrom, right_lower_bound, right_upper_bound):
            if read.is_unmapped or read.is_duplicate:
                continue
            if read.query_name in right_reads:
                right_reads[read.query_name].append(read)
            else:
                right_reads[read.query_name] = []
                right_reads[read.query_name].append(read)

        reads = {}
        for read_list in left_reads.values():
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
        for read_list in right_reads.values():
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

        ref_support = 0
        alt_support = 0
        # Instead of ifs, get probabilities of support 
        for read_list in reads.values():
            for read in read_list:
                #read_count += 1
                if spans_breakpoint(read, left_pos, min_aligned, min_pct_aligned):
                    ref_support += 0.5
                if spans_breakpoint(read, right_pos, min_aligned, min_pct_aligned):
                    ref_support += 0.5
                if clipped_by_breakpoint(read, left_pos, right_pos, split_slop):
                    alt_support += 0.01
                # TODO: This currently isn't helping us do better but probably has the potential to do so
                #coverage_diff = get_coverage_diff_skip(read, left_pos, right_pos, chrom_length, read_depth_skip, read_depth_interval)
                #coverage_diff = get_coverage_diff(read, left_pos, right_pos, chrom_length, read_depth_interval)
                #if coverage_diff >= 0.8:
                #    ref_support += min(coverage_diff - 1, 2)
                #elif coverage_diff < 0.1:
                #    alt_support += (1 - coverage_diff) * 0.20
        #if read_count > 0:
        #    coverage_support = float(coverage_diff) / float(read_count)
        #    if coverage_support > 1.5:
        #        ref_support += read_count * 0.1
        #    elif coverage_support < 0.5:
        #        alt_support += read_count * 0.1

        for read_list in reads.values():
            if split_by_breakpoint(read_list, left_pos, right_pos, split_slop):
                alt_support += 1

        total_support = alt_support + ref_support
        if total_support > 0:
            min_pct_het = 0.15
            is_het = False
            if ref_support > total_support * min_pct_het and alt_support > total_support * min_pct_het:
                is_het = True
            if is_het:
                genotype = '0/1'
            elif alt_support > ref_support:
                genotype = '1/1'
            else:
                genotype = '0/0'
        else:
            genotype = '1/1'

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


def is_splitter (read):
    """

    """
    return read.has_tag('SA')

def enters_ref_region (read, pos, split_slop):
    """

    """
    reference_region = pos + split_slop
    return read.reference_end > reference_region

def begins_before_breakpoint (read, pos):
    """

    """
    return read.reference_start < pos

if __name__ == '__main__':
    """

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
