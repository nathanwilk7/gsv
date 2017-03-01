"""
TODO: How to test funcs that depend on stuff, test directory with sample stuff?
TODO: VCF output pretty print and file match up and take care of fields with values like True (IMPRECISE)
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
def cigarstring_to_tuple(cigarstring):
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

# TODO: Consider conf_int
def spans_breakpoint (read, chrom, pos, conf_int, min_aligned, min_pct_aligned):
    """

    """
    #pdb.set_trace()
    return read.get_overlap(max(0, pos - min_aligned), pos + min_aligned) >= 2 * min_aligned * min_pct_aligned

# TODO: Consider conf_int
def split_by_breakpoint (read, chrom, left_pos, left_conf_int, right_pos, right_conf_int, split_slop):
    """

    """
    #pdb.set_trace()
    if not read.has_tag('SA'):
        return False
    
    split_left = read.reference_end >= left_pos - split_slop and read.reference_end <= left_pos + split_slop
    split_right = read.reference_start >= right_pos - split_slop and read.reference_start <= right_pos + split_slop
    if not split_left and not split_right:
        return False
    
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

    return False

def clipped_by_breakpoint (read, chrom, left_pos, left_conf_int, right_pos, right_conf_int, split_slop):
    """

    """
    

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
    #pdb.set_trace()
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

        fetch_flank = 200

        # Get range from which to get BAM reads
        chrom_length = input_bam.lengths[input_bam.gettid(chrom)]
        left_lower_bound, left_upper_bound = get_query_bounds(svtype, left_pos, left_conf_int, fetch_flank, chrom_length)
        right_lower_bound, right_upper_bound = get_query_bounds(svtype, right_pos, right_conf_int, fetch_flank, chrom_length)

        # Get the reads
        # TODO: What does it mean when 2 reads have the same query_name?
        reads = {}
        left_reads = {}
        for read in input_bam.fetch(chrom, left_lower_bound, left_upper_bound):
            if read.is_unmapped or read.is_duplicate:
                continue
            if read.query_name in left_reads:
                left_reads[read.query_name].append(read)
            else:
                left_reads[read.query_name] = []
                left_reads[read.query_name].append(read)
            if read.query_name in reads:
                reads[read.query_name].append(read)
            else:
                reads[read.query_name] = []
                reads[read.query_name].append(read)

        right_reads = {}
        for read in input_bam.fetch(chrom, right_lower_bound, right_upper_bound):
            if read.is_unmapped or read.is_duplicate:
                continue
            if read.query_name in right_reads:
                right_reads[read.query_name].append(read)
            else:
                right_reads[read.query_name] = []
                right_reads[read.query_name].append(read)
            if read.query_name in reads:
                reads[read.query_name].append(read)
            else:
                reads[read.query_name] = []
                reads[read.query_name].append(read)

        min_aligned = 200
        min_pct_aligned = 0.80
        split_slop = 150
        total_count = 0
        ref_count = 0
        alt_count = 0
        for read_list in reads.values():
            for read in read_list:
                total_count += 1
                if spans_breakpoint(read, chrom, left_pos, left_conf_int, min_aligned, min_pct_aligned):
                    ref_count += 0.5
                if spans_breakpoint(read, chrom, right_pos, right_conf_int, min_aligned, min_pct_aligned):
                    ref_count += 0.5
                if split_by_breakpoint(read, chrom, left_pos, left_conf_int, right_pos, right_conf_int, split_slop):
                    alt_count += 0.5

        if total_count > 0:
            prior_ref = .15
            prior_alt = 1
            pct_ref = (float(ref_count) / float(total_count)) * prior_ref
            pct_alt = (float(alt_count) / float(total_count)) * prior_alt
            abs_diff = abs(pct_ref - pct_alt)
            #pdb.set_trace()
            if abs_diff < 0.05:
                genotype = '0/1'
            elif pct_alt > pct_ref:
                genotype = '1/1'
            else:
                genotype = '0/0'
        else:
            genotype = '1/1'
        #if total_count > 0:
        #    prob_ref = float(ref_count) / float(total_count)
        #    prob_alt = float(alt_count) / float(total_count)
        #else:
        #    prob_ref = 0
        #    prob_alt = 1

        #prob_hom_ref = prob_ref * (1 - prob_alt) * 0.01
        #prob_hom_alt = prob_alt * (1 - prob_ref)
        #if prob_hom_alt >= prob_hom_ref:
        #    genotype = '1/1'
        #else:
        #    genotype = '0/0'
        #prob_het = ((prob_hom_ref + prob_hom_alt) / 2)
        #if prob_hom_alt >= prob_het and prob_hom_alt >= prob_hom_ref:
        #    genotype = '1/1'
        #elif prob_het >= prob_hom_alt and prob_het >= prob_hom_ref:
            #genotype = '0/1'
        #else:
        #    genotype = '0/0'

        #if prob_ref > 0.50:
        #    genotype = '0/0'
        #elif prob_ref > 0.25:
        #    genotype = '0/1'
        #else:
        #    genotype = '1/1'

        #left_ref_count = 0
        #left_alt_count = 0
        #for read_list in left_reads.values():
        #    for read in read_list:
        #        if enters_ref_region(read, left_pos, split_slop) and begins_before_breakpoint(read, left_pos):
        #            left_ref_count += 1
        #        else:
        #            left_alt_count += 1
        
        #left_total_read_count = left_ref_count + left_alt_count
        #left_pct_ref = float(left_ref_count) / left_total_count
        #left_pct_alt = float(left_alt_count) / left_total_count

        #if left_pct_ref > .75:
        #    genotype = '0/0'
        #elif left_pct_ref > .25:
        #    genotype = '0/1'
        #else:
        #    genotype = '1/1'

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
            #pdb.set_trace()
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
