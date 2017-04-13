import argparse
import cyvcf2

def convert_cyvcf_gt_to_str (gts):
    if gts[0][:2] == [0,0]:
        return '0/0'
    elif gts[0][:2] == [0,1]:
        return '0/1'
    elif gts[0][:2] == [1,1]:
        return '1/1'
    else:
        return './.'

parser = argparse.ArgumentParser()
parser.add_argument('-t', '--truth_vcf')
parser.add_argument('-c', '--csv')
parser.add_argument('-o', '--output')
args = parser.parse_args()

truth_gts = {}
for variant in cyvcf2.VCF(args.truth_vcf):
    truth_gts[(variant.CHROM, variant.start)] = convert_cyvcf_gt_to_str(variant.genotypes)

is_header = True
with open(args.output, 'w') as w:
    with open(args.csv, 'r') as r:
        for line in r:
            if is_header:
                w.write(line.strip() + ',' + 'true_genotype\n')
                is_header = False
                continue
            
            chrom = line.strip().split(',')[0]
            start = int(line.strip().split(',')[2])
            try:
                truth_gt = truth_gts[(chrom, start)]
            except KeyError:
                truth_gt = './.'
            w.write(line.strip() + ',' + truth_gt + '\n')
