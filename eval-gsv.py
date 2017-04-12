import argparse
import cyvcf2
import pdb

# parse args
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--truth', help='Truth VCF')
parser.add_argument('-p', '--predicted', help='Predicted VCF')
parser.add_argument('-d', '--debug', action='store_true', help='Debug')
args = parser.parse_args()

debug = args.debug

# read variants from truth vcf
truth_chrom_pos_to_svtype_gt = {}
for truth_variant in cyvcf2.VCF(args.truth):
    truth_chrom_pos_to_svtype_gt[(truth_variant.CHROM, truth_variant.start)] = (truth_variant.INFO['SVTYPE'], truth_variant.genotypes)

# read variant from predicted vcf
pred_chrom_pos_to_svtype_gt = {}
for pred_variant in cyvcf2.VCF(args.predicted):
    pred_chrom_pos_to_svtype_gt[(pred_variant.CHROM, pred_variant.start)] = (pred_variant.INFO['SVTYPE'], pred_variant.genotypes)

REF = [0,0]
HET = [0,1]
ALT = [1,1]

# results by category
ref_ref = 0
ref_het = 0
ref_alt = 0
het_ref = 0
het_het = 0
het_alt = 0
alt_ref = 0
alt_het = 0
alt_alt = 0

# compare truth and pred for all DELs by iterating through keys of truth_dict and then massive if bidness from other...
for truth_chrom_pos in truth_chrom_pos_to_svtype_gt.keys():
    truth_svtype_gt = truth_chrom_pos_to_svtype_gt[truth_chrom_pos]
    pred_svtype_gt = pred_chrom_pos_to_svtype_gt[truth_chrom_pos]
    if truth_svtype_gt[0] == 'DUP':
        continue
    t = truth_svtype_gt[1][0][:2] # drop off the True/False field of genotype cause I don't know what it is
    p = pred_svtype_gt[1][0][:2]
    if t == REF:
        if p == REF:
            ref_ref += 1
            if debug:
                print('ref_ref:', i)
        elif p== HET:
            ref_het += 1
            if debug:
                print('ref_het:', i)
        elif p == ALT:
            ref_alt += 1
            if debug:
                print('ref_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == HET:
        if p == REF:
            het_ref += 1
            if debug:
                print('het_ref:', i)
        elif p == HET:
            het_het += 1
            if debug:
                print('het_het:', i)
        elif p == ALT:
            het_alt += 1
            if debug:
                print('het_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == ALT:
        if p == REF:
            alt_ref += 1
            if debug:
                print('alt_ref:', i)
        elif p == HET:
            alt_het += 1
            if debug:
                print('alt_het:', i)
        elif p == ALT:
            alt_alt += 1
            if debug:
                print('alt_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    else:
        print("Truth not supported: ", truth_chrom_pos, t)
        continue

num_correct = ref_ref + het_het + alt_alt
num_wrong = ref_het + ref_alt + het_ref + het_alt + alt_ref + alt_het
pct_correct = float(num_correct) / (float(num_correct) + float(num_wrong))
print('% Correct:', pct_correct)
print('False Positives:', ref_het + ref_alt)
print()
print('ref_ref:', ref_ref)
print('ref_het:', ref_het)
print('ref_alt:', ref_alt)
print('het_ref:', het_ref)
print('het_het:', het_het)
print('het_alt:', het_alt)
print('alt_ref:', alt_ref)
print('alt_het:', alt_het)
print('alt_alt:', alt_alt)
