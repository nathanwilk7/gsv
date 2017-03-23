from cyvcf2 import VCF
import pdb

HOM_REF = 0
HET = 1
UNKNOWN = 2
HOM_ALT = 3

chr_pos_id_to_truth_gts = {}

#pdb.set_trace()
for variant in VCF('ill.gst-symlink.vcf'):
    chr_pos_id_to_truth_gts[variant.CHROM, variant.start, variant.ID] = variant.gt_types[0]

print('Chrom,Left_Pos,ID,SV_Type,Called_Genotype,Spans_Left_Breakpoint,Spans_Right_Breakpoint,Clipped_By_Breakpoint,Mismatched_Over_SV,Coverage_Diff,Splitters,Reference_Splitters,True_Genotype')

#pdb.set_trace()
natee = 0
with open('full-to-plot-2.csv', 'r') as f:
    for line in f:
        if natee == 0:
            natee += 1
            continue
        features = line.strip().split(',')
        chr = features[0]
        pos = features[1]
        ID = features[2]
        try:
            true_gt = chr_pos_id_to_truth_gts[str(chr), int(pos) - 1, str(ID)]
        except KeyError:
            true_gt = UNKNOWN
        if true_gt == HOM_REF:
            true_gt = '0/0'
        elif true_gt == HET:
            true_gt = '0/1'
        elif true_gt == HOM_ALT:
            true_gt = '1/1'
        elif true_gt == UNKNOWN:
            true_gt = './.'
            
        for i in range(len(features)):
            print(features[i], end=',')
        print(true_gt)
"""
chr_pos_id_to_called_gts = {}
for variant in VCF('gsv-results2.vcf'):
    chr_pos_id_to_called_gts[variant.CHROM, variant.start, variant.ID] = variant.gt_types[0]

ref_ref = 0
ref_het = 0
ref_alt = 0
het_ref = 0
het_het = 0
het_alt = 0
alt_ref = 0
alt_het = 0
alt_alt = 0

debug = False # TODO: debug

for k in chr_pos_id_to_truth_gts.keys():
    try:
        t = chr_pos_id_to_truth_gts[k]
        p = chr_pos_id_to_called_gts[k]
    except KeyError:
        continue

    if t == HOM_REF:
        if p == HOM_REF:
            ref_ref += 1
            if debug:
                print('ref_ref:', i)
        elif p== HET:
            ref_het += 1
            if debug:
                print('ref_het:', i)
        elif p == HOM_ALT:
            ref_alt += 1
            if debug:
                print('ref_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == HET:
        if p == HOM_REF:
            het_ref += 1
            if debug:
                print('het_ref:', i)
        elif p == HET:
            het_het += 1
            if debug:
                print('het_het:', i)
        elif p == HOM_ALT:
            het_alt += 1
            if debug:
                print('het_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == HOM_ALT:
        if p == HOM_REF:
            alt_ref += 1
            if debug:
                print('alt_ref:', i)
        elif p == HET:
            alt_het += 1
            if debug:
                print('alt_het:', i)
        elif p == HOM_ALT:
            alt_alt += 1
            if debug:
                print('alt_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    else:
        print("Truth not supported: ", t)
        exit()

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
"""
