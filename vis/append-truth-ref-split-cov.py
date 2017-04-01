output_stuff = []
dups = set()
truth_gts = []
with open('../results/truth-gts.txt', 'r') as f:
    for line in f:
        truth_gts.append(line.strip())

predicted_gts = []
natee = -1
with open('ref-split-cov-no-truth.csv', 'r') as f:
    for line in f:
        if natee == -1:
            natee += 1
            continue
        features = line.strip().split(',')
        svtype = features[3]
        genotype = features[4]
        predicted_gts.append(genotype)
        is_dup = svtype == 'DUP'
        if is_dup:
            dups.add(natee)
        natee += 1
        output_stuff.append((features[0], features[1], features[2], features[3], features[4], float(features[5]), 
                             float(features[6]), float(features[7]), float(features[8]), float(features[9]), float(features[10]), 
                             float(features[11]), float(features[12]), float(features[13])))

# Chrom,Left_Pos,ID,SV_Type,Genotype,Spans_Left_Breakpoint,Spans_Right_Breakpoint,Clipped_By_Breakpoint,Mismatched_Over_SV,Coverage_Diff,Splitters,
# Reference_Split\ters,Ref_Avg_Base_Coverage,Alt_Avg_Base_Coverage

if len(truth_gts) != len(predicted_gts):
    print("Length of truth set and predictions is different")
    exit()

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

def print_as_csv(mylist, mystr):
    for el in mylist:
        print(el, end=',')
    print(mystr)

print('SV_Type,Called_Genotype,Spans_Left_Breakpoint,Spans_Right_Breakpoint,Clipped_By_Breakpoint,Mismatched_Over_SV,Coverage_Diff,Splitters,Reference_Splitters,Truth_Genotype')

for i in range(200):
    if i in dups:
        continue
    t = truth_gts[i]
    p = predicted_gts[i]
    if t == '0/0':
        if p == '0/0':
            print_as_csv(output_stuff[i], truth_gts[i])
            ref_ref += 1
            if debug:
                print('ref_ref:', i)
        elif p== '0/1':
            ref_het += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('ref_het:', i)
        elif p == '1/1':
            ref_alt += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('ref_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '0/1':
        if p == '0/0':
            het_ref += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('het_ref:', i)
        elif p == '0/1':
            het_het += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('het_het:', i)
        elif p == '1/1':
            het_alt += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('het_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '1/1':
        if p == '0/0':
            alt_ref += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('alt_ref:', i)
        elif p == '0/1':
            alt_het += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('alt_het:', i)
        elif p == '1/1':
            alt_alt += 1
            print_as_csv(output_stuff[i], truth_gts[i])
            if debug:
                print('alt_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    else:
        print("Truth not supported: ", t)
        exit()
"""
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
