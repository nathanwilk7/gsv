truth_gts = []
with open('truth-gts.txt', 'r') as f:
    for line in f:
        truth_gts.append(line.strip())

predicted_gts = []
with open('svtyper-quick-test.vcf', 'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        gt = line.strip().split('\t')[10].split(':')[0]
        predicted_gts.append(gt)

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

for i in range(200):
    t = truth_gts[i]
    p = predicted_gts[i]
    if t == '0/0':
        if p == '0/0':
            ref_ref += 1
        elif p== '0/1':
            ref_het += 1
        elif p == '1/1':
            ref_alt += 1
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '0/1':
        if p == '0/0':
            het_ref += 1
        elif p == '0/1':
            het_het += 1
        elif p == '1/1':
            het_alt += 1
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '1/1':
        if p == '0/0':
            alt_ref += 1
        elif p == '0/1':
            alt_het += 1
        elif p == '1/1':
            alt_alt += 1
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

