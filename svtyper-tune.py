truth_gts = []
with open('truth-gts.txt', 'r') as f:
    for line in f:
        truth_gts.append(line.strip())

dups = set()
predicted_gts = []
natee = 0
with open('svtyper-quick-test.vcf', 'r') as f:
    for line in f:
        if line[0] == '#':
            continue
        gt = line.strip().split('\t')[10].split(':')[0]
        predicted_gts.append(gt)
        is_dup = line.split('\t')[4].strip() == '<DUP>'
        if is_dup:
            dups.add(natee)
        natee += 1

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

debug = True # TODO: debug

for i in range(200):
    if i in dups:
        continue
    t = truth_gts[i]
    p = predicted_gts[i]
    if t == '0/0':
        if p == '0/0':
            ref_ref += 1
            if debug:
                print('ref_ref:', i)
        elif p== '0/1':
            ref_het += 1
            if debug:
                print('ref_het:', i)
        elif p == '1/1':
            ref_alt += 1
            if debug:
                print('ref_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '0/1':
        if p == '0/0':
            het_ref += 1
            if debug:
                print('het_ref:', i)
        elif p == '0/1':
            het_het += 1
            if debug:
                print('het_het:', i)
        elif p == '1/1':
            het_alt += 1
            if debug:
                print('het_alt:', i)
        else:
            print("Predicted not supported: ", p)
            exit()
    elif t == '1/1':
        if p == '0/0':
            alt_ref += 1
            if debug:
                print('alt_ref:', i)
        elif p == '0/1':
            alt_het += 1
            if debug:
                print('alt_het:', i)
        elif p == '1/1':
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

