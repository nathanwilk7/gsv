import sys
N = 100
n = 0

#bad = open('bad.vcf', 'w')
#bad = sys.stdout

#tps, fps, tns, fns = 0, 0, 0, 0

ref_called_ref = 0
ref_called_het = 0
ref_called_alt = 0
het_called_ref = 0
het_called_het = 0
het_called_alt = 0
alt_called_ref = 0
alt_called_het = 0
alt_called_alt = 0

for line in sys.stdin:
    if line[0] == "#": 
        #bad.write(line)
        continue
    n+=1
    gt = line.split("\t")[9].split(":")[0]
    if "1" in gt and n <= N:
       tps+=1
    elif n <= N:
       fns+=1
       #bad.write(line)
    else: # into the true negative sets
        if "1" in gt:
           fps+=1
        else:
           tns+=1

print("tn:", tns)
print("fn:", fns)
print("tp:", tps)
print("fp:", fps)
