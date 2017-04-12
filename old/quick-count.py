from cyvcf2 import VCF
import argparse
p = argparse.ArgumentParser()
p.add_argument('--vcf')
p.add_argument('-vout', default=False)
p.add_argument('-p', default=False)
args = p.parse_args()

HOM_REF = [0]
HET = [1]
HOM_ALT = [3]
UNKNOWN = [2]
EITHER = [1,3]

if args.vout:
    vcf = VCF(args.vcf)
    print vcf.raw_header

ha_counts = [0,0,0]
het_counts = [0,0,0]
hr_counts = [0,0,0]

def count_ont(ont, li):
    if ont == 0:
        li[0] += 1
    elif ont == 1:
        li[1] += 1
    elif ont == 3:
        li[2] += 1

either_y, either_n = 0, 0

for v in VCF(args.vcf):
    ad, rd = v.gt_ref_depths, v.gt_alt_depths
    ont_ad, ont_rd = ad[1], rd[1]
    ill_ad, ont_ad = ad[0], rd[0]
    if ont_ad + ont_rd < 20:
        continue
    if v.INFO.get('SVTYPE') == 'DUP':
        continue
    gts = v.gt_types
    ill, ont = gts[0], gts[1]
    if ill == 0:
        count_ont(ont, hr_counts)
        if ont != 0 and args.vout:
            print str(v),
    elif ill == 1:
        count_ont(ont, het_counts)
    elif ill == 3:
        count_ont(ont, ha_counts)
    if (ill == 1 or ill == 3) and (ont == 1 or ont == 3):
        either_y += 1
    elif (ill == 1 or ill == 3) and (ont != 1 and ont != 3):
        either_n += 1
print 'either: ', either_y / float(either_y + either_n)
if args.p:
    print 'I_HA\t', '\t'.join(['%.3f' % (n / float(sum(ha_counts))) for n in ha_counts])
    print 'I_H\t', '\t'.join(['%.3f' % (n / float(sum(het_counts))) for n in het_counts])
    print 'I_HR\t', '\t'.join(['%.3f' % (n / float(sum(hr_counts))) for n in hr_counts])
else:
    print 'I_HA\t', '\t'.join([str(n) for n in ha_counts])
    print 'I_H\t', '\t'.join([str(n) for n in het_counts])
    print 'I_HR\t', '\t'.join([str(n) for n in hr_counts])
print '\tO_HR\tO_H\tO_HA'

