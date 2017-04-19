setwd('/Users/nathanwilkinson/quinlan-lab/gsv/vis')
dats = read.csv('features-train.csv')
dats.dels = subset(dats, svtype == 'DEL')
# summary(dats)
refs = subset(dats.dels, true_genotype == '0/0')
hets = subset(dats.dels, true_genotype == '0/1')
alts = subset(dats.dels, true_genotype == '1/1')
refs[['spanners']] = refs[['left_breakpoint_spanners']] + refs[['right_breakpoint_spanners']]
hets[['spanners']] = hets[['left_breakpoint_spanners']] + hets[['right_breakpoint_spanners']]
alts[['spanners']] = alts[['left_breakpoint_spanners']] + alts[['right_breakpoint_spanners']]

refs[['splitters']] = refs[['left_splitters']] + refs[['right_splitters']]
hets[['splitters']] = hets[['left_splitters']] + hets[['right_splitters']]
alts[['splitters']] = alts[['left_splitters']] + alts[['right_splitters']]

# pdf('classify-full-genotypes-dels.pdf')
plot(hets[['ref_avg_base_coverage']], hets[['alt_avg_base_coverage']], col='blue', xlim=c(0,50), ylim=c(0,50), 
     main='Genotypes (ref_avg_base_coverage, alt_avg_base_coverage)', xlab='ref_avg_base_coverage', ylab='alt_avg_base_coverage')
points(alts[['ref_avg_base_coverage']], alts[['alt_avg_base_coverage']], col='green', pch=0)
points(refs[['ref_avg_base_coverage']], refs[['alt_avg_base_coverage']], col='red', pch=2)
legend(40, 10, pch=c(0,1,2), col=c('green', 'blue', 'red'), 
       c('Hom Alt', 'Het', 'Hom Ref'), text.width=21, y.intersp=1.3)
# dev.off()