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

# histogram analysis
## Prepare data for input to barplot
field.name = 'ref_avg_base_coverage'
breaks <- pretty(range(c(refs[field.name], 
                         hets[field.name], 
                         alts[field.name])), 
                 n=20)
refs.hist <- hist(data.matrix(refs[field.name]), breaks=breaks, plot=FALSE)$counts
hets.hist <- hist(data.matrix(hets[field.name]), breaks=breaks, plot=FALSE)$counts
alts.hist <- hist(data.matrix(alts[field.name]), breaks=breaks, plot=FALSE)$counts
dat <- rbind(refs.hist, hets.hist, alts.hist)
colnames(dat) <- paste(breaks[-length(breaks)], breaks[-1], sep="-")

## Plot it
barplot(dat, beside=TRUE, space=c(0, 0.1), las=2, 
        col = c('red', 'blue', 'green'), main=field.name)

# Remove alt/het calls with <= 1 splitter
#alts = subset(alts, Splitters > 1)
#hets = subset(hets, Splitters > 1)

# summary(refs[["Reference_Splitters"]])
# summary(hets[["Reference_Splitters"]])
# summary(alts[["Reference_Splitters"]])

# pdf('classify-full-genotypes-dels.pdf')
#plot(hets[['spanners']], hets[['spanners']], col='blue')
#     #main='Genotypes (Spanners, Splitters)', xlab='spanners', ylab='Splitters')
#points(alts[['spanners']], alts[['spanners']], col='green', pch=0)
#points(refs[['spanners']], refs[['spanners']], col='red', pch=2)
#legend(100, 60, pch=c(0,1,2), col=c('green', 'blue', 'red'), 
       c('Hom Alt', 'Het', 'Hom Ref'), text.width=21, y.intersp=1.3)
# dev.off()
