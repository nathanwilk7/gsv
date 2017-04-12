# for debugging, run gsv.py but don't redirect stdout so we can see the output on the screen if we hit a breakpoint
python ../gsv.py -b ../data/nanopore-NA12878.params.merged.bam -v ../data/ill.gts.vcf -o debug.vcf 
# add this if we want to get actual mismatches: -f ../data/g1k_v37_decoy.fa
