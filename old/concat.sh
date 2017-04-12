ILL=$1
ONT=$2

sed -i '/^#/ d' $ONT
cut -f 10 $ONT > $ONT.temp.gts.vcf
echo 'ont' > $ONT.name.vcf
cat $ONT.name.vcf $ONT.temp.gts.vcf > ont.gts.vcf
rm $ONT.name.vcf
rm $ONT.temp.gts.vcf

tail -n +50 $ILL > $ILL.temp.vcf
head -n 49 $ILL > head.vcf

paste $ILL.temp.vcf ont.gts.vcf > gts.vcf
cat head.vcf gts.vcf > plat-ont.vcf

rm $ILL.temp.vcf
rm head.vcf
rm gts.vcf
rm ont.gts.vcf

