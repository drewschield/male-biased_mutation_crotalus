chromlist=$1
for chrom in `cat $chromlist`; do
	echo parsing $chrom VCF
	bcftools view --threads 16 -r $chrom -O z -o ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz
	tabix -p vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz
done
