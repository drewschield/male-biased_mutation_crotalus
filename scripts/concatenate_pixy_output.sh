# concatenate pi output
head -n 1 ./pixy_results/pilot_analysis_v2.male.intergenic.scaffold-ma1_pi.txt > ./pixy_results/pilot_analysis_v2.male.intergenic.all_pi.txt
for chrom in `cat ./processing_files/chrom.list`; do
	tail -n +2 ./pixy_results/pilot_analysis_v2.male.intergenic.${chrom}_pi.txt >> ./pixy_results/pilot_analysis_v2.male.intergenic.all_pi.txt
done
# concatenate dxy output
head -n 1 ./pixy_results/pilot_analysis_v2.male.intergenic.scaffold-ma1_dxy.txt > ./pixy_results/pilot_analysis_v2.male.intergenic.all_dxy.txt
for chrom in `cat ./processing_files/chrom.list`; do
	tail -n +2 ./pixy_results/pilot_analysis_v2.male.intergenic.${chrom}_dxy.txt >> ./pixy_results/pilot_analysis_v2.male.intergenic.all_dxy.txt
done
