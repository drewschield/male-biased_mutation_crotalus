for chrom in `cat chrom.list`; do
	pixy --stats pi dxy --vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz --zarr_path ./pixy_zarr --window_size 100000 --populations pilot_v2.male.popmap --variant_filter_expression 'DP>=5' --invariant_filter_expression 'DP>=5' --outfile_prefix ./pixy_results/pilot_analysis_v2.male.intergenic.$chrom
	rm -r pixy_zarr/$chrom/*
done
