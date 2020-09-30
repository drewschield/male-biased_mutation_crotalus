# Pipeline for male-biased mutation rate analysis in rattlesnakes

This file contains details about the data processing and analysis for estimation of the ratio of male-to-female mutation rates in rattlesnakes based on comparisons of Z-linked and autosomal sequence divergence. 
Created 09/30/2020 by Drew Schield

The steps described here require the following software:

* trimmomatic
* bwa
* samtools/bcftools/htslib/bgzip/tabix
* vcftools
* GATK (v3.8-1-0 and v4.0.8.1)
* conda
* Pixy

List and miscellaneous files are in the `processing_files` directory.
Example shell scripts are in the `scripts` directory. These will need to be moved accordingly.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.

Note that you may need to adjust the organization of your environment (e.g., script locations) to suite your workflow.

__*If you want to go straight to mutation rate ratio calculations, skip ahead to the [Analysis in R](#analysis-in-r) section and go nuts!*__

## Contents

* [Read filtering (update will link to raw reads)](#read-filtering)
* [Read mapping](#read-mapping)
* [Variant calling](#variant-calling)
* [Variant filtering](#variant-filtering)
* [Pixy analysis](#pixy-analysis)
* [Analysis in R](#analysis-in-r)

### Read filtering

We impose these filters to trim reads:

* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

#### Set up environment

Get raw fastq data into `fastq` directory. <br /> Make a `fastq_filtered` directory for trimmed output.

```
mkdir fastq
mkdir fastq_filtered
```

#### Filter reads with `trimmomatic`

The script below will loop through samples in `processing_files/sample.list`.

trimmomatic.sh:

```
list=$1
for line in `cat $list`; do
	name=$line
	echo processing and filtering ${name}.
	trimmomatic PE -phred33 -threads 16 ./fastq/${name}_R1_001.fastq.gz ./fastq/${name}_R2_001.fastq.gz ./fastq_filtered/${name}_R1_P.trim.fastq.gz ./fastq_filtered/${name}_R1_U.trim.fastq.gz ./fastq_filtered/${name}_R2_P.trim.fastq.gz ./fastq_filtered/${name}_R2_U.trim.fastq.gz LEADING:20 TRAILING:20 MINLEN:32 AVGQUAL:30
done
```

`sh trimmomatic.sh ./processing_files/sample.list`

### Read mapping

We will map filtered reads to the C. viridis reference genome using `bwa`.

You can download the genome [here](https://ndownloader.figshare.com/files/16522091).

#### Set up environment, map reads with `bwa`, sort with `samtools`

`mkdir bam`

bwa_mem.sh:

```
list=$1
for line in `cat $list`; do
	name=$line
	echo "Mapping filtered $name data to reference"
	bwa mem -t 16 -R "@RG\tID:$name\tLB:CVOS\tPL:illumina\tPU:NovaSeq6000\tSM:$name" CroVir_genome_L77pg_16Aug2017.final_rename.fasta ./fastq_filtered/${name}_*_R1_P.trim.fastq.gz ./fastq_filtered/${name}_*_R2_P.trim.fastq.gz | samtools sort -@ 16 -O bam -T temp -o ./bam/$name.bam -
done
```

`sh bwa_mem.sh .processing_files/sample.list`

### Variant calling

We will call variants using `GATK`, specifying 'all-sites' output for downstream calculation of summary statistics using variant and invariant sites.

#### Set up environment

Make `gvcf` and `vcf` directories.

```
mkdir gvcf
mkdir vcf
```

#### Create sequence dictionary for reference genome

`./gatk-4.0.8.1/gatk CreateSequenceDictionary -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta`

#### Call individual variants using `GATK HaplotypeCaller`

This step will also compress and index output using `bgzip` and `tabix`.

GATK_HaplotypeCaller.sh

```
list=$1
for i in `cat $list`; do
	./gatk-4.0.8.1/gatk HaplotypeCaller -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta --ERC GVCF -I ./bam/$i.bam -O ./gvcf/$i.raw.snps.indels.g.vcf
	bgzip ./gvcf/$i.raw.snps.indels.g.vcf
	tabix -p vcf ./gvcf/$i.raw.snps.indels.g.vcf.gz
done
```

`sh GATK_HaplotypeCaller.sh ./processing_files/sample.list`

#### Call cohort variants using `GATK GenotypeGVCFs`, specifying 'all-sites' output

This command will use the gvcf files in `./processing_files/gvcf.pilot_analysis_v2.male.list`.

```
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T GenotypeGVCFs -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta -V ./processing_files/gvcf.pilot_analysis_v2.male.list -allSites -o ./vcf/pilot_analysis_v2.male.raw.vcf.gz
```

### Variant filtering

We will filter to set sites meeting these criteria as missing genotypes:

* Low depth or low quality genotypes
* Sites overlapping with repeat or gene annotations

Note: zipped copies of the gene and repeat BED annotations are in `./processing_files`.

#### Index feature files for masking using `GATK IndexFeatureFile`

```
./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file CroVir_genome_L77pg_16Aug2017.final.reformat.repeat.masked.sort.chrom.bed

./gatk-4.0.8.1/gatk IndexFeatureFile --feature-file CroVir_genome_gene.chrom.sort.bed
```

#### Mask repeat bases using `GATK VariantFiltration`

```
java -jar ./gatk-3.8-1-0/GenomeAnalysisTK.jar -T VariantFiltration -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta --mask CroVir_genome_L77pg_16Aug2017.final.reformat.repeat.masked.sort.chrom.bed --maskName REP --setFilteredGtToNocall --variant ./vcf/pilot_analysis_v2.male.raw.vcf.gz --out ./vcf/pilot_analysis_v2.male.mask.vcf.gz
```

#### Filter low quality and repeat bases using `bcftools` and index using `tabix` 

```
bcftools filter --threads 24 -e 'FORMAT/DP<5 | FORMAT/GQ<30 || TYPE="indel" || FILTER="REP"' --set-GTs . -O z -o ./vcf/pilot_analysis_v2.male.mask.HardFilter.vcf.gz ./vcf/pilot_analysis_v2.male.mask.vcf.gz

tabix -p vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.vcf.gz
```

#### Mask genes using `GATK VariantFiltration`

We want to analyze intergenic regions only to best approximate neutral mutation rates across the genome.

You're probably thinking, "we could have masked genes in the same step as repeats". You're right! I was lazy and returned to the hard-filtered VCF after deciding to analyze intergenic regions only. Fix as you please!

```
./gatk-3.8-1-0/GenomeAnalysisTK.jar -T VariantFiltration -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta --mask CroVir_genome_gene.chrom.sort.bed --maskName GENE --setFilteredGtToNocall --variant ./vcf/pilot_analysis_v2.male.mask.HardFilter.vcf.gz --out ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.vcf.gz
```

#### Filter masked gene bases using `bcftools`

Again, this is a totally needless step if you mask genes and repeats at the same time. If you're smarter than me and do so, adding `|| FILTER="GENE"` to the `bcftools filter` "-e" flag will do the trick.

```
bcftools filter --threads 24 -e 'FILTER="GENE"' --set-GTs . -O z -o ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.vcf.gz

tabix -p vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz
```

#### Parse all-sites VCFs for each chromosome

We'll perform analysis in `pixy` one chromosome at a time, calling commands in a loop on the set of chromosome-specific VCFs.

We'll parse the all-sites VCF using `bcftools` and `tabix`.

parse_chrom_vcf.sh:

```
chromlist=$1
for chrom in `cat $chromlist`; do
	echo parsing $chrom VCF
	bcftools view --threads 16 -r $chrom -O z -o ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz
	tabix -p vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz
done
```

`sh parse_chrom_vcf.sh ./processing_files/chrom.list`

### Pixy analysis

Check out [pixy](https://pixy.readthedocs.io/en/latest/), the new method for unbiased estimation of diversity and divergence statistics from all-sites VCFs.

#### Set up environment

Pixy is straightforward to install via conda and by following the authors' steps [here](https://pixy.readthedocs.io/en/latest/installation.html).
In my case, I made a new conda environment for `pixy`, so prior to running analyses I needed to run:

```
conda deactivate

conda activate pixy
```

Make `pixy_results` and `pixy_zarr` directories.

```
mkdir pixy_results

mkdir pixy_zarr
```

#### Run `pixy`

We'll run `pixy` on each chromosome in 100 Kb sliding windows to estimate pi and dxy.

Analyses will use the sample 'population map' in `./processing_files/pilot_v2.male.popmap` and the `chrom.list` above.

pixyloop.sh:

```
for chrom in `cat chrom.list`; do
	pixy --stats pi dxy --vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.$chrom.vcf.gz --zarr_path ./pixy_zarr --window_size 100000 --populations pilot_v2.male.popmap --variant_filter_expression 'DP>=5' --invariant_filter_expression 'DP>=5' --outfile_prefix ./pixy_results/pilot_analysis_v2.male.intergenic.$chrom
	rm -r pixy_zarr/$chrom/*
done
```

`sh pixyloop ./processing_files/chrom.list`

Note: `pixyloop.sh` removes the intermediate zarr files used in `pixy` analysis once the analysis is complete.
These files can get quite large, which is why I remove them to save HD space, but they are useful for re-running analyses.
Adjust the script to keep your zarr files if you wish!

#### Concatenate `pixy` results

Bring the chromosome-specific results together for each statistic.

concatenate_pixy_output.sh:

```
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
```

`sh concatenate_pixy_output.sh`

### Analysis in R

The results from `pixy` can now be read into R to calculate the ratio of male-to-female mutation rate and the ratio of Z chromosome-to-autosome mutation rate.

Summary statistics results files are in the `./pixy_results` directory.

Perform analyses using `male-biased_mutation_calculations_crotalus.R`











































 