# Pipeline for male-biased mutation rate analysis in rattlesnakes

This file contains details about the data processing and analysis for estimation of the ratio of male-to-female mutation rates in rattlesnakes based on comparisons of Z-linked and autosomal sequence divergence. 
Created 09/30/2020 by Drew Schield

The steps described here require the following software:

* trimmomatic
* bwa
* samtools/bcftools/htslib/bgzip/tabix
* vcftools
* GATK (v3.8-1-0 and v4.0.8.1)
* Pixy

List and miscellaneous files are in the `processing_files` directory.
Example shell scripts are in the `scripts` directory. These will need to be moved accordingly.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.

Note that you may need to adjust the organization of your environment (e.g., script locations) to suite your workflow.

## Contents

* Read filtering (update will link to raw reads)
* Read mapping (update will link to reference genome & annotation)
* Variant calling
* Variant filtering
* Pixy analysis

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

We will map filtered reads to the C. viridis reference genome (add link) using `bwa`.

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

Again, this is a totally needless step if you mask genes and repeats at the same time. If you're smarter than me and do so, just adding `FILTER="GENE"` to the `bcftools filter` -e flag will do the trick.

```
bcftools filter --threads 24 -e 'FILTER="GENE"' --set-GTs . -O z -o ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.vcf.gz

tabix -p vcf ./vcf/pilot_analysis_v2.male.mask.HardFilter.intergenic.filter.vcf.gz
```

#### Parse chromosome all-sites VCFs






































 