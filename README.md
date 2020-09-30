# Pipeline for male-biased mutation rate analysis in rattlesnakes

This file contains details about the data processing and analysis for estimation of the ratio of male-to-female mutation rates in rattlesnakes based on comparisons of Z-linked and autosomal sequence divergence. 
Created 09/30/2020 by Drew Schield

The steps described here require the following software:

* trimmomatic
* bwa
* samtools/bcftools/htslib/bgzip
* vcftools
* GATK (v3 and v4.0.8.1)
* Pixy

List files and shell script examples are in the `processing_files` directory.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.

Note that you may need to adjust the organization of your environment to suite your workflow.

## Contents

* Read filtering (update will link to raw reads)
* Read mapping (update will link to reference genome & annotation)
* Variant calling
* Variant filtering
* Pixy analysis

### Read Filtering

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

We will call variants using GATK, specifying 'all-sites' output for downstream calculation of summary statistics using variant and invariant sites.

#### Set up environment

Make `gvcf` and `vcf` directories.

```
mkdir gvcf
mkdir vcf
```

#### Create sequence dictionary for reference genome

`./gatk-4.0.8.1/gatk CreateSequenceDictionary -R CroVir_genome_L77pg_16Aug2017.final_rename.fasta`

#### Call individual variants using GATK HaplotypeCaller

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














































 