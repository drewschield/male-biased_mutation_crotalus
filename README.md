# Pipeline for male-biased mutation rate analysis in rattlesnakes

This file contains details about the data processing and analysis for estimation of the ratio of male-to-female mutation rates in rattlesnakes based on comparisons of Z-linked and autosomal sequence divergence. 
Created 09/30/2020 by Drew Schield

The steps described here require the following software:

* trimmomatic
* bwa
* samtools/bcftools/htslib
* vcftools
* GATK (v3 and v4)
* Pixy

## Contents

* Read filtering (update will link to raw reads)
* Read mapping
* Variant calling
* Variant filtering
* Pixy analysis

### Read Filtering

#### Set up environment

Get raw fastq data into `fastq` directory.