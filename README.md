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

List files and shell script examples are in the `processing_files` directory.
Population genetic summary statistics output from `pixy` are in the `pixy_results` directory.

Note, you may need to adjust the organization of your environment to suite your workflow.

## Contents

* Read filtering (update will link to raw reads)
* Read mapping
* Variant calling
* Variant filtering
* Pixy analysis

### Read Filtering

I will impose these filters to trim reads:

* Remove 5' end bases if quality is below 20
* Remove 3' end bases if quality is below 20
* Minimum read length = 32
* Remove reads if average quality is < 30

#### Set up environment

Get raw fastq data into `fastq` directory.