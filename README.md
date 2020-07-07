# Singletons-snakemake
Snakemake pipeline for singleton scores and stats calculation

This pipeline is based on the work presented in:

Mezzavilla M, Cocca M, Guidolin F, Gasparini P. A population-based approach for gene prioritization in understanding complex traits. Hum Genet. 2020;139(5):647-655. doi:10.1007/s00439-020-02152-4

---
##Requirements

+ Python <= 3.6
+ snakemake
+ bedtools
+ vcftools
+ bcftools
+ ...

---

##Pipeline steps:

1. Provide genomic region input
	* Get gene info from GENCODE data (if requested)
	* Provide a genomic region using chr:start-end informations

2. Calculate singletons for the study population
	* Computationally intensive step
	* This step will be performed in any case, so it should be a stand-alone procedure, to give the option to reuse already generated data

3. Calculate singleton stats
	+ This step will allow two option:
		1. Genomic region stats
			+ Singleton density for each gene across all samples
			+ Singleton count for each gene across all samples
		2. Sample level stats
			+ Genome wide singleton density for each sample
			+ Singleton density for each sample in the provided genomic region
			+ Singleton count for each sample in the provided genomic region

4. Calculate singleton scores
	+ This step will be performed only in case of genomic region pipeline

---
##Input files

+ Multisample vcf files for the study population
+ ...

---
##Sample command

