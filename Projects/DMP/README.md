By: Weersma Group, Fu Group and Zhernakova Group, UMCG, The Netherlands 

Last update: 17/12/2020

# Dutch Microbiome Project

This github repo describes workflow and codes using in The Dutch Microbiome Project (DMP) study:

## Contents:

- Microbiome profiling
- Heritability analysis
- Identification of core and keystone microbes 
- Microbiome clustering
- Microbiome-association analyses
#

## Microbiome profiling

Metagenomes were profiled consistent with previous data analysis of 1000IBD and Lifelines-DEEP10 cohorts, as follows. KneadData tools (v0.5.1) were used to process metagenomic reads (in fastq format) by trimming the reads to PHRED quality 30 and removing Illumina adapters. Following trimming, the KneadData integrated Bowtie2 tool (v2.3.4.1) was used to remove reads that aligned to the human genome (GRCh37/hg19).
Taxonomic composition of metagenomes was profiled by MetaPhlAn2 tool (v2.7.2) using the MetaPhlAn database of marker genes mpa_v20_m200. Profiling of genes encoding microbial biochemical pathways was performed using the HUMAnN2 pipeline (v0.11.1) integrated with the DIAMOND alignment tool (v0.8.22), UniRef90 protein database (v0.1.1) and ChocoPhlAn pan-genome database (v0.1.1). As a final quality control step, samples with unrealistic microbiome composition (eukaryotic or viral abundance > 25% of total microbiome content or total read depth < 10 million) were excluded, leaving 8,208 samples for further analyses. Analyses were performed using locally installed tools and databases on CentOS (release 6.9) on the high-performance computing infrastructure available at our institution and using the MOLGENIS data platform2

An example of metagenome processing jobs is provided in *microbiome_profiling* folder

