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
- Miscellaneous scripts
- Supporting R scripts
#

## Microbiome profiling

Metagenomes were profiled consistent with previous data analysis of 1000IBD and Lifelines-DEEP10 cohorts, as follows. KneadData tools (v0.5.1) were used to process metagenomic reads (in fastq format) by trimming the reads to PHRED quality 30 and removing Illumina adapters. Following trimming, the KneadData integrated Bowtie2 tool (v2.3.4.1) was used to remove reads that aligned to the human genome (GRCh37/hg19).
Taxonomic composition of metagenomes was profiled by MetaPhlAn2 tool (v2.7.2) using the MetaPhlAn database of marker genes mpa_v20_m200. Profiling of genes encoding microbial biochemical pathways was performed using the HUMAnN2 pipeline (v0.11.1) integrated with the DIAMOND alignment tool (v0.8.22), UniRef90 protein database (v0.1.1) and ChocoPhlAn pan-genome database (v0.1.1). As a final quality control step, samples with unrealistic microbiome composition (eukaryotic or viral abundance > 25% of total microbiome content or total read depth < 10 million) were excluded, leaving 8,208 samples for further analyses. Analyses were performed using locally installed tools and databases on CentOS (release 6.9) on the high-performance computing infrastructure available at our institution and using the MOLGENIS data platform2

An example of metagenome processing jobs is provided in *microbiome_profiling* folder

## Heritability analysis

We estimated the heritability of microbiome features using a variance components model implemented in the software POLY (v.0.5.1). We first considered a base model in which variance is partitioned into a polygenic component, Vg, shared between individuals that is proportional to their kinship coefficient, and an environmental component, Ve, that is unique to each individual. Thus, if Y is the measured trait, the variance of Y is Var(Y) =Vg+Ve, and its broad heritability is H2=Vg/(Vg+Ve). Of note, this measure reflects the overall impact of genes on the phenotype, thus including all potential models of action, as opposed to narrow heritability, which only includes additive effects. After fitting this base model, we also considered two refined models that included i) an additional variance component to model a shared environment for people belonging to the same family, Vh, and ii) a covariate with values 0/1 that distinguishes family members currently living in the same house from those who do not. We compared the significance of these two models to the basic model using a likelihood ratio test and found the base model to be the most appropriate fit for all microbiome traits.
We restricted the analysis of heritability to the relative abundances of the 242 microbial taxa present in at least 250 individuals and focused on 4,745 individuals in 2,756 families in which at least two individuals had available microbiome data. In total, the analysis included 2,756 parent-child pairs, 530 sibling pairs and 815 pairs with second-degree or more distant relationships. All models were adjusted for age, sex, BMI, read depth and stool frequency, and values of relative abundances of taxa were transformed using the centred log-ratio (clr) transformation36. Benjamini-Hochberg correction was used to control the multiple testing false discovery rate (FDR), and results with FDR < 0.1 were considered significant. 

POLY workflow and models are included in *heritability_analysis* folder

## Analysis of core and keystone microbiome features

To identify core microbial species and pathways, we used a bootstrapping-based selection approach. We randomly sampled 1% to 100% of the samples of the cohort a hundred times and calculated the standard deviation of the presence rate of each microbial species/pathway at different sampling percentages. Microbial features with a presence rate of more than 95% of samples were defined as the core microbiome. 
To analyse microbiome community structure, we constructed microbial species and pathway co-abundance networks using SparCC tool for network inference. Relative abundances of taxa were converted to estimated read counts by multiplying abundance percentages by total sequenced reads per sample after quality control. For pathway analysis, the read counts (RPKM) from HUMAnN2 were directly used for SparCC. Significant co-abundance was controlled at FDR 0.05 level using 100 permutations. In each permutation, the abundance of each microbial feature was randomly shuffled across samples.

Scripts for identification of core microbiome and keystone feature are in *core_keystone_microbes* folder


## Microbiome clustering

To identify microbial clusters and assess the presence of gut enterotypes in our cohort, we performed the partitioning around the medoid method on the relative abundances of microbial species and used the Calinski-Harabasz index to select the optimal number of clusters, as previously published in a study of gut enterotypes. Enrichment of phenotypes in each cluster was assessed by logistic regression in R

Codes used for clustering, plotting of clusters and enrichment analysis are in *microbiome_clustering* folder

## Miscellaneous scripts

*misc_scripts* folder contains miscellaneous used in data analysis not covered by other groups. For example, scripts for comparison of diet questionnaires over time and assigment of gastrointestinal disorders using ROMEIII criteria

## Supporting R scripts

*r_scripts_library* folder contains R functions used for various tasks (such as generation of supplementary data, plots, data parsing...). These functions were developed by Weersma group for internal use and are provided as-is, without comprehensive documentation
