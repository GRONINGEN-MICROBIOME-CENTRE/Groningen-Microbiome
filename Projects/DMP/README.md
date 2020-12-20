By: Weersma Group, Fu Group and Zhernakova Group, UMCG, The Netherlands 

Last update: 17/12/2020

# Dutch Microbiome Project

This github repo describes workflow and codes using in The Dutch Microbiome Project (DMP) study:

## Contents:

- Microbiome profiling
- Heritability analysis
- Identification of core and keystone microbes 
- Microbiome clustering
- Calculation of microbiome variance explained by phenotypes
- Microbiome-phenotypes association analyses
- Calculation of General Microbiome Health Index (GMHI)
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

## Calculation of microbiome variance explained by phenotypes

The microbiome composition variance explained by phenotypes was calculated by permutational multivariate analysis of variance using distance matrices, implemented in the adonis function for R package vegan (v.2.4-6), using 20,000 permutations and a Bray-Curtis distance matrix calculated using relative abundances of microbial species. A separate analysis was performed to calculate the microbiome functional potential explained by phenotypes using equivalent methodology. The functional dissimilarity matrix was calculated using the Bray-Curtis dissimilarity index calculated on the relative abundances of MetaCyc microbial biochemical pathways.

Scripts used for calculation are in *microbiome_variance*

## Microbiome-phenotypes association analyses

Prior to the association analysis of phenotypes and microbiome features, the microbiome data was transformed using the clr transformation. The geometric mean for clr transformation of relative abundances of taxa was calculated on species-level and applied to higher levels. The associations between phenotypes and microbial features (microbial taxa, MetaCyc functional pathways, CARD and VFDB entities) were calculated using linear regression, adjusting for age, sex and BMI of the individual along with Bristol stool scale of the faecal sample and technical factors (DNA concentration, sequencing read depth, sequencing batch and sampling season). Benjamini-Hochberg correction was used to control for multiple testing with the number of tests equal to the number of tested feature‒phenotype pairs. 

*association_analysis* folder contains scripts used for association analysis

## Calculation of microbiome signatures predictive of diseases and health

We calculated the microbial signatures predictive of the 36 most common (Ncases > 100) diseases in our dataset. In addition, we defined a “healthy” phenotype as an absence of any self-reported disease. Using this definition, 2,937 (36%) out of 8,208 individuals were defined as “healthy”. To build prediction models for common diseases, the dataset was randomly split into training (90%) and test (10%) sets. Next, we performed elastic net L1/L2 regularized regression (R package glmnet v.4.0) on the training set, using Shannon diversity, clr-transformed microbial taxa, clr-transformed MetaCyc bacterial pathways and age, sex and BMI as fixed covariates (not penalized in the models). The model for each disease was calculated independently using five-fold cross-validation to select the optimal lambda penalization factor (at L1/L2 mixing parameter alpha fixed at 0.5). The lambda with minimal cross-validation error was used in the downstream analysis. In total, we defined three probabilistic models: a “null” signature that only includes effects of general covariates (age, sex and BMI), a “microbiome” signature that includes all selected microbiome features and a “combined” signature that includes both the effects of microbiome features and general covariates.

*health_disease_prediction* lists scripts used for training and testing of prediction models

## Calculation of General Microbiome Health Index (GMHI)

We calculated the recently developed Gut Microbiome Health Index (GMHI) for DMP data, using the parameters identified in GMHI study (Gupta et al., Nature Comm. 2020) on our data. 

*ghmi* folder contains scripts used for calcuation of GMHI

## Miscellaneous scripts

*misc_scripts* folder contains miscellaneous used in data analysis not covered by other groups. For example, scripts for comparison of diet questionnaires over time and assigment of gastrointestinal disorders using ROMEIII criteria

## Supporting R scripts

*r_scripts_library* folder contains R functions used for various tasks (such as generation of supplementary data, plots, data parsing...). These functions were developed by Weersma group for internal use and are provided as-is, without comprehensive documentation
