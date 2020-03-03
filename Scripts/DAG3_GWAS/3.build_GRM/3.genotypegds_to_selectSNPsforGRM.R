###################################
### step 3. Select SNPs for GRM from UGLI genotype
### version 1.0
### date: 29-01-2020
### Author:EALM
###################################
##New


### load necessary packages and create functions

library(SNPRelate)
library(GWASTools)
library(data.table)
library(tidyverse)
library(optparse)

####################################
#output paths and/or input files####

#opt<-list()
#opt$genotypes_VCF<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/genotypes/genotype.full.gds"
#opt$genotype_for_grm<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/genotypes/UGLI_SNPs_for_GRM.gds"


#############################################
######################## arguments ##########

option_list = list(
  
  make_option(c("-g", "--genotypes_VCF"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-r", "--genotype_for_grm"), type="character", default=NULL, 
              help="", metavar="character")
) 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#############################################
########## main #############################


###TAKE NOTE: as this is for the GRM, we need only the snps that are informative of relationship between the individuals. COnsidering most of the 
# variability and IBD is already contained in the genotype it would be unnecesary and redundanr to use the imputed data. therefore the opt$genotypes_VCF comes from the genotyping and not the imputed vcf files. Also note in spite of the name opt$genotypes_VCF  is a GDS file!!

## open the created GDS file
gds_grm <- seqOpen(opt$genotypes_VCF)

## Prune SNPs for GRM
set.seed(100)## make the random choosing reproducible
snpset <- snpgdsLDpruning(gds_grm,maf=0.05,ld.threshold = 0.2, slide.max.bp=500000)
snpset.id <- unlist(snpset, use.names=FALSE) 

## select pruned snps from the genotype_gds
seqSetFilter(gds_grm, variant.id=snpset.id)

## export selected genotypes without all the annotation information 
seqExport(gds_grm, opt$genotype_for_grm, info.var=character(), fmt.var=character(), samp.var=character())

