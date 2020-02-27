###################################
### SAIGE gwas per microbe per chromosome
### version 2.0
### date: 29-01-2020
### Author:EALM
###################################
## New
## 06-02-2020
## change MAF filter for 1% in the association calculation
## changed test genotype files for the quality-imputation-filtered ones
## 13-02-2020
## added BMI to model
## commented the lines to equal phenotypes file to the list (given that the list this time is directly taken from this file no change is needed)
## renamed covariates to adjust to full features

library(SNPRelate)
library(SAIGEgds)
library(GWASTools)
library(data.table)
library(tidyverse)
library(optparse)

####################################
#output paths and/or input options####

## List of features or outcomes to be evaluated
#opt$outcome<-"k_Bacteria.p_Actinobacteria" 

## snps to use for the relationship matrix and fitting the null model
#opt$snps_for_grm<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/genotypes/UGLI_SNPs_for_GRM.gds"

## folder with imputed genotypes  of DAG3 samples
#opt$testgenotypes<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_gds/"

##output folder
#opt$out<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/outputs"

## microbiomefeatures & phenotypes files
#opt$features<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/features_and_phenotypes/DAG3_metaphlan_bac_log_normalized.txt"

## chromosome files
#opt$chr<-"2"

###test lauch
#Rscript 

#############################################
######################## arguments ##########

option_list = list(
  
  make_option(c("-m", "--outcome"), type="character", default=NULL, 
              help="List of features or outcomes to be evaluated", metavar="character"),
  
  make_option(c("-s", "--snps_for_grm"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-t", "--testgenotypes"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-f", "--features"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-c", "--chromosome"), type="character", default=NULL, 
              help="", metavar="character")
  
) 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

#############################################
########## main #############################
###  1. Bring phenotypes
# sample IDs in the phenotype file has to be the same as in genetic file (samples to use in the association analysis) 
mbfeatures <- fread(opt$features,data.table = F,header=T)
###format to match the list (remove change | for .) (and change  _ for __)
#colnames(mbfeatures)<-gsub(colnames(mbfeatures),pattern = "\\|",replacement = "\\.")
#colnames(mbfeatures)<-gsub(colnames(mbfeatures),pattern = "__",replacement = "_")

### 2. fit null model
    mb<-opt$outcome
    col.filter<-which(colnames(mbfeatures)==mb)
    select.rows<- which(mbfeatures[,col.filter]!=0)
    ecuation<-paste0(opt$outcome," ~ "," age +"," gender +"," conc.ng.ul +"," vol.ul +", " postclean.gc +"," ANTHRO.BMI "   )
    glmm2 <- seqFitNullGLMM_SPA( as.formula(ecuation) , 
                                 data=mbfeatures[select.rows,], 
                                 opt$snps_for_grm, 
                                 trait.type="quantitative", 
                                 sample.col="ID" ) #Thu Jan 30 21:50:48 2020 ## terminated Thu Jan 30 21:58:00 2020 with whole genome GRM
    ### correlate each covariate to see if they are correlated btween them.
### 3. calculate association
    genotype.file<-file.path(opt$testgenotypes,paste0("filtered.genotype.",opt$chr,".gds"))
    assoc <- seqAssocGLMM_SPA(genotype.file, glmm2, maf=0.01, parallel=2) ###modified from maf=0.05
    ### 4. save p-values
    outdir<-file.path(opt$out,paste0(mb,"_gwas"))
    dir.create(outdir,recursive = T)
    assoc.out<-file.path(outdir,paste0("p-values.","chr.",opt$chr))
    write.table(assoc,assoc.out,quote=F,row.names = F,col.names = T,sep='\t')
    
######




