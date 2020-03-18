###################################
### File convertion from VCF to GDS
### version 1.0
### date: 01-20-2020
### Author:EALM
###################################
##New

### load necessary packages and create functions

library(SNPRelate)
library(SAIGEgds)
library(GWASTools)
library(data.table)
library(tidyverse)
library(optparse)
#################
#test paths
#opt<-list()
#opt$chr=22
#opt$workdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG"
#opt$root=".DAG3_imputed_DAG_ID.vcf.gz"

###test run
#Rscript vcftoGDS.R -w "/groups/umcg-ugli/tmp01/umcg-elopera/testDAG" \
#                   -c "22" \
#                   -r ".DAG3_imputed_DAG_ID.vcf.gz"


#############################################
######################## arguments ##########

option_list = list(
  
  make_option(c("-w", "--workdir"), type="character", default=NULL, 
              help="working directory of both input output files", metavar="character"),
  
  make_option(c("-c", "--chr"), type="character", default=NULL, 
              help="", metavar="character"),
  
  make_option(c("-r", "--root"), type="character", default=NULL, 
              help="root (fixed part) file name of the file(s) to convert", metavar="character")
  
  
) 

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
################################################################################

## turn the vcf GENOTYPE file into a seqGDS file (it has to be seq file bcause SAIGEgds does not accept snpGDS calss)

## chromosome list files
working.dir<-file.path(opt$workdir)
setwd(opt$workdir)

  file<-paste0(opt$chr,opt$root) 
  seqVCF2GDS(file, paste0("genotype.",opt$chr,".gds"), fmt.import="DS") 

##done
  
