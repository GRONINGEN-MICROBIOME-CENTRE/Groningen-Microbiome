##version 1.0.0 from 7  July 2020.
##Written by Alex Kur

# Instruction manual (simple)
# 1. go to Gearshift server
# 2. load RPlus package:
#   $ module load RPlus
# 3. open R  
#   $ R
# 4. Load script as a source 
#   > source("run_multiple_phenos.R")
#   > run.associations(phenotypes, trait, trait_group = "taxa")
# In this function, 'phenotypes' is a vector with phenotypes of your interest, 'trait' is microbial trait of corresponding type,
# and trait_group is a group of microbiome data, should be 'taxa','pathways','KOterms','CARD' , 'VFDB' or 'alpha'
# for "trait_group = alpha", trait should be: trait='shannon'  
# Phenotypes could be vector of one or few phenotype names, 
# for example:
#   > c("MED.BLOOD.Lymphocytes_10E9perL","MED.BLOOD.Bpressure.systolic.mm.Hg")
# OUTPUT:
# result = run_associations(parameters)
# result$anova - overall significance of all phenotypes (similar to anova() output)
# result$summary - conditional significance of each phenotypes (similar to summary(lm) output)



# load R libraries --------------------------------------------------------
library(data.table)
library(vegan)



# load and prepare all input data -----------------------------------------------------

do_clr_externalWeighting = function(interest_matrix, core_matrix){
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}
load("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/pheno_release26.RData")
covar = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/covariates.txt")[,-8]

taxa = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/taxa.txt")

alpha = data.frame(shannon=diversity(taxa[,grep("s__",colnames(taxa))],index = "shannon"))

pathways = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/pathways_metacyc.txt")

taxa_transformed = do_clr_externalWeighting(taxa,taxa[,grep("[.]s__",colnames(taxa))])
taxa_transformed = taxa_transformed[,colSums(taxa>0)>nrow(taxa) * 0.05]

pathways_transformed = do_clr_externalWeighting(pathways,pathways)
pathways_transformed = pathways_transformed[,colSums(pathways>0)>nrow(pathways)*0.05]

KOterms = fread("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/processed/dag3_kegg_unfiltered_lvl3.csv",header=T,sep=",")
KOterms = as.data.frame(KOterms)
rownames(KOterms) = KOterms[,ncol(KOterms)]
KOterms = KOterms[,-ncol(KOterms)]
KOterms = KOterms[intersect(rownames(KOterms),rownames(taxa)),]
KOterms_transformed = do_clr_externalWeighting(KOterms,KOterms)
KOterms_transformed = KOterms_transformed[,colSums(KOterms>0)>0.05 * nrow(KOterms)]
CARD = fread("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/processed/DAG3_CARD_nofiltering.txt",header=T)
CARD = as.data.frame(CARD,stringsAsFactors=F)
rownames(CARD) = CARD[,1]
CARD = CARD[,-1]
CARD = CARD[intersect(rownames(CARD),rownames(taxa)),]
CARD_transformed = do_clr_externalWeighting(CARD,CARD)
CARD_transformed = CARD_transformed[,colSums(CARD>0)>0.05 * nrow(CARD)]

VFDB = fread("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/processed/DAG3_VFDB_VFs_nofiltering.txt",header=T)
VFDB = as.data.frame(VFDB,stringsAsFactors=F)
rownames(VFDB) = VFDB[,1]
VFDB = VFDB[,-1]
VFDB = VFDB[intersect(rownames(taxa),rownames(VFDB)),]
VFDB_transformed = do_clr_externalWeighting(VFDB,VFDB)
VFDB_transformed = VFDB_transformed[,colSums(VFDB>0)>0.05 * nrow(VFDB)]



# function for associations -----------------------------------------------

run_associations = function(phenotypes, trait, trait_group = "taxa"){
  if (trait_group == "taxa") {trait_interest = taxa_transformed[,trait,drop = F]
  } else if ((trait_group == "pathways")) {trait_interest = pathways_transformed[,trait,drop = F]
  } else if ((trait_group == "KOterms")) {trait_interest = KOterms_transformed[,trait,drop = F]
  } else if ((trait_group == "CARD")) {trait_interest = CARD_transformed[,trait,drop = F]
  } else if ((trait_group == "alpha")) {trait_interest = alpha[,trait,drop = F]
  } else if ((trait_group == "VFDB")) {trait_interest = VFDB_transformed[,trait,drop = F]
  }else stop("trait_group or trait name is wrong! trait_group should be one of 'taxa','pathways','KOterms','CARD' or 'VFDB'. Check trait names by using traitnames() function)")  
  
  covar_local = covar[rownames(trait_interest),]
  pheno_local = pheno26[rownames(trait_interest),,drop = F]
  predictors = try(data.frame(covar_local,pheno_local[,phenotypes,drop = F]))
  if (class(predictors) == "try-error") stop ("wrong phenotype names!")
  lm_null = lm(trait_interest[complete.cases(predictors),1] ~ .,data = covar_local[complete.cases(predictors),])
  lm_alter = lm(trait_interest[complete.cases(predictors),1] ~ .,data = predictors[complete.cases(predictors),])
  anova.overall = anova(lm_alter,lm_null)
  summary_object = summary(lm_alter)$coef[-c(1:(1+ncol(covar_local))),,drop =F]
  
  result = list()
  result$anova = anova.overall
  result$summary = summary_object
  rownames(result$summary) = names(lm_alter$coef[-c(1:(ncol(covar_local)+1))])
  result
}
