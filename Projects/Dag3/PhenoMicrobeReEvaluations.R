##version 1.0.1 from 7 May 2020.
##Written by Alex Kur

# Instruction manual (simple)
# 1. go to Gearshift server
# 2. load RPlus package:
#   $ module load RPlus
# 3. open R  
#   $ R
# 4. Load script as a source 
#   > source("PhenoMicrobeReEvaluation_v.1.0.0.R")
# 5. two main functions:
#   > run.association.pathway(pathway_name,phenotypes,model.type = "complete")
#   > run.association.taxa(taxon_name,phenotypes,model.type = "complete")
# Pathway_name and taxon_name should be single element characters of your taxon/pathway of interest.
#
# Phenotypes could be vector of one or few phenotype names, 
# for example:
#   > c("MED.BLOOD.Lymphocytes_10E9perL","MED.BLOOD.Bpressure.systolic.mm.Hg")
#
# model.type should be "complete", "mu", "nu" and "sigma". Mainly, complete and mu should be in frequent use. 




#loading data for analysis
message("Start loading data")
library(gamlss)

covar_path = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/covariates_path.txt")
pathways = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/pathways_metacyc.txt")
pathways = pathways[,colSums(pathways>0)>nrow(pathways) * 0.05]

covar_tax = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/covariates.txt")
taxa = read.table("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/taxa.txt")
taxa = taxa[,colSums(taxa>0)>nrow(taxa) * 0.05]

load("/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov/Association.analysis.input/phenos.RData")

if(all(unlist(lapply(c("taxa","pathways","pheno","covar_path","covar_tax"),exists)))) message("All data succesfully loaded")

#functions for public use
show.pathways = function() print(colnames(pathways))
show.taxa = function() print(colnames(taxa))
show.phenotypes = function() print(colnames(pheno))

#main function for phenotypes
run.association.pathway = function(pathway_name,phenotypes,model.type = "complete"){
  if(!(pathway_name %in% colnames(pathways))) stop ("check your pathway names. Stopped")
  if(!all(phenotypes %in% colnames(pheno))) stop ("check your phenotype names. Stopped")
  
  #preparing data for analysis
  predictors = data.frame(covar_path,pheno[,phenotypes,drop = F])
  good.samples = complete.cases(predictors)
  pathway = pathways[good.samples,pathway_name]
  cleaned_data = predictors[good.samples,]
  formula.str = as.formula(paste("pathway ~ ",
                      paste(collapse="+",c("Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7",phenotypes))))
  formula.shrinked = as.formula(paste(collapse="+",c("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7",phenotypes)))
  
  if(model.type == "complete"){
    bezi.model = gamlss(formula.str,
                        sigma.formula = formula.shrinked,
                        nu.formula = formula.shrinked,
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = c(rep("mu",length(phenotypes)),rep("sigma",length(phenotypes)),rep("nu",length(phenotypes))),
                                      summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = sub("[.][0-9]$","",rownames(output)),pathway = pathway_name,output,row.names=NULL)
  } else if (model.type == "mu"){
    bezi.model = gamlss(formula.str,
                        sigma.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        nu.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = "mu",summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = rownames(output),pathway = pathway_name,output, row.names=NULL)
  } else if (model.type == "nu"){
    bezi.model = gamlss(as.formula("pathway ~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                                     DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        sigma.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        nu.formula = formula.shrinked,
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = "nu",summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = rownames(output),pathway = pathway_name,output, row.names=NULL)
  }
    else if (model.type == "sigma"){
      bezi.model = gamlss(as.formula("pathway ~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                                     DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                          sigma.formula = formula.shrinked,
                          nu.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                          data = cleaned_data,
                          family = "BEZI",
                          method = RS(50))
      summary.object = summary(bezi.model)
      output = data.frame(parameter = "sigma",summary.object[rownames(summary.object) %in% phenotypes,])
      output = data.frame(phenotype = rownames(output),pathway = pathway_name,output, row.names=NULL)
    } else {stop("incorrect parameter selected. Should be mu,nu or sigma")}
  output
}
  
  
run.association.taxa = function(taxon_name,phenotypes,model.type = "complete"){
  if(!(taxon_name %in% colnames(taxa))) stop ("check your taxon name. Stopped")
  if(!all(phenotypes %in% colnames(pheno))) stop ("check your phenotype names. Stopped")
  
  #preparing data for analysis
  predictors = data.frame(covar_tax,pheno[,phenotypes,drop = FALSE])

  good.samples = complete.cases(predictors)
  taxon = taxa[good.samples,taxon_name]
  cleaned_data = predictors[good.samples,]
  formula.str = as.formula(paste("taxon ~ ",
                                 paste(collapse="+",c("Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7",phenotypes))))
  formula.shrinked = as.formula(paste(collapse="+",c("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7",phenotypes)))
  
  if(model.type == "complete"){
    bezi.model = gamlss(formula.str,
                        sigma.formula = formula.shrinked,
                        nu.formula = formula.shrinked,
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = c(rep("mu",length(phenotypes)),rep("sigma",length(phenotypes)),rep("nu",length(phenotypes))),
                        summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = sub("[.][0-9]$","",rownames(output)),taxon = taxon_name,output,row.names=NULL)
  } else if (model.type == "mu"){
    bezi.model = gamlss(formula.str,
                        sigma.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        nu.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = "mu",summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = rownames(output),taxon = taxon_name,output, row.names=NULL)
  } else if (model.type == "nu"){
    bezi.model = gamlss(as.formula("taxon ~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                                     DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        sigma.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        nu.formula = formula.shrinked,
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = "nu",summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = rownames(output),taxon = taxon_name,output, row.names=NULL)
  }
  else if (model.type == "sigma"){
    bezi.model = gamlss(as.formula("taxon ~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                                     DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        sigma.formula = formula.shrinked,
                        nu.formula = as.formula("~ Age+Sex+Season.spring+Season.summer+Season.unknown+Season.winter+
                              DNA.con+GeoMean+batch.2+batch.3+batch.4+batch.5+batch.6+batch.7"),
                        data = cleaned_data,
                        family = "BEZI",
                        method = RS(50))
    summary.object = summary(bezi.model)
    output = data.frame(parameter = "sigma",summary.object[rownames(summary.object) %in% phenotypes,])
    output = data.frame(phenotype = rownames(output),taxon = taxon_name,output, row.names=NULL)
  } else {stop("incorrect parameter selected. Should be mu,nu or sigma")}
  output
}


