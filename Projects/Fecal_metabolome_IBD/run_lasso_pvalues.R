library (hdi)
library(tidyverse)
library (dplyr)
library(glmnet)
library(glmnetUtils)


Fit_lasso_pvalues_v2 = function(Dependent, Regressors){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors2
  
  #Perform regular LASSO
  Regressors2=Regressors2[,-1]
  cv.glmnet(Dependent ~ ., Regressors2, alpha = 1, nfolds = 10, type.measure="mse",standardize=T) -> cvfit
  Param_best <- coef(cvfit, s = "lambda.min")
  Explained_variance = cvfit$glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
  as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
  colnames(Param_best) = c("Variable","Beta_LASSO")
  Param_best=Param_best[-1,]
  
  #Boostrapped LASSO to obtain p-values
  my_models <- model.matrix(Input_model[,Metabolite]~. , Input_model)[,-1]
  my_var=Input_model[,Metabolite]
  my_models2=my_models[,-1]
  #Pvalue estimation using a bootstrap approach
  #Adjust number of boostraps and cores
  #fit.bootlasso <- boot.lasso.proj(my_models2, my_var, parallel = T, multiplecorr.method = "WY", B=1000, robust = T, ncores = 5)
  
  #Pvalue estimation using asymptotic gaussian
  #This method is faster
  fit.bootlasso <- lasso.proj(my_models2, my_var, parallel = T, multiplecorr.method = "BH", robust = T, ncores = 5)
  my_boot=data.frame("Beta_proj"=fit.bootlasso$betahat, "Beta_proj_pval"=fit.bootlasso$bhat,"P-value_proj"=fit.bootlasso$pval, "FDR_proj"=fit.bootlasso$pval.corr)
  
  #LM after LASSO (just for curiosity)
  select=Param_best$Variable[Param_best$Beta_LASSO!=0]
  Input2=Input_model[,colnames(Input_model) %in% select]
  my_lm=summary(lm(Dependent~.,Input2))
  
  my_lm=data.frame("Beta_LM"=my_lm$coefficients[,1], "P-value_LM"=my_lm$coefficients[,4])
  my_lm=my_lm[-1,]
  #Merge results
  m1=merge(Param_best,my_boot, by.x="Variable", by.y="row.names", all.x = T)
  lasso_merged=merge(m1,my_lm, by.x="Variable", by.y="row.names",all.x=T)
  lasso_merged$Metabolite=Metabolite
  return(lasso_merged)
}

#For the cluster
#library(glmnetUtils, lib.loc = "./Rlib/")
#library(scalreg, lib.loc = "./Rlib/")
#library(hdi, lib.loc = "./Rlib/")

Summary_IBD=read.table("~/Desktop/Metabolomics_v2/3.Preliminary_results/3.Regressions/within_IBD_associations_quantitative_corrected_location.txt", header = T, sep="\t", stringsAsFactors = F)
Input_IBD=read.table("~/Desktop/Metabolomics_v2/2.Input/IBD_phenos_recoded.txt", sep = "\t", header = T, row.names = 1)
Input_IBD$blood_ACPA.IgA_pos_neg=NULL
Input_IBD$blood_ANCAposneg=NULL
Input_IBD$blood_ASCAposneg=NULL
Input_IBD$blood_IgARF_pos_neg=NULL
Input_IBD$blood_IgMRF_pos_neg=NULL
Input_IBD$Calprot200=as.character(Input_IBD$Calprot200)
Input_IBD$Calprot200[Input_IBD$Calprot200=="no"]=1
Input_IBD$Calprot200[Input_IBD$Calprot200=="yes"]=2
Input_IBD$Calprot200=as.numeric(as.character(Input_IBD$Calprot200))

#Remove phenotypes with small groups or redundant (e.g diet elements vs diet groups)
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="ibd_FecalCalprotectinOver200yesno")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="meds_Biologicals")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="med_alpha_blockers")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="med_anti_androgen_oral_contraceptive")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="med_melatonine")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="med_methylphenidate")
Summary_IBD=subset(Summary_IBD, Summary_IBD$phenotype!="med_parasympathicolytic_inhaler")
Summary_IBD2=Summary_IBD[!grepl("diet_group_", Summary_IBD$phenotype),]

#Rename phenotypes
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="Shannon_Index"]="sp_Shannon_Index"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_DiseaseLocation"]="clinical_DiseaseLocation"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_EverHadStomaOrPouch"]="clinical_EverHadStomaOrPouch"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_IleocecalValveInSitu"]="clinical_IleocecalValveInSitu"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_NumberOfResectionColonic"]="clinical_NumberOfResectionColonic"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_NumberOfResectionsAny"]="clinical_NumberOfResectionsAny"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_NumberOfResectionsIleal"]="clinical_NumberOfResectionsIleal"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_NumberOfResetionsIleoCecal"]="clinical_NumberOfResetionsIleoCecal"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_ResectionAny"]="clinical_ResectionAny"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_ResectionColonicAny"]="clinical_ResectionColonicAny"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_ResectionIlealAny"]="clinical_ResectionIlealAny"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_ResectionIlealCecalAny"]="clinical_ResectionIlealCecalAny"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_ActiveDisease"]="clinical_ActiveDisease"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_Diagnosis"]="clinical_Diagnosis"
Summary_IBD2$phenotype[Summary_IBD2$phenotype=="ibd_DiseaseDurationYears.x"]="clinical_DiseaseDurationYears.x"
Summary_IBD2$FDR=p.adjust(Summary_IBD2$Pr...t.., method = "BH")
Summary_IBD2$Bonferroni=p.adjust(Summary_IBD2$Pr...t.., method = "bonferroni")



mcount=1
flag=1
Metabolites <- unique(Summary_IBD2$metabolite)
#Metabolitess=Metabolites[685:length(Metabolites)]
#Test metabolites
#Metabolites= c("cholate","glucose", "fructose")

for (Metabolite in Metabolites){
  print(paste("Testing:",Metabolite, sep=" "))
  print(paste(mcount,length(Metabolites), sep = "/"))
  Summary_IBD2 %>% filter(metabolite == Metabolite & FDR<0.05)  -> Selected_phenotypes
  if (dim(Selected_phenotypes)[[1]] < 2){
    Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex") 
    Input_IBD %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model
    Dependent <- as.vector(as_vector(Input_model[,1]))
    Input_model %>% mutate(Dependent = Dependent) -> Regressors2
    Regressors2=Regressors2[,-1]
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    lasso_merged=data.frame("Variable"= names(summary(Fitted)$coefficients[,1])[2],"Beta_LASSO"=NA,"Beta_proj"=NA,"Beta_proj_pval"=NA,"P.value_proj"=NA,"FDR_proj"=NA,"Beta_LM"=Beta,"P.value_LM"=summary(Fitted)$coefficients[2,4],"Metabolite"=Metabolite)
  } else{
    Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex") 
    Input_IBD %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model
    Dependent <- as.vector(as_vector(Input_model[,1]))
    lasso_merged=Fit_lasso_pvalues_v2(Dependent = Dependent, Regressors = Input_model)
  }
  if (flag==1){
    LASSO_models=lasso_merged
    flag=50
  }else{
    LASSO_models=rbind(LASSO_models,lasso_merged)
  }
  mcount=mcount+1
}

write.table(Summary_IBD2,"~/Desktop/Metabolomics_v2/3.Preliminary_results/4.Lasso/Univariate_results_pre_LASSO.txt", sep = "\t", quote = F)
write.table(LASSO_models,"~/Desktop/Metabolomics_v2/3.Preliminary_results/4.Lasso/IBD_LASSO.txt", sep = "\t", quote = F)


##Controls


Summary_controls=read.table("~/Desktop/Metabolomics_v2/3.Preliminary_results/3.Regressions/Controls_associations_quantitative.txt", header = T, sep="\t", stringsAsFactors = F)
Input_controls=read.table("~/Desktop/Metabolomics_v2/2.Input/Controls_phenos_recoded.txt", sep = "\t", header = T, row.names = 1)

Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="Shannon_Index.1")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="BioMK_BetaDefensin2")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="BioMK_ChromograninA")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="antrop_height")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_tricyclic_antidepressant")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_methylphenidate")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_bisphosphonates")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_calcium")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_folic_acid")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_anti_androgen_oral_contraceptive")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_alpha_blockers")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_ca_channel_blocker")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_ferrum")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_insulin")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_K_saving_diuretic")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_melatonine")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_metformin")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_opiat")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_oral_anti_diabetics")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_oral_steroid")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_other_antidepressant")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_paracetamol")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_parasympathicolytic_inhaler")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_thyrax")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_vitamin_D")
Summary_controls=subset(Summary_controls, Summary_controls$phenotype!="med_vitamin_K_antagonist")
Summary_controls2=Summary_controls[!grepl("diet_group_", Summary_controls$phenotype),]


Summary_controls2$phenotype[Summary_controls2$phenotype=="Shannon_Index"]="sp_Shannon_Index"
Summary_controls2$phenotype[Summary_controls2$phenotype=="ever_asthma"]="clinical_ever_asthma"
Summary_controls2$phenotype[Summary_controls2$phenotype=="heart_rhythm_problems"]="clinical_heart_rhythm_problems"
Summary_controls2$phenotype[Summary_controls2$phenotype=="bloodpressure_ever_high"]="clinical_bloodpressure_ever_high"
Summary_controls2$phenotype[Summary_controls2$phenotype=="migraine"]="clinical_migraine"
Summary_controls2$phenotype[Summary_controls2$phenotype=="IBS"]="clinical_IBS"
Summary_controls2$phenotype[Summary_controls2$phenotype=="anemia"]="clinical_anemia"
Summary_controls2$phenotype[Summary_controls2$phenotype=="arthrosis"]="clinical_arthrosis"
Summary_controls2$phenotype[Summary_controls2$phenotype=="neck_or_back_hernia"]="clinical_neck_or_back_hernia"
Summary_controls2$phenotype[Summary_controls2$phenotype=="eczema"]="clinical_eczema"
Summary_controls2$phenotype[Summary_controls2$phenotype=="food_allergy"]="clinical_food_allergy"
Summary_controls2$phenotype[Summary_controls2$phenotype=="dust_allergy"]="clinical_dust_allergy"
Summary_controls2$phenotype[Summary_controls2$phenotype=="animals_allergy"]="clinical_animals_allergy"
Summary_controls2$phenotype[Summary_controls2$phenotype=="pollen_allergy"]="clinical_pollen_allergy"
Summary_controls2$phenotype[Summary_controls2$phenotype=="medication_allergy"]="clinical_medication_allergy"
Summary_controls2$phenotype[Summary_controls2$phenotype=="contact_allergy"]="clinical_contact_allergy"



mcount=1
flag=1
Metabolites <- unique(Summary_controls2$metabolite)
#Metabolitess=Metabolites[408:length(Metabolites)]
#Test metabolites
#Metabolites= c("cholate","glucose", "fructose")

for (Metabolite in Metabolites){
  print(paste("Testing:",Metabolite, sep=" "))
  print(paste(mcount,length(Metabolites), sep = "/"))
  Summary_controls2 %>% filter(metabolite == Metabolite & FDR<0.05)  -> Selected_phenotypes
  if (dim(Selected_phenotypes)[[1]] < 2){
    Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex") 
    Input_controls %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model
    Dependent <- as.vector(as_vector(Input_model[,1]))
    Input_model %>% mutate(Dependent = Dependent) -> Regressors2
    Regressors2=Regressors2[,-1]
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    lasso_merged=data.frame("Variable"= names(summary(Fitted)$coefficients[,1])[2],"Beta_LASSO"=NA,"Beta_proj"=NA,"Beta_proj_pval"=NA,"P.value_proj"=NA,"FDR_proj"=NA,"Beta_LM"=Beta,"P.value_LM"=summary(Fitted)$coefficients[2,4],"Metabolite"=Metabolite)
  } else{
    Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex") 
    Input_controls %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model
    Dependent <- as.vector(as_vector(Input_model[,1]))
    lasso_merged=Fit_lasso_pvalues_v2(Dependent = Dependent, Regressors = Input_model)
  }
  if (flag==1){
    LASSO_models=lasso_merged
    flag=50
  }else{
    LASSO_models=rbind(LASSO_models,lasso_merged)
  }
  mcount=mcount+1
}




write.table(Summary_controls2,"~/Desktop/Metabolomics_v2/3.Preliminary_results/4.Lasso/Univariate_results_pre_LASSO_controls.txt", sep = "\t", quote = F)
write.table(LASSO_models,"~/Desktop/Metabolomics_v2/3.Preliminary_results/4.Lasso/controls_LASSO.txt", sep = "\t", quote = F)

























cv.glmnet(Dependent ~ ., Regressors2, alpha = 1, nfolds = 10, type.measure="mse",standardize=T, type.multinomial="grouped") -> cvfit2

Param_best <- coef(cvfit2, s = 0.014/369)
as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
colnames(Param_best) = c("Variable","Beta")

betas=Param_best$Beta[-1]

out=fixedLassoInf(as.matrix(Regressors3),x,betas, lambda = 0.014)

output <- cbind(Param_best[-1,][order(out$vars),],out$pv)

fit.lasso <- boot.lasso.proj(Regressors3, x)

Input_model %>% mutate(Dependent = Dependent) -> Regressors2

Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex")

Metabolite="cholate"
Summary_IBD2 %>% filter(metabolite == Metabolite & FDR<0.05)  -> Selected_phenotypes
Input_IBD %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model

#if (Logit == F){ Fit_lasso(Dependent = Dependent, Regressors = Input_model) -> Lasso_results

Dependent <- as.vector(as_vector(Input_model[,1]))

my_models <- model.matrix(glucose~. , Input_model)[,-1]
my_var=Input_model$glucose

fit.lasso <- boot.lasso.proj(my_models, my_var, parallel = T, multiplecorr.method = "fdr", B=10000, robust = T)







