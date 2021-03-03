
Scripts for LASSO regression with cross-validation steps
====


*Creator: Sergio Andreu Sánchez*

*Adapted by: Arnau Vich*

*Date: 2020-201*



Description
---

This script uses L1-regularization (LASSO) (as implemented in the *gmlnet* library) to estimate the proportion of metabolite variability that can be explained by different models

For each metabolites and each of the models it performs a 5-fold cross-validation step to determine the best set of predictiors with an extra 10-fold cv to fine-tune the lambda parameter.

Input files

Input: Dataframe with samples in rows and phenotypes/data layer in columns as well all metabolites to be tested

Summary: List of metabolites that should match the column names in the Input files

Output files

Output[[1]]: percentage of variation explained per each model & metabolite

Output[[2]]: betas of each phenotype in each model-metabolites relation   



Load libraries
---

```{r}
library(tidyverse)
library (dplyr)
library(glmnet)
library(glmnetUtils)
library(doParallel)
registerDoParallel(5)

```


Load data and run functions (load functions below)
---

```{r}
my_all <- read.delim("./input_enzymes.txt", row.names=1)
my_ibd=read.table("./IBD_with_enzyme.txt", header=T)
my_cnt=read.table("./CNT_with_enzyme.txt", header=T)

metabolites=colnames(input_cc)[597:ncol(input_cc)]

Metabolite_iteration_v3(Input=my_all, Summary=metabolites) -> CC_model_enzymes
Metabolite_iteration_v3(Input=my_ibd, Summary=metabolites) -> IBD_model_enzymes
Metabolite_iteration_v3(Input=my_cnt, Summary=metabolites) -> CNT_model_enzymes

```

Functions
---

```{r}

R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}

Fit_lasso_cv = function(Dependent, Regressors,my_folds=5){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors2
  if (dim(Regressors)[[2]] < 2){
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    Param_best = tibble(Variable = colnames(Regressors), Beta= Beta)
  } else{
    #Create folds
    set.seed(2021)
    row_folds=createFolds(1:nrow(Regressors2),k=my_folds)
    Best_error = NULL
    for (i in names(row_folds)){
      #print (paste("Running fold" ,i))
      #Split training / testing
      select_rows=row_folds[[i]]
      test_reg=Regressors2[select_rows,]
      train_reg=Regressors2[-select_rows,]
      #Lasso with 10 fold 
      cvfit=cv.glmnet(Dependent ~ ., train_reg, alpha = 1, nfolds = 10, type.measure="mse",standardize=T, parallel=TRUE)
      #Predict in the test dataset and estimate the R2
      lasso_predicted = predict(cvfit, s=cvfit$lambda.min, newdata=data.matrix(test_reg[,-ncol(test_reg)]))
      tmp_cv=data.frame(real=test_reg$Dependent,pred=lasso_predicted)
      colnames(tmp_cv)=c("real", "pred")
      RSSM = sum((tmp_cv$real - tmp_cv$pred)^2)/dim(tmp_cv)[1]
      if (length(Best_error) == 0){ 
          Best_error = RSSM
          Best_model = cvfit
      }else if  (RSSM < Best_error){
        Best_error = RSSM
        Best_model = cvfit
      }
    }
    cvfit = Best_model
    Param_best <- coef(cvfit, s = "lambda.min")
    as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
    colnames(Param_best) = c("Variable","Beta")
    lasso_predicted = predict(cvfit, s=cvfit$lambda.min, newdata=data.matrix(Regressors2))
    Explained_variance = R2_calc(Regressors2$Dependent, lasso_predicted)
    if (Explained_variance < 0){ Explained_variance = 0 }
    #print(Explained_variance)
  }
  return(list(Param_best, Explained_variance))
}


Metabolite_iteration_v3 = function(Input, Summary){

  #List of all metabolites to iterate 
  Metabolites <- Summary
  #Divide phenotypes in categories
  Microbes <- colnames(Input)[grepl("sp_",colnames(Input))]
  Diet <- colnames(Input)[grepl("diet_",colnames(Input))]
  Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex", "Amount_sample_gram", "host_SmokeCurrentSmoker") #Possibly include others?
  Medication <- colnames(Input)[grepl("med_",colnames(Input))]
  #IBD <- c("IBD")
  #IBD=colnames(Input)[grepl("clinical_",colnames(Input))]
  Shannon=c("richness_Shannon_Index")
  genetics=colnames(Input)[grepl("enzyme_",colnames(Input))]
  #biomarkers=c("clinical_HBD2", "clinical_ChrA", "clinical_Calprot200")
  biomarkers=c("HBD2", "ChrA", "Calprot200")
  IBD=setdiff(colnames(Input), c(Metabolites,Microbes,Diet,Covariates, Medication, Shannon,genetics,biomarkers)) 
  #Other <- setdiff(colnames(Input), c(Diet,Microbes,Covariates, Medication, IBD, Metabolites)) #whatever was tested and is not already in diet/microbes
  
  #Output dataframe
  Variability_explained = tibble()
  #Name of the models
  Name_models <- c("Null","Complete", "Microbes", "Diet", "Medication", "Disease", "Shannon", "Genetics","Biomarkers")
  #This are the column names that are used to save the Betas, so if you need the beta of a specific varaible, should be in Variables_col vector
  Input %>% select(c(Microbes, Diet, Covariates, Medication, IBD, Shannon, genetics, biomarkers)) -> Variables_col 
  colnames(Variables_col) -> Variables_col
  hmm=0
  All_model_info ={}
    for (Metabolite in Metabolites){
    #foreach(i=1:length(Metabolites), .combine=rbind) %dopar% {
    #Metabolite=Metabolites[i]
    Logit = F
    hmm=hmm+1
    tested=paste(hmm,length(Metabolites), sep="/")
    #Get the metabolite of interest and their associations 
    #Variables_col  -> Selected_phenotypes
    Selected_phenotypes=c(Microbes, Diet, Medication, IBD, Shannon, genetics, biomarkers)
    print(paste(tested,Metabolite, sep=" "))
    #From the Input after transforming it to numeric select the Metabolite (dependent), phenotypes assocaited and Covariates. Remove all records with NA.
    Input %>% select(one_of(c(Metabolite, Selected_phenotypes,Covariates ))) %>% drop_na() -> Input_model
    #Make a vector out of dependent
    Dependent <- as.vector(as_vector(Input_model[,1]))
    #If dependent is a character, then do logistic
    if (class(Dependent[0]) == "character"){ Logit = T}
    #Prepare the different inputs for each model
    Variables_complete <- Input_model[,2:dim(Input_model)[2]]
    #Variables_clinical  <- select(Input_model, one_of(c(Clinical, Covariates)))
    Variables_IBD  <- select(Input_model, one_of(c(Covariates, IBD)))
    Variables_microbiome <- select(Input_model, one_of(c(Microbes, Covariates)))
    Variables_diet <- select(Input_model, one_of(c(Diet, Covariates)))
    Variables_medication <- select(Input_model, one_of(c(Medication, Covariates)))
    Variables_null <- select(Input_model, one_of(Covariates))
    
    Variables_Shannon <- select(Input_model, one_of(c(Shannon, Covariates)))
    Variables_genetics <- select(Input_model, one_of(c(genetics, Covariates)))
    Variables_biomarkers <- select(Input_model, one_of(c(biomarkers, Covariates)))
    
    Models <- list( Variables_null, Variables_complete,Variables_microbiome, Variables_diet,Variables_medication, Variables_IBD, Variables_Shannon, Variables_genetics, Variables_biomarkers)
    #Vector of R2s for output 1
    Variability_model = c()
    #Data.frame of variables for output 2
    tibble(Variable = Variables_col) -> Variables
    #For each model, fit lasso (normal or logistic) and save R2 and variables
    for (N in seq(1:length(Name_models)) ){
      Input_model <- Models[[N]] ; Name <- Name_models[[N]]
      #If 0 significant features add as 0 all beta and R2 and go to next
      if (dim(Input_model)[2] < 1){ 
          Variability_model = c(Variability_model,0)
          Model_summary = tibble(Variable=Variables, Beta=NA) %>% t() %>% as_tibble() %>% `colnames<-`(Metabolites)
          Model_summary %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
          next 
        }

      if (Logit == F){ Fit_lasso_cv(Dependent = Dependent, Regressors = Input_model) -> Lasso_results
      }else{ Fit_logistic_lasso(Dependent = Dependent, Regressors = Input_model) -> Lasso_results }
      #Add the betas to the Data.frame of features (only features included in that data.frame are going to get the Beta saved)   
      left_join(Variables,Lasso_results[[1]],by = "Variable") %>% t() %>% as_tibble() %>% `colnames<-`(Variables$Variable) -> Model_summary
      Model_summary[2,] %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
      All_model_info = rbind(All_model_info, Model_summary)
      #Save R2
      Variability_model = c(Variability_model,as.numeric(Lasso_results[[2]]))
    }
    #Make R2s into a data.frame with each model per column
    as_tibble(matrix(Variability_model,nrow = 1,ncol = 9)) %>% mutate(V10 = Metabolite) -> Variability_model
    colnames(Variability_model) = c(Name_models, "Metabolite")
    rbind(Variability_explained, Variability_model) -> Variability_explained
    #Variability_explained

  }
  return(list(Variability_explained, All_model_info))
}

```
