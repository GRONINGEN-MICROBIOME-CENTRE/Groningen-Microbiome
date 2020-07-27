library(tidyverse)
library(glmnet)
library(glmnetUtils) #this one for an easier syntax on glmnet, not using matrices and makes the dummy coding (the matrices of glmnet only accept numbers, not factors)
setwd("/Users/sergio/Documents/PhD/IBD_metabolomics_Arnau2020/6.Multivariate")

###Comments
# I am not sure if the logistic is calculating R2 as it should. If there is only one variable I do logistic instead of lasso. And the R2 calculation compares 1/0 to the Risk.
# The main function (Metabolite iteration) right now is iterating only thorugh metabolites included in the Summary file. It only takes features under FDR 0.05 and covariates
# The covariates are not including 'run_day_cat'. Basically because it is multi-factorial and lasso recodes it with dummy variables and I would have to hardcode the variable name first into
#Variables_col in order to save the betas
# Similarly the other multi-factorials phenotypes probably dont have the betas recorded. The reason is, that if you have a variable (Origin) and two levels: Oral/Gut, Lasso calls the beta OriginOral and OriginGut (gives two).
#as I am appending the betas to a data.frame with the variable names saved in Variables_col, as the names do not match they wont be in the final dataset. My idea was to hardcode them in line 88: c(colnames(Variables_col), OriginOral, OriginGut) -> Variables_col







###############
###Functions###
###############
R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}
Fit_lasso = function(Dependent, Regressors){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors2
  if (dim(Regressors)[[2]] < 2){
    lm(Dependent ~ ., Regressors2) -> Fitted
    Explained_variance = summary(Fitted)$r.squared
    Beta = summary(Fitted)$coefficients[2]
    Param_best = tibble(Variable = colnames(Regressors), Beta= Beta)
  } else{
    cv.glmnet(Dependent ~ ., Regressors2, alpha = 1, nfolds = 10, type.measure="mse",standardize=T) -> cvfit
    Param_best <- coef(cvfit, s = "lambda.min")
    Explained_variance = cvfit$glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
    as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
    colnames(Param_best) = c("Variable","Beta")
  }
  return(list(Param_best, Explained_variance))
  
}
Fit_logistic_lasso = function(Dependent, Regressors){
  Regressors %>% mutate(Dependent = Dependent) -> Regressors2
  if (dim(Regressors)[[2]] < 2){
    glm(as.factor(Dependent) ~ ., Regressors2, family=binomial(link="logit")) -> Fitted
    Explained_variance = R2_calc(as.factor(Dependent),predict(Fitted,type="response") )
    Beta = summary(Fitted)$coefficients[2]
    Param_best = tibble(Varaible=colnames(Regressors), Beta=Beta)
  } else{
    cv.glmnet(as.factor(Dependent) ~ ., Regressors2, alpha = 1, nfolds = 10, family = "binomial", type.measure = "class",standardize=T,) -> cvfit
    Param_best <- coef(cvfit, s = "lambda.min")
    Explained_variance = cvfit$glmnet.fit$dev.ratio[which(cvfit$glmnet.fit$lambda == cvfit$lambda.min)]
    as.data.frame(as.matrix(Param_best)) %>% rownames_to_column() %>% as_tibble() -> Param_best
    colnames(Param_best) = c("Variable","Beta")
  }
  return(list(Param_best, Explained_variance))
  
  
}
  


Metabolite_iteration = function(Input, Summary){
  ##there are some duplicated records with a _1 at the end, remove them from Input and from summary (in summary instead of _ ther is a .)
  colnames(Input)[grepl("[a-z]_1$",colnames(Input))] -> Remove
  Input %>% select(-Remove) -> Input ; Summary %>% select(-one_of(str_replace(Remove,"_","."))) -> Summary
  #List of all metabolites to iterate 
  Metabolites <- unique(Summary$metabolite)
  #Divide phenotypes in categories
  Microbes <- c(unique(Summary$phenotype)[grepl("sp_",unique(Summary$phenotype))],  "Shannon_Index")
  Diet <- unique(Summary$phenotype)[grepl("diet_",unique(Summary$phenotype))]
  Covariates <- c("LC.COLUMN","clinical_BowelMovementADayDef","metabolon_Month_in_freezer","host_Age","host_BMI","host_Sex") #Possibly include others?
  Other <- setdiff(unique(Summary$phenotype), c(Diet,Microbes)) #whatever was tested and is not already in diet/microbes
  
  #Output dataframe
  Variability_explained = tibble()
  #Name of the models
  Name_models <- c("Complete", "Clinical", "Microbes", "Diet", "Null")
  #Make all variables in numberic/character, currently the function does not accept highly multifactorial variables. 2 level factors become numeric.
  select(Input, -one_of(c("Row.names", "run_day_cat"))) -> Variables_Input
  apply(Variables_Input,2, FUN=Make_numeric) %>% as_tibble() -> Variables_Input #Check Make_numeric function
  #This are the column names that are used to save the Betas, so if you need the beta of a specific varaible, should be in Variables_col vector
  Variables_Input %>% select(c(Microbes, Diet, Covariates, Other)) -> Variables_col 
  colnames(Variables_col) -> Variables_col
  
  All_model_info ={}
    for (Metabolite in Metabolites){
    Logit = F
    #Get the metabolite of interest and their associations 
    Summary %>% filter(metabolite == Metabolite & FDR<0.05)  -> Selected_phenotypes
    if (dim(Selected_phenotypes)[1] == 0 ){ next }  #if no associations, go to next metabolite
    print(Metabolite)
    #From the Input after transforming it to numeric select the Metabolite (dependent), phenotypes assocaited and Covariates. Remove all records with NA.
    Variables_Input %>% select(one_of(c(Metabolite, unique(Selected_phenotypes$phenotype),Covariates ))) %>% drop_na() -> Input_model
    #Make a vector out of dependent
    Dependent <- as.vector(as_vector(Input_model[,1]))
    #If dependent is a character, then do logistic
    if (class(Dependent[0]) == "character"){ Logit = T}
    #Prepare the different inputs for each model
    Variables_complete <- Input_model[,2:dim(Input_model)[2]]
    Variables_clinical  <- select(Input_model, one_of(c(Other, Covariates)))
    Variables_microbiome <- select(Input_model, one_of(c(Microbes, Covariates)))
    Variables_diet <- select(Input_model, one_of(c(Diet, Covariates)))
    Variables_null <- select(Input_model, one_of(Covariates))
    
    Models <- list( Variables_complete, Variables_clinical, Variables_microbiome, Variables_diet,Variables_null)
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
          Model_summary = tibble(Variable=Variables, Beta=NA) %>% t() %>% as_tibble() %>% `colnames<-`(Summary$metabolite)
          Model_summary %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
          next 
        }
      if (Logit == F){ Fit_lasso(Dependent = Dependent, Regressors = Input_model) -> Lasso_results
      }else{ Fit_logistic_lasso(Dependent = Dependent, Regressors = Input_model) -> Lasso_results }
      #Add the betas to the Data.frame of features (only features included in that data.frame are going to get the Beta saved)   
      left_join(Variables,Lasso_results[[1]],by = "Variable") %>% t() %>% as_tibble() %>% `colnames<-`(Variables$Variable) -> Model_summary
      Model_summary[2,] %>% mutate(Model = Name, Metabolite = Metabolite) -> Model_summary
      All_model_info = rbind(All_model_info, Model_summary)
      #Save R2
      Variability_model = c(Variability_model,as.numeric(Lasso_results[[2]]))
    }
    #Make R2s into a data.frame with each model per column
    as_tibble(matrix(Variability_model,nrow = 1,ncol = 5)) %>% mutate(V6 = Metabolite) -> Variability_model
    colnames(Variability_model) = c(Name_models, "Metabolite")
    rbind(Variability_explained, Variability_model) -> Variability_explained
  }
  return(list(Variability_explained, All_model_info))
}

Make_numeric = function(x){
  if(length(unique(x)) == 2 ){
    x = as.numeric(as.factor(x))-1
  } else if (length(unique(x)) < 5){
    x = as.factor(x)
  }else{
    x=as.numeric(x)
  }
  return(x)
}



###############
####Data#######
###############
read_tsv("Input/Controls/controls_phenos2.txt") -> Input_controls
read_tsv("Input/IBD/ibd_phenos2.txt") -> Input_IBD

##Open summary statistics used for filtering
read_tsv("Output/Controls_associations_quantitative.txt") -> Summary_controls
Summary_controls %>% filter(!grepl("\\.1",phenotype)) ->Summary_controls
read_tsv("Output/within_IBD_associations_quantitative.txt") -> Summary_IBD


#Calculate R2 per metabolite based on Lasso
#Output: 1. R2 per metabolute  ; 2. Table of variables chosen per metabolite and per model

Metabolite_iteration(Input = Input_IBD, Summary= Summary_IBD) -> Output_model
Output_model[[1]] %>%  arrange(desc(Complete)) -> IBD_variance_explained


Metabolite_iteration(Input = Input_controls, Summary= Summary_controls) -> Output_model2
Output_model2[[1]] %>%  arrange(desc(Complete)) -> Controls_variance_explained

#Plot of distributions, use patchwork to add the ggplots together

library(patchwork)
IBD_variance_explained %>% gather(Model,R2, c("Complete","Clinical", "Microbes", "Diet", "Null")) %>% mutate(Model = fct_relevel(Model,c("Complete","Microbes", "Clinical", "Diet", "Null"))) %>% 
  ggplot(aes(x=R2, fill=Model)) + geom_histogram(position="identity") + theme_bw() + ggtitle("IBD R2") -> IBD_distribution
Controls_variance_explained %>% gather(Model,R2, c("Complete","Clinical", "Microbes", "Diet", "Null")) %>% mutate(Model = fct_relevel(Model,c("Complete","Microbes", "Clinical", "Diet", "Null"))) %>% 
  ggplot(aes(x=R2, fill=Model)) + geom_histogram(position="identity") + theme_bw() + ggtitle("Control R2") -> Control_distribution
IBD_distribution + Control_distribution


#Compare output 2, only Complete models by doing a correlation ; option 2, make 1 and 0 and correlate that 


Output_model[[2]] %>% filter(Model=="Complete") %>% select(one_of(colnames(Output_model2[[2]]))) %>% filter(Metabolite %in%  Output_model2[[2]]$Metabolite) -> Complete_params_IBD
Output_model2[[2]] %>% filter(Model=="Complete") %>% select(one_of(colnames(Complete_params_IBD))) %>% filter(Metabolite %in% Complete_params_IBD$Metabolite) -> Complete_params_controls
tibble_correlations = tibble()
for(Meta in unique(Complete_params_IBD$Metabolite)){
  Complete_params_IBD %>% filter(Metabolite == Meta) %>% select(-c(Model, Metabolite)) %>% as_vector() %>% as.vector() -> IBD_meta
  IBD_meta[is.na(IBD_meta)] = 0 ; as.numeric(IBD_meta) -> IBD_meta
  
  Complete_params_controls %>% filter(Metabolite == Meta) %>% select(-c(Model, Metabolite)) %>% as_vector() %>% as.vector() -> Controls_meta
  Controls_meta[is.na(Controls_meta)] = 0 ; as.numeric(Controls_meta) -> Controls_meta
  
  cor.test(Controls_meta,IBD_meta) -> Correlations
  tibble(Metabolite=Meta , Correlation=Correlations$estimate) -> Meta_cor
  rbind(tibble_correlations,Meta_cor) -> tibble_correlations
}
tibble_correlations %>% arrange(desc(Correlation))
tibble_correlations %>% ggplot(aes(x=Correlation)) + geom_histogram() + theme_bw() 

#Check the strongest Correlations

head(arrange(tibble_correlations,desc(abs(Correlation))), 9) -> Strongest_corr
for (Meta in Strongest_corr$Metabolite){
  IBD_variance_explained %>% filter(Metabolite==Meta) %>% print()
  Controls_variance_explained %>% filter(Metabolite==Meta) %>% print()
}
