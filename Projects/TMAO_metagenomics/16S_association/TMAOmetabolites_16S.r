#Association of 16S data with TMAO-related metabolits

###Package load
library(tidyverse) #Standard packages for data load (readr), visualization (ggplot) and processing (dplyr)
library(microbiome) #Library for clr

######For Cluster version, get file path from command line
#Input = commandArgs(trailingOnly = TRUE)
#Count_table = options[1] #Output Step 1.3. Generate summary statistics from https://github.com/alexa-kur/miQTL_cookbook
#Linking_table = options[2] #Column 1: Cohort ID, Column2: Sequencing ID ; No headers
#Covariates = options[3] #Metadata should contain the following columns: ID, Sex, BMI, Age.
#Phenos = option[4] #Should contain a column called ID, the rest will be used as dependent variables in the model

######Running from desktop
Count_table =  "LLD_rarefDefault.taxonomyTable.txt.gz"
Linking_table = "coupling_pheno_16s.txt"
Covariates = "Cov_file.tsv"
Phenos = "Pheno_file.tsv"
##########

Count_table = read_tsv(Count_table) #Rows are individuals, columns are Bacteria
Linking_table = read_tsv(Linking_table, col_names=F) 
Covariates = read_tsv(Covariates)
Phenos = read_tsv(Phenos)



#Functions
Choose_replicate = function(Linking_table){
  #Random choice if Same individual's microbiome has been sequenced more than once
  to_remove = vector()
  for (Entry in unique(Linking_table$X1[duplicated(Linking_table$X1)])){
    Linking_table %>% filter(X1 == Entry) %>% sample_n(1) -> Choice
    Linking_table %>% filter(X1 == Entry) %>% filter(!X2 == Choice$X2) -> Out
    to_remove = c(to_remove, Out$X2)
  }
  Linking_table %>% filter(!X2 %in% to_remove) -> Filtered_table
  return(Filtered_table)
}

Filter_unclassified = function(Count_table){
  #NOTAX is the notation for the classifier unclassified samples
  Remove_columns = colnames(Count_table)[grepl("NOTAX", colnames(Count_table))]
  Count_table %>% select(-Remove_columns) %>% select(-rootrank.Root) -> Count_table
  return(Count_table)
}

Transformation = function(Count_table){
  SampleID = Count_table$SampleID
  Count_table %>% select(-SampleID) -> Counts
  #Transform counts
  
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
  ###
  Counts_transformed %>% mutate(SampleID = SampleID) -> Counts_transformed
  return(Counts_transformed)
}
Filter_by_abundance = function(Count_table, threshold=10){
  SampleID = Count_table$SampleID
  Count_table %>% select(-SampleID) -> abundance_matrix
  abundance_matrix = abundance_matrix[,((colSums(abundance_matrix !=0) / nrow(abundance_matrix)) *100 )>threshold]
  abundance_matrix %>% mutate(SampleID = SampleID) -> abundance_matrix
  return(abundance_matrix)
}

Fit_model = function(Count_table, Phenotype, Covariates, BMI_model=T){
    Total_results = tibble()
    
    for (Number in seq(1:dim(Count_table)[2])){
      Bug = colnames(Count_table)[Number]
      Bug_abundance = as.vector(as_vector(Count_table[,Number]))
      for (Number_metabolite in seq(1:(dim(Phenotype)[2]-1))){
          Number_metabolite = Number_metabolite + 1
          Metabolite = colnames(Phenotype)[Number_metabolite]
          Metabolite_measurement = as.vector(as_vector(Phenotype[,Number_metabolite]))
          Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement)))
          ###Fit model
          Covariates %>% mutate(Metabolite_measurement = Metabolite_measurement, Bug_abundance = Bug_abundance ) -> Model_data
          if (length(unique(Model_data$Gender) > 1)){
              if (BMI_model == T){
                MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + Gender + BMI, data = Model_data)
              }else{ MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + Gender, data = Model_data) }
              Summary = as_tibble(summary(MODEL_FIT)$coefficients)
              Beta = Summary$Estimate[2]
              pvalue = Summary$`Pr(>|t|)`[2]
              stand = Summary$`Std. Error`[2]
          }else{
            if (BMI_model == T){
              MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age + BMI, data = Model_data)
            }else{ MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance + Age, data = Model_data) }
            Summary = as_tibble(summary(MODEL_FIT)$coefficients)
            Beta = Summary$Estimate[2]
            pvalue = Summary$`Pr(>|t|)`[2]
            stand = Summary$`Std. Error`[2]
          }
          Results = as_tibble(t(c(Bug, Metabolite, Beta, stand, pvalue)))
          colnames(Results) = c("Bug", "Metabolite", "Beta", "SE", "Pvalue")
          Total_results = rbind(Total_results,Results)    
      }
   
    }
    return(Total_results)  
}

Calculate_fdr = function(Count_table, Phenotype, Covariates, Real_pvalues, BMI_model=T){
  FDR_vector = vector()
  #1st, Number of permutations
  Random_distribution  = vector()
  for (i in seq(100)){
    #Per each permutation rearrange data and fit model
    #rearrange
    sample_n(Count_table, size= nrow(Count_table)) -> Rearrange
    #Fitmodel
    Fit_model(Rearrange, Phenotype, Covariates,BMI_model) -> H0
    #Save Pvalues in the random distribution of pvalues
    Random_distribution = c(Random_distribution,as.numeric(H0$Pvalue))
  }
  #2nd, calculate FDR per each threshold
  Real_pvalues = Real_pvalues$Pvalue 
  for (Threshold in Real_pvalues){
    Threshold = as.numeric(Threshold)
    P = sum((Real_pvalues <= Threshold) *1)
    FP = sum((Random_distribution <= Threshold)*1)
    FDR = (FP/100)/P
    if (FDR >1){ FDR = 1}
    #If the FDR value of a higher Pvalue is lower, make the FDR of the lower Pvalue at least as low
    if (length(FDR_vector) > 0){ if (FDR > FDR_vector[length(FDR_vector)]){ FDR = FDR_vector[length(FDR_vector)]} }
    FDR_vector = c(FDR_vector,FDR)
  }
  return(FDR_vector)
}


###QTL Analysis
Count_table = Filter_unclassified(Count_table)
Linking_table = Choose_replicate(Linking_table)
taxonomic_levels = c("genus","order", "phylum","class","family")
set.seed(111)

for (Level in taxonomic_levels){
  print(paste0(c("Working in",Level, "taxonomic level"),collapse=" "))
  
  Count_table$SampleID -> Sample_ID
  colnames(Count_table)[grepl(Level, colnames(Count_table))] -> Level_filter
  Count_table %>% select(Level_filter) %>% mutate(SampleID =  Sample_ID) -> Taxa_table
  
  #Merging SeqID and cohort ID
  Linking_table %>% arrange(X2) %>% filter(X2 %in% Taxa_table$SampleID) -> Linking_table_taxa
  Taxa_table %>% arrange(SampleID) %>% filter(SampleID %in% Linking_table$X2) %>% 
    mutate(SampleID = Linking_table_taxa$X1) %>% arrange(SampleID) -> Taxa_table
  
  #Get same number of entries in Covariates and Phenotyes than in the 16S data
  Covariates %>% arrange(ID) %>% filter(ID %in% Taxa_table$SampleID) -> Covariates
  Phenos %>% arrange(ID) %>% filter(ID %in% Taxa_table$SampleID) -> Phenos
  
  Taxa_table %>% filter(SampleID %in% Phenos$ID) -> Taxa_table
  
  
  #Fit each of the models: Females, Males  or All
  for (S in c(unique(Covariates$Gender), 3) ){
    if (is.na(S)){next}
    print(paste0(c("Working in",S, "Gender group"),collapse=" "))
    
    if (S == 3){ Covariates -> Cov_model
    } else{ Covariates %>% filter(Gender == S) -> Cov_model }
    
    Phenos %>% filter(ID %in% Cov_model$ID) -> Pheno_model
    Taxa_table %>% filter(SampleID %in% Cov_model$ID) -> Taxa_model
    
    Transformation(Taxa_model) -> Taxa_model
    Taxa_model %>% select(-SampleID) -> Taxa_model
    
    for (BMI_model in c(T,F)){
      if (BMI_model == T){
        Output_name = paste(c("16S",Level,S, "tsv"), collapse= ".")
      }else{Output_name = paste(c("16S",Level,S,"noBMI","tsv"), collapse= ".")}
      
      Result_table = Fit_model(Taxa_model, Pheno_model, Cov_model, BMI_model)
  
      Result_table %>% arrange(desc(as.numeric(Pvalue))) -> Result_table
      Calculate_fdr(Taxa_model, Pheno_model, Cov_model, Result_table, BMI_model) -> FDR_vector
      Result_table %>% mutate(FDR = FDR_vector) -> Result_table
      Result_table %>% arrange(FDR) -> Result_table
      write_tsv(Result_table,path = Output_name)
    }

  }

  
}

