###Package load
library(tidyverse) #Standard packages for data load (readr), visualization (ggplot) and processing (dplyr)
library(patchwork)

######Running from desktop
setwd("~/Resilio Sync/Transfer/PhD/TMAO_16S")
Covariates = "Cov_file.tsv"
Phenos = "Pheno_LLD.tsv"
Diet = "20150722_Diet__1135patients.txt"
##########

Covariates = read_tsv(Covariates)
Phenos = read_tsv(Phenos)
colnames(Phenos) = c("ID","TMAO", "Choline", "Betaine", "Deoxycarnitine", "Carnitine")
Phenos %>% mutate(TMAO.Choline = TMAO/Choline , TMAO.Betaine = TMAO/Betaine , 
                  TMAO.Butyrobetaine = TMAO/Deoxycarnitine, TMAO.Carnitine = TMAO/Carnitine,
                  Butyrobetain.Carnitine =  Deoxycarnitine/Carnitine ) -> Phenos
Diet = read_tsv(Diet)
Diet %>% mutate(ID = LL_ID_corrected) %>% dplyr::select(-LL_ID_corrected) -> Diet

read_csv("Diet/LLS_IOP1-FFQ-Codebook_Foodgroups_20210628.csv") -> Diet_LLS
#Diet %>% dplyr::select( colnames(Diet)[grepl("wine", colnames(Diet))] )
#Diet %>% filter(LLD_present == "TRUE") -> Items
Diet %>% dplyr::select(c("ID", filter(Diet_LLS, LLD_present=="TRUE")$X6)) -> Diet



#Normalization. Metabolite measurement and Diet measurements are both Rank normal transformed
Phenos %>% gather(Metabolite, Concentration, 2:11)  %>% ggplot(aes(x=Concentration,fill=Metabolite)) +geom_density(alpha=0.3) + facet_wrap(~Metabolite, scales = "free")+ theme_bw()
normalization = function(Metabolite_measurement){ Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement))) }
apply(dplyr::select(Phenos, -ID), 2, FUN=normalization) %>% as_tibble() %>% mutate(ID=Phenos$ID) -> Phenos
Phenos %>% gather(Metabolite, Concentration, 1:10)  %>% ggplot(aes(x=Concentration,fill=Metabolite)) + geom_density(alpha=0.3) + facet_wrap(~Metabolite, scales = "free")+ theme_bw()

apply(dplyr::select(Diet, -ID), 2, FUN=normalization) %>% as_tibble() %>% mutate(ID=Diet$ID) -> Diet


#####Partial correlation on rank transformed food groups
library(ppcor)
result_partial = tibble()

for (d in colnames(Diet)){
    if(d ==  "ID"){ next }
  for (m in colnames(Phenos)){
    if(m ==  "ID"){ next }
      D = dplyr::select(Diet, c("ID", d)) ; M = dplyr::select(Phenos, c("ID", m))
      left_join(left_join(Covariates, D, by="ID"), M, by="ID") %>% drop_na() %>% as.data.frame() -> Model_input
      res <- pcor.test(Model_input[5], Model_input[6], Model_input[2:4]) 
      
      as.data.frame(res) %>% as_tibble() %>% mutate(food = d, .before=1) %>% mutate(metabolite = m, .before=1) -> res
      rbind(result_partial, res) -> result_partial
  }
}

write_tsv(result_partial,"Diet/Diet_Stats_partialCorr_LLD.tsv")





Previous_analyses = function(Diet, Phenos){
  #Diet should not be normalized prior to this
    Summary_results = tibble()
    
    scale(select(Diet, -ID)) %>% as_tibble() %>% mutate(ID = Diet$ID) -> Diet_scaled
    
    for (Interest in colnames(Diet_scaled) ){
      if (Interest == "ID"){ next }
    
      left_join( Phenos, Covariates, by= "ID") -> D_test
      left_join( D_test, select(Diet_scaled, c("ID", Interest)), by= "ID") -> D_test
    
      for (Metabolite in  colnames(Phenos)){
        if (Metabolite == "ID"){ next }
        Formula = paste(c(Metabolite, " ~ ", Interest, " + Age + Gender + BMI"), collapse="")
        lm(as.formula(Formula), D_test)  -> M1
        as.data.frame(summary(M1)$coefficients)[Interest,] %>% rownames_to_column("Item") %>% as_tibble() %>% mutate(Metabolite = Metabolite) -> R_test
        Summary_results = rbind(Summary_results, R_test)
      }
    }
    Summary_results %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>% arrange(FDR) -> Summary_results
    
    write_tsv(x = Summary_results,path =  "Diet/LLD_diet_stats.tsv")
    Summary_results %>% ggplot(aes(x=Estimate, y=-log10(`Pr(>|t|)`))) + geom_point(aes(col= FDR<0.05)) + theme_bw() + facet_wrap(~Metabolite, scales = "free")
    Summary_results %>% mutate(P=`Pr(>|t|)`) %>% select(Item, Metabolite, Estimate, P) %>% gt::gt()
    
    tibble_diet = tibble()
    for (Product in colnames(Diet)){
      if(Product == "ID"){ next }
      Diet %>% select(Product) %>% as.vector() %>% as_vector() -> V
      NA_proportion = sum(is.na(V))
      V = V[! is.na(V)]
      Not_NA = length(V)
      Mean = mean(V)
      Median = median(V)
      SD = sd(V)
      Max = max(V)
      Min = min(V)
      SU = tibble(Fodd_item = Product, NAs = NA_proportion, Total_noNA = Not_NA, Mean =Mean, Median=Median, SD=SD, Max=Max, Min =Min )
      rbind(tibble_diet, SU) -> tibble_diet
    }
    write_tsv(x = tibble_diet, path = "Diet/Diet_stats.tsv")
    
    
    ##########################
    ###Test interactions #####
    ##########################
    
    ###Prepare microbial data
    Results_meta = read_tsv("Metanalysis/Summary_stats_All.tsv")
    Microbes = read_tsv("LLD_rarefDefault.taxonomyTable.txt.gz")
    Linking_table = read_tsv("coupling_pheno_16s.txt", col_names=F) 
    Transformation = function(Count_table){
      SampleID = Count_table$ID
      Count_table %>% select(-ID) -> Counts
      #Transform counts
      
      colnames(Counts) -> N
      Counts %>% t() -> Check
      Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
      ###
      Counts_transformed %>% mutate(ID = SampleID) -> Counts_transformed
      return(Counts_transformed)
    }
    Filter_unclassified = function(Count_table){
      #NOTAX is the notation for the classifier unclassified samples
      Remove_columns = colnames(Count_table)[grepl("NOTAX", colnames(Count_table))]
      Count_table %>% select(-Remove_columns) %>% select(-rootrank.Root) -> Count_table
      return(Count_table)
    }
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


    Microbes = Filter_unclassified(Microbes)
    Linking_table = Choose_replicate(Linking_table)
    Results_meta %>% filter(meta_FDR < 0.05) -> Filter_bugs
    
    
    Linking_table %>% arrange(X2) %>% filter(X2 %in% Microbes$SampleID) -> Linking_table_taxa
    Microbes %>% arrange(SampleID) %>% filter(SampleID %in% Linking_table$X2) %>% 
      mutate(SampleID = Linking_table_taxa$X1) %>% arrange(SampleID) -> Microbes
    ####
    Diet_results = read_tsv("Diet/LLD_diet_stats.tsv")
    Diet_results %>% filter(`Pr(>|t|)` < 0.05) -> Diet_results
    ####
    
    taxonomic_levels = c("genus","order", "phylum","class","family")
    set.seed(111)
    Total_results = tibble()
    for (Level in taxonomic_levels){
      colnames(Microbes)[grepl(Level, colnames(Microbes))] -> Level_filter
      Microbes %>% select(Level_filter) %>% mutate(ID = Microbes$SampleID) %>% arrange(ID)  -> Microbes_level
      
      Covariates %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Cov
      Phenos %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Phen
      Diet %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) ->  Di
      
      Transformation(Microbes_level) -> Microbes_level
      
      Microbes_level %>% select(c(ID, one_of(Filter_bugs$Bug))) -> Microbes_level 
      Reduce(intersect, list(Cov$ID, Phen$ID, Di$ID, Microbes_level$ID)) -> Keep
      Microbes_level %>% filter(ID %in% Keep) -> Microbes_level
      Cov %>% filter(ID %in% Keep) -> Cov ; Phen %>% filter(ID %in% Keep) -> Phen ; Di %>% filter(ID %in% Keep) -> Di
      
      for (Number in seq(1:dim(Microbes_level)[2])){
        Bug = colnames(Microbes_level)[Number]
        if (Bug == "ID"){ next }
        Bug_abundance = as.vector(as_vector(Microbes_level[,Number]))
        Results_meta %>% filter(meta_FDR <0.05) %>% filter(Bug == colnames(Microbes_level)[Number]) -> Associated
        for (Metabolite_n in unique(Associated$Metabolite) ){
          Phen %>% select(Metabolite_n) -> Phen_2
          Metabolite_name = colnames(Phen_2)[1]
          Metabolite_measurement = as.vector(as_vector(Phen_2))
          Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement)))
          
          Diet_results %>% filter(Metabolite == Metabolite_name) -> food_items
          
          for (Food in unique(food_items$Item)){
            Di %>% select(Food) -> Diet_item
            Diet_name = colnames(Diet_item)[1]
            Diet_ab = as.vector(as_vector(Diet_item))
            ###Fit model
            Cov %>% mutate(Metabolite_measurement = Metabolite_measurement, Bug_abundance = Bug_abundance, Food_item=Diet_ab ) -> Model_data
          
            MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance*Food_item + Age + Gender + BMI, data = Model_data)
            Summary = as.data.frame(summary(MODEL_FIT)$coefficients)
          
            Beta = Summary$Estimate[7]
            pvalue = Summary$`Pr(>|t|)`[7]
            stand = Summary$`Std. Error`[7]
            N_samples = length(MODEL_FIT$residuals)
            
            Results = as_tibble(t(c(Bug, Diet_name, Metabolite_name, Beta, stand, pvalue, N_samples)))
            colnames(Results) = c("Bug","Food_item", "Metabolite", "Beta", "SE", "Pvalue", "N_samples")
            Total_results = rbind(Total_results,Results)    
          }  
        }
      }
    }  
    
    
    Total_results %>% mutate(Model = 1) %>% mutate(Beta = as.numeric(Beta), Pvalue=as.numeric(Pvalue)) -> Total_results
    Total_results %>% ggplot(aes(x=Beta, y=-log10(Pvalue)))+geom_point()+theme_bw()
    ###Prevalence analysis
    Total_results_prev = tibble()
    for (Level in taxonomic_levels){
      colnames(Microbes)[grepl(Level, colnames(Microbes))] -> Level_filter
      Microbes %>% select(Level_filter) %>% mutate(ID = Microbes$SampleID) %>% arrange(ID)  -> Microbes_level
      
      Covariates %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Cov
      Phenos %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Phen
      Diet %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) ->  Di
      
      Microbes_level$ID -> ID_save
      Microbes_level %>% select(-ID) -> Microbes_level ; Microbes_level[Microbes_level != 0 ] = 1 ; Microbes_level %>% mutate(ID=ID_save) -> Microbes_level
      
      Microbes_level %>% select(c(ID, one_of(Filter_bugs$Bug))) -> Microbes_level 
      Reduce(intersect, list(Cov$ID, Phen$ID, Di$ID, Microbes_level$ID)) -> Keep
      Microbes_level %>% filter(ID %in% Keep) -> Microbes_level
      Cov %>% filter(ID %in% Keep) -> Cov ; Phen %>% filter(ID %in% Keep) -> Phen ; Di %>% filter(ID %in% Keep) -> Di
      
      for (Number in seq(1:dim(Microbes_level)[2])){
        Bug = colnames(Microbes_level)[Number]
        if (Bug == "ID"){ next }
        Bug_abundance = as.vector(as_vector(Microbes_level[,Number]))
        Results_meta %>% filter(meta_FDR <0.05) %>% filter(Bug == colnames(Microbes_level)[Number]) -> Associated
        for (Metabolite_n in unique(Associated$Metabolite) ){
          Phen %>% select(Metabolite_n) -> Phen_2
          Metabolite_name = colnames(Phen_2)[1]
          Metabolite_measurement = as.vector(as_vector(Phen_2))
          Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement)))
          
          Diet_results %>% filter(Metabolite == Metabolite_name) -> food_items
          
          for (Food in unique(food_items$Item)){
            Di %>% select(Food) -> Diet_item
            Diet_name = colnames(Diet_item)[1]
            Diet_ab = as.vector(as_vector(Diet_item))
            ###Fit model
            Cov %>% mutate(Metabolite_measurement = Metabolite_measurement, Bug_abundance = Bug_abundance, Food_item=Diet_ab ) -> Model_data
            
            MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance*Food_item + Age + Gender + BMI, data = Model_data)
            Summary = as.data.frame(summary(MODEL_FIT)$coefficients)
            
            Beta = Summary$Estimate[7]
            pvalue = Summary$`Pr(>|t|)`[7]
            stand = Summary$`Std. Error`[7]
            N_samples = length(MODEL_FIT$residuals)
            
            Results = as_tibble(t(c(Bug, Diet_name, Metabolite_name, Beta, stand, pvalue, N_samples)))
            colnames(Results) = c("Bug","Food_item", "Metabolite", "Beta", "SE", "Pvalue", "N_samples")
            Total_results_prev = rbind(Total_results_prev,Results)    
          }  
        }
      }
    }  
    
    
    Total_results_prev %>% mutate(Model = 2) %>% mutate(Beta = as.numeric(Beta), Pvalue=as.numeric(Pvalue)) -> Total_results_prev
    Total_results_prev %>% ggplot(aes(x=Beta, y=-log10(Pvalue)))+geom_point()+theme_bw()
    
    
    ###COmbined-species analysis - we can also cluster in quantiles, but results should not be too far off.
    
    Total_results_cluster = tibble()
    for (Level in taxonomic_levels){
      colnames(Microbes)[grepl(Level, colnames(Microbes))] -> Level_filter
      Microbes %>% select(Level_filter) %>% mutate(ID = Microbes$SampleID) %>% arrange(ID)  -> Microbes_level
      
      Covariates %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Cov
      Phenos %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) -> Phen
      Diet %>% filter(ID %in% Microbes_level$ID) %>% arrange(ID) ->  Di
      
      #Transformation(Microbes_level) -> Microbes_level
      
      Microbes_level %>% select(c(ID, one_of(Filter_bugs$Bug))) -> Microbes_level 
      Reduce(intersect, list(Cov$ID, Phen$ID, Di$ID, Microbes_level$ID)) -> Keep
      Microbes_level %>% filter(ID %in% Keep) -> Microbes_level
      Cov %>% filter(ID %in% Keep) -> Cov ; Phen %>% filter(ID %in% Keep) -> Phen ; Di %>% filter(ID %in% Keep) -> Di
      
      for (Metabolite_n in colnames(Phen) ){
          if (Metabolite_n== "ID"){ next }
          
         Results_meta %>% filter(meta_FDR < 0.05) %>% filter(Metabolite == Metabolite_n ) -> Filter_bugs_n
         Microbes_level %>% select(one_of(Filter_bugs_n$Bug)) -> Sample_Data
         
         Filter_bugs_n %>% filter(Bug %in% colnames(Microbes_level)) -> Interest    
         DIs = c()
         for (x in seq(1, dim(Sample_Data)[1])){
           x = Sample_Data[x,]
           Interest %>% filter(Treatment_effect < 0) -> Negative
           Interest %>% filter(Treatment_effect > 0) -> Positive
           x  %>% select(Negative$Bug) %>% as_vector() %>% as.vector() -> x_p
           Denominator = sqrt(prod(x_p+1))
           x  %>% select(Positive$Bug) %>% as_vector() %>% as.vector() -> x_n
           Numerator = sqrt(prod(x_n+1))
           DI = log10(Numerator/Denominator)
           DIs = c(DIs, DI)
         } 
    
          Phen %>% select(Metabolite_n) -> Phen_2    
          Metabolite_name = colnames(Phen_2)[1]
          Metabolite_measurement = as.vector(as_vector(Phen_2))
          Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement)))
    
          
          Diet_results %>% filter(Metabolite == Metabolite_name) -> food_items
          
          for (Food in unique(food_items$Item)){
            Di %>% select(Food) -> Diet_item
            Diet_name = colnames(Diet_item)[1]
            Diet_ab = as.vector(as_vector(Diet_item))
            ###Fit model
            Cov %>% mutate(Metabolite_measurement = Metabolite_measurement, Bug_abundance = DIs, Food_item=Diet_ab ) -> Model_data
            
            MODEL_FIT = lm(Metabolite_measurement ~ Bug_abundance*Food_item + Age + Gender + BMI, data = Model_data)
            Summary = as.data.frame(summary(MODEL_FIT)$coefficients)
            
            Beta = Summary$Estimate[7]
            pvalue = Summary$`Pr(>|t|)`[7]
            stand = Summary$`Std. Error`[7]
            N_samples = length(MODEL_FIT$residuals)
            
            Results = as_tibble(t(c(Level, Diet_name, Metabolite_name, Beta, stand, pvalue, N_samples)))
            colnames(Results) = c("Level","Food_item", "Metabolite", "Beta", "SE", "Pvalue", "N_samples")
            Total_results_cluster = rbind(Total_results_cluster,Results)    
          }  
        }
      }
    
    
    
    Total_results_cluster %>% mutate(Model=3) ->Total_results_cluster
    
    
    Total_results %>% arrange(Pvalue)  %>% select(-SE,N_samples, Model) %>% gt::gt()
    Total_results_prev %>% arrange(Pvalue)  %>% select(-SE,N_samples, Model) %>% gt::gt()
    Total_results_cluster %>% arrange(Pvalue)  %>% select(-SE,N_samples, Model) %>% gt::gt()



}

#######################################
##Dietary patterns#####################
#######################################

read_tsv("/Users/sergio/Resilio Sync/Transfer/Diet_patterns_LLD.tsv") -> Diet_patterns
Result = tibble()
for (Pheno in colnames(Phenos)){
  if (Pheno == "ID"){ next }
  Phenos %>% dplyr::select(c("ID", Pheno)) -> D
  left_join(left_join(Diet_patterns, D), Covariates) -> Input_model
  Formula = paste(c("`" , Pheno, "`", " ~ RC1 + RC2 + RC4 + RC5 + RC3 + Age+Gender+BMI" ), collapse= "")
  lm(Formula, Input_model) -> Model
  as.data.frame(summary(Model)$coefficients) %>% rownames_to_column("Parameter") %>% mutate(Metabolite=Pheno) -> DF
  Result = rbind(Result, DF)
}
Result %>% as_tibble() %>%  filter(! Parameter %in% c("(Intercept)", "Age", "Gender", "BMI")) -> Result
Result %>% dplyr::select(Estimate, Metabolite, Parameter) %>% spread(Parameter, Estimate) %>% as.data.frame() %>% column_to_rownames("Metabolite") -> wide_results
Result %>% dplyr::select(`Pr(>|t|)`, Metabolite, Parameter) %>% mutate(Significance = ifelse(`Pr(>|t|)`<0.0005, "***", ifelse(`Pr(>|t|)`<0.005, "**", ifelse(`Pr(>|t|)` < 0.05, "*", "")))) %>% dplyr::select(-`Pr(>|t|)`)  %>% spread(Parameter,Significance) %>% as.data.frame() %>% column_to_rownames("Metabolite") -> wide_annotation
paletteLength <- 50
myColor = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(paletteLength)
myBreaks <- c(seq(min(wide_results), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(wide_results)/paletteLength, max(wide_results), length.out=floor(paletteLength/2)))
Plot = pheatmap::pheatmap(wide_results, display_numbers = wide_annotation,fontsize_number = 15, color=myColor, breaks=myBreaks)
#3 and 5 are kind healthy mediterranian habits. 5 is the beverages habit. 3 is the healthy Dutch.
Result %>% mutate(FDR = p.adjust(`Pr(>|t|)`, "fdr")) %>% filter(FDR<0.05)






read_tsv("~/Resilio Sync/Transfer/PhD/TMAO_project/MSS/Data/LLD/Cov_file.tsv") %>% filter(BMI > 30) -> Remove
Data = gather(Data, Metabolite, Concentration, TMAO:`L-Carnitine`, factor_key=TRUE)
Data %>% filter(!ID %in% Remove$ID)  %>% ggplot(aes(x=log10(Concentration), col=Cohort)) + geom_density() + theme_bw() + facet_wrap(~Metabolite, scales = "free") + scale_color_manual(values = c("black","grey"))
for (metabolite in unique(Data$Metabolite) ){
  print(metabolite)
  Data %>% filter(!ID %in% Remove$ID) %>% filter(Metabolite == metabolite) -> Filt
  lm( Concentration ~ Cohort, Filt) %>% summary() %>% print()
}
