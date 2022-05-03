###Estimation of variability explained by covariates, microbiome and genetics
set.seed(111)
library(glmnet)
library(caret)
library(tidyverse)

options(collapse_mask = "manip") 
library(collapse)

select = dplyr::select
mutate = dplyr::mutate 
setwd("~/Resilio Sync/Transfer/PhD/TMAO_16S/")
  Count_table =  "LLD_rarefDefault.taxonomyTable.txt.gz"
  Linking_table = "coupling_pheno_16s.txt"
  Covariates = "Cov_file.tsv"
  Phenos = "Pheno_LLD.tsv"
  ##########
  
  Count_table = read_tsv(Count_table) #Rows are individuals, columns are Bacteria
  Linking_table = read_tsv(Linking_table, col_names=FALSE) 
  Covariates = read_tsv(Covariates)
  Phenos = read_tsv(Phenos)
  Phenos %>% mutate(TMAO.Choline = TMAO/Choline , TMAO.Betaine = TMAO/Betaine , 
                    TMAO.Butyrobetaine = TMAO/`y-butyrobetaine`, TMAO.Carnitine = TMAO/`L-Carnitine`,
                    Butyrobetain.Carnitine =  `y-butyrobetaine`/`L-Carnitine` ) -> Phenos
  Diet = read_tsv("20150722_Diet__1135patients.txt")
  colnames(Diet)[1] = "ID"
  
read_tsv("Manuscript/Additional_material/Summary_stats_All.tsv") ->Stats
readxl::read_excel("Stats_diet_all.xlsx",sheet=2) -> Stats_diet_all
read_csv("Diet/LLS_IOP1-FFQ-Codebook_Foodgroups_20210628.csv") -> Stats_diet
Stats_diet %>% 
  filter(X6 %in% c(casefold(Stats_diet_all$...2), "dairy", "fruits" )) -> Stats_diet2

Diet  %>% dplyr::select(c("ID", Stats_diet2$X6)) %>%
  mutate(dairy_products = dairy, fruit = fruits  ) %>% select(-c("dairy",  "fruits")) -> Diet


#Diet %>% dplyr::select(c("ID", casefold(unique(as_vector(Stats_diet[,2])) ))) -> Diet
#colnames(Diet)


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
  Count_table %>% dplyr::select(-Remove_columns) %>% dplyr::select(-rootrank.Root) -> Count_table
  return(Count_table)
}
normalization_inversranknorm = function(Metabolite_measurement){ Metabolite_measurement = qnorm((rank(Metabolite_measurement,na.last="keep")-0.5)/sum(!is.na(Metabolite_measurement))) }



#Format 16S table
Count_table = Filter_unclassified(Count_table)
Linking_table = Choose_replicate(Linking_table)
colnames(Linking_table) = c("ID", "SampleID")
left_join(Count_table, Linking_table) %>% drop_na() %>% select(-SampleID) -> Count_table2
Count_table2 %>% select(ID, Stats$Bug) -> Count_table2 #Keep only bugs in meta-analysis (present in both cohorts)


#Format Genetics table
Get_genetic_table = function(Metabolite){
  Path_file = paste(c("/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_16S/Variance_explained/Data_genetics/", Metabolite, ".tsv"), collapse = "")
  ID_info = read_tsv("/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_16S/Variance_explained/Data_genetics/ID_genetics.tsv", col_names=FALSE)
  read_tsv(Path_file, col_names=FALSE) -> Matrix_genetics
  colnames(Matrix_genetics) = c("ID", ID_info$X1)
  Matrix_genetics[,2:dim(Matrix_genetics)[2]] %>% apply(2,  function(x){ as.numeric(as.factor(x)) })  %>% t() %>% as.data.frame() %>% rownames_to_column("ID") %>% as_tibble()  -> Matrix_genetics2
  colnames(Matrix_genetics2) = c("ID", Matrix_genetics$ID )
  Matrix_genetics2 %>% mutate(ID = as.character(ID)) -> Matrix_genetics2
  return(Matrix_genetics2)
}



#Match tables
Covariates %>% drop_na() -> Covariates2
Phenos %>% drop_na() -> Phenos2
Count_table2 %>% arrange(ID) %>% filter(ID %in% Covariates2$ID) %>% filter(ID %in% Phenos2$ID) -> Count_table2
Covariates2 %>% filter(ID %in% Count_table2$ID) %>% arrange(ID) %>% mutate(ID = as.character(ID)) -> Covariates2
Phenos2 %>% filter(ID %in% Count_table2$ID) %>% select(-Butyrobetain.Carnitine) %>% mutate(TMAO.Deoxycarnitine = TMAO.Butyrobetaine, Deoxycarnitine = `y-butyrobetaine`) %>% 
  select(! c(TMAO.Butyrobetaine, `y-butyrobetaine`) ) %>% arrange(ID) %>% mutate(ID = as.character(ID)) -> Phenos2


#Transform data
#Microbiome transformation
Geom_mean = function(x){
  exp(mean(log(x)))
}
CLR = function(D){
  log(D / Geom_mean(D) )
  
}
Compute_CLR_taxonomy = function(Data, Taxonomy = "genus"){
  paste( c(Taxonomy, "."), collapse="" ) -> Ta
  colnames(select(Data, -c(ID))) -> Taxa
  Taxa[grepl(Ta, Taxa)] -> Taxa
  Data %>% select( c("ID", Taxa) ) -> Data_taxonomy
  
  #computation of a different pseudocoount per taxa.
  CLR( apply( select(Data_taxonomy, -ID), 2, as.numeric)  + 1 ) %>% as_tibble()  %>% mutate(ID = Data_taxonomy$ID, .before=1) -> Data_taxonomy
  return(Data_taxonomy)
}  
All_taxonomy = Compute_CLR_taxonomy(Count_table2, Taxonomy = "genus")
for (i in c("phylum", "class", "order", "family")){
  Compute_CLR_taxonomy(Count_table2, Taxonomy = i) -> Transformed_data
  cbind(All_taxonomy, select(Transformed_data, -ID)) %>% as_tibble() -> All_taxonomy
} 
#Metabolite transformation
apply(select(Phenos2, -ID), 2, normalization_inversranknorm) %>% as_tibble() %>% mutate(ID = Phenos2$ID, .before=1) -> Phenos3
apply(dplyr::select(Diet, -ID), 2, FUN=normalization_inversranknorm) %>% as_tibble() %>% mutate(ID=Diet$ID) -> Diet



Model_metabolite = function( Data,Metabolite_n = "TMAO" ){
  Data %>% select(-ID) -> Data
  Data %>% drop_na() -> Data
  Data  %>% select(-Metabolite_n) %>% apply(2,  scale) %>% as_tibble() %>% mutate(Metabolite = as_vector(select(Data, Metabolite_n)) ) -> Data
  
  Formula = as.formula("Metabolite ~ .")
  #Parallel modelling
  library(doParallel)
  cl <- makePSOCKcluster(5)
  registerDoParallel(cl)
  model = train( Formula, data = Data , method = "glmnet", trControl = trainControl("repeatedcv", number = 10, repeats = 5), tuneLength = 10, allowParallel=TRUE)  
  stopCluster(cl)
  #If warning messages, check: model$resample %>% as_tibble()
  Model_name = paste(c("~/Documents/GitHub/Groningen-Microbiome/Projects/TMAO_metagenomics/16S_association/Variability_explained/Models/Model_fitted_", Metabolite_n, ".rds"), collapse="" )
  saveRDS(model, Model_name)
}

Predict_by_layers = function(Model, Covariates = Covariates2, Genetics, Microbiome = All_taxonomy, Diet_i = Diet, Metabolite = Phenos3$TMAO){
  #All categories should have an "ID" column
  left_join(left_join(left_join(Covariates, Genetics), Microbiome), Diet) -> Regressors
  layers = tibble()
  for (i in  c("Cov", "Gene", "Micr", "Diet")){
    Regressors_0s = Regressors
      if (i == "Cov"){
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in% colnames(Covariates))  == FALSE]
      }else if (i == "Gene" ){
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Genetics),colnames(Covariates)) ) == FALSE]
      } else if (i == "Mic" ){
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Diet_i) )) == TRUE]
      } else { 
        To_0 = c() 
      }
      Regressors_0s[, To_0 ] = 0
      Model %>% predict(select(Regressors_0s, -ID)) -> y_hat
      Variance = R2( y_hat, Metabolite )
      rbind(layers, tibble(Layer = i, R2 = Variance ) ) -> layers
  }
  return(layers)
}
Predict_by_layers2 = function(Model, Covariates = Covariates2, Genetics, Microbiome = All_taxonomy, Diet_i = Diet, Metabolite = Phenos3$TMAO){
  #All categories should have an "ID" column
  left_join(left_join(left_join(Covariates, Genetics), Microbiome), Diet) -> Regressors
  select(Regressors, -ID) %>% apply(2,  scale) %>% as_tibble() %>% mutate(ID = Regressors$ID ) -> Regressors
  layers = tibble()
  for (i in  c("Cov", "Gene", "Micr", "Diet")){
    Regressors_0s = Regressors
    if (i == "Cov"){
      To_0 = colnames(Regressors)[ (colnames(Regressors) %in% colnames(Covariates))  == FALSE]
    }else if (i == "Gene" ){
      To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Genetics),colnames(Covariates)) ) == FALSE]
    } else if (i == "Diet" ){
      To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Microbiome) )) == TRUE]
    } else { 
      To_0 = c() 
    }
    Regressors_0s[, To_0 ] = 0
    Model %>% predict(select(Regressors_0s, -ID)) -> y_hat
    Variance = R2( y_hat, Metabolite )
    rbind(layers, tibble(Layer = i, R2 = Variance ) ) -> layers
  }
  return(layers)
}





All_taxonomy = All_taxonomy[!duplicated(All_taxonomy),]
All_variants = tibble()
for (Met in c("TMAO","Betaine" ,"Choline", "Carnitine", "Deoxycarnitine", "TMAO.Choline","TMAO.Betaine","TMAO.Carnitine", "TMAO.Deoxycarnitine") ){
  print(Met)
  Genetics_table = Get_genetic_table(Met)
  Genetics = tibble(Variants = colnames(select(Genetics_table, -ID) ), Metabolite=Met)
  All_variants = rbind(All_variants, Genetics)
  }  
All_variants %>% dplyr::group_by(Metabolite) %>% dplyr::summarise(N = n())
write_tsv(All_variants, "~/Documents/Variants_per_metabolite.tsv")
################TRAINING##########
###Just run with the training cohorts, or with a train/test split if no training cohort available#########
set.seed(888)
colnames(Phenos3)[5] = "Carnitine"
for (Met in c("TMAO","Betaine" ,"Choline", "Carnitine", "Deoxycarnitine", "TMAO.Choline","TMAO.Betaine","TMAO.Carnitine", "TMAO.Deoxycarnitine") ){
  print(Met)
  Genetics_table = Get_genetic_table(Met)
  #Check if left_join is doing it by ID
  D = left_join(left_join(left_join(left_join(All_taxonomy, Covariates2, by="ID" ), select(Phenos3, c("ID", Met)) , by="ID"), Genetics_table, by="ID"), Diet, by="ID")
  Model_metabolite(Data= D , Metabolite_n = Met  )
}
#Phenos3 %>% mutate(Carnitine = `L-Carnitine`) -> Phenos3
###########R2 estimation in test cohort###################
Test_model = function( Covariates2, All_taxonomy, Phenos3, Diet, Cohort_train = "LLD" ){
  Results_all_metabolites = tibble()
  Coefficients_models = tibble()
  All_coefficients = tibble()
  for (Met in c("TMAO", "Choline", "Carnitine", "Deoxycarnitine", "Betaine", "TMAO.Choline", "TMAO.Carnitine", "TMAO.Deoxycarnitine", "TMAO.Betaine" ) ){
    print(Met)
    if (Cohort_train == "LLD"){
    readRDS(file = paste(c("~/Documents/GitHub/Groningen-Microbiome/Projects/TMAO_metagenomics/16S_association/Variability_explained/Models/Model_fitted_", Met,".rds"), collapse= "")) -> Model
    } else { 
    readRDS(file = paste0("/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_16S/Variance_explained/Models_Rotterdam/Model_fitted_Rotterdam_Study_",Met,".rds" ) ) -> Model
    }
    coef(Model$finalModel, Model$bestTune$lambda) %>% as.matrix() %>% as.data.frame() %>% as.data.frame() %>% rownames_to_column("Feature") %>% as_tibble() -> Model_coeff
    Model_coeff %>% mutate(Metabolite = Met) -> Model_coeff
    rbind(All_coefficients, Model_coeff) -> All_coefficients
    Genetics_table = Get_genetic_table(Met)
    left_join(left_join(left_join(left_join(Covariates2, Genetics_table), All_taxonomy), Phenos3), Diet) %>% drop_na() -> For_prediction
  
    Predict_by_layers2(Model = Model, Covariates = filter(Covariates2, ID %in% For_prediction$ID), Genetics = filter(Genetics_table, ID %in% For_prediction$ID)  , Microbiome = filter(All_taxonomy, ID %in% For_prediction$ID) , Metabolite= as_vector(select( For_prediction, Met)), Diet_i = filter(Diet, ID %in% For_prediction$ID)) -> Results_metabolite
    Results_metabolite %>% mutate(Metabolite = Met) -> Results_metabolite
    rbind(Results_all_metabolites, Results_metabolite) -> Results_all_metabolites
  }
  #Model coeff
  All_coefficients %>% spread(Metabolite, `1`) -> All_coefficients
  write_tsv(All_coefficients, paste0("Variance_explained/Coefficients_",Cohort_train,".tsv"))

  Colors = c("purple", "light blue", "salmon", "grey")
  Results_all_metabolites %>% mutate(R2 = ifelse(is.na(R2), 0, R2 )) %>%
    spread(Layer, R2) %>% mutate(Micr = ifelse(Micr - Diet > 0,Micr - Diet, 0), Diet = ifelse(Diet - Gene > 0,Diet - Gene, 0)  , Gene = ifelse(Gene-Cov > 0, Gene-Cov, 0)  ) %>% gather(Layer, R2, 2:5, factor_key=TRUE) %>%
    mutate(Layer = factor(Layer, levels = c("Micr","Diet", "Gene", "Cov") )) %>%
    ggplot(aes(x=R2, y=Metabolite, fill=Layer)) + geom_bar(position = "stack", stat = "identity") + theme_bw() + scale_fill_manual(values = Colors) + 
    theme(text = element_text(size=20),legend.position="none") + ylab("") -> PLOT
  ggsave(paste0("Variance_explained/PlotR2_",Cohort_train,".pdf" ) ,PLOT, width = 12, height = 15, units = "cm" )
  print(PLOT)
  return(Results_all_metabolites)
}

Test_model(Covariates2, All_taxonomy, Phenos3, Diet, Cohort_train = "LLD") -> LLD_R2
Test_model(Covariates2, All_taxonomy, Phenos3, Diet, Cohort_train = "RS") -> RS_R2


#Check if boostrap changes the findings
Layered_variability_call = function( Covariates2, All_taxonomy, Phenos3, Diet ){
  Results_all_metabolites = tibble()
  for (Met in c("TMAO", "Choline", "Carnitine", "Deoxycarnitine", "Betaine", "TMAO.Choline", "TMAO.Carnitine", "TMAO.Deoxycarnitine", "TMAO.Betaine" ) ){
    readRDS(file = paste(c("~/Documents/GitHub/Groningen-Microbiome/Projects/TMAO_metagenomics/16S_association/Variability_explained/Models/Model_fitted_", Met,".rds"), collapse= "")) -> Model
    Genetics_table = Get_genetic_table(Met)
    #Check if left_join is doing it by ID
    left_join(left_join(left_join(left_join(Covariates2, Genetics_table), All_taxonomy), Phenos3), Diet) %>% drop_na() -> For_prediction
  
    Predict_by_layers2(Model = Model, Covariates = filter(Covariates2, ID %in% For_prediction$ID), Genetics = filter(Genetics_table, ID %in% For_prediction$ID)  , Microbiome = filter(All_taxonomy, ID %in% For_prediction$ID) , Metabolite= as_vector(select( For_prediction, Met)), Diet_i = filter(Diet, ID %in% For_prediction$ID)) -> Results_metabolite
    Results_metabolite %>% mutate(Metabolite = Met) -> Results_metabolite
    rbind(Results_all_metabolites, Results_metabolite) -> Results_all_metabolites
  }
  Results_all_metabolites %>% 
    spread(Layer, R2) %>% mutate(Micr = ifelse(Micr - Diet > 0,Micr - Diet, 0), Diet = ifelse(Diet - Gene > 0,Diet - Gene, 0)  , Gene = ifelse(Gene-Cov > 0, Gene-Cov, 0)  ) %>% gather(Layer, R2, 2:5, factor_key=TRUE) %>%
    mutate(Layer = factor(Layer, levels = c("Micr","Diet", "Gene", "Cov") )) %>%
    ggplot(aes(x=R2, y=Metabolite, fill=Layer)) + geom_bar(position = "stack", stat = "identity") + theme_bw() + scale_fill_manual(values = Colors) -> P
  print(P)
  return(Results_all_metabolites)
}
Boostraps = tibble()
for (i in seq(10)){
  sample(Covariates2$ID, length(Covariates2$ID)*0.8) -> Bootstrap
  Covariates2 %>% filter(ID %in% Bootstrap) -> Covs_boot
  Layered_variability_call(Covs_boot, All_taxonomy, Phenos3, Diet) ->  Results_all_metabolites2
  rbind(Boostraps, Results_all_metabolites2) -> Boostraps 

}
Boostraps %>% dplyr::group_by(Layer, Metabolite) %>% dplyr::summarise( R2 = mean(R2) ) -> Boostraps_all
Boostraps_all %>%
  spread(Layer, R2) %>% mutate(Micr = ifelse(Micr - Diet > 0,Micr - Diet, 0), Diet = ifelse(Diet - Gene > 0,Diet - Gene, 0)  , Gene = ifelse(Gene-Cov > 0, Gene-Cov, 0)  ) %>% gather(Layer, R2, 2:5, factor_key=TRUE) %>%
  mutate(Layer = factor(Layer, levels = c("Micr","Diet", "Gene", "Cov") )) %>%
  ggplot(aes(x=R2, y=Metabolite, fill=Layer)) + geom_bar(position = "stack", stat = "identity") + theme_bw() + scale_fill_manual(values = Colors)



##Comparing regularized coefficients with meta-analysis

read_tsv("Variance_explained/Coefficients_RS.tsv") %>%gather(Metabolite, s1, 2:10) %>% drop_na() -> Coefficients_RS
read_tsv("Variance_explained/Coefficients_LLD.tsv") %>%gather(Metabolite, s1, 2:10) %>% drop_na() -> Coefficients_LLD
read_tsv("Manuscript/Additional_material/Summary_stats_All.tsv") -> Bug_MA
files <- list.files(path="/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_16S/Variance_explained/Signiciant_Genetics/", pattern="*.txt", full.names=TRUE, recursive=FALSE)
Genetic_variants = tibble()
for (file in files){ Name = str_split(file, "/")[[1]] -> y ; N =y[length(y)] ; read.table(file,header = T) %>% as_tibble() %>% mutate(Metabolite = N) -> y2 ; rbind(Genetic_variants, y2)->Genetic_variants }
Genetic_variants$Metabolite =  sapply(Genetic_variants$Metabolite,function(x){ str_replace(str_replace(x,"_LLS_300OB_sugg10e-05_forSergio.txt", ""), "METAANALYSIS_", "") } )

left_join(Coefficients_RS, Coefficients_LLD, by=c("Feature", "Metabolite"), suffix=c("_RS", "_LLD")) -> Combined
Combined %>% mutate(Coefficient_RS.train = s1_RS, Coefficient_LLD.train = s1_LLD ) %>% select(Feature, Metabolite, Coefficient_RS.train, Coefficient_LLD.train ) -> Combined_save
Combined_save %>% mutate(Layer = ifelse( Feature %in% Bug_MA$Bug, "Microbiome", ifelse(grepl(":", Feature), "Genetics", ifelse(Feature %in% c("Age", "Gender", "BMI"), "Anthropometrics", ifelse(Feature == "(Intercept)", "Intercept", "Diet" ) ) )) ) -> Combined_save
write_tsv(Combined_save, "Manuscript/Additional_material/ElasticNet_coefficients.tsv")

Combined_save %>% mutate( Selection = ifelse( Coefficient_RS.train == 0 & Coefficient_LLD.train == 0, "Not_selected", ifelse( Coefficient_RS.train != 0 & Coefficient_LLD.train == 0, "RS", ifelse( Coefficient_RS.train == 0 & Coefficient_LLD.train != 0, "LLD", "Both")))) -> Combined_save
Combined_save %>% mutate(Selection = factor(Selection, c("Not_selected","RS", "LLD", "Both") )) %>% filter(! Layer=="Intercept") %>% ggplot(aes(x=Layer, fill=Selection)) + geom_bar(position = "fill") + theme_bw() + coord_flip() + facet_wrap(~Metabolite, ncol = 3) + scale_fill_manual(values=c("grey", "light blue", "dark blue", "purple" ))
###Microbial Similarity analysis###

Combined %>% filter(Feature %in% Bug_MA$Bug) -> Bugs

Bug_MA %>% select(Bug, Metabolite, MetaP, Treatment_effect, Heterogeneity_Q,Heterogeneity_Pvalue, Beta_LLD, Beta_Rott, Pvalue_LLD, Pvalue_Rott, meta_FDR ) -> Bug_MA
colnames(Bug_MA)[1] = "Feature"

left_join(Bugs, Bug_MA) -> Bacteria_comparison
Bacteria_comparison %>% ggplot(aes(x=s1_RS, y=s1_LLD, col =-log10(MetaP), shape=Heterogeneity_Pvalue<0.05 )) +
geom_point() + theme_bw() + facet_wrap(~Metabolite)
#It is more likely you are picked up if is significantly associated to each cohort, it is heterogenous and there is no effect of the meta-analysis stats
glm( Present ~ log10(Pvalue_Rott) + Heterogeneity_Q + Treatment_effect  ,mutate(Bacteria_comparison, Present= ifelse( s1_RS == 0, 0,1)), family=binomial()) %>% summary()
glm( Present ~ log10(Pvalue_LLD) + Heterogeneity_Q + Treatment_effect  ,mutate(Bacteria_comparison, Present= ifelse( s1_LLD == 0, 0,1)), family=binomial()) %>% summary()
#Check what happened with the significant Meta-analysis findings
Bacteria_comparison %>% dplyr::group_by(s1_LLD ==0,meta_FDR < 0.1 ) %>%
  dplyr::summarise(n())
Bacteria_comparison %>% dplyr::group_by(s1_RS ==0, meta_FDR < 0.1 ) %>%
  dplyr::summarise(n())
#Relation between Coefficients in regularizaiton and not-regularized models
Bacteria_comparison %>% ggplot(aes(x=Beta_LLD, y=Beta_Rott, shape=meta_FDR<0.1, col=(s1_LLD != 0 |  s1_RS != 0)  )) + geom_point() + theme_bw()
Bacteria_comparison %>% ggplot(aes(x=Beta_LLD, y=s1_LLD, shape=meta_FDR<0.1)) + geom_point() + theme_bw()
Bacteria_comparison %>% ggplot(aes(x=Beta_Rott, y=s1_RS, shape=meta_FDR<0.1)) + geom_point() + theme_bw()
#Check selected microbes
Bacteria_comparison %>% dplyr::group_by(Metabolite, s1_LLD != 0, s1_RS != 0 ) %>%
  dplyr::summarise(N = n()) %>% mutate( Match_group = paste0("LLD.",`s1_LLD != 0`,"-","RS.",`s1_RS != 0`) )  %>%
  filter(! Match_group == "LLD.FALSE-RS.FALSE" ) %>%
  ggplot(aes(x=Metabolite, y=N, fill=Match_group)) + geom_bar(stat = "identity") +coord_flip()
#Phylogenetic analysis?
Bacteria_comparison$Feature 


###Genetics ###
Combined$Metabolite = tolower(Combined$Metabolite)
Genetic_variants$Metabolite = tolower(Genetic_variants$Metabolite)
Combined %>% filter(! Feature %in% Bugs$Feature ) %>% filter(grepl(":", Feature)) -> Combined_genetics
Combined_genetics$Feature = sapply(Combined_genetics$Feature, function(x){ str_replace(str_replace_all(str_replace_all(x, "`",""), ":", "_" ), "_", ":") } ) 
left_join()
Genetic_variants %>% filter(! MarkerName %in% Combined_genetics$Feature) -> Change
Change$MarkerName %>% sapply(function(x){ y = str_split(x, "_")[[1]] ; paste0(y[1],"_",y[3],"_",y[2]) }) -> Changed_marker
Genetic_variants$MarkerName[Genetic_variants$MarkerName %in% Change$MarkerName] = Changed_marker
Genetic_variants %>% filter( MarkerName %in% Combined_genetics$Feature) %>% mutate(Feature = MarkerName) %>% select(-MarkerName) -> Variants_matching
Variants_matching$Metabolite %>% sapply(function(x){ str_replace(x, "_", "." )} ) ->Variants_matching$Metabolite

left_join(Combined_genetics, Variants_matching, by=c("Feature", "Metabolite")) -> Genetics_comparison
  
Combined_genetics %>% dplyr::group_by(Metabolite, s1_LLD != 0, s1_RS != 0 ) %>%
  dplyr::summarise(N = n()) %>% mutate( Match_group = paste0("LLD.",`s1_LLD != 0`,"-","RS.",`s1_RS != 0`) )  %>%
  filter(! Match_group == "LLD.FALSE-RS.FALSE" ) %>%
  ggplot(aes(x=Metabolite, y=N, fill=Match_group)) + geom_bar(stat = "identity") +coord_flip()

Genetics_comparison %>% ggplot(aes(x=s1_RS, y=s1_LLD, col =-log10(P.value), shape=HetPVal<0.05 )) +
  geom_point() + theme_bw() + facet_wrap(~Metabolite, scales = "free")
