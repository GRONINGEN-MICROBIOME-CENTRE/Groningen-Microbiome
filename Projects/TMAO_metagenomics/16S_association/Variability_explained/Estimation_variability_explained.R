###Estimation of variability explained by covariates, microbiome and genetics
set.seed(111)
library(glmnet)
library(caret)
library(tidyverse)

setwd("~/Resilio Sync/Transfer/PhD/TMAO_16S/")
  Count_table =  "LLD_rarefDefault.taxonomyTable.txt.gz"
  Linking_table = "coupling_pheno_16s.txt"
  Covariates = "Cov_file.tsv"
  Phenos = "Pheno_LLD.tsv"
  ##########
  
  Count_table = read_tsv(Count_table) #Rows are individuals, columns are Bacteria
  Linking_table = read_tsv(Linking_table, col_names=F) 
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
Stats_diet %>% filter(X6 %in% c(casefold(Stats_diet_all$...2), "dairy","how_often_yoghurt_milk_based_puddings", "fruits" )) -> Stats_diet2

casefold(unique(as_vector(Stats_diet_all[,2])) )
Stats_diet2

Diet

#Diet %>% dplyr::select(c("ID", filter(Stats_diet, LLD_present=="TRUE")$X6)) -> Diet
Diet %>% dplyr::select(c("ID", casefold(unique(as_vector(Stats_diet[,2])) ))) -> Diet
colnames(Diet)


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
  ID_info = read_tsv("/Users/sergio/Resilio Sync/Transfer/PhD/TMAO_16S/Variance_explained/Data_genetics/ID_genetics.tsv", col_names=F)
  read_tsv(Path_file, col_names=F) -> Matrix_genetics
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
  Formula = as.formula(paste(c(Metabolite_n, " ~ ."), collapse=""))
  model = train( Formula, data = Data , method = "glmnet", trControl = trainControl("cv", number = 10), tuneLength = 10)  
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
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in% colnames(Covariates))  == F]
      }else if (i == "Gene" ){
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Genetics),colnames(Covariates)) ) == F]
      } else if (i == "Mic" ){
        To_0 = colnames(Regressors)[ (colnames(Regressors) %in%  c(colnames(Diet_i) )) == T]
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
Results_all_metabolites = tibble()
for (Met in c("TMAO", "Choline", "Carnitine", "Deoxycarnitine", "Betaine", "TMAO.Choline", "TMAO.Carnitine", "TMAO.Deoxycarnitine", "TMAO.Betaine" ) ){
  print(Met)
  readRDS(file = paste(c("~/Documents/GitHub/Groningen-Microbiome/Projects/TMAO_metagenomics/16S_association/Variability_explained/Models/Model_fitted_", Met,".rds"), collapse= "")) -> Model
  Genetics_table = Get_genetic_table(Met)
  #Check if left_join is doing it by ID
  left_join(left_join(left_join(left_join(Covariates2, Genetics_table), All_taxonomy), Phenos3), Diet) %>% drop_na() -> For_prediction
  
  Predict_by_layers(Model = Model, Covariates = filter(Covariates2, ID %in% For_prediction$ID), Genetics = filter(Genetics_table, ID %in% For_prediction$ID)  , Microbiome = filter(All_taxonomy, ID %in% For_prediction$ID) , Metabolite= as_vector(select( For_prediction, Met)), Diet_i = filter(Diet, ID %in% For_prediction$ID)) -> Results_metabolite
  Results_metabolite %>% mutate(Metabolite = Met) -> Results_metabolite
  rbind(Results_all_metabolites, Results_metabolite) -> Results_all_metabolites
}
Colors = c("light blue", "purple", "salmon", "grey")
Results_all_metabolites %>% mutate(Layer = factor(Layer, levels = c("Diet","Micr", "Gene", "Cov") )) -> Results_all_metabolites
Results_all_metabolites %>% spread(Layer, R2) %>% mutate(Diet = ifelse(Diet - Micr > 0,Diet - Micr, 0) ,Micr = ifelse(Micr - Gene > 0,Micr - Gene, 0) , Gene = ifelse(Gene-Cov > 0, Gene-Cov, 0)  ) %>% gather(Layer, R2, Diet:Cov, factor_key=TRUE) %>%
ggplot(aes(x=R2, y=Metabolite, fill=Layer)) + geom_bar(position = "stack", stat = "identity", col="black") + theme_bw() + scale_fill_manual(values = Colors)






