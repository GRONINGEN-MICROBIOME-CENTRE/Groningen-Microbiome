### Telomere Analysis Script ###
#Author :Sergio Andreu-Sanchez
#This script includes exploration and statistical analysis of telomere lengths of 6 cell types in the Lifelines cohort. 

#Data needed:
#Bioage markers: Include telomere lengts, sj-TRECS and Methylation age
#Lifelines phenotypes: Over 90 phenotypes, including blood biochemical markers, lifestyle questionaires, blood cell counts, cytokine measurements, etc

###WORKING ENV###

#Packages needed for this script
library(tidyverse)

library(RColorBrewer)
library(wesanderson)
library(viridis)

library(pheatmap)
library(glmnet)
library(car)
library("survival")


###Prepare working path
Prefix_path = "" #Prefix of the path where to work
Make_path = function(P){
  paste(c(Prefix_path,P), collapse="/")
}
setwd(Make_path(""))
#########



###########################
#########DATA##############
###########################



######Phenotypes and selectiong of covariates to control for in statistical modelling

#Merged_phenotyopes_LLD is a file containing a selection of phenotypes from Lifelines Deep that were specifically selected for this project. This data should be requested to Lifelines if wants to be used
Phenotypes <- read_delim(Make_path("Merged_Phenotypes_LLD.csv"), delim="|")
#adding a new covariate on whether either father or mather smokes
Phenotypes %>% mutate(Parental_smoking = ifelse(smk13==1 | smk14==1, 1, 0)) -> Phenotypes
#We control for several of those phenotypes in most analysis, including chronological age, sex and cell proportions
#Columns with cell measurements
column_cells = c("Lymph", "Eryth", "Mono", "Neutro", "Thrombo", "Eosino","Blood_baso")
#Columns used for covariates
column_covariates = c("Sex", "AgeIncludingMonth", column_cells)

#It seems taht  a couple of samples had a wrong parental age, we correct it. 
Phenotypes %>% mutate(age_m= ifelse(ID == "LLDeep_0968", 24, age_m)) %>% mutate(age_f = ifelse(ID == "LLDeep_1127", 30, age_f)) -> Phenotypes


######Read biological aging markers
Age_markers <- read_tsv("Combined_data_uncorrected.tsv")
##Selection of telomere markers
Telomere_markers <- colnames(Age_markers)[grepl("MTL",colnames(Age_markers))]
Age_markers %>% dplyr::select(c("Sample",Telomere_markers)) %>% mutate(MTL_gran = as.numeric(MTL_gran), `MTL_CD20+`=as.numeric(`MTL_CD20+`), `MTL_CD45+_CD20-`=as.numeric(`MTL_CD45+_CD20-`),`MTL_CD57+`=as.numeric(`MTL_CD57+`)) %>%
  drop_na() -> Telomere_markers
colnames(Telomere_markers) = c("Sample", "Lymphocytes", "Granulocytes", "Naïve T-cells", "Memory T-cells","B-cells","NK-cells")

### Generate a file with all biological aging markers after removing all linear variability attributed to Age. We do this by retrieving the residuals of a model Bioaging ~ Age , residuals are named as Delta
#Calculate Deltas for each age marker
Get_deltas = function(){
  select(Age_markers, Sample) -> Outcome_delta
  for (i in c( "Methyl.Age.Hannum.", "dCT", "MTL_lymp","MTL_gran","MTL_CD45+_CD20-","MTL_CD45-","MTL_CD20+","MTL_CD57+")){
  #Y --> Bio Age vector
  Age_markers %>% select(i) %>% as_vector() %>% as.vector() -> Marker
  #X1 AgincludingMonth, X2 Sex, X3... Cell lines
  Phenotypes %>% select(c("ID", "AgeIncludingMonth", "Sex", column_cells)) %>% mutate(Sample = ID) %>% select(-ID) -> To_add
  #Add Y and Xs using the common IDs, remove all columns missing data
  left_join(select(Age_markers, Sample), To_add, by="Sample")%>% mutate(Marker = as.numeric(Marker))  %>% drop_na() -> Age_marker
  #Get residuals of the fit
  as.vector(as_vector(lm(Marker ~ ., select(Age_marker, - Sample))$residuals)) -> Delta
  #Add the delta to the Age_marker dataset (after removing NA) where we know the Sample
  Age_marker %>% mutate(Delta = Delta) %>% select(Sample, Delta) -> Output_delta
  #Add column to the output
  colnames(Output_delta) = c("Sample", paste(c("Delta_",i),collapse="")) 
  #Add data to the output
  left_join(Outcome_delta, Output_delta, by="Sample") -> Outcome_delta
}
  write_tsv(Make_path("LLD phenotypes for Sergio/Phenotypes/Delta_BioAge.tsv"),x = Outcome_delta)
}



#################
###FUNCTIONS#####
#################



Linear_regression2 = function(Dependent, Regressor,Covariates, Regressor_name, Dependent_name, Plot = F,cells=F, BMI=F){
  Regressor = as.numeric(Regressor)
  if (BMI == T){  Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Covariates$ID, BMI= Covariates$`Body Mass Index (kg/M^2)`)
  }else{Model_data = tibble(Dependent=Dependent, Regressor=Regressor, Sex= factor(Covariates$Sex), Age= as.numeric(Covariates$Age), ID = Covariates$ID)}
  Model_data %>% drop_na() -> Model_data
  if (BMI == T){
    left_join(Model_data, select(Covariates, c("ID",cell_types)),by="ID") -> Model_data
    Model_data %>% drop_na() -> Model_data
    paste(c("Regressor ~ Sex + Age + BMI",cell_types, "Dependent"), collapse = "+") -> Model_f
  }
  else if (cells==T){
    left_join(Model_data, select(Covariates, c("ID",cell_types)),by="ID") -> Model_data
    Model_data %>% drop_na() -> Model_data
    paste(c("Regressor ~ Sex + Age",cell_types, "Dependent"), collapse = "+") -> Model_f
  }else{
    c("Regressor ~Sex + Age + Dependent") -> Model_f
  }
  Model = lm(Model_f , Model_data )
  Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Dependent",])
  
  Beta = Summary$Estimate
  pvalue = Summary$`Pr(>|t|)`
  stand = Summary$`Std. Error`
  N = length(summary(Model)$residuals)
  Results = as_tibble(t(c(Dependent_name, Regressor_name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
  colnames(Results) = c("Regressor", "Dependent", "Beta", "SE", "Pvalue", "N")
  
  if (Plot == T){
    RES = lm(Regressor ~ Age ,Model_data)$residuals
    
    Model_data %>% mutate(RES = RES) -> Model_data
    if (length(unique(Dependent)) > 8){
      paste(c(Regressor_name, Dependent_name), collapse= " ") -> NAME
      Figure1 = ggplot(data=Model_data,aes(x=RES, y=Dependent)) + geom_point() + theme_bw() + geom_smooth(method ="lm") + ggtitle(NAME)
    } else{ Figure1 = ggplot(data=Model_data,aes(x=factor(Dependent), y=RES)) + geom_boxplot() + theme_bw() + geom_smooth(method ="lm") }
    print(Figure1)
    #try(suppressMessages(ggsave(plot=Figure1,filename=paste(c("C:/Users/Sergio/Documents/PhD/WORK/ImmunoAging/Plots/",Regressor_name,"vs",Dependent_name,".png"),collapse=""))))
  }
  
  
  return(Results)
}



clean_name = function(Name){
  Name= str_trim(Name)
  Name = gsub(" ", "_",Name, fixed= T)
  Name = gsub("(", "",Name, fixed= T)
  Name = gsub(")", "",Name, fixed= T)
  Name = gsub("%", "",Name, fixed= T)
  Name = gsub("/", "",Name, fixed= T)
  #Name = gsub("`", "",Name, fixed= T)
  return(Name)
}

Model_phenotypes = function(Phenotypes, Data, Covariates, cells=NULL, BMI=F, Proteomics = F ){ 
  Overall = tibble()
  for (Col_n in 1:ncol(Phenotypes)){
    Name = colnames(Phenotypes)[Col_n]
    print(Name)
    Name= clean_name(Name)
    Dependent = as.vector(as_vector((Phenotypes[Col_n])))
    #If it is a binary outcome, perform a logistic regression. Otherwise, OLS
    if (length(unique(Dependent[!is.na(Dependent)])) == 2){ 
      Categorical = T
      Dependent = as.factor(Dependent)
    }else if (length(unique(Dependent[!is.na(Dependent)])) >= 2){ 
      Categorical = F
      Dependent = as.numeric(Dependent)
    }else{ next }
    
    for (Col_m in 1:ncol(Data)){
      Explanatory  = as.vector(as_vector(Data[,Col_m]))
      Name_2 = colnames(Data)[Col_m]
      Name_2= clean_name(Name_2)
      if (Name %in% cell_types){ Dependent = scale(as.numeric(Dependent))  ;  Explanatory = as.numeric(Explanatory)}
      To_plot = c()#c("Thrombo", "Neutro")#c("age_f", "age_m", "BMI", "Parental_smoking","smk13","Waist_circum" )
      if (Name %in% To_plot){ PL = T
      } else{ PL = F}
      if (Categorical == T){
        if (is.null(cells)){
          Results = Linear_regression2(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL)
        } else{ Results = Linear_regression2(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL, cells =T) }
      }else{
        if (Proteomics == T){ 
          Results = lm_proteomics(Dependent,Explanatory, Covariates, Name_2, Name)
        }else if (is.null(cells)){
          Results = Linear_regression2(Dependent, Explanatory, Covariates,Name_2, Name, Plot = PL)
        }else{ Results = Linear_regression2(Dependent, Explanatory, Covariates, Name_2, Name, Plot= PL, cells =T, BMI=BMI) }  
      }
      
      #Count how many entries per level
      if (length(unique(Dependent[!is.na(Dependent)])) < 10){ 
        LEVEL =  paste(unique(Dependent[!is.na(Dependent)]),collapse=",")
        N_LEVEL = c()
        for (L in str_split(LEVEL, ",")[[1]]){
          N_LEVEL = c(N_LEVEL,length(which(as.character(Dependent) == L)))
        }
        N_LEVEL = paste(N_LEVEL, collapse=",")
      }else{ 
        LEVEL =  "Numeric" 
        N_LEVEL = length(Dependent[!is.na(Dependent)])
      }
      
      Results %>% mutate(Levels = LEVEL, Levels_n = N_LEVEL) -> Results
      Overall = rbind(Overall, Results)
    }
  }
  return(Overall)     
}

Compare_outlaiers_bioage2 = function(DF = Age_markers, Phenos = Phenotypes, Covariates=Cov, cells=cell_types ){
  DF_all_shared = tibble()
  left_join(Phenos , select(Covariates,c("ID","Age","Sex"))) -> Phenos
  apply(select(DF, -Sample),2,FUN=as.numeric) %>% as_tibble() %>% mutate(Sample = DF$Sample) -> DF
  for (BioAge in colnames(DF)[grepl("MTL",colnames(DF))]){
    if (BioAge %in% c("ID","Age")){ next }
    DF %>% select(c(Sample,BioAge)) %>% drop_na() -> temp_df
    Phenos %>% filter(ID %in% temp_df$Sample) %>% arrange(ID) %>% drop_na(Age) -> AGE  
    temp_df %>% filter(Sample %in% AGE$ID) %>% arrange(Sample) %>% select(BioAge) %>% as_vector() %>% as.vector() -> bioage
    lm(bioage ~ AGE$Age) -> Model
    
    residual_sd = summary(Model)$sigma
    predict(Model) -> Mean_value
    Upper_threshold = Mean_value + 1*residual_sd #1.2sd
    Lower_threshold = Mean_value - 1*residual_sd
    
    AGE %>% mutate(bioage) %>% filter(bioage > Upper_threshold) %>% select(ID) %>% as_vector() %>% as.vector() ->outliers_upper
    AGE %>% mutate(bioage) %>% filter(bioage < Lower_threshold) %>% select(ID) %>% as_vector() %>% as.vector() ->outliers_lower
    
    DF %>% select(Sample) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Up_outlier = ifelse(ID %in% outliers_upper, T, F), Down_outlier = ifelse(ID %in% outliers_lower, T, F)) %>%
      mutate(Marker = BioAge) -> DF_all
    rbind(DF_all_shared, DF_all) -> DF_all_shared
  }
  DF_all_shared %>% filter(Up_outlier == T | Down_outlier == T) %>%  group_by(Marker) %>% summarise(n())
  DF_all_shared %>% mutate(Status = Up_outlier+Down_outlier) %>% group_by(ID) %>% summarise(N=sum(Status)) %>% arrange(desc(N)) %>% filter(N>=2) ->
    Outliers
  DF_all_shared  %>% filter(ID %in% Outliers$ID) %>% mutate(Up_outlier = ifelse(Up_outlier == TRUE, 1, 0), Down_outlier = ifelse(Down_outlier == TRUE, -1, 0)) %>% filter(!Down_outlier+Up_outlier == 0) %>% 
    mutate(Outlier = Up_outlier+Down_outlier) %>% group_by(ID) %>% summarise(N= sum(Outlier), C=n() ) %>% filter(! abs(N) != C ) ->
    Groups 
  
  Down_group = filter(Groups, N<0 )$ID
  Up_group =  filter(Groups, N>0 )$ID
  
  select(DF, Sample, Methyl.Age.Hannum.) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Group = ifelse(ID %in% Up_group, "1-Long",ifelse(ID %in% Down_group, "0-Short",NA))) %>% drop_na() -> Methylation
  left_join(Methylation,Phenos) -> Methylation
  summary(glm(as.factor(Group) ~ Methyl.Age.Hannum. + Age + Sex, Methylation, family = binomial(link="logit")))
  Methylation %>% select(Group,  Methyl.Age.Hannum.,Age,Sex, cells ) -> Methylation2
  summary(glm(as.factor(Group) ~ . , Methylation2, family = binomial(link="logit")))
  
  select(DF, Sample, dCT) %>% mutate(ID=Sample) %>% select(-Sample) %>% mutate(Group = ifelse(ID %in% Up_group, "1-Long",ifelse(ID %in% Down_group, "0-Short",NA))) %>% drop_na() -> TRECS
  left_join(TRECS,Phenos) -> TRECS
  summary(glm(as.factor(Group) ~ dCT + Age + Sex, TRECS, family = binomial(link="logit")))
  TRECS %>% select(Group,  dCT,Age,Sex, cells ) -> TRECS2
  summary(glm(as.factor(Group) ~ . , TRECS2, family = binomial(link="logit"))) 
}



SAVE = function(NAME, PLOT, Width = 6.3, Height = 6.3, Dpi = 300, Scale = 1){
  ggsave(filename = NAME, PLOT,
         width = Width, height = Height, dpi = Dpi, units = "in", scale = Scale)
}


#################################
##Prepare input for Analysis#####
#################################

###Phenotypes
Phenotypes %>% filter(ID %in% Telomere_markers$Sample) %>% arrange(ID) %>%
  filter(! is.na(AgeIncludingMonth)) %>% filter(! is.na(Sex)) -> Phenotypes
Phenotypes %>% mutate(smk_now = ifelse(smk_now == "2", "0", smk_now) ) %>% select(-one_of("HR")) -> Phenotypes
Phenotypes %>% mutate(Duration_breastfeeding = ifelse(Duration_breastfeeding == 7, NA, Duration_breastfeeding))-> Phenotypes
Phenotypes %>% mutate(Type_delivery = ifelse(Type_delivery == 3, 2, Type_delivery))-> Phenotypes

###Telomeres
Telomere_markers %>% filter(Sample %in% Phenotypes$ID) %>% arrange(Sample) -> Telomere_markers
Telomere_markers %>% filter(Sample %in% Proteomics$ID) %>% arrange(Sample) -> Telomere_markers_protein
###Prepare covaraites for Phenotypes
Phenotypes %>% dplyr::select(c("ID",column_covariates)) %>% mutate(Age = AgeIncludingMonth) %>% dplyr::select(-AgeIncludingMonth) -> Cov
Phenotypes %>% dplyr::select(- any_of(c("AgeIncludingMonth", "AgeInDays", "Sex"))) -> Phenotypes

cell_types = column_cells

cell_types_PC %>% arrange(`-`) %>% filter(`-` %in% Cov$ID) -> cell_types_PC
all_cell_types %>% arrange(ID) %>% filter(ID %in% Cov$ID) -> all_cell_types

###########################
##Exploratory analysis#####
###########################

#Create dataframe with telomeres and covariates
left_join(mutate(Telomere_markers, ID = Sample), select(Cov,c("ID", "Age" , "Sex")), by="ID") %>% gather(Cell_line, Telomere_length, 2:7) %>% select(-Sample) %>% mutate(Telomere_length = as.numeric(Telomere_length)) -> Exploration_telomeres

#distribution of telomere length by gender
Exploration_telomeres %>% ggplot(aes(x=Telomere_length, col=as.factor(Sex) )) + geom_density() + theme_bw() + facet_wrap(~ Cell_line) +
  scale_color_manual(values = wesanderson::wes_palette("Royal1"),name="Sex",labels = c("Male", "Female")) -> Distribution_plot
SAVE(NAME =  Make_path("Figures/Distribution_Telomeres.pdf"),PLOT =Distribution_plot)
#Distribution telomere length
Exploration_telomeres %>% mutate(Cell_line=factor(Cell_line, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>% #c("Naïve T-cells","B-cells","Granulocytes", "Lymphocytes","Memory T-cells", "NK-cells"))) %>%
  ggplot(aes(x=Age, y=Telomere_length, col=Cell_line ))  + theme_bw(base_size = 12)  + geom_smooth(method="lm", formula=y~x, size=1.5) + 
  scale_color_manual(values = viridis::viridis(6)) + ylab("Telomere length (Kb)") -> telomere_trends
SAVE(NAME = Make_path("Figures/Telomeres_comparison_trend.pdf"),PLOT = telomere_trends)

Exploration_telomeres %>% mutate(Cell_line=factor(Cell_line, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>% #c("Naïve T-cells","B-cells","Granulocytes", "Lymphocytes","Memory T-cells", "NK-cells"))) %>% 
ggplot(aes(y=Telomere_length, x=Cell_line, fill=Cell_line )) +  geom_violin() + geom_boxplot(width=0.1, outlier.alpha = 0) + scale_fill_manual(values=viridis::viridis(6)) + 
 theme_bw(base_size = 12) + theme(axis.title.x=element_blank(),axis.text.x=element_blank() ,axis.ticks.x=element_blank()) +ylab("Telomere length (Kb)") -> distribution_length_by_Cell
                                                                                                                                                                             
SAVE(NAME = Make_path("Figures/Telomeres_Distribution_lengths.pdf"),PLOT = distribution_length_by_Cell)

pairwise.wilcox.test(Exploration_telomeres$Telomere_length,Exploration_telomeres$Cell_line,paired = T)

#Number of Males/Females
Exploration_telomeres %>% group_by(Sex,Cell_line) %>% summarise(n())
#correlation between telomere lengths
Exploration_telomeres %>% spread(Cell_line, Telomere_length) %>% drop_na() -> Exploration_telomeres_wide
cor(Exploration_telomeres_wide[4:9],method = "pearson") -> Correlation_matrix
pheatmap(Correlation_matrix, color=RColorBrewer::brewer.pal(name ="YlOrRd", n=9)) -> Pheatmap_correlation_Telomeres
SAVE(NAME = Make_path("Figures/Heatmap_Correlation_Telomeres.pdf"),PLOT = Pheatmap_correlation_Telomeres)
#In absolute values
pheatmap(abs(Correlation_matrix), color=viridis::viridis(10)[10:20]) -> Pheatmap_correlation_Telomeres_abs
SAVE(NAME = Make_path("Figures/Heatmap_Correlation_Telomeres_abs.pdf"),PLOT = Pheatmap_correlation_Telomeres_abs)



#Correlation with Age and with Gender
Exploration_telomeres %>% mutate(Cell_line=factor(Cell_line, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>% 
  ggplot(aes(x=Age, y=Telomere_length, col=as.factor(Sex) )) + geom_point(size=1, alpha=0.3) + theme_bw(base_size = 12) + facet_wrap(~ Cell_line, ncol = 2) +
  scale_color_manual(values = viridis::viridis(2),name="Sex",labels = c("Male", "Female")) + geom_smooth(method="lm", formula=y~x) + scale_x_continuous(name="Age") + ylab("Telomere length (Kb)") -> Age_vs_Telomere
SAVE(NAME = Make_path("Figures/Age_vs_Telomere.pdf"),PLOT = Age_vs_Telomere, Width=12.6, Height = 12.6)
#Just lymphocytes
Exploration_telomeres %>% filter(Cell_line == "Lymphocytes") %>%  ggplot(aes(x=Age, y=Telomere_length, col=as.factor(Sex) )) + geom_point(size=1, alpha=0.3) + theme_bw(base_size = 12) + 
  scale_color_manual(values = viridis::viridis(2),name="Sex",labels = c("Male", "Female")) + geom_smooth(method="lm", formula=y~x) + scale_x_continuous(name="Age") +ylab("Telomere length (Kb)") -> Age_vs_Telomere
SAVE(NAME = Make_path("Figures/Age_vs_Telomere_lymph.pdf"),PLOT = Age_vs_Telomere)



#Relation of Telomeres with other Biological aging values - uncorrected for chronological Age
Age_markers %>% select(c(Sample, dCT, Methyl.Age.Hannum.)) %>% mutate(sjTREC = as.numeric(dCT), ID=Sample) %>% select(-c(Sample,dCT)) -> Other_markers
left_join(Exploration_telomeres_wide, Other_markers, by="ID") %>% drop_na() %>% select(-Sex) -> Subset_all
cor(select(Subset_all, -ID),method = "pearson") -> Correlation_matrix
pheatmap(Correlation_matrix, color=rev(RColorBrewer::brewer.pal(name ="RdYlBu", n=9))) -> Corre_bioages
SAVE(NAME = Make_path("Figures/Correlation_BioAges.pdf"),PLOT = Corre_bioages)
pheatmap(abs(Correlation_matrix), color=viridis::viridis(10)[2:10]) -> Pheatmap_correlation_bioages_abs
SAVE(NAME = Make_path("Figures/Heatmap_Correlation_Telomeres_abs.pdf"),PLOT = Pheatmap_correlation_bioages_abs)
#Relation of Telomeres with other Biological aging values - corrected for chronological Age
Corrected_all = tibble(Age = Subset_all$Age)
for (Age_1 in seq(dim(Subset_all)[2]) ){
  n_Age_1 = colnames(Subset_all)[Age_1]
  if (n_Age_1 %in% c("ID", "Age") ){ next }
  New = tibble(Age_m1 = as_vector(Subset_all[n_Age_1]), Age = Subset_all$Age )
  R1 = lm(scale(Age_m1) ~ scale(Age), New)$residuals
  NT = as_tibble(R1) ; colnames(NT) =  n_Age_1
  Corrected_all = cbind(Corrected_all, NT)
}
cor(Corrected_all,method = "pearson") -> Correlation_matrix_corrected
pheatmap(abs(Correlation_matrix_corrected), color=viridis::viridis(10)[2:10] ) -> Corre_bioages_corrected
SAVE(NAME = Make_path("Figures/Heatmap_Correlation_Telomeres_abs_corrected.pdf"),PLOT = Corre_bioages_corrected)


#Are samples aging too-fast or too-slow the same in the different bio-aging measurements?
Compare_outlaiers_bioage2(Phenos = left_join(Phenotypes, mutate(Cov, AgeIncludingMonth=Age)))

#Exploration of phenotypes availables, basic metrics
Summary_Stats = tibble()
for (Phenotype in colnames(Phenotypes)){ 
  if (Phenotype == "ID"){ next }
  Phenotypes %>% select(Phenotype) %>% as_vector() %>% as.vector() -> Info_pheno
  Number_NA = sum(is.na(Info_pheno))
  Info_pheno = Info_pheno[!is.na(Info_pheno)]
  Not_NA = length(Info_pheno)
  if (length(unique(Info_pheno)) < 7){
    All_info = ""
    Info_pheno = as.factor(Info_pheno)
    for (i in levels(Info_pheno)){
      Number_level = length(Info_pheno[Info_pheno == i])
      paste(c("Level:",as.character(i),",N:",as.character(Number_level)), collapse="") -> Level_info
      All_info = paste(All_info, Level_info, collapse=" ")
    }
  } else{ 
    Info_pheno = as.numeric(Info_pheno)
    All_info = paste(c("Mean:",as.character(round(mean(Info_pheno),2))," Sd:", as.character(round(sd(Info_pheno),2))," Max:", as.character(max(Info_pheno)), " Min:",as.character(min(Info_pheno))),collapse="") 
  }
  Pheno_stats= tibble(Phenotype = Phenotype, NAs= Number_NA, Not_NA= Not_NA, Info=All_info)
  Summary_Stats = rbind(Summary_Stats, Pheno_stats)
}
write_tsv(Summary_Stats, Make_path("Results/Summary_phenotypes.tsv"))


#############################
#####Phenotype association###
#############################


##Association without controlling for cell number

Model_phenotypes(select(Phenotypes,-"ID"), select(Telomere_markers, -"Sample"), Cov) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results
Model_results %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results1
write_tsv(Model_results1, Make_path("Results/Phenotype_Telomere_assoc_Nocell.tsv"))

##Association while controlling for cell number

Model_phenotypes(select(Phenotypes,-"ID"), select(Telomere_markers, -"Sample"), Cov, cells=T) %>% mutate(Pvalue=as.numeric(Pvalue)) -> Model_results
Model_results %>% mutate(FDR = p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue) -> Model_results

#Add associations in a table
left_join(select(Model_results, c(Regressor, Dependent, ,Beta, Pvalue, FDR, Levels_n)),select(Model_results1, c(Regressor, Dependent,Beta, Pvalue, FDR, Levels_n)), by=c("Regressor", "Dependent", "Levels_n"),suffix=c(".cellcorrec",".notcorrec") ) -> Table_results
Table_results %>% filter(FDR.cellcorrec < 0.05 | FDR.notcorrec<0.05 ) %>% arrange(Pvalue.cellcorrec) %>% mutate(Beta.cellcorrec=as.numeric(Beta.cellcorrec), Beta.notcorrec=as.numeric(Beta.notcorrec)) -> Table_results_sig
write_tsv(Table_results_sig,path = "Results/Phenos_vs_Telomere_allSig.tsv")
write_tsv(Model_results,path = "Results/Phenotype_Telomere_assoc.tsv")

Model_results %>% filter(FDR<0.05) %>% group_by(Regressor) %>% summarise(n())
Model_results %>% filter(Regressor %in% filter(Model_results,FDR<0.05)$Regressor ) %>% print()

Model_results %>% filter(Regressor %in% filter(Model_results, FDR < 0.05)$Regressor) %>% select(c(Regressor,Dependent,Beta,SE,Pvalue)) -> SR
Final = NULL
for (TL in unique(Model_results$Dependent)){
  paste(c("Estimate",TL), collapse="\n" ) -> N_1 ; paste(c("SE",TL), collapse="\n" ) -> N_2 ; paste(c("Pvalue",TL), collapse="\n" ) -> N_3
  SR %>% filter(Dependent == TL) %>% select(-Dependent) -> SR2
  colnames(SR2) = c("Regressor", N_1, N_2, N_3)
  if (is.null(Final)){ Final = SR2
  }else{ Final = left_join(Final, SR2, by="Regressor" ) }
}

  
#################Check the signiciant hits.
######BMI
# Are associations to BMI independent to parental smoking and parental age?

Phenotypes %>% select(Waist_circum,BMI) -> Phenotype_BMI
left_join(Cov, select(Phenotypes, c(ID,age_f, age_m, Parental_smoking))) -> Cov_BMI
select(Telomere_markers, -"Sample") -> Data_BMI

BMI_associations = tibble()

for (Col_n in 1:ncol(Phenotype_BMI)){
  Dependent = as.vector(as_vector((Phenotype_BMI[,Col_n])))
  Name = colnames(Phenotype_BMI)[Col_n]
  Name= clean_name(Name)
  for (Col_m in 1:ncol(Data_BMI)){
    Explanatory  = as.vector(as_vector(Data_BMI[,Col_m]))
    Name_2 = colnames(Data_BMI)[Col_m]
    Name_2= clean_name(Name_2)
    Cov_BMI %>% select(-ID) %>% mutate(Dependent = Dependent, Explanatory = as.numeric(Explanatory)) -> Model_input
    #lm(Dependent ~ Explanatory + ., Model_input) -> Model
    lm(Explanatory ~ Dependent + ., Model_input) -> Model
    
    #Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Explanatory",])
    Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Dependent",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|t|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, Name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    colnames(Results) = c("Dependent", "Regressor", "Beta", "SE", "Pvalue", "N")
    BMI_associations = rbind(BMI_associations, Results)
  }
}
library(patchwork)
#Plotting to compare Beta before and after corrections and Pvalues before and after corrections
inner_join(Model_results, BMI_associations, c("Regressor", "Dependent"), suffix=c("Regular", "Parentsmk_corrected")) -> Check_BMI
Check_BMI %>%
  ggplot(aes(x=as.numeric(BetaRegular), y=as.numeric(BetaParentsmk_corrected), col=Dependent)) + geom_point() + theme_bw() -> BMI_1

Check_BMI %>%
  ggplot(aes(x=-log10(as.numeric(PvalueRegular)), y=-log10(as.numeric(PvalueParentsmk_corrected)),col=Dependent)) + geom_point() + theme_bw() + geom_hline(yintercept = -log10(0.05)) -> BMI_2
BMI_1 + BMI_2 + plot_layout(guides = "collect")


####Smoking
#Effect of different smoking phenotypes
Model_results %>% filter(grepl("smk",Regressor) ) %>% filter(!Regressor %in% c("smk2","Age_smk_stop")) %>% mutate(Beta = as.numeric(Beta), Phenotype = Regressor) %>%
  ggplot(aes(x= Phenotype, y=Beta, fill=-log10(Pvalue), col=Pvalue<0.05 )) + geom_bar(size=1, stat = "identity") +
  theme_bw(base_size = 12) + facet_wrap(~Dependent) + coord_flip() + scale_color_manual(values = c("white", "black")) + scale_fill_viridis_c()+ #scale_fill_distiller(palette="YlOrRd" ,direction=1) +
  scale_x_discrete(labels=c("smk1" = "Have_you_smoked_a_year?", "smk2" ="How_old_when_started_smoking?","smk13"= "father_smk", "smk14" = "mother_smk", "smk15"="mother_smk_pregnancy","smk10"="Number_people_smoke_at_home","smk11"= "Do_people_smoke_work","smk_now"="Current_smoker")) +
  geom_errorbar(aes(ymin=Beta-as.numeric(SE) , ymax=Beta + as.numeric(SE)),width=.2, col="black")-> Smoking_Associations
SAVE(NAME = "Figures/Smoking associations.pdf" ,PLOT = Smoking_Associations)
#compare effect size of Father and  Mother smoking
Model_results %>% filter(grepl("smk",Regressor) ) %>% filter(Regressor %in% c("smk13", "smk14")) -> Man_vs_woman
Man_vs_woman %>% group_by(Regressor) %>% summarise(mean(Beta))
summary(lm(Beta ~ Regressor, Man_vs_woman)) #Significantly higher Beta in men
left_join(filter(Man_vs_woman,Regressor=="smk13"),filter(Man_vs_woman,Regressor=="smk14"), by="Dependent") %>% ggplot(aes(x=as.numeric(Beta.x), y=as.numeric(Beta.y), col=Dependent)) + geom_point()+theme_bw() #+geom_abline(intercept= -0.28317, slope=0.08748)

#Age association - combine mother-parent age in model an calculate inflation.
#Compute age gap and calculate association
Phenotypes %>% select(c(ID, age_f, age_m)) %>% mutate(gap = age_f - age_m) -> Phenotype_age
Age_associations = tibble()

for (Col_m in 2:ncol(Telomere_markers)){
    Dependent  = as.vector(as_vector(Telomere_markers[,Col_m]))
    Name_2 = colnames(Telomere_markers)[Col_m]
    Name_2= clean_name(Name_2)
    Cov %>% select(-ID) %>% mutate(Dependent = Dependent) -> Model_input
    cbind(Model_input,Phenotype_age) %>% as_tibble() -> Model_input
    Model_input %>% drop_na() -> Model_input
    lm(Dependent ~ age_m + age_f + Sex+ Lymph+ Eryth + Mono+ Neutro+ Thrombo +Eosino + Blood_baso+   Age  , Model_input) -> Model
    car::vif(Model) -> VIF
    lm(Dependent ~ gap + Sex+ Lymph+ Eryth + Mono+ Neutro+ Thrombo +Eosino + Blood_baso+   Age  , Model_input) -> Model2
    
    Summary = as_tibble(as.data.frame(summary(Model2)$coefficients)["gap",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|t|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, "Gap_phenotype", as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    Results %>% mutate(VIF_mother = VIF["age_m"], VIF_father = VIF["age_f"]) -> Results
    colnames(Results) = c("Dependent", "Regressor", "Beta", "SE", "Pvalue", "N", "VIF_mother", "VIF_father")
    Age_associations = rbind(Age_associations, Results)
  }

Age_associations %>% mutate(FDR=p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue)



##CVD-related phenotypes with BMI correction
Phenotypes %>% select(c(ID, Trigly, Heart_attack, HDL_cholesterol, LDL_cholesterol, Pulse_rate)) -> Phenotype_CVD

left_join(Cov, select(Phenotypes, c(ID, BMI))) -> Cov_CVD
CVD_associations = tibble()
select(Telomere_markers, -"Sample") -> Data_CVD

for (Col_n in 2:ncol(Phenotype_CVD)){
  Dependent = as.vector(as_vector((Phenotype_CVD[,Col_n])))
  Name = colnames(Phenotype_CVD)[Col_n]
  Name= clean_name(Name)
  for (Col_m in 1:ncol(Data_CVD)){
    Explanatory  = as.vector(as_vector(Data_CVD[,Col_m]))
    Name_2 = colnames(Data_CVD)[Col_m]
    Name_2= clean_name(Name_2)
    Cov_CVD %>% select(-ID) %>% mutate(Dependent = Dependent, Explanatory = as.numeric(Explanatory)) -> Model_input
    Model_input %>% drop_na() -> Model_input
    lm(Explanatory ~ Dependent + ., Model_input) -> Model
    
    Summary = as_tibble(as.data.frame(summary(Model)$coefficients)["Dependent",])
    Beta = Summary$Estimate
    pvalue = Summary$`Pr(>|t|)`
    stand = Summary$`Std. Error`
    N = length(summary(Model)$residuals)
    Results = as_tibble(t(c(Name_2, Name, as.numeric(Beta), as.numeric(stand), as.numeric(pvalue), N)))
    colnames(Results) = c("Dependent", "Regressor", "Beta", "SE", "Pvalue", "N")
    CVD_associations = rbind(CVD_associations, Results)
  }
}
CVD_associations %>% mutate(FDR=p.adjust(Pvalue, "fdr")) %>% arrange(Pvalue)




#####################
######GENETICS#######
#####################

#This part of the script uses Polygenic Risk Scores

#We compare the association founds with PRS computeted from to different studies, using weighted and unweighted betas

###PRS association
##Load PRS##
Codd_beta = as_tibble(read.table("PRS/codd_prs.profile",header = T))
Codd_nobeta = as_tibble(read.table("PRS/codd_prs_nobeta.profile", header=T))
Li_beta = as_tibble(read.table("PRS/li_prs.profile",header=T))
Li_nobeta = as_tibble(read.table("PRS/li_prs_nobeta.profile",header=T))
Results_PRS = tibble()
for (N_PRS in seq(1,4)){
  PRS = list(Codd_beta,Codd_nobeta,Li_beta,Li_nobeta)[[N_PRS]]
  PRS_name = c("Codd_beta","Codd_nobeta","Li_beta","Li_nobeta")[N_PRS]
  PRS %>% select(IID,SCORESUM ) -> PRS
  lapply(PRS$IID, FUN= function(x){ str_split(x,"1_")[[1]][2] }) %>% unlist() -> IDs
  PRS %>% mutate(ID = IDs, PRS_value = SCORESUM) -> PRS
  
  for (Cell_line in colnames(Telomere_markers)){
    if (Cell_line == "Sample"){ next }
    Telomere_markers %>% filter(Sample %in% PRS$ID) %>% arrange(Sample) %>% select(Cell_line) %>% as_vector() %>% as.vector() %>% as.numeric() -> Dependent
    PRS %>% arrange(ID) %>% filter(ID %in% Telomere_markers$Sample ) %>% mutate(Dependent =Dependent) -> Entry_model
    Cov -> Cov_to_add
    left_join(Entry_model, Cov_to_add, by="ID") -> Entry_model
    Entry_model %>% drop_na() -> Entry_model
    MODEL = lm(Dependent ~ PRS_value + Age, Entry_model)
    as.data.frame(confint(MODEL))["PRS_value",] -> Confidence
    as.data.frame(summary(MODEL)$coefficients)["PRS_value",][,c(1,4)] %>% as_tibble() %>% mutate(Line = Cell_line, Input = PRS_name) -> RR
    RR %>% mutate(Low = Confidence[,1], Up=Confidence[,2] ) -> RR
    Results_PRS = rbind(Results_PRS, RR)
  }
  
}
Results_PRS %>% filter(Input == "Li_beta") %>% mutate(Pvalue = `Pr(>|t|)`, Telomere=Line) %>% mutate(Telomere=factor(Telomere, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>%
  ggplot(aes(y=Estimate,x= Telomere)) + geom_point(aes(col = -log10(Pvalue)), size=4, shape=15) +
  coord_flip() + theme_bw(base_size = 12) + geom_hline(yintercept=0, linetype="dotted", size=1.5) +
  geom_errorbar(aes(ymin=as.numeric(Low), ymax=as.numeric(Up)),width=.2) + scale_colour_viridis_c() -> PRS_effect
  

SAVE(PLOT = PRS_effect,NAME ="~/Figures/PRS_associations.pdf", Height = 2)

PRS = Li_beta 


#####Plot significant results

clean_name(colnames(Phenotypes)) -> colnames(Phenotypes)
Phenotypes %>% select(c(unique(filter(Model_results,FDR<0.05)$Regressor),ID)) -> Phenotypes_significant


##PLOT1
Model_results %>% filter(Regressor %in% colnames(Phenotypes_significant)) -> Model_results2
Model_results2$Regressor[Model_results2$Regressor == "age_f"] = "Age Father(years)" ;Model_results2$Regressor[Model_results2$Regressor == "age_m"] = "Age Mother(years)" ;Model_results2$Regressor[Model_results2$Regressor == "smk15"] = "Mother smoke while pregnant(categories)" ; Model_results2$Regressor[Model_results2$Regressor == "smk13"] = "Father smoke(yes/no)" ;Model_results2$Regressor[Model_results2$Regressor == "Parental_smoking"] = "Any parent smoke(yes/no)" ; Model_results2$Regressor[Model_results2$Regressor == "smk14"] = "Mother smoke(yes/no)"
Model_results2$Regressor[Model_results2$Regressor == "Weight_kg"] = "Weight(Kg)" ; Model_results2$Regressor[Model_results2$Regressor == "Waist_circum"] = "Waist circumference(cm)" ; Model_results2$Regressor[Model_results2$Regressor == "Pulse_rate"] = "Pulse rate(per minute)" ; Model_results2$Regressor[Model_results2$Regressor == "Poorly_healing_wounds_feet"] = "Poorly healing wounds in feet(yes/no)"


Model_results2  %>% select(Regressor, Dependent, Beta) %>% mutate(Direction = ifelse(Beta>0, 1, -1)) -> Plot_direction
Plot_direction %>% filter(! Regressor == "Smk_mother_childhood") -> Plot_direction

Plot_direction$Regressor[Plot_direction$Regressor == "age_f"] = "Age Father" ;Plot_direction$Regressor[Plot_direction$Regressor == "age_m"] = "Age Mother" ;Plot_direction$Regressor[Plot_direction$Regressor == "smk15"] = "Mother smoke while pregnant" ; Plot_direction$Regressor[Plot_direction$Regressor == "smk13"] = "Father smoke" ;Plot_direction$Regressor[Plot_direction$Regressor == "Parental_smoking"] = "Any parent smoke"
Plot_direction %>% mutate(Direction = as.numeric(Beta)) %>% select(Dependent, Regressor, Direction) %>% spread(Dependent, Direction) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> wide_results2


Model_results2 %>%filter(! Regressor == "Smk_mother_childhood") %>% select(Regressor, Dependent, Pvalue) %>% mutate(Significance = ifelse(Pvalue<0.0005, "***", ifelse(Pvalue<0.005, "**", ifelse(Pvalue<0.05, "*", "")))) -> Plot_significance
Plot_significance$Regressor[Plot_significance$Regressor == "age_f"] = "Age Father" ;Plot_significance$Regressor[Plot_significance$Regressor == "age_m"] = "Age Mother" ;Plot_significance$Regressor[Plot_significance$Regressor == "smk15"] = "Mother smoke while pregnant" ; Plot_significance$Regressor[Plot_significance$Regressor == "smk13"] = "Father smoke" ;Plot_significance$Regressor[Plot_significance$Regressor == "Parental_smoking"] = "Any parent smoke" ; Plot_significance$Regressor[Plot_significance$Regressor == "smk14"] = "Mother smoke"
Plot_significance %>% select(Dependent, Regressor, Significance) %>% spread(Dependent, Significance) %>% as.data.frame() %>% column_to_rownames(var = "Regressor") -> annotation_heat


COl = viridis::viridis(10)[4:20]
pheatmap::pheatmap(wide_results2,color = COL,clustering_distance_rows="correlation",clustering_distance_cols="correlation",display_numbers = annotation_heat ) -> Heatmap_significant_ass
SAVE(PLOT = Heatmap_significant_ass,NAME ="Figures/Heatmap_assocaitions_FDR.pdf",Height = 5)

##PLOT2
symlog = function(x, C){    
  #Where x is your data with positive and negative values
  #where the scaling constant C determines the resolution of the data around zero.  
  
  t = sign(x)*log(1+abs(x)/10^C)
  return(t)
}
Model_results2%>%filter(! Regressor == "Smk_mother_childhood") %>% mutate(Up =Beta+2*SE, Down = Beta-2*SE ) -> Model_results2
N= -2

max(Model_results2$Up) ; min(Model_results2$Down)
Plot = c(-0.5, -0.25 ,-0.1, -0.05, -0.01, 0 , 0.01, 0.05, 0.1, 0.25 ,0.5, 1)
plot_ticks<-symlog(Plot, N)

B =symlog(Model_results2$Beta, N) ; S =symlog(Model_results2$SE, N) ; U= symlog(Model_results2$Up, N) ; D = symlog(Model_results2$Down, N)

Model_results2%>% mutate(Beta = B, SE= S,Up = U, Down= D ) %>% ggplot(aes(x=Regressor, y= Beta,col=Dependent)) + geom_point(position=position_dodge(width=0.8), size=2) +  
  coord_flip() + geom_hline(yintercept = symlog(0, N) ,linetype="dashed")  + theme_bw() +
  geom_errorbar( aes(ymin = Down, ymax = Up),width = 0.7,linetype = "dotted",position=position_dodge(width=0.8)) +  scale_color_manual(values = viridis::viridis(6)) + ylab("Estimated effect on telomere length") + 
  scale_y_continuous(breaks = plot_ticks, labels = Plot) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) -> Phenotype_associations
SAVE(PLOT = Phenotype_associations,NAME ="Results/Manuscript/Figures/Phenotype_associations.pdf", Height = 7)

#########################################################################
#####Combine all phenotypes and calculate total variance explained#######
#########################################################################

set.seed(2020)
R2_calc = function(Real,pred){
  rss = sum((pred - Real)^2)
  tss = sum((Real - mean(Real))^2)
  rsq = 1 - rss/tss
  return(rsq)
}


Variance_explained = tibble()

lapply(PRS$IID, FUN= function(x){ str_split(x,"1_")[[1]][2] }) %>% unlist() -> IDs ; PRS %>% mutate(ID = IDs, PRS_value = SCORESUM) %>% select(ID, PRS_value) %>% arrange(ID) -> PRS


Variance_explained2 = tibble()
for (Cell_line in colnames(Telomere_markers)){
  if (Cell_line == "Sample"){ next }
  Telomere_markers %>% select(Cell_line) %>% as_vector() %>% as.vector() %>% as.numeric() -> Dependent
  
  Cov -> Cov_to_add
  Phenotypes_significant %>% mutate(Dependent = Dependent) -> Entry_model
  
  
  left_join(Entry_model, Cov_to_add, by="ID") -> Entry_model
  Entry_model %>% drop_na() -> Entry_model
  
  #####
  #Age model
  lm(Dependent ~ Age, Entry_model) ->  Null_0
  Entry_model %>% mutate(Dependent  = Null_0$residuals) %>% select(-Age) -> Entry_model
  
  #Host factors model
  Entry_model %>% select(-c(Parental_smoking, smk13,age_f, age_m)) -> Entry_model2
  #1. Sex
  lm(Dependent ~ Sex, Entry_model) ->  r2_sex
  r2_sex = R2_calc(as.vector(Entry_model$Dependent), predict(r2_sex))
  #2. Cells
  cvfit = cv.glmnet(y=as.vector(Entry_model2$Dependent), x=as.matrix(select(Entry_model2,c("Sex", cell_types))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Estimation = predict(cvfit, as.matrix(select(Entry_model2,c("Sex", cell_types))), s="lambda.min")
  r2_cells = R2_calc(as.vector(Entry_model2$Dependent), Estimation)
  #3. BMI
  cvfit = cv.glmnet(y=as.vector(Entry_model2$Dependent), x=as.matrix(select(Entry_model2,-c(ID,Dependent))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Estimation = predict(cvfit, as.matrix(select(Entry_model2,-c(ID,Dependent))), s="lambda.min")
  r2_BMI = R2_calc(as.vector(Entry_model2$Dependent), Estimation)
  
  #Model complete (+Parental factors)
  cvfit = cv.glmnet(y=as.vector(Entry_model$Dependent), x=as.matrix(select(Entry_model,-c(Dependent,ID))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Param_best = coef(cvfit, s = "lambda.min")
  Estimation = predict(cvfit, as.matrix(select(Entry_model,-c(Dependent,ID))), s="lambda.min")
  r2 = R2_calc(as.vector(Entry_model$Dependent), Estimation)
  
  #Model with Genetics 
  left_join(Entry_model, PRS, by="ID") -> Entry_model
  Entry_model %>% mutate(PRS_value = as.numeric(PRS_value)) %>% drop_na() -> Entry_model
  cvfit = cv.glmnet(y=as.vector(Entry_model$Dependent), x=as.matrix(select(Entry_model,-c(Dependent,ID))), nfolds = 10,type.measure="mse",standardize=T, alpha=1)
  Param_best = coef(cvfit, s = "lambda.min")
  Estimation = predict(cvfit, as.matrix(select(Entry_model,-c(Dependent,ID))), s="lambda.min")
  r2_genetics = R2_calc(as.vector(Entry_model$Dependent), Estimation)
  #Output
  rbind(Variance_explained2,tibble(Telomere=Cell_line, R2_parents= r2, R2_Sex = r2_sex, R2_Cells = r2_cells, R2_BMI = r2_BMI, R2_genetics = r2_genetics)) -> Variance_explained2
}

PALETTE  =  c("#000004FF", "#56106EFF", "#BB3754FF", "#F98C0AFF", "#FCFFA4FF")
Variance_explained2 %>% mutate(Genetics = R2_genetics - R2_parents, Parental = R2_parents-R2_BMI, BMI= R2_BMI-R2_Cells, Cells = R2_Cells - R2_Sex, Sex=R2_Sex) %>% select(-c("R2_parents","R2_Sex","R2_BMI","R2_Cells","R2_genetics")) -> MODEL
MODEL %>% mutate(Intrinsic = BMI+Cells+Sex ) %>% summarise(mean(Parental), mean(Genetics), mean(Intrinsic))

MODEL %>%
  gather("Model","R2", c("Parental","Sex","BMI","Cells","Genetics")) %>% mutate(Model = fct_relevel(Model,c("Genetics","Parental", "BMI","Cells","Sex"))) %>%
  mutate(R2 = ifelse(R2 < 0,0,R2)) %>%
  mutate(Telomere=factor(Telomere, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>%
  ggplot(aes(x=Telomere, y=R2, fill=Model)) + geom_bar(position="stack", stat="identity",width = 0.5) +theme_bw(base_size=12) +
  coord_flip() + theme_bw() + scale_fill_manual(values = viridis::viridis(5)) -> partition_R2
SAVE(NAME = "Figures/partition_R2_noage.pdf" ,PLOT=partition_R2,Height = 2)



###################
##Heritability####
##################

##Heritability plot. Requires heritability estimations

heritability_summary = read_tsv("Heritability_total.tsv") 
heritability_summary$Cell_type = c("Naïve T-cells", "Memory T-cells","Lymphocytes","B-cells", "Granulocytes","NK-cells")

heritability_summary %>% mutate(Error = ifelse((h2+Sd) > 1, 1, h2+Sd )) %>% 
  #mutate(Cell_line=factor(Cell_line, levels= c("Naïve T-cells","B-cells","Granulocytes", "Lymphocytes","Memory T-cells", "NK-cells"))) %>%
  mutate(Cell_type=factor(Cell_type, levels= c("Granulocytes","Lymphocytes","B-cells","Naïve T-cells","Memory T-cells","NK-cells"))) %>%
  ggplot(aes(x = Cell_type, y = h2, fill=Cell_type)) + geom_bar(stat="identity",position = "dodge",col="black", width = 0.5) + coord_flip() +
  geom_errorbar(aes(ymin=h2, ymax=Error), width=.2,position=position_dodge(.9)) + theme_bw(base_size = 12) +
  scale_fill_manual(values = viridis::viridis(6)) + xlab(label = "Telomere") + theme(legend.position = "none")  -> Heritability_plot
SAVE(NAME = Make_path("Figures/Heritability_GREML.pdf") ,PLOT = Heritability_plot, Height = 2)


#########################
#Prediction of mortality#
#########################
#Prediction (association) of death. Requires extra information.
dead = as_tibble(read.table(Make_path("LLD_passed.txt"),header = T))
Info = read_tsv(Make_path("baseline.txt"))
End_date = as.Date("2021-06-08")
dead   %>% mutate (DATE = as.Date(paste(deathyear,deathmonth,"01",sep="-") )) -> dead
colnames(dead)[2] = "LLDeep_SampleID"

Info %>% mutate(First_time = as.Date(paste(bl1year,"01-01",sep="-")) ) -> Info

Second_time = c()
for (ID in Info$LLDEEPID){
  Info[Info$LLDEEPID==ID,]$First_time -> First_time
  if(ID %in% dead$LLDeep_SampleID){
    dead[dead$LLDeep_SampleID==ID,]$DATE -First_time -> Datet
  } else { End_date-First_time -> Datet }
  Second_time = c(Second_time, Datet)
}
Info %>% mutate(Second_time = Second_time, Status = ifelse(LLDEEPID %in% dead$LLDeep_SampleID, 1, 0)) -> Info



Predictors = Telomere_markers

Cov %>% filter(ID %in% Info$LLDEEPID) %>% arrange(ID) -> Cov_s
Predictors %>% filter(Sample %in% Info$LLDEEPID) %>% filter(Sample %in% Cov$ID) %>% arrange(Sample) -> Predictors

Results = tibble()
for( Phenotype in colnames(select(Predictors,-c("Sample")))){
  print(Phenotype)
  Regressor = Predictors %>% select(Phenotype) %>% as_vector() %>% as.numeric()
  Info %>% filter(LLDEEPID %in% Predictors$Sample) -> Info2 ;  Cov_s %>% filter(ID  %in% Info2$LLDEEPID) -> Cov2
  cbind(Info2 %>% select(Second_time,Status), select(Cov2, - ID)) %>% as_tibble() %>% mutate(Pheno=Regressor) ->  All_info
  All_info %>% group_by(Status) %>% summarise(n()) -> Number_samples
  coxph(Surv(Second_time, Status) ~ ., data=All_info) -> res.cox
  test.ph <- cox.zph(res.cox)
  confint(res.cox)["Pheno",] -> CI

  as.data.frame(summary(res.cox)$coefficients)["Pheno",] %>% as_tibble() %>% mutate(Regressor = Phenotype) -> res.cox
  mutate(res.cox, CI_low = exp(CI[[1]]), CI_high = exp(CI[[2]]), N_cases= Number_samples[2,2]$`n()`, N_controls = Number_samples[1,2]$`n()` ) -> res.cox
  res.cox %>% mutate(P_ProportionalHazard_Telomere = test.ph$table["Pheno", "p"], P_ProportionalHazard_Global = test.ph$table["GLOBAL", "p"]) -> res.cox
  rbind(res.cox, Results) -> Results
}
arrange(Results, `Pr(>|z|)`)


  Results %>% mutate(`Hazard ratio` = as.numeric(`exp(coef)`), Pvalue= as.numeric(`Pr(>|z|)`), Telomere=Regressor) %>% 
  ggplot(aes(y=`Hazard ratio`,x= Telomere)) + geom_point(aes(col = -log10(Pvalue)), size=4, shape=15) +
  coord_flip() + theme_bw(base_size = 12) + geom_hline(yintercept=1, linetype="dotted", size=1.5) +
  geom_errorbar(aes(ymin=as.numeric(CI_low), ymax=as.numeric(CI_high)),width=.2) +
  scale_colour_viridis_c() -> Forest_mortality


SAVE(PLOT = Forest_mortality,NAME ="Figures/mortality_associations.pdf", Height = 2)
SAVE(PLOT = Forest_mortality,NAME ="Figures/mortality_associations.png", Height = 2)
write_tsv(Results, "Results/Cox_regression.tsv")
  
                                            


##########################
##Plot data availability #
##########################
#Requires info about data availability and UpsetR package
Info_avail = read.table(Make_path("Genetic_and_methylation_availability.tsv"),header = T)
Telomere_markers <- colnames(Age_markers)[grepl("MTL",colnames(Age_markers))]
Age_markers %>% dplyr::select(c("Sample",Telomere_markers)) %>% mutate(MTL_gran = as.numeric(MTL_gran), `MTL_CD20+`=as.numeric(`MTL_CD20+`), `MTL_CD45+_CD20-`=as.numeric(`MTL_CD45+_CD20-`),`MTL_CD57+`=as.numeric(`MTL_CD57+`)) %>%
  filter(! (is.na(MTL_lymp) & is.na(MTL_gran) & is.na(`MTL_CD45+_CD20-`) & is.na(`MTL_CD45-`) & is.na(`MTL_CD20+`) & is.na(`MTL_CD57+`))) -> Info_tel
Info_tel2 = drop_na(Info_tel)
Info_avail %>% mutate(Telomere_ID = ifelse(DEEPID %in% Info_tel2$Sample, 1, 0)   ,Telomere_any = ifelse(DEEPID %in% Info_tel$Sample, 1, 0)) -> Info_avail

Age_markers %>% mutate(DEEPID= Sample) %>%  select(DEEPID,Methyl.Age.Hannum., dCT) -> Add_markers
left_join(Info_avail, Add_markers) -> Info_avail

SC_available = read_tsv(Make_path("Participants_with_transcriptomics.txt"), col_names = F)
Info_avail %>% mutate(Single_cell = ifelse(DEEPID %in% SC_available$X1, 1, NA) ) -> Info_avail

apply(Info_avail,2, FUN= as.character) %>% as.matrix() -> Info_avail
Info_avail[! Info_avail==0] = 1 ; Info_avail[Info_avail==0] = 0 ; Info_avail[is.na(Info_avail)] = 0

apply(Info_avail, 2, as.numeric) -> Info_avail
colnames(Info_avail) = c("LLDeep", "Methylation", "Genetics", "Olink", "Telomere all", "Telomere any", "Methylation age","sjTREC", "Single cell")

Up_input <- select(filter(as.data.frame(Info_avail), `Telomere any`==1), -LLDeep) %>% select(-Olink)
Up_input %>% summarise_all(sum) %>% t()
UpSetR::upset(Up_input ,nsets = 10,text.scale = 2, order.by = "freq") -> Data_avail

pdf(file="Figures/Data_availability_UpSet.pdf",width = 12)
Data_avail
dev.off()
