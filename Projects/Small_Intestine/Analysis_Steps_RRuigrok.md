Small Intestine Microbiota Project
==================================

Creator: Renate Ruigrok
Year: 2020

### 1. Cohort Summary Statistics

```{r}
# Functions Required:

library(forcats)
summary_statistics_metadata   #see 'Functions' file, point 1.
Remove_FactCols #see 'Functions' file, point 2.
kruskal.test_function. #see 'Functions' file, point 3.


# Data Import & Processing

LLD_IBD_MetTaxS <- read.table(file = "../RScripts/LLD_IBD_MetaSpecies150420.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) 
LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS[-c(212,218,379,1554,1671),]
LLD_IBD_MetTaxS <- subset(LLD_IBD_MetTaxS, DiagnosisCurrent != "MicroscopicColitis" & DiagnosisCurrent != "IBDI")
LLD_IBD_MetTaxS$DiagnosisCurrent <- fct_drop(LLD_IBD_MetTaxS$DiagnosisCurrent)

LLD_IBD_MetTaxS$CurrentStomaOrPouchType <- as.factor(gsub("Resection NoStoma", "Resections", LLD_IBD_MetTaxS$CurrentStomaOrPouchType))
LLD_IBD_MetTaxS$CurrentStomaOrPouchType <- as.factor(gsub("colostoma", "Resections", LLD_IBD_MetTaxS$CurrentStomaOrPouchType))

#Subsetting
LLD.IBD_MetTaxS.F <- LLD_IBD_MetTaxS[c(2,4:6,8,9,14:17,20,22,24:28,30:37,42:48,50:69,71:86,88:181,183:197,202,203,408:416,418:425,427:678,681:1332)]

x <- gsub('CU','UC', LLD.IBD_MetTaxS.F$DiagnosisFirst) 
LLD.IBD_MetTaxS.F$DiagnosisFirst <- as.factor(x)

LLD.IBD_MetTaxS.F$DiagnosisCurrent <- as.factor(gsub("ReconsideringDiagnosis", "IBDU", LLD.IBD_MetTaxS.F$DiagnosisCurrent))

LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation <- as.numeric(as.character(LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation))


LLD.IBD_MetTaxS.F$ReasonUnclearSurgeryResectionProcedures <- NULL


LLD.IBD_Final <- Remove_FactCols(LLD.IBD_MetTaxS.F) 

LLD.IBD_Final <- LLD.IBD_Final[c(1:114,142:1080)] #remove logical variables, not informative in this case 

#library(forcats)
LLD.IBD_Final[408,128] <- NA #change 'NA' to NA 
LLD.IBD_Final$NumberIndicationabcessOrfistula <- fct_drop(LLD.IBD_Final$NumberIndicationabcessOrfistula) #remove level 'NA'

LLD.IBD_Final$IncludedSamples <- NULL

#include BinaryBSS and bowel movement variable to LLD.IBD_Final dataframe
BinaryBSS <- read.table(file = "../Data/SI_project/intestinal_info.txt", header = TRUE, sep = "\t", quote = "\\")
LLD.IBD_Final <- merge(BinaryBSS, LLD.IBD_Final, by.x = 'UMCGIBDResearchIDorLLDeepID', by.y = 'row.names')
rownames(LLD.IBD_Final) <-LLD.IBD_Final$UMCGIBDResearchIDorLLDeepID
LLD.IBD_Final$UMCGIBDResearchIDorLLDeepID <- NULL
LLD.IBD_Final <- LLD.IBD_Final[c(3:388,1,2,389:1054)] #change order

#include updated PFReads variable to LLD.IBD_Final dataframe 
PFReads <- read.table(file = "../Data/SI_project/PFReads31.03.20.txt", header = TRUE, sep = "\t", quote = "\\")
LLD.IBD_Final$PFReads <- NULL
LLD.IBD_Final <- merge(PFReads, LLD.IBD_Final, by.x = 'SampleID', by.y = 'row.names')
rownames(LLD.IBD_Final) <-LLD.IBD_Final$SampleID
LLD.IBD_Final$SampleID <- NULL
LLD.IBD_Final <- LLD.IBD_Final[c(2,1,3:1054)] #change order
colnames(LLD.IBD_Final)[2] <- "PFReads"

#Final data with ileostoma and pouch individuals as one group
LLD.IBD_Finalx <-  LLD.IBD_Final
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("No Resections", "Colon", LLD.IBD_Finalx$CurrentStomaOrPouchType))
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("Resections", "Colon", LLD.IBD_Finalx$CurrentStomaOrPouchType))

LLD.IBD_Finalx$NumberIndicationabcessOrfistula <- as.numeric(as.character(LLD.IBD_Finalx$NumberIndicationabcessOrfistula))

Final_phenotypes1 <- read.table(file = "UnivariateAssociations_asin.txt",header = TRUE, sep = "\t", row.names = 1)
Final_phenotypes1 <- as.character(unique(Final_phenotypes1$Phenotype))
LLD.IBD.Stats <- LLD.IBD_Finalx[c(Final_phenotypes1)]

# Run Summary statistics function

StomaCat_table <- LLD.IBD.Stats[c("CurrentStomaOrPouchType")] #Stoma type category table
LLD.IBD_Meta_Input <- LLD.IBD.Stats[,-c(45)]

summary_statistics_metadata(metadata_input = LLD.IBD_Meta_Input, category_table = StomaCat_table)


# Kruskal-Wallis test

#IBD vs SI
LLD.IBD.Stats1 <- subset(LLD.IBD.Stats, CurrentStomaOrPouchType %in% c('Colon', 'IleostomaPouch'))
LLD.IBD.Stats1$CurrentStomaOrPouchType <- fct_drop(LLD.IBD.Stats1$CurrentStomaOrPouchType)

LLD.IBD.Stats1a <- LLD.IBD.Stats1[,c(45,1:44,46:107,112:119)]

Test1 <- kruskal.test_function(LLD.IBD.Stats1a)

write.table(Test1, "Pval_DescrpStats_SIvsIBD.txt", sep = "\t",col.names = T)

#LLD vs SI
LLD.IBD.Stats2 <- subset(LLD.IBD.Stats, CurrentStomaOrPouchType %in% c("IleostomaPouch",  'HealthyControl'))
LLD.IBD.Stats2$CurrentStomaOrPouchType <- fct_drop(LLD.IBD.Stats2$CurrentStomaOrPouchType)

LLD.IBD.Stats2a <- LLD.IBD.Stats2[,c(45,1:3,9,14,24,36,41,43,44,62:69,74:107,119)]

Test2 <- kruskal.test_function(LLD.IBD.Stats2a)

write.table(Test2, "Pval_DescrpStats_SIvsLLD.txt", sep = "\t",col.names = T)
```

### 2. Genus Composition  

```{r}
# Packages & Functions Required:

library(dplyr)
library(ggplot2)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(forcats)
TopNTaxa  #see 'Functions' file, point 4.



# Data Import & Process

#metadata (LLD and IBD faecal ID's combined in one column to allow for merge (see following steps))
LLD.IBD_Meta <- read.table(file = "../Data/SI_project/Metadata27022019_new.txt", header = TRUE, sep = "\t", quote = "\\") #Metadata 

LLD.IBD_Meta <- LLD.IBD_Meta[-c(1783,1793,2000,1349,1839),]

LLD.IBD_Meta <- subset(LLD.IBD_Meta, DiagnosisCurrent != "MicroscopicColitis" & DiagnosisCurrent != "IBDI")

LLD.IBD_Meta$DiagnosisCurrent <- fct_drop(LLD.IBD_Meta$DiagnosisCurrent)

#To distinguish patients with no stoma and no resections from patients with resections 
LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="yes", gsub('none', 'Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="no", gsub('none', 'No Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(LLD.IBD_Meta$CurrentStomaOrPouchType)

#Add additional level to 'CurrentStomaOrPouchType' column to represent LLD individuals 
levels <- levels(LLD.IBD_Meta$CurrentStomaOrPouchType)
levels[length(levels) + 1] <- "HealthyControl"

# refactor CurrentStomaOrPouchType to include "HealthyControl" as a level
# and replace NAs with "HealthyControl"
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = levels)
LLD.IBD_Meta$CurrentStomaOrPouchType[is.na(LLD.IBD_Meta$CurrentStomaOrPouchType)] <- "HealthyControl"

#Final data with ileostoma and pouch individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Final data with colotoma and resection individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("colostoma", "Resections", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Change order of the levels for variable CurrentStomaOrPouchType
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = c('HealthyControl', 'No Resections', 'Resections', 'IleostomaPouch'))

#IBD Taxonomy 
IBD_Taxa <- read.table(file = "../Data/SI_project/IBD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") 
rownames(IBD_Taxa) <- IBD_Taxa$ID
IBD_Taxa$ID <- NULL
IBD_Taxa <- as.data.frame(t(IBD_Taxa)) #transform

#LLD Taxonomy
LLD_Taxa <- read.table(file = "../Data/SI_project/LLD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #LLD Taxonomy 
rownames(LLD_Taxa) <- LLD_Taxa$ID
LLD_Taxa$ID <- NULL

LLD <- gsub("^X", "", colnames(LLD_Taxa)) #Remove 'X' at the beginning of ID numbers
colnames(LLD_Taxa) <- LLD

#Merge the two cohort taxa's
IBD_Taxa1 <- as.data.frame(t(IBD_Taxa))
LLD_IBD_Taxa <- merge(LLD_Taxa, IBD_Taxa1, by = 'row.names', all = T)
rownames(LLD_IBD_Taxa) <- LLD_IBD_Taxa$Row.names
LLD_IBD_Taxa$Row.names <- NULL
LLD_IBD_Taxa <- as.data.frame(t(LLD_IBD_Taxa))

#Merge metadata with taxa
LLD_IBD_MetTax <- merge(LLD.IBD_Meta,LLD_IBD_Taxa, by.x = 'CombinedIDs', by.y = 'row.names')


****



# Top Genus composition according to Healthy control group
  
#Extract genus only from Taxa table
LLD_IBD_TaxaG <- LLD_IBD_Taxa[,grepl('g__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaG <- LLD_IBD_TaxaG[,!grepl('s__',colnames(LLD_IBD_TaxaG),ignore.case = T)]
  
  
#Merge genus Taxo with LLD_IBD_MetaTaxa
  LLD_IBD_MetTaxG <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaG, by.x = 'CombinedIDs', by.y = 'row.names')
  rownames(LLD_IBD_MetTaxG) <- LLD_IBD_MetTaxG$CombinedIDs
  LLD_IBD_MetTaxG$CombinedIDs <- NULL
  
  remove_cols <- vector()
  for (i in 2:ncol(LLD_IBD_MetTaxG)) {
    cname <- colnames(LLD_IBD_MetTaxG)[i]
    if(colSums(LLD_IBD_MetTaxG[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxG) *100 < 15){
      remove_cols <- c(remove_cols, cname) 
    }
  }
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG %>% select(-remove_cols)
  print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
  
  #subset 'HealthyControl' individuals only
  HC <- subset(LLD_IBD_MetTaxG, CurrentStomaOrPouchType == 'HealthyControl')
  
  #Vector with mean relative abundances in descending order
  Means <- colMeans(HC[,c(2:ncol(HC))]) #mean
  Means <- sort(Means, decreasing = T) #order
  
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]
  
  LLD_IBD_MetTaxG$g__Other <- 100 - rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) # minus numbers?
  
  #Plot graph by stoma type
  #column mean by group
  LLD_IBD_MetTaxG <- aggregate(LLD_IBD_MetTaxG[,2:12], list(LLD_IBD_MetTaxG$CurrentStomaOrPouchType), mean)
  colnames(LLD_IBD_MetTaxG)[1] <- "CurrentStomaOrPouchType"
  
  #Convert each mean relative abundance to a percentage of total means 
  LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))] <- LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]/
    rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))])*100
  
  rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage
  
  LLD <- gsub('.*\\|g__','',colnames(LLD_IBD_MetTaxG)) #Keep taxonomy level name only i.e. 'class name' 
  colnames(LLD_IBD_MetTaxG) <- LLD
  
  #Graph
  LLD_IBD_MetTaxG<- melt(LLD_IBD_MetTaxG, 'CurrentStomaOrPouchType')
  
  LLD_IBD_MetTaxG %>%
    group_by(CurrentStomaOrPouchType) %>%
    summarise(Sum = sum(value)) #Check each group sums to 100
  
cols <- c("Eubacterium" = "#66C2A5", "Bifidobacterium" = "#FC8D62", "Ruminococcus" = "#8DA0CB", "Blautia" = "#FFFF33", "Subdoligranulum" = "#E78AC3", "Faecalibacterium" = "#A6D854", "g__Other" = "#D9D9D9", "Coprococcus" = "#377EB8", "Bacteroides" = "#E5C494", "Enterococcus" = "#483D8B", "Dorea" = "#EF3A3A", "Collinsella" = "#9846AC" )
  
  ggplot(data = LLD_IBD_MetTaxG, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
    geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +
    labs(y = "Mean Relative Abundance", x = 'Current Phenotype', title = 'Microbial Genus Composition Top HealthyControl') + scale_x_discrete(labels=c('HealthyControl' = 'General population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine')) + theme_minimal() 
  #+ theme(axis.text.x = element_text(angle = 45)) 
  
  write.table(LLD_IBD_MetTaxG, file = 'LLD_IBD_GenusComposition_top10HC.txt', sep = '\t', col.names = T)
  

  
#Top Genus composition according to IleostomaPouch group
  
  #Merge Genera with LLD_IBD_MetaTaxa
  LLD_IBD_MetTaxG <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaG, by.x = 'CombinedIDs', by.y = 'row.names')
  rownames(LLD_IBD_MetTaxG) <- LLD_IBD_MetTaxG$CombinedIDs
  LLD_IBD_MetTaxG$CombinedIDs <- NULL
  
  remove_cols <- vector()
  for (i in 2:ncol(LLD_IBD_MetTaxG)) {
    cname <- colnames(LLD_IBD_MetTaxG)[i]
    if(colSums(LLD_IBD_MetTaxG[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxG) *100 < 15){
      remove_cols <- c(remove_cols, cname) 
    }
  }
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG %>% select(-remove_cols)
  print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
  #subset ileostomaPouch individuals only
  Ileo <- subset(LLD_IBD_MetTaxG, CurrentStomaOrPouchType == 'IleostomaPouch')
  
  #Vector with mean relative abundances in descending order
  Means <- colMeans(Ileo[,c(2:ncol(Ileo))]) #mean
  Means <- sort(Means, decreasing = T) #order
  
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]
  
  LLD_IBD_MetTaxG$g__Other <- 100 - rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) 
  
  #Plot graph by stoma type
  
  #column mean by group
  LLD_IBD_MetTaxG <- aggregate(LLD_IBD_MetTaxG[,2:12], list(LLD_IBD_MetTaxG$CurrentStomaOrPouchType), mean)
  colnames(LLD_IBD_MetTaxG)[1] <- "CurrentStomaOrPouchType"
  
  #Convert each mean relative abundance to a percentage of total means 
  LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))] <- LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]/
    rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))])*100
  
  rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage
  
  LLD <- gsub('.*\\|g__','',colnames(LLD_IBD_MetTaxG)) #Keep taxonomy level name only i.e. 'class name' 
  colnames(LLD_IBD_MetTaxG) <- LLD
  
  #Graph
  LLD_IBD_MetTaxG<- melt(LLD_IBD_MetTaxG, 'CurrentStomaOrPouchType')
  
  LLD_IBD_MetTaxG %>%
    group_by(CurrentStomaOrPouchType) %>%
    summarise(Sum = sum(value)) #Check each group sums to 100
  
  cols <- c("Eubacterium" = "#66C2A5", "Bifidobacterium" = "#FC8D62", "Ruminococcus" = "#8DA0CB", "Blautia" = "#FFFF33", "Subdoligranulum" = "#E78AC3", "Faecalibacterium" = "#A6D854", "g__Other" = "#D9D9D9", "Coprococcus" = "#377EB8", "Bacteroides" = "#E5C494", "Streptococcus" = "#E41A1C", "Escherichia" = "#C2A5CF", "Peptostreptococcaceae_noname" = "#BF5B17", "Lactobacillus" = "#CCEBC5", "Veillonella" = "#A6CEE3", "Clostridium" = "#F0027F", "Enterococcus" = "#483D8B")
  
  ggplot(data = LLD_IBD_MetTaxG, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
    geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +
    labs(y = "Mean Relative Abundance", x = 'Current Phenotype', title = 'Microbial Genus Composition Top SI') + scale_x_discrete(labels=c('HealthyControl' = 'General population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine')) + theme_minimal() 
  #+ theme(axis.text.x = element_text(angle = 45)) 
  
  
write.table(LLD_IBD_MetTaxG, file = 'LLD_IBD_GenusComposition_top10SI.txt', sep = '\t', col.names = T)
 
```


### 3. Genus Composition: Wilcoxon test

```{r}
Top10HCGenera <- read.table(file = "LLD_IBD_GenusComposition_top10HC.txt", header = TRUE, sep = "\t",row.names = 1)
Top10HCGenera <- Top10HCGenera[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10HCGenera[c(2:ncol(Top10HCGenera))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10HCGenera)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10HCGenera[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10HCGenera[,j], Top10HCGenera[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10HCGenera)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10HCGenera)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10HCGenera[,k])
      my_matrix[a,6] <- sum(Top10HCGenera[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10HCGenera[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10HCGenera[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10HCGenera, Top10HCGenera[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10HCGenera, Top10HCGenera[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
  }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_HC_Genera.txt", sep = "\t", col.names = T)


#Top SI Genera
Top10SIGenera <- read.table(file = "LLD_IBD_GenusComposition_top10SI.txt", header = TRUE, sep = "\t",row.names = 1)
Top10SIGenera <- Top10SIGenera[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10SIGenera[c(2:ncol(Top10SIGenera))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10SIGenera)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10SIGenera[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10SIGenera[,j], Top10SIGenera[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10SIGenera)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10SIGenera)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10SIGenera[,k])
      my_matrix[a,6] <- sum(Top10SIGenera[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10SIGenera[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10SIGenera[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10SIGenera, Top10SIGenera[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10SIGenera, Top10SIGenera[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_SI_Genera.txt", sep = "\t", col.names = T)

```

### 4. Alpha & Beta Diversity 

```{r}

# Packages & Functions Required:
library(vegan)
library(dplyr)
library(ggplot2)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(ggpubr)
library(forcats)
plot_pcoa #see 'Functions' file, point 5.

# Data Import & Process

#_Metadata
#metadata (LLD and IBD faecal ID's combined in one column to allow for merge (see following steps))
LLD.IBD_Meta <- read.table(file = "../Data/SI_project/Metadata27022019_new.txt", header = TRUE, sep = "\t", quote = "\\") #Metadata 

LLD.IBD_Meta <- LLD.IBD_Meta[-c(1783,1793,2000,1349,1839),]

LLD.IBD_Meta <- subset(LLD.IBD_Meta, DiagnosisCurrent != "MicroscopicColitis" & DiagnosisCurrent != "IBDI")

LLD.IBD_Meta$DiagnosisCurrent <- fct_drop(LLD.IBD_Meta$DiagnosisCurrent)

#To distinguish patients with no stoma and no resections from patients with resections 
LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="yes", gsub('none', 'Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="no", gsub('none', 'No Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(LLD.IBD_Meta$CurrentStomaOrPouchType)

#Add additional level to 'CurrentStomaOrPouchType' column to represent LLD individuals 
levels <- levels(LLD.IBD_Meta$CurrentStomaOrPouchType)
levels[length(levels) + 1] <- "HealthyControl"

# refactor CurrentStomaOrPouchType to include "HealthyControl" as a level
# and replace NAs with "HealthyControl"
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = levels)
LLD.IBD_Meta$CurrentStomaOrPouchType[is.na(LLD.IBD_Meta$CurrentStomaOrPouchType)] <- "HealthyControl"

#Final data with ileostoma and pouch individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Final data with colotoma and resection individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("colostoma", "Resections", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Change order of the levels for variable CurrentStomaOrPouchType
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = c('HealthyControl', 'No Resections', 'Resections', 'IleostomaPouch'))


#_IBD_Taxonomy
IBD_Taxa <- read.table(file = "../Data/SI_project/IBD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #IBD Taxonomy 
rownames(IBD_Taxa) <- IBD_Taxa$ID
IBD_Taxa$ID <- NULL
IBD_Taxa <- as.data.frame(t(IBD_Taxa)) #transform


#_LLD_Taxonomy
LLD_Taxa <- read.table(file = "../Data/SI_project/LLD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #LLD Taxonomy 
rownames(LLD_Taxa) <- LLD_Taxa$ID
LLD_Taxa$ID <- NULL

LLD <- gsub("^X", "", colnames(LLD_Taxa)) #Remove 'X' at the beginning of ID numbers
colnames(LLD_Taxa) <- LLD

#Merge the cohort taxa's
IBD_Taxa1 <- as.data.frame(t(IBD_Taxa))
LLD_IBD_Taxa <- merge(LLD_Taxa, IBD_Taxa1, by = 'row.names', all = T)
rownames(LLD_IBD_Taxa) <- LLD_IBD_Taxa$Row.names
LLD_IBD_Taxa$Row.names <- NULL
LLD_IBD_Taxa <- as.data.frame(t(LLD_IBD_Taxa))

#Species Taxa
LLD_IBD_TaxaS <- LLD_IBD_Taxa[,grepl('s__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaS<- LLD_IBD_TaxaS[,!grepl('t__',colnames(LLD_IBD_TaxaS),ignore.case = T)]

#Merge metadata with taxa
LLD_IBD_MetTax <- merge(LLD.IBD_Meta,LLD_IBD_TaxaS, by.x = 'CombinedIDs', by.y = 'row.names')

# Remove bacteria from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 748:ncol(LLD_IBD_MetTax)) {
  cname <- colnames(LLD_IBD_MetTax)[i]
  if(colSums(LLD_IBD_MetTax[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTax) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTax15 <- LLD_IBD_MetTax %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

rownames(LLD_IBD_MetTax15) <- LLD_IBD_MetTax15$CombinedIDs
LLD_IBD_MetTax15$CombinedIDs <- NULL

#Category table
Stoma1 <- LLD_IBD_MetTax15[c(169)]



# Alpha diversity (Shannon Index)

Alpha1 <- as.data.frame(diversity(LLD_IBD_MetTax15[,c(747:880)], index="shannon"))
colnames(Alpha1)[1] <- 'AlphaDiversity' 
  
Alpha_Stoma1 <- merge(Stoma1, Alpha1, by = 'row.names')
rownames(Alpha_Stoma1) <- Alpha_Stoma1$Row.names
Alpha_Stoma1$Row.names <- NULL

 my_comparisions=list(c("HealthyControl","IleostomaPouch"),c("Resections","IleostomaPouch"),c("No Resections","IleostomaPouch"))

m <- ggplot(Alpha_Stoma1, aes (x=CurrentStomaOrPouchType, y=AlphaDiversity, fill=CurrentStomaOrPouchType)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white") + theme_minimal() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 0.5,vjust = 1, size=8,color="black")) + scale_fill_manual(values=c('grey54','slateblue2','gold1','red2')) + stat_compare_means(comparisons = my_comparisions, method = "wilcox.test") + ylab ("Shannon Index") + xlab ("Bowel Phenotype") + scale_x_discrete(labels=c('HealthyControl' = 'general population\n N = 1178', 'No Resections' = 'IBD non-resected bowel\n N = 309', 'Resections' = 'IBD resected bowel\n N = 169', 'IleostomaPouch' = 'IBD small intestine\n N = 57'))

ggplotly(m)

table(Alpha_Stoma1$CurrentStomaOrPouchType)

write.table(Alpha_Stoma1, file = "LLD_IBD_ShannonDiversity_ileopouch.txt", sep = "\t", col.names = T)

Alpha_Stoma1x<- melt(Alpha_Stoma1, 'CurrentStomaOrPouchType')
Alpha_Stoma1x <- aggregate(Alpha_Stoma1[, 2], list(Alpha_Stoma1$CurrentStomaOrPouchType), mean)
colnames(Alpha_Stoma1x) <- c("CurrentStomaOrPouchType", "MeanShannonIndex")



#Beta diversity (PCoA, bray-curtis dissimilarity)

#PCoA function (plot_pcoa) required
#Taxa variables from LLD_IBD_MetTax required
#Category table: Stoma1 required

#Remove samples/observations whose taxa sums to 0 (function doesn't work otherwise)
#LLD_IBD_TaxaS1 <- LLD_IBD_TaxaS[rowSums(LLD_IBD_TaxaS) != 0,]

LLD_IBD_TaxaS1 <- LLD_IBD_MetTax15[,c(747:880)]

#Normalise data
LLD_IBD_TaxaS1 <- LLD_IBD_TaxaS1/100
LLD_IBD_TaxaS1 <- asin(sqrt(LLD_IBD_TaxaS1))

#Run PCoA function
plot_pcoa(LLD_IBD_TaxaS1,5,Stoma1) #with 5 PCoA elements - table renamed to pcoa_table5LLD

#Import table generated from 'plot_pcoa' function
PCoA_values <- read.table(file = "../RScripts/pcoa_tableALL.txt", header = TRUE, sep = "\t", quote = "\\") 
rownames(PCoA_values) <- PCoA_values$X
PCoA_values$X <- NULL

PCoA_values$CurrentStomaOrPouchType <- factor(PCoA_values$CurrentStomaOrPouchType, levels=c("HealthyControl", "No Resections", "Resections", "IleostomaPouch"))

PCoA_values1<- melt(PCoA_values, 'CurrentStomaOrPouchType')
PCoA_values2 <- aggregate(PCoA_values1[, 3], list(PCoA_values1$CurrentStomaOrPouchType,PCoA_values1$variable), mean)
colnames(PCoA_values2) <- c("CurrentStomaOrPouchType","PCoA", "MeanPCoA")

write.table(PCoA_values2, file = "PCoA_AllGroups_Means.txt", sep = "\t", col.names = T)

#Plot PCoA graphs
#PCoA1 vs PCoA2
h <- ggplot(PCoA_values, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + labs(y="PCoA2", x="PCoA1", color = "Bowel Phenotype") + geom_point(size = 2) + scale_color_manual(labels = c('HealthyControl' = 'general population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + theme_classic()  

#PCoA2 vs PCoA3
i <- ggplot(PCoA_values, aes(x=V2, y=V3, color= CurrentStomaOrPouchType)) + labs(y="PCoA3", x="PCoA2", color = "Bowel Phenotype") + geom_point(size = 2) + scale_color_manual(labels = c('HealthyControl' = 'general population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + theme_classic() 

#PCoA1 vs PCoA3
j <- ggplot(PCoA_values, aes(x=V1, y=V3, color= CurrentStomaOrPouchType)) + labs(y="PCoA3", x="PCoA1", color = "Bowel Phenotype") + geom_point(size = 2) + scale_color_manual(labels = c('HealthyControl' = 'general population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + theme_classic() 

#PCoA1 vs PCoA4
k <- ggplot(PCoA_values, aes(x=V1, y=V4, color= CurrentStomaOrPouchType)) + labs(y="PCoA4", x="PCoA1", color = "Bowel Phenotype") + geom_point(size = 2) + scale_color_manual(labels = c('HealthyControl' = 'general population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + theme_classic() 


ggplotly(h)
ggplotly(i)
ggplotly(j)
ggplotly(k)


#Wilcoxon test - comparing mean bray-curtis distances per group
??
 ```


### 5. Bacterial species abundance vs. number of resections

```{r}
#Functions required:
transform_and_filter_taxa #see 'Functions' file, point 6.

#__ Data import

#_Metadata
#metadata (LLD and IBD faecal ID's combined in one column to allow for merge (see following steps))
LLD.IBD_Meta <- read.table(file = "../Data/SI_project/Metadata27022019_new.txt", header = TRUE, sep = "\t", quote = "\\") #Metadata 

LLD.IBD_Meta <- LLD.IBD_Meta[-c(1783,1793,2000,1349,1839),]

LLD.IBD_Meta <- subset(LLD.IBD_Meta, DiagnosisCurrent != "MicroscopicColitis" & DiagnosisCurrent != "IBDI")

LLD.IBD_Meta$DiagnosisCurrent <- fct_drop(LLD.IBD_Meta$DiagnosisCurrent)

#To distinguish patients with no stoma and no resections from patients with resections 
LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="yes", gsub('none', 'Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="no", gsub('none', 'No Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(LLD.IBD_Meta$CurrentStomaOrPouchType)

#Add additional level to 'CurrentStomaOrPouchType' column to represent LLD individuals 
levels <- levels(LLD.IBD_Meta$CurrentStomaOrPouchType)
levels[length(levels) + 1] <- "HealthyControl"

# refactor CurrentStomaOrPouchType to include "HealthyControl" as a level
# and replace NAs with "HealthyControl"
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = levels)
LLD.IBD_Meta$CurrentStomaOrPouchType[is.na(LLD.IBD_Meta$CurrentStomaOrPouchType)] <- "HealthyControl"

#Final data with ileostoma and pouch individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Final data with colotoma and resection individuals as one group
LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(gsub("colostoma", "Resections", LLD.IBD_Meta$CurrentStomaOrPouchType))

#Change order of the levels for variable CurrentStomaOrPouchType
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = c('HealthyControl', 'No Resections', 'Resections', 'IleostomaPouch'))


#_IBD_Taxonomy
IBD_Taxa <- read.table(file = "../Data/SI_project/IBD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #IBD Taxonomy 
rownames(IBD_Taxa) <- IBD_Taxa$ID
IBD_Taxa$ID <- NULL
IBD_Taxa <- as.data.frame(t(IBD_Taxa)) #transform


#_LLD_Taxonomy
LLD_Taxa <- read.table(file = "../Data/SI_project/LLD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #LLD Taxonomy 
rownames(LLD_Taxa) <- LLD_Taxa$ID
LLD_Taxa$ID <- NULL

LLD <- gsub("^X", "", colnames(LLD_Taxa)) #Remove 'X' at the beginning of ID numbers
colnames(LLD_Taxa) <- LLD

#Merge the cohort taxa's
IBD_Taxa1 <- as.data.frame(t(IBD_Taxa))
LLD_IBD_Taxa <- merge(LLD_Taxa, IBD_Taxa1, by = 'row.names', all = T)
rownames(LLD_IBD_Taxa) <- LLD_IBD_Taxa$Row.names
LLD_IBD_Taxa$Row.names <- NULL
LLD_IBD_Taxa <- as.data.frame(t(LLD_IBD_Taxa))

#Species Taxa
LLD_IBD_TaxaS <- LLD_IBD_Taxa[,grepl('s__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaS<- LLD_IBD_TaxaS[,!grepl('t__',colnames(LLD_IBD_TaxaS),ignore.case = T)]
#LLD_IBD_TaxaS[is.na(LLD_IBD_TaxaS)] <- 0

#Merge metadata with taxa
LLD_IBD_MetTaxS <- merge(LLD.IBD_Meta,LLD_IBD_TaxaS, by.x = 'CombinedIDs', by.y = 'row.names')

x <- gsub('.*\\|s__','',colnames(LLD_IBD_MetTaxS)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxS) <- x

rownames(LLD_IBD_MetTaxS) <- LLD_IBD_MetTaxS$CombinedIDs
LLD_IBD_MetTaxS$CombinedIDs <- NULL

# Remove bacteria from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 747:ncol(LLD_IBD_MetTaxS)) {
  cname <- colnames(LLD_IBD_MetTaxS)[i]
  if(colSums(LLD_IBD_MetTaxS[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxS) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTax15 <- LLD_IBD_MetTaxS %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

#Asin Transform
LLD_IBD_MetTax.asin <- transform_and_filter_taxa(LLD_IBD_MetTax15[,c(747:880)])
LLD_IBD_MetTax.asin <- cbind(LLD_IBD_MetTax15[c(1:746)], LLD_IBD_MetTax.asin)


#subset Resection Individuals only
IBD_Resection_MetTaxS <- subset(LLD_IBD_MetTax.asin, CurrentStomaOrPouchType == 'Resections')

#For loop to correlate Number of Resections with each bacterial species
my_matrix <-  matrix(nrow = ncol(IBD_Resection_MetTaxS[c(747:ncol(IBD_Resection_MetTaxS))]), ncol = 6)
colnames(my_matrix) <- c('Taxonomy', 'CorrelationCoeff', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'NSamples')
a=1
for (i in 747:ncol(IBD_Resection_MetTaxS)) {
  my_matrix[a,1] <- colnames(IBD_Resection_MetTaxS)[i]
  my_matrix[a,2] <- (cor.test(IBD_Resection_MetTaxS[,i], IBD_Resection_MetTaxS[,228], method = 'spearman', exact = F))$estimate #Spearman correlation coefficient values 
  my_matrix[a,3] <- (cor.test(IBD_Resection_MetTaxS[,i], IBD_Resection_MetTaxS[,228], method = 'spearman', exact = F))$p.value #Spearman correlation p-values
  a=a+1
  }
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,4] <- p.adjust(c(my_matrix[,3]), method = 'bonferroni')
my_matrix[,5] <- p.adjust(c(my_matrix[,3]), method = 'fdr')
my_matrix[,6] <- nrow(IBD_Resection_MetTaxS)

NumResectionsVsBacterialSpecies <- my_matrix

write.table(NumResectionsVsBacterialSpecies, file = 'Within_IBDResection_NumOfResectionsVsBacterialSpecies.asin.txt', sep = "\t", col.names = TRUE)

```

### 6. Bacterial species abundance vs. intestinal resection type

```{r}
# Correlate Ileo vs colonic resection

ResectionType <- subset(LLD_IBD_MetTax.asin, CurrentStomaOrPouchType == 'Resections' & ResectionIlealAny == 'yes' | CurrentStomaOrPouchType == 'Resections' & ResectionColonicAny == 'yes')

ResectionType <- ResectionType[-c(53),]

a=1
for (i in 1:nrow(ResectionType)) {
  if(ResectionType$ResectionIlealAny[i] == "yes" & ResectionType$ResectionColonicAny[i] == "yes") {
    ResectionType$ResectionIlealColonicOrBoth[a] <- 'Both'
    a=a+1
  }
  else{
   if(ResectionType$ResectionIlealAny[i] == "yes"){
     ResectionType$ResectionIlealColonicOrBoth[a] <- 'Ileal'
     a=a+1
   }
    else{
      ResectionType$ResectionIlealColonicOrBoth[a] <- 'Colonic'
      a=a+1
    }
  }
}


ResectionType <- ResectionType[c(232,234,881,747:880)]

ResectionType$ResectionIlealColonicOrBoth <- as.factor(ResectionType$ResectionIlealColonicOrBoth)

#For loop to correlate Resections type with each bacterial species

my_matrix <-  matrix(nrow = ncol(ResectionType[c(4:ncol(ResectionType))]) * 3, ncol = 10)
colnames(my_matrix) <- c('Taxonomy', 'Group1', 'Group2', 'NGroup1', 'NGroup2', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'MeanGroup1', 'MeanGroup2')

combos <- combn(levels(ResectionType[,3]),2)
a=1
for (i in 4:ncol(ResectionType)) {
  wilcox <- (pairwise.wilcox.test(ResectionType[,i], ResectionType[,3], p.adjust.method = 'none'))$p.value
  for (z in 1:ncol(combos)) {
    my_matrix[a,1] <- colnames(ResectionType)[i]
    my_matrix[a,2] <- combos[1,z]
    my_matrix[a,3] <- combos[2,z]
    my_matrix[a,4] <- sum(ResectionType[,3]==combos[1,z], na.rm = T) #Number of samples in group1
    my_matrix[a,5] <- sum(ResectionType[,3]==combos[2,z], na.rm = T)
    my_matrix[a,6] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test coefficient values 
    my_matrix[a,9] <- mean(subset(ResectionType, ResectionType[c(3)] == combos[1,z], select = i)[,1]) #Mean taxonomy relative abundances group 3
    my_matrix[a,10] <- mean(subset(ResectionType, ResectionType[c(3)] == combos[2,z], select = i)[,1]) #Mean taxonomy relative abundances group 2
    a=a+1
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')


ResectionTypeVsBacterialSpecies <- my_matrix
write.table(ResectionTypeVsBacterialSpecies, file = 'IlealvsColonicResectionsSpeciesAbun.asin.txt', sep = "\t",col.names = TRUE)





#Including IleocecalValve Resection

ResectionType <- subset(LLD_IBD_MetTax.asin, CurrentStomaOrPouchType == 'Resections' & ResectionIlealAny == 'yes' | CurrentStomaOrPouchType == 'Resections' & ResectionColonicAny == 'yes' | CurrentStomaOrPouchType == 'Resections' & ResectionIlealCecalAny == 'yes')

ResectionType <- ResectionType[-c(33,93,111),]

a=1
for (i in 1:nrow(ResectionType)) {
  if(ResectionType$ResectionIlealAny[i] == "yes" & ResectionType$ResectionColonicAny[i] == "yes" & ResectionType$ResectionIlealCecalAny[i] == "yes") {
    ResectionType$ResectionType1[a] <- 'All'
    a=a+1
  }
  else{
    if(ResectionType$ResectionIlealAny[i] == "yes" & ResectionType$ResectionColonicAny[i] == "yes"){
      ResectionType$ResectionType1[a] <- 'IlealColonic'
      a=a+1
    }
    else{
      if(ResectionType$ResectionIlealAny[i] == "yes" & ResectionType$ResectionIlealCecalAny[i] == "yes"){
        ResectionType$ResectionType1[a] <- 'IlealIlealcecal'
        a=a+1
      }
      else{
        if(ResectionType$ResectionColonicAny[i] == "yes" & ResectionType$ResectionIlealCecalAny[i] == "yes"){
          ResectionType$ResectionType1[a] <- 'ColonicIlealCecal'
          a=a+1
        }
      else{
        if(ResectionType$ResectionIlealAny[i] == "yes"){
          ResectionType$ResectionType1[a] <- 'Ileal'
          a=a+1
        }
        else{
          if(ResectionType$ResectionColonicAny[i] == "yes"){
            ResectionType$ResectionType1[a] <- 'Colonic'
            a=a+1
          }
          else {
            ResectionType$ResectionType1[a] <- 'IlealCecal'
            a=a+1
          }
        }
      }
     }
    }
  }
}



ResectionType <- ResectionType[c(229,232,234,881,747:880)]

ResectionType$ResectionType1 <- as.factor(ResectionType$ResectionType1)



my_matrix <-  matrix(nrow = ncol(ResectionType[c(4:ncol(ResectionType))]) * 21, ncol = 10)
colnames(my_matrix) <- c('Taxonomy', 'Group1', 'Group2', 'NGroup1', 'NGroup2', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'MeanGroup1', 'MeanGroup2')

combos <- combn(levels(ResectionType[,4]),2)
a=1
for (i in 5:ncol(ResectionType)) {
  wilcox <- (pairwise.wilcox.test(ResectionType[,i], ResectionType[,4], p.adjust.method = 'none'))$p.value
  for (z in 1:ncol(combos)) {
    my_matrix[a,1] <- colnames(ResectionType)[i]
    my_matrix[a,2] <- combos[1,z]
    my_matrix[a,3] <- combos[2,z]
    my_matrix[a,4] <- sum(ResectionType[,4]==combos[1,z], na.rm = T) #Number of samples in group1
    my_matrix[a,5] <- sum(ResectionType[,4]==combos[2,z], na.rm = T)
    my_matrix[a,6] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test coefficient values 
    my_matrix[a,9] <- mean(subset(ResectionType, ResectionType[c(4)] == combos[1,z], select = i)[,1]) #Mean taxonomy relative abundances group 3
    my_matrix[a,10] <- mean(subset(ResectionType, ResectionType[c(4)] == combos[2,z], select = i)[,1]) #Mean taxonomy relative abundances group 2
    a=a+1
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')


ResectionTypeVsBacterialSpecies1 <- my_matrix
write.table(ResectionTypeVsBacterialSpecies1, file = 'IlealColonic&IlealCecalResectionsSpeciesAbunCorrelation.AsinTranform1.txt', sep = "\t",col.names = TRUE)

```


### 6. Univariate Analysis

```{r}

# Packages & Functions Required:
library(ggplot2)
library(forcats)
library(knitr)
library(dplyr)
library(plotly)
library(rmarkdown)
library(UpSetR)
library(stringr)

Remove_FactCols #see 'Functions' file, point 2.
transform_and_filter_taxa #see 'Functions' file, point 6.
univariate.analysis2 #see 'Functions' file, point 7.
remove.rows #see 'Functions' file, point 8.
Prevalence.filter #see 'Functions' file, point 9.


#Data import & process:
#Import updated Meta/Species... See 'pathways' R script for origin 
LLD_IBD_MetTaxS <- read.table(file = "../RScripts/LLD_IBD_MetaSpecies150420.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) 

LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS[-c(212,218,379,1554,1671),]
LLD_IBD_MetTaxS <- subset(LLD_IBD_MetTaxS, DiagnosisCurrent != "MicroscopicColitis" & DiagnosisCurrent != "IBDI")
LLD_IBD_MetTaxS$DiagnosisCurrent <- fct_drop(LLD_IBD_MetTaxS$DiagnosisCurrent)

LLD_IBD_MetTaxS$CurrentStomaOrPouchType <- as.factor(gsub("Resection NoStoma", "Resections", LLD_IBD_MetTaxS$CurrentStomaOrPouchType))
LLD_IBD_MetTaxS$CurrentStomaOrPouchType <- as.factor(gsub("colostoma", "Resections", LLD_IBD_MetTaxS$CurrentStomaOrPouchType))

#Final subsetting .. 
#LLD.IBD_MetTaxS.F <- LLD_IBD_MetTaxS[c(2,4:6,8,9,14:17,20,22,24:28,30:37,42:48,50:69,71:86,88:181,183:197,414:422,424:431,433:684,687:1338)]
LLD.IBD_MetTaxS.F <- LLD_IBD_MetTaxS[c(2,4:6,8,9,14:17,20,22,24:28,30:37,42:48,50:69,71:86,88:181,183:197,202,203,408:416,418:425,427:678,681:1332)]

x <- gsub('CU','UC', LLD.IBD_MetTaxS.F$DiagnosisFirst) 
LLD.IBD_MetTaxS.F$DiagnosisFirst <- as.factor(x)

LLD.IBD_MetTaxS.F$DiagnosisCurrent <- as.factor(gsub("ReconsideringDiagnosis", "IBDU", LLD.IBD_MetTaxS.F$DiagnosisCurrent))

LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation <- as.numeric(as.character(LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation))

y <- lapply(LLD.IBD_MetTaxS.F[,1:448], class) #check class of each phenotype variable

#check how many levels per group, remove variables with too many levels i.e. too few samples per group for statistical testing
w <- lapply(LLD.IBD_MetTaxS.F[,1:448], nlevels) 
LLD.IBD_MetTaxS.F$ReasonUnclearSurgeryResectionProcedures <- NULL
#LLD.IBD_MetTaxS.F$NumberIndicationabcessOrfistula <- NULL #this group contains an NA level 

#write.table(LLD.IBD_MetTaxS.F, file = "LLD.IBD_MetTaxS.F.txt", sep = "\t", col.names = TRUE)

#**Remove factor columns function : removes factor variables with 'nlevels' (default = 1)
Remove_FactCols <- function(x, nlevels = 1){
  Remove_cols <- vector()
  for (i in 1:ncol(x)) {
    Name <- colnames(x)[i]
    if(nlevels(x[,i]) == 1){
      Remove_cols <- c(Remove_cols, Name)
    }
  }
  x_new <- x %>% select(-Remove_cols)
  print(paste(c("Function removed a total of ",length(Remove_cols), "variables"), collapse= ""))
  return(x_new)
}

LLD.IBD_Final <- Remove_FactCols(LLD.IBD_MetTaxS.F) 
y <- lapply(LLD.IBD_Final[,1:428], class)

LLD.IBD_Final <- LLD.IBD_Final[c(1:114,142:1080)] #remove logical variables, not informative and cannot be applied to the 'univariate.analysis' function
y <- lapply(LLD.IBD_Final[,1:401], class)      


#library(forcats)
LLD.IBD_Final[408,128] <- NA #change 'NA' to NA 
LLD.IBD_Final$NumberIndicationabcessOrfistula <- fct_drop(LLD.IBD_Final$NumberIndicationabcessOrfistula) #remove level 'NA'

LLD.IBD_Final$IncludedSamples <- NULL

#include BinaryBSS and bowel movement variable to LLD.IBD_Final dataframe
BinaryBSS <- read.table(file = "../Data/SI_project/intestinal_info.txt", header = TRUE, sep = "\t", quote = "\\")
LLD.IBD_Final <- merge(BinaryBSS, LLD.IBD_Final, by.x = 'UMCGIBDResearchIDorLLDeepID', by.y = 'row.names')
rownames(LLD.IBD_Final) <-LLD.IBD_Final$UMCGIBDResearchIDorLLDeepID
LLD.IBD_Final$UMCGIBDResearchIDorLLDeepID <- NULL
LLD.IBD_Final <- LLD.IBD_Final[c(3:388,1,2,389:1054)] #change order

#include updated PFReads variable to LLD.IBD_Final dataframe 
PFReads <- read.table(file = "../Data/SI_project/PFReads31.03.20.txt", header = TRUE, sep = "\t", quote = "\\")
LLD.IBD_Final$PFReads <- NULL
LLD.IBD_Final <- merge(PFReads, LLD.IBD_Final, by.x = 'SampleID', by.y = 'row.names')
rownames(LLD.IBD_Final) <-LLD.IBD_Final$SampleID
LLD.IBD_Final$SampleID <- NULL
LLD.IBD_Final <- LLD.IBD_Final[c(2,1,3:1054)] #change order
colnames(LLD.IBD_Final)[2] <- "PFReads"

#remove individuals with <10 million PFReads
#LLD.IBD_Final1 <- subset(LLD.IBD_Final, PFReads >= 10000000) #pouch = 8 (instead of 9), ileostoma = 45 (instead of 49)

#Final data with ileostoma and pouch individuals as one group
LLD.IBD_Finalx <-  LLD.IBD_Final
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))

LLD.IBD_Finalx$NumberIndicationabcessOrfistula <- as.numeric(as.character(LLD.IBD_Finalx$NumberIndicationabcessOrfistula))

#Remove 'correlated/related' variables
LLD.IBD_Finalxx <- LLD.IBD_Finalx[c(1,2,5:13,15,18,21:31,33,36:43,45,47,49,52:108,110,111,114:134,302,310:317,319,321,323:370,373,377,379,381:385,387,389,403:1054)]
LLD.IBD_Finalxx <- LLD.IBD_Finalxx[-c(74,76:78,80:93,96,98,113)]

# Remove bacteria from dataframe with a prevalence < 15%
LLD.IBD_MetaSp15 <- Prevalence.filter(LLD.IBD_Finalxx, RowSpecies = 165)

#Transform species abundance values
LLD.IBD_MetaSp15.asin <-transform_and_filter_taxa(LLD.IBD_MetaSp15[,c(165:298)])
LLD.IBD_MetaSp15.asin <-cbind(LLD.IBD_MetaSp15[c(1:164)],LLD.IBD_MetaSp15.asin)


#Analyses
# LLD.IBD_Finalx... ileostoma and pouch individuals as one group

LLD.IBD.UniA.asin <- univariate.analysis2(LLD.IBD_MetaSp15.asin,164,165,134) 

LLD.IBD.UniA.asin$PrevalenceGroup1 <- as.numeric(as.character(LLD.IBD.UniA.asin$PrevalenceGroup1)) 
LLD.IBD.UniA.asin$PrevalenceGroup2 <- as.numeric(as.character(LLD.IBD.UniA.asin$PrevalenceGroup2))
LLD.IBD.UniA.asin$`NA` <- as.numeric(as.character(LLD.IBD.UniA.asin$`NA`))
LLD.IBD.UniA.asin$`NSamplesG1/Nnonzero` <- as.numeric(as.character(LLD.IBD.UniA.asin$`NSamplesG1/Nnonzero`))
LLD.IBD.UniA.asin$`NSampleG2/Nzero` <- as.numeric(as.character(LLD.IBD.UniA.asin$`NSampleG2/Nzero`))

LLD.IBD.UniA.asin1 <- subset(LLD.IBD.UniA.asin, `NA` < 1720*0.95) #remove rows with NA > 95%

# Filter out tests with a group sample number < 20
LLD.IBD.UniA.asin2 <- remove.rows(LLD.IBD.UniA.asin1)

library(plotly)

#____ Top 25 phenotypes with most bacterial associations ___

Significant <- subset(LLD.IBD.UniA.asin2, BelowAdjustPval == 'yes') #subset associations with p-values below adjusted p-value

Significant <- transform(Significant, Combinations=paste(Phenotype, Group1, Group2, sep=":")) #merge Phenotype, Group1, Group2 variables into a new variable

paged_table(Significant, options = list(rows.print = 15))

#plot top 25 most frequent variables 
sorted_TopPhe.asin <- Significant %>% #create dataframe with 'combinations' variable in descending order of frequency
  count(Combinations) %>%
  arrange(-n) %>% 
  mutate(Combinations = factor(Combinations, levels = unique(Combinations)))

paged_table(sorted_TopPhe.asin, options = list(rows.print = 15))

a <- ggplot(sorted_TopPhe.asin[1:25,], aes(x = Combinations, y = n)) +
  geom_bar(stat="identity", fill = "grey") + xlab('Phenotype') + 
  ylab('Significant Associations (n)') + theme_classic() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size=5)) + scale_x_discrete(labels=c('CurrentStomaOrPouchType:HealthyControl:IleostomaPouch' = 'Bowel Phenotype\n (LLD General Population vs IBD Small Intestine)','CurrentStomaOrPouchType:IleostomaPouch:No Resections' = 'Bowel Phenotype\n (IBD Small Intestine vs IBD Intact Colon)','EverHadStomaOrPouch:no:yes' = 'Ever Had A Stoma Or Pouch\n (No/Yes)','NumberOfResectionColonic:NumberOfResectionColonic:NumberOfResectionColonic' = 'Number Of Colonic Resections','ResectionColonicAny:no:yes' = 'Colonic Resection Any\n (No/Yes)','DiagnosisCurrent:CD:generalpopulation' = 'Current IBD Diagnosis\n (CD vs General Population)','IleocecalValveInSitu:no:yes' = 'Ileocecal Valve In Situ\n (No/Yes)','CurrentStomaOrPouchType:HealthyControl:Resections' = 'Bowel Phenotype\n (LLD General Population vs IBD Partially Intact Colon)','ResectionAny:no:yes' = 'Resection Any\n (No/Yes)','CurrentStomaOrPouchType:IleostomaPouch:Resections' ='Bowel Phenotype\n (IBD Small Intestine vs IBD Partially Intact Colon)', 'NumberOfResectionsAny:NumberOfResectionsAny:NumberOfResectionsAny' = 'Total Number Of Resections','DiagnosisCurrent:generalpopulation:UC' = 'Current IBD Diagnosis\n (General Population vs UC)','PFReads:PFReads:PFReads' = 'Reads Depth', 'vitamin_B12:no:yes' = 'Vitamin B12\n (No/Yes)', 'Alc_en:Alc_en:Alc_en' = 'Alcohol\n (% Total Energy Intake)','CurrentStomaOrPouchType:HealthyControl:No Resections' = 'Bowel Phenotype\n (LLD General Population vs IBD Intact Colon)','vitamin_D:no:yes' = 'Vitamin D\n (No/Yes)','calcium:no:yes' = 'Calcium\n (No/Yes)', 'CurrentStomaOrPouchType:No Resections:Resections' = 'Bowel Phenotype\n (IBD Intact Colon vs IBD Partially Intact Colon)','FecalCalprotectinOver200yesno:no:yes' = 'Fecal Calprotectin Over 200\n (No/Yes)','oral_steroid:no:yes'= 'Oral Steroids\n (No/Yes)','PPI:no:yes' = 'PPI\n (No/Yes)','ResectionIlealCecalAny:no:yes' = 'IlealCecal Resection Any\n (No/Yes)','NumberOfResetionsIleoCecal:NumberOfResetionsIleoCecal:NumberOfResetionsIleoCecal' = 'Number Of IleoCecal Resections','BowelMovementADayDef:BowelMovementADayDef:BowelMovementADayDef' = 'Bowel Movement A Day Def'))

ggplotly(a)

## Top 25 species according to number of significant phenotype associations  
sorted_TopTax.asin <- Significant %>% #create dataframe with 'Taxonomy' variable in descending order of frequency
  count(Taxonomy) %>%
  arrange(-n) %>% 
  mutate(Taxonomy = factor(Taxonomy, levels = unique(Taxonomy)))

paged_table(sorted_TopTax.asin, options = list(rows.print = 15))

sorted_TopTax.asin$Taxonomy <- as.character(gsub("_", " ", sorted_TopTax.asin$Taxonomy))

b <- ggplot(sorted_TopTax.asin[1:25,], aes(x = reorder(Taxonomy, -n), y = n)) +
  geom_bar(stat="identity", fill = "darkgrey") +  theme_classic() + geom_text(aes(label=n), vjust=-0.3, size=2) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size=5)) + ylab('Significant Associations (n)') + xlab('Species') #plot top 25 most frequent variables

ggplotly(b)

```



