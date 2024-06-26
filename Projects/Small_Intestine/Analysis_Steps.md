Small Intestine Microbiota Project
==================================

Creator: Renate Ruigrok
Year: 2020

### Key

```{r}

# HealthyControl = General population (LifelinesDEEP samples)
# No Resections = IBD-NoRes (samples from 1000IBD patients with no intestinal resections)
# Resections = IBD-Res (samples from 1000IBD patients with segmental intestinal resections)
# Ileostoma + pouch (IleostomaPouch) = IBD-SI (samples from 1000IBD patients with an ileostomy or ileoanal pouch)

```

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

#renaming diagnosis from dutch to english abbreviation
x <- gsub('CU','UC', LLD.IBD_MetTaxS.F$DiagnosisFirst) 
LLD.IBD_MetTaxS.F$DiagnosisFirst <- as.factor(x)

LLD.IBD_MetTaxS.F$DiagnosisCurrent <- as.factor(gsub("ReconsideringDiagnosis", "IBDU", LLD.IBD_MetTaxS.F$DiagnosisCurrent))

LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation <- as.numeric(as.character(LLD.IBD_MetTaxS.F$TimeEndPreviousExacerbation))

LLD.IBD_MetTaxS.F$ReasonUnclearSurgeryResectionProcedures <- NULL

LLD.IBD_Final <- Remove_FactCols(LLD.IBD_MetTaxS.F) #remove factor variables with only 1 level

LLD.IBD_Final <- LLD.IBD_Final[c(1:114,142:1080)] #remove logical variables, not informative in this case 

#library(forcats)
LLD.IBD_Final[408,128] <- NA #change 'NA' to NA 
LLD.IBD_Final$NumberIndicationabcessOrfistula <- fct_drop(LLD.IBD_Final$NumberIndicationabcessOrfistula) #remove the level 'NA'

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

===================================================================================================================================
Revised Manuscript - 17.12.20
===================================================================================================================================

# Read depth range per group
aggregate(x = LLD.IBD_Finalx$PFReads,  # Specify data column
          by = list(LLD.IBD_Finalx$CurrentStomaOrPouchType), # Specify group indicator
          FUN = range, na.rm = T)  

# Histogram - sequencing read depth per group

LLD.IBD_Finalx <-  LLD.IBD_Final
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("ileostoma", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))
LLD.IBD_Finalx$CurrentStomaOrPouchType <- as.factor(gsub("pouch", "IleostomaPouch", LLD.IBD_Finalx$CurrentStomaOrPouchType))

LLD.IBD_Finalx$CurrentStomaOrPouchType <- factor(LLD.IBD_Finalx$CurrentStomaOrPouchType, levels=c("HealthyControl", "No Resections", "Resections", "IleostomaPouch"))
ReadDepthDiv1m <- LLD.IBD_Finalx[,c("CurrentStomaOrPouchType","PFReads")]
ReadDepthDiv1m$PFReads <- ReadDepthDiv1m$PFReads/1000000
ReadDepthDiv1m$ID <- rownames(ReadDepthDiv1m)
ReadDepthDiv1m <- ReadDepthDiv1m[order(ReadDepthDiv1m$CurrentStomaOrPouchType),]

names <- c(
  "HealthyControl" = 'General Population',
  'No Resections' = 'IBD non-resected intestine',
  'Resections' = 'IBD resected intestine',
  "IleostomaPouch" = 'IBD small intestine'
)

# plot histogram
ggplot(ReadDepthDiv1m) + geom_histogram(aes(x=PFReads,fill=CurrentStomaOrPouchType),color="black") +
  facet_grid(scales = "free", CurrentStomaOrPouchType~., labeller = as_labeller(names)) + 
  theme_minimal() + scale_fill_manual(values=c('grey54','slateblue2', 'gold1', 'red2'))+
  labs(x=bquote('Sequencing Read Depth' ~(x10^6)), color = "Bowel Phenotype") + theme (legend.position = "none")

```

### 2. Taxonomic Composition  

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
#IBD_Taxa <- as.data.frame(t(IBD_Taxa)) #transform

#LLD Taxonomy
LLD_Taxa <- read.table(file = "../Data/SI_project/LLD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #LLD Taxonomy 
rownames(LLD_Taxa) <- LLD_Taxa$ID
LLD_Taxa$ID <- NULL

LLD <- gsub("^X", "", colnames(LLD_Taxa)) #Remove 'X' at the beginning of ID numbers
colnames(LLD_Taxa) <- LLD

#Merge the two cohort taxa's
#IBD_Taxa1 <- as.data.frame(t(LLD_Taxa))
LLD_IBD_Taxa <- merge(LLD_Taxa, IBD_Taxa, by = 'row.names', all = T)
rownames(LLD_IBD_Taxa) <- LLD_IBD_Taxa$Row.names
LLD_IBD_Taxa$Row.names <- NULL
LLD_IBD_Taxa <- as.data.frame(t(LLD_IBD_Taxa))


#Merge metadata with taxa
LLD_IBD_MetTax <- merge(LLD.IBD_Meta,LLD_IBD_Taxa, by.x = 'CombinedIDs', by.y = 'row.names')

*************
GENUS Profile
*************

# Top Genera within the Healthy control group
  
#Extract only genera from Taxa table
LLD_IBD_TaxaG <- LLD_IBD_Taxa[,grepl('g__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaG <- LLD_IBD_TaxaG[,!grepl('s__',colnames(LLD_IBD_TaxaG),ignore.case = T)]
  
  
#Merge genus Taxa with LLD_IBD_MetaTaxa
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
  
  LLD_IBD_MetTaxG$g__Other <- 100 - rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))])
  
  write.table(LLD_IBD_MetTaxG, file = 'LLD_IBD_GenusComposition_top10HC.txt', sep = '\t', col.names = T)
 
  #Plot graph by stoma type
  #column mean by group
  LLD_IBD_MetTaxG <- aggregate(LLD_IBD_MetTaxG[,2:12], list(LLD_IBD_MetTaxG$CurrentStomaOrPouchType), mean)
  colnames(LLD_IBD_MetTaxG)[1] <- "CurrentStomaOrPouchType"
  
  #Convert each mean relative abundance to a percentage of total means 
  LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))] <- LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]/
    rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))])*100
  
  rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage
  
  LLD <- gsub('.*\\|g__','',colnames(LLD_IBD_MetTaxG)) #Keep taxonomy level name only  
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
  
  
  
#Top Genera within the IleostomaPouch group
  
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
  
  write.table(LLD_IBD_MetTaxG, file = 'LLD_IBD_GenusComposition_top10SI.txt', sep = '\t', col.names = T)

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
  
======================================================================================================================================================         
REVISED Manuscript - 17.12.20
======================================================================================================================================================         

*******************************
GENUS Profile cont'd - Boxplots
*******************************

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

#Select top 12 abundant genera in General population
#subset 'HealthyControl' individuals only
  HC <- subset(LLD_IBD_MetTaxG, CurrentStomaOrPouchType == 'HealthyControl')
  
  #Vector with mean relative abundances in descending order
  Means <- colMeans(HC[,c(2:ncol(HC))]) #mean
  Means <- sort(Means, decreasing = T) #order
  
LLD_IBD_MetTaxG1 <- LLD_IBD_MetTaxG[c("CurrentStomaOrPouchType",TopNTaxa(Means,12))]

#Select top 12 abundant genera in SI
SI <- subset(LLD_IBD_MetTaxG, CurrentStomaOrPouchType == 'IleostomaPouch')
  
  #Vector with mean relative abundances in descending order
  Means <- colMeans(SI[,c(2:ncol(SI))]) #mean
  Means <- sort(Means, decreasing = T) #order

Dataframe top 12 genera in SI and General population group
LLD_IBD_MetTaxG2 <- LLD_IBD_MetTaxG[c(colnames(LLD_IBD_MetTaxG1),TopNTaxa(Means,12))]
  
#Transform genera abundances
MetTaxG_transformed <- transform_and_filter_taxa(LLD_IBD_MetTaxG2[,c(2:ncol(LLD_IBD_MetTaxG2))])
MetTaxG_transformed <-cbind(LLD_IBD_MetTaxG[1],MetTaxG_transformed)

colnames(MetTaxG_transformed) <- gsub('.*g__','',colnames(MetTaxG_transformed)) 

#Plot graph
MetTaxG_transformed1<- melt(MetTaxG_transformed, 'CurrentStomaOrPouchType')

#Facet boxplots - Top 12 General population
ggplot(data = MetTaxG_transformed1[1:20556,], aes(x=variable, y=value,fill=CurrentStomaOrPouchType),color = "black") + scale_fill_manual(values=c('grey54','slateblue2', 'gold1', 'red2')) +  geom_point(position=position_jitterdodge(),alpha=0.2) + theme_minimal()+ geom_boxplot(outlier.shape = NA) + labs(x="", y="Relative Abundance (log10)", color = "Bowel Phenotype") +
     coord_cartesian(ylim = c(0, 1.25)) + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,1,1.25)) + facet_wrap( ~ variable, scales = "free_x") 

#Facet boxplots - Top 12 SI
ggplot(data = MetTaxG_transformed1[20557:41112,], aes(x=variable, y=value,fill=CurrentStomaOrPouchType),color = "black") + scale_fill_manual(values=c('grey54','slateblue2', 'gold1', 'red2')) +  geom_point(position=position_jitterdodge(),alpha=0.2) + theme_minimal()+ geom_boxplot(outlier.shape = NA) + labs(x="", y="Relative Abundance (log10)", color = "Bowel Phenotype")+ facet_wrap( ~ variable, scales = "free_x") 


**************
PHYLUM Profile
**************

#Extract only phyla from LLD_IBD_Taxa table
LLD_IBD_TaxaP <- LLD_IBD_Taxa[,grepl('p__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaP<- LLD_IBD_TaxaP[,!grepl('c__',colnames(LLD_IBD_TaxaP),ignore.case = T)]


#Merge Phyla Taxa with LLD_MetaTaxa IDs & CurrentStomaOrPouchType variable
LLD_IBD_MetTaxP <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaP, by.x = 'CombinedIDs', by.y = 'row.names')
rownames(LLD_IBD_MetTaxP) <- LLD_IBD_MetTaxP$CombinedIDs
LLD_IBD_MetTaxP$CombinedIDs <- NULL

# Remove phyla from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 2:ncol(LLD_IBD_MetTaxP)) {
  cname <- colnames(LLD_IBD_MetTaxP)[i]
  if(colSums(LLD_IBD_MetTaxP[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxP) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTaxP1 <- LLD_IBD_MetTaxP %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

LLD_IBD_MetTaxP1 <- LLD_IBD_MetTaxP1[,-c(2,8)]

#Vector with mean relative abundances in descending order
#PMeans <- colMeans(LLD_IBD_MetTaxP1[,c(2:ncol(LLD_IBD_MetTaxP1))]) #mean
#PMeans <- sort(PMeans, decreasing = T) #order

#LLD_IBD_MetTaxP1 <- LLD_IBD_MetTaxP1[c("CurrentStomaOrPouchType",TopNTaxa(PMeans,7))]

LLD_IBD_MetTaxP1$p__Other <- 100 - rowSums(LLD_IBD_MetTaxP1[,grep('__',colnames(LLD_IBD_MetTaxP1))]) 

write.table(LLD_IBD_MetTaxP1, file = 'LLD_IBD_PhylumComposition.txt', sep = '\t', col.names = T)
  
#Plot graph by stoma type

#column mean by group
LLD_IBD_MetTaxP1 <- aggregate(LLD_IBD_MetTaxP1[,2:7], list(LLD_IBD_MetTaxP1$CurrentStomaOrPouchType), mean)
colnames(LLD_IBD_MetTaxP1)[1] <- "CurrentStomaOrPouchType"

#Convert each mean relative abundance to a percentage of total means 
LLD_IBD_MetTaxP1[,grep('__',colnames(LLD_IBD_MetTaxP1))] <- LLD_IBD_MetTaxP1[,grep('__',colnames(LLD_IBD_MetTaxP1))]/
  rowSums(LLD_IBD_MetTaxP1[,grep('__',colnames(LLD_IBD_MetTaxP1))])*100

rowSums(LLD_IBD_MetTaxP1[,grep('__',colnames(LLD_IBD_MetTaxP1))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage

colnames(LLD_IBD_MetTaxP1) <- gsub('.*p__','',colnames(LLD_IBD_MetTaxP1))

#Plot graph
LLD_IBD_MetTaxP1<- melt(LLD_IBD_MetTaxP1, 'CurrentStomaOrPouchType')

LLD_IBD_MetTaxP1 %>%
  group_by(CurrentStomaOrPouchType) %>%
  summarise(Sum = sum(value)) #Check each group sums to 100

ggplot(data = LLD_IBD_MetTaxP1, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
  geom_bar(stat = "identity",color = "black") + scale_fill_brewer(palette = 'Set2') +
  labs(y = "Relative Abundance", x = 'Bowel Phenotype', title = 'Microbial Phylum Composition') + scale_x_discrete(labels=c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD Small Intestine')) + theme_minimal() + theme(axis.text.x = element_text(angle = 45)) 


**************
FAMILY Profile
**************

# Top Families within the Healthy control group

#Extract Family only from Taxa table
LLD_IBD_TaxaF <- LLD_IBD_Taxa[,grepl('f__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaF <- LLD_IBD_TaxaF[,!grepl('g__',colnames(LLD_IBD_TaxaF),ignore.case = T)]


#Merge Family Taxo with LLD_IBD_MetaTaxa
LLD_IBD_MetTaxF <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaF, by.x = 'CombinedIDs', by.y = 'row.names')
rownames(LLD_IBD_MetTaxF) <- LLD_IBD_MetTaxF$CombinedIDs
LLD_IBD_MetTaxF$CombinedIDs <- NULL

remove_cols <- vector()
for (i in 2:ncol(LLD_IBD_MetTaxF)) {
  cname <- colnames(LLD_IBD_MetTaxF)[i]
  if(colSums(LLD_IBD_MetTaxF[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxF) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTaxF <- LLD_IBD_MetTaxF %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

#subset 'HealthyControl' individuals only
HC <- subset(LLD_IBD_MetTaxF, CurrentStomaOrPouchType == 'HealthyControl')

#Vector with mean relative abundances in descending order
Means <- colMeans(HC[,c(2:ncol(HC))]) #mean
Means <- sort(Means, decreasing = T) #order

LLD_IBD_MetTaxF <- LLD_IBD_MetTaxF[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]

LLD_IBD_MetTaxF$f__Other <- 100 - rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))]) 

write.table(LLD_IBD_MetTaxF, file = 'LLD_IBD_FamilyComposition_top10HC.txt', sep = '\t', col.names = T)

#Plot graph by stoma type
#column mean by group
LLD_IBD_MetTaxF <- aggregate(LLD_IBD_MetTaxF[,2:12], list(LLD_IBD_MetTaxF$CurrentStomaOrPouchType), mean)
colnames(LLD_IBD_MetTaxF)[1] <- "CurrentStomaOrPouchType"

#Convert each mean relative abundance to a percentage of total means 
LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))] <- LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))]/
  rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))])*100

rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage

LLD <- gsub('.*f__','',colnames(LLD_IBD_MetTaxF)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxF) <- LLD

#Graph
LLD_IBD_MetTaxF<- melt(LLD_IBD_MetTaxF, 'CurrentStomaOrPouchType')

LLD_IBD_MetTaxF %>%
  group_by(CurrentStomaOrPouchType) %>%
  summarise(Sum = sum(value)) #Check each group sums to 100

cols <- c("Eubacteriaceae" = "#66C2A5", "Bifidobacteriaceae" = "#FC8D62", "Ruminococcaceae" = "#8DA0CB", "Lachnospiraceae" = "#FFFF33", "Methanobacteriaceae" = "#E78AC3", "Verrucomicrobiaceae" = "#BC8F8F", "Other" = "#D9D9D9", "Veillonellaceae" = "#A6CEE3", "Bacteroidaceae" = "#E5C494", "Streptococcaceae" = "#E41A1C", "Coriobacteriaceae" = "#FF6347")

ggplot(data = LLD_IBD_MetTaxF, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
  geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +
  labs(y = "Relative Abundance", x = 'Bowel Phenotype', title = 'Most Abundant Families: GP') + scale_x_discrete(labels=c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD Small Intestine')) + theme_minimal() + theme(axis.text.x = element_text(angle = 45)) 



# Top Families within the IleostomaPouch group

#Merge Family Taxo with LLD_IBD_MetaTaxa
LLD_IBD_MetTaxF <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaF, by.x = 'CombinedIDs', by.y = 'row.names')
rownames(LLD_IBD_MetTaxF) <- LLD_IBD_MetTaxF$CombinedIDs
LLD_IBD_MetTaxF$CombinedIDs <- NULL

remove_cols <- vector()
for (i in 2:ncol(LLD_IBD_MetTaxF)) {
  cname <- colnames(LLD_IBD_MetTaxF)[i]
  if(colSums(LLD_IBD_MetTaxF[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxF) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTaxF <- LLD_IBD_MetTaxF %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

#subset 'SI' individuals only
SI <- subset(LLD_IBD_MetTaxF, CurrentStomaOrPouchType == 'IleostomaPouch')

#Vector with mean relative abundances in descending order
Means <- colMeans(SI[,c(2:ncol(SI))]) #mean
Means <- sort(Means, decreasing = T) #order

LLD_IBD_MetTaxF <- LLD_IBD_MetTaxF[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]

LLD_IBD_MetTaxF$f__Other <- 100 - rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))])

write.table(LLD_IBD_MetTaxF, file = 'LLD_IBD_FamilyComposition_top10SI.txt', sep = '\t', col.names = T)

#Plot graph by stoma type
#column mean by group
LLD_IBD_MetTaxF <- aggregate(LLD_IBD_MetTaxF[,2:12], list(LLD_IBD_MetTaxF$CurrentStomaOrPouchType), mean)
colnames(LLD_IBD_MetTaxF)[1] <- "CurrentStomaOrPouchType"

#Convert each mean relative abundance to a percentage of total means 
LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))] <- LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))]/
  rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))])*100

rowSums(LLD_IBD_MetTaxF[,grep('__',colnames(LLD_IBD_MetTaxF))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage

LLD <- gsub('.*f__','',colnames(LLD_IBD_MetTaxF)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxF) <- LLD

#Graph
LLD_IBD_MetTaxF<- melt(LLD_IBD_MetTaxF, 'CurrentStomaOrPouchType')

LLD_IBD_MetTaxF %>%
  group_by(CurrentStomaOrPouchType) %>%
  summarise(Sum = sum(value)) #Check each group sums to 100


cols <- c("Eubacteriaceae" = "#66C2A5", "Bifidobacteriaceae" = "#FC8D62", "Ruminococcaceae" = "#8DA0CB", "Lachnospiraceae" = "#FFFF33", "Methanobacteriaceae" = "#E78AC3", "Enterococcaceae" = "#A6D854", "Other" = "#D9D9D9", "Coprococcus" = "#377EB8", "Bacteroidaceae" = "#E5C494", "Streptococcaceae" = "#E41A1C", "Enterobacteriaceae" = "#C2A5CF", "Peptostreptococcaceae" = "#BF5B17", "Lactobacillaceae" = "#CCEBC5", "Veillonellaceae" = "#A6CEE3", "Clostridiaceae" = "#F0027F")


ggplot(data = LLD_IBD_MetTaxF, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
  geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +
  labs(y = "Relative Abundance", x = 'Bowel Phenotype', title = 'Most Abundant Families: SI') + scale_x_discrete(labels=c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD Small Intestine')) + theme_minimal() + theme(axis.text.x = element_text(angle = 45)) 


***************
SPECIES Profile
***************

#Extract Species only
LLD_IBD_TaxaS <- LLD_IBD_Taxa[,grepl('s__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaS<- LLD_IBD_TaxaS[,!grepl('t__',colnames(LLD_IBD_TaxaS),ignore.case = T)]

# Top Species within the General population group

#Merge Species Taxo with stoma/pouch type and ID's
LLD_IBD_MetTaxS <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaS, by.x = 'CombinedIDs', by.y = 'row.names')
rownames(LLD_IBD_MetTaxS) <- LLD_IBD_MetTaxS$CombinedIDs
LLD_IBD_MetTaxS$CombinedIDs <- NULL

remove_cols <- vector()
for (i in 2:ncol(LLD_IBD_MetTaxS)) {
  cname <- colnames(LLD_IBD_MetTaxS)[i]
  if(colSums(LLD_IBD_MetTaxS[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxS) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

#subset ileostomy individuals only
HC <- subset(LLD_IBD_MetTaxS, CurrentStomaOrPouchType == 'HealthyControl')

#Vector with mean relative abundances in descending order
Means <- colMeans(HC[,c(2:ncol(HC))]) 
Means <- sort(Means, decreasing = T) #order from most abundant species to least 

LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]
LLD_IBD_MetTaxS$s__Other <- 100 - rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]) 

write.table(LLD_IBD_MetTaxS, file = 'LLD_IBD_SpeciesComposition_top10HC.txt', sep = '\t', col.names = T)
  

#Plot graph by stoma type

#column mean by group
LLD_IBD_MetTaxS <- aggregate(LLD_IBD_MetTaxS[,2:12], list(LLD_IBD_MetTaxS$CurrentStomaOrPouchType), mean)
colnames(LLD_IBD_MetTaxS)[1] <- "CurrentStomaOrPouchType"

#Convert each mean relative abundance to a percentage of total means 
LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))] <- LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]/
  rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))])*100

rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage

IBD <- gsub('.*s__','',colnames(LLD_IBD_MetTaxS)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxS) <- IBD


#Plot graph
LLD_IBD_MetTaxS<- melt(LLD_IBD_MetTaxS, 'CurrentStomaOrPouchType')

LLD_IBD_MetTaxS %>%
  group_by(CurrentStomaOrPouchType) %>%
  summarise(Sum = sum(value)) #Check each group sums to 100

cols <- c("Eubacterium_rectale" = "#66C2A5", "Bifidobacterium_adolescentis" = "#FC8D62", 
          "Ruminococcus_sp_5_1_39BFAA" = "#8DA0CB", "Ruminococcus_bromii" = "#E5C494", 
          "Subdoligranulum_unclassified" = "#E78AC3", "Faecalibacterium_prausnitzii" = 
            "#A6D854", "Bifidobacterium_longum" = "#FFD92F", "Collinsella_aerofaciens" ="#9846AC", 
          "Dorea_longicatena" = "#EF3A3A", "Eubacterium_hallii" = "#66C2A5", "Other" = "#D9D9D9" )

r <- ggplot(data = LLD_IBD_MetTaxS, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +labs(y = "Mean Relative Abundance", x = 'Bowel Phenotype', title = 'Most Abundant Species in General Population') + scale_x_discrete(labels=c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD Small Intestine')) + theme_minimal() 



# Top Species within the IleostomaPouch group

#Merge Species Taxo with stoma/pouch type and ID's
LLD_IBD_MetTaxS <- merge(LLD_IBD_MetTax[c(1,170)], LLD_IBD_TaxaS, by.x = 'CombinedIDs', by.y = 'row.names')
rownames(LLD_IBD_MetTaxS) <- LLD_IBD_MetTaxS$CombinedIDs
LLD_IBD_MetTaxS$CombinedIDs <- NULL

remove_cols <- vector()
for (i in 2:ncol(LLD_IBD_MetTaxS)) {
  cname <- colnames(LLD_IBD_MetTaxS)[i]
  if(colSums(LLD_IBD_MetTaxS[i] != 0, na.rm = T) / nrow(LLD_IBD_MetTaxS) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

#subset IleostomaPouch individuals only
Ileo <- subset(LLD_IBD_MetTaxS, CurrentStomaOrPouchType == 'IleostomaPouch')

#Vector with mean relative abundances in descending order
Means <- colMeans(Ileo[,c(2:ncol(Ileo))]) 
Means <- sort(Means, decreasing = T) #order from most abundant species to least 

LLD_IBD_MetTaxS <- LLD_IBD_MetTaxS[c("CurrentStomaOrPouchType",TopNTaxa(Means,10))]
LLD_IBD_MetTaxS$s__Other <- 100 - rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]) 

write.table(LLD_IBD_MetTaxS, file = 'LLD_IBD_SpeciesComposition_top10SI.txt', sep = '\t', col.names = T)

#Plot graph by stoma type

#column mean by group
LLD_IBD_MetTaxS <- aggregate(LLD_IBD_MetTaxS[,2:12], list(LLD_IBD_MetTaxS$CurrentStomaOrPouchType), mean)
colnames(LLD_IBD_MetTaxS)[1] <- "CurrentStomaOrPouchType"

#Convert each mean relative abundance to a percentage of total means 
LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))] <- LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]/
  rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))])*100

rowSums(LLD_IBD_MetTaxS[,grep('__',colnames(LLD_IBD_MetTaxS))]) #to check the sum of each row adds to 100... i.e correctly calculated mean as a percentage

IBD <- gsub('.*s__','',colnames(LLD_IBD_MetTaxS)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxS) <- IBD

#Plot graph
LLD_IBD_MetTaxS<- melt(LLD_IBD_MetTaxS, 'CurrentStomaOrPouchType')

LLD_IBD_MetTaxS %>%
  group_by(CurrentStomaOrPouchType) %>%
  summarise(Sum = sum(value)) #Check each group sums to 100

cols <- c("Escherichia_coli" = "#C2A5CF","Streptococcus_salivarius" = "#E41A1C", "Ruminococcus_gnavus" = 
            "#CCEBC5","Peptostreptococcaceae_noname_unclassified" = "#BF5B17","Streptococcus_parasanguinis" = "#377EB8",
          "Streptococcus_thermophilus" = "#F0027F","Veillonella_unclassified" = "#A6CEE3","Eubacterium_rectale" = "#66C2A5",
          "Bifidobacterium_longum" = "#FFD92F", "Prevotella_copri" = "#FF6347", "Other" = "#D9D9D9")
          

ggplot(data = LLD_IBD_MetTaxS, aes(x = CurrentStomaOrPouchType, y = value, fill = variable)) +
  geom_bar(stat = "identity",color = "black") + scale_fill_manual(values = cols) +
  labs(y = "Mean Relative Abundance", x = 'Bowel Phenotype', title = 'Most Abundant Species in IBD Small Intestine') + 
  scale_x_discrete(labels=c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 
                            'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD Small Intestine')) + theme_minimal() 

```


### 3. Taxonomic Composition: Wilcoxon test

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

======================================================================================================================================================         
REVISED Manuscript - 17.12.20
======================================================================================================================================================         

******
GENERA
******

#Top genera within the General population group
  
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
  
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG[c("CurrentStomaOrPouchType",TopNTaxa(Means,12))]
  
  LLD_IBD_MetTaxG$Other <- 100 - rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) 
  
  LLD <- gsub('.*\\|g__','',colnames(LLD_IBD_MetTaxG)) #Keep taxonomy level name only i.e. 'class name' 
  colnames(LLD_IBD_MetTaxG) <- LLD
  
  Top12HCGenera <- LLD_IBD_MetTaxG
  Top12HCGenera <- Top12HCGenera[-c(14)]
  
  my_matrix<-  matrix(nrow = 6 * ncol(Top12HCGenera[c(2:ncol(Top12HCGenera))]), ncol = 13)
  colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')
  
  
  #__Add correlation analysis data to the matrix 
  a=1 
  #Iterate through taxonomy columns
  for (j in 2:ncol(Top12HCGenera)) {
    #iterate through phenotype columns
    for (k in 1:1) {
      #following conditions applied to all factor variables
      combos <- combn(levels(Top12HCGenera[,k]),2)
      #pairwise.wilcox.test 
      wilcox <- (pairwise.wilcox.test(Top12HCGenera[,j], Top12HCGenera[,k], p.adjust.method = 'none'))$p.value
      #condition applied to all factor variables with more than 2 groups
      for (z in 1:ncol(combos)) {
        #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
        my_matrix[a,1] <- colnames(Top12HCGenera)[j]
        #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
        my_matrix[a,2] <- colnames(Top12HCGenera)[k]
        #adding the group combinations to my_matrix
        my_matrix[a,3] <- combos[1,z]
        my_matrix[a,4] <- combos[2,z]
        my_matrix[a,5] <- nlevels(Top12HCGenera[,k])
        my_matrix[a,6] <- sum(Top12HCGenera[,k]==combos[1,z], na.rm = T) #Number of samples in group1
        my_matrix[a,7] <- sum(Top12HCGenera[,k]==combos[2,z], na.rm = T) #Number of samples in group2
        my_matrix[a,8] <- colSums(is.na(Top12HCGenera[1]))
        my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
        my_matrix[a,12] <- mean(subset(Top12HCGenera, Top12HCGenera[c(k)] == combos[1,z], select = j)[,k])
        my_matrix[a,13] <- mean(subset(Top12HCGenera, Top12HCGenera[c(k)] == combos[2,z], select = j)[,k])
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
  my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
  my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')
  
  
  write.table(my_matrix, file = "Top12_HC_Genera.txt", sep = "\t", col.names = T)
  
  

  
#Top genera within the IleostomaPouch group
  
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
  
  LLD_IBD_MetTaxG <- LLD_IBD_MetTaxG[c("CurrentStomaOrPouchType",TopNTaxa(Means,12))]
  
  LLD_IBD_MetTaxG$Other <- 100 - rowSums(LLD_IBD_MetTaxG[,grep('__',colnames(LLD_IBD_MetTaxG))]) 
  
  LLD <- gsub('.*\\|g__','',colnames(LLD_IBD_MetTaxG)) #Keep taxonomy level name only i.e. 'class name' 
  colnames(LLD_IBD_MetTaxG) <- LLD

  

  #Top SI Genera
  
  Top12SIGenera <- LLD_IBD_MetTaxG
  Top12SIGenera <- Top12SIGenera[-c(14)]
  
  my_matrix<-  matrix(nrow = 6 * ncol(Top12SIGenera[c(2:ncol(Top12SIGenera))]), ncol = 13)
  colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')
  
  
  #__Add correlation analysis data to the matrix 
  a=1 
  #Iterate through taxonomy columns
  for (j in 2:ncol(Top12SIGenera)) {
    #iterate through phenotype columns
    for (k in 1:1) {
      #following conditions applied to all factor variables
      combos <- combn(levels(Top12SIGenera[,k]),2)
      #pairwise.wilcox.test 
      wilcox <- (pairwise.wilcox.test(Top12SIGenera[,j], Top12SIGenera[,k], p.adjust.method = 'none'))$p.value
      #condition applied to all factor variables with more than 2 groups
      for (z in 1:ncol(combos)) {
        #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
        my_matrix[a,1] <- colnames(Top12SIGenera)[j]
        #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
        my_matrix[a,2] <- colnames(Top12SIGenera)[k]
        #adding the group combinations to my_matrix
        my_matrix[a,3] <- combos[1,z]
        my_matrix[a,4] <- combos[2,z]
        my_matrix[a,5] <- nlevels(Top12SIGenera[,k])
        my_matrix[a,6] <- sum(Top12SIGenera[,k]==combos[1,z], na.rm = T) #Number of samples in group1
        my_matrix[a,7] <- sum(Top12SIGenera[,k]==combos[2,z], na.rm = T) #Number of samples in group2
        my_matrix[a,8] <- colSums(is.na(Top12SIGenera[1]))
        my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
        my_matrix[a,12] <- mean(subset(Top12SIGenera, Top12SIGenera[c(k)] == combos[1,z], select = j)[,k])
        my_matrix[a,13] <- mean(subset(Top12SIGenera, Top12SIGenera[c(k)] == combos[2,z], select = j)[,k])
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
  my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
  my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')
  
  
  write.table(my_matrix, file = "Top12_SI_Genera.txt", sep = "\t", col.names = T)

*****
PHYLA
*****

Phyla <- read.table(file = "LLD_IBD_PhylumComposition.txt", header = TRUE, sep = "\t",row.names = 1)
Phyla <- Phyla[-c(7)]

my_matrix<-  matrix(nrow = 6 * ncol(Phyla[c(2:ncol(Phyla))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Phyla)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Phyla[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Phyla[,j], Phyla[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Phyla)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Phyla)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Phyla[,k])
      my_matrix[a,6] <- sum(Phyla[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Phyla[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Phyla[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Phyla, Phyla[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Phyla, Phyla[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "TopPhlya.txt", sep = "\t", col.names = T)


********
FAMILIES
********

####______Top HC Families_________

Top10HCFamilies <- read.table(file = "LLD_IBD_FamilyComposition_top10HC.txt", header = TRUE, sep = "\t",row.names = 1)
Top10HCFamilies <- Top10HCFamilies[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10HCFamilies[c(2:ncol(Top10HCFamilies))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10HCFamilies)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10HCFamilies[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10HCFamilies[,j], Top10HCFamilies[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10HCFamilies)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10HCFamilies)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10HCFamilies[,k])
      my_matrix[a,6] <- sum(Top10HCFamilies[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10HCFamilies[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10HCFamilies[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10HCFamilies, Top10HCFamilies[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10HCFamilies, Top10HCFamilies[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_HC_Families.txt", sep = "\t", col.names = T)


####______Top SI Families_________

Top10SIFamilies <- read.table(file = "LLD_IBD_FamilyComposition_top10SI.txt", header = TRUE, sep = "\t",row.names = 1)
Top10SIFamilies <- Top10SIFamilies[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10SIFamilies[c(2:ncol(Top10SIFamilies))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10SIFamilies)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10SIFamilies[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10SIFamilies[,j], Top10SIFamilies[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10SIFamilies)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10SIFamilies)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10SIFamilies[,k])
      my_matrix[a,6] <- sum(Top10SIFamilies[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10SIFamilies[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10SIFamilies[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10SIFamilies, Top10SIFamilies[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10SIFamilies, Top10SIFamilies[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_SI_Families.txt", sep = "\t", col.names = T)


*******
SPECIES
*******

####### _______ Top HC Species _______

Top10HCSpecies <- read.table(file = "LLD_IBD_SpeciesComposition_top10HC.txt", header = TRUE, sep = "\t",row.names = 1)
Top10HCSpecies <- Top10HCSpecies[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10HCSpecies[c(2:ncol(Top10HCSpecies))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10HCSpecies)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10HCSpecies[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10HCSpecies[,j], Top10HCSpecies[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10HCSpecies)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10HCSpecies)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10HCSpecies[,k])
      my_matrix[a,6] <- sum(Top10HCSpecies[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10HCSpecies[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10HCSpecies[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10HCSpecies, Top10HCSpecies[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10HCSpecies, Top10HCSpecies[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_HC_Species.txt", sep = "\t", col.names = T)


####### _______ Top SI Species _______

Top10SISpecies <- read.table(file = "LLD_IBD_SpeciesComposition_top10SI.txt", header = TRUE, sep = "\t",row.names = 1)
Top10SISpecies <- Top10SISpecies[-c(12)]

my_matrix<-  matrix(nrow = 6 * ncol(Top10SISpecies[c(2:ncol(Top10SISpecies))]), ncol = 13)
colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'NGroups', 'NSamplesG1', 'NSampleG2', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'Mean1', 'Mean2')


#__Add correlation analysis data to the matrix 
a=1 
#Iterate through taxonomy columns
for (j in 2:ncol(Top10SISpecies)) {
  #iterate through phenotype columns
  for (k in 1:1) {
    #following conditions applied to all factor variables
    combos <- combn(levels(Top10SISpecies[,k]),2)
    #pairwise.wilcox.test 
    wilcox <- (pairwise.wilcox.test(Top10SISpecies[,j], Top10SISpecies[,k], p.adjust.method = 'none'))$p.value
    #condition applied to all factor variables with more than 2 groups
    for (z in 1:ncol(combos)) {
      #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
      my_matrix[a,1] <- colnames(Top10SISpecies)[j]
      #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
      my_matrix[a,2] <- colnames(Top10SISpecies)[k]
      #adding the group combinations to my_matrix
      my_matrix[a,3] <- combos[1,z]
      my_matrix[a,4] <- combos[2,z]
      my_matrix[a,5] <- nlevels(Top10SISpecies[,k])
      my_matrix[a,6] <- sum(Top10SISpecies[,k]==combos[1,z], na.rm = T) #Number of samples in group1
      my_matrix[a,7] <- sum(Top10SISpecies[,k]==combos[2,z], na.rm = T) #Number of samples in group2
      my_matrix[a,8] <- colSums(is.na(Top10SISpecies[1]))
      my_matrix[a,9] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
      my_matrix[a,12] <- mean(subset(Top10SISpecies, Top10SISpecies[c(k)] == combos[1,z], select = j)[,k])
      my_matrix[a,13] <- mean(subset(Top10SISpecies, Top10SISpecies[c(k)] == combos[2,z], select = j)[,k])
      a=a+1
    }
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix[,9] <- as.numeric(as.character(my_matrix[,9]))
my_matrix[,10] <- p.adjust(c(my_matrix[,9]), method = 'bonferroni')
my_matrix[,11] <- p.adjust(c(my_matrix[,9]), method = 'fdr')


write.table(my_matrix, file = "Top10_SI_Species.txt", sep = "\t", col.names = T)

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


***************
ALPHA DIVERSITY
***************

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


**************
BETA DIVERSITY
**************

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


***********************************
IBD Only (i.e. IBD-NoRes & IBD-Res)
***********************************

StomaIBD <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = CurrentStomaOrPouchType)

StomaIBD$CurrentStomaOrPouchType <- as.numeric(StomaIBD$CurrentStomaOrPouchType)

IBD_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))

#Run PCoA function
plot_pcoa(IBD_taxa,5,StomaIBD) 

#Import table generated from 'plot_pcoa' function
PCoA_valuesIBD <- read.table(file = "../RScripts/pcoa_tableIBD.txt", header = TRUE, sep = "\t",row.names = 1)

PCoA_valuesIBD$CurrentStomaOrPouchType <- factor(PCoA_valuesIBD$CurrentStomaOrPouchType, levels=c("No Resections", "Resections"))

#Plot PCoA graphs
#PCoA1 vs PCoA2

l <- ggplot(PCoA_valuesIBD, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + labs(y="PCoA2", x="PCoA1", color = "Bowel Phenotype") + geom_point(size = 2) + scale_color_manual(labels = c('No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel'), values=c('slateblue2', 'gold1')) + theme_classic() 

ggplotly(l)

******************************************************
Small Intestine Only (i.e. Ileostomy & Ileoanal Pouch)
******************************************************

StomaSI <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = CurrentStomaOrPouchType)

StomaSI$CurrentStomaOrPouchType <- as.numeric(StomaSI$CurrentStomaOrPouchType)

SmallIntestine_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))

#Run PCoA function
plot_pcoa(SmallIntestine_taxa,3,StomaSI) #with 5 PCoA elements - table renamed to pcoa_table5LLD

#Import table generated from 'plot_pcoa' function
PCoA_valuesSI <- read.table(file = "../RScripts/pcoa_tableSI.txt", header = TRUE, sep = "\t",row.names = 1)

#Plot PCoA graphs
#PCoA1 vs PCoA2

ggplot(PCoA_valuesSI, aes(x=V1, y=V2)) + labs(y="PCoA2", x="PCoA1") + geom_point(size = 2, colour = "red2") + theme_classic() 


================================================================================================================================================
REVISED Manuscript - 17.12.20
================================================================================================================================================

#Functions required
Linear_model_2, point 14

***************
ALPHA DIVERSITY
***************

# --- Add all p-values for all possible comparisons ---

my_comparisons=list(c("Resections","IleostomaPouch"),c("No Resections","Resections"),c("HealthyControl","No Resections"),c("No Resections","IleostomaPouch"),c("HealthyControl","Resections"),c("HealthyControl","IleostomaPouch"))

ggplot(Alpha_Stoma1, aes (x=CurrentStomaOrPouchType, y=AlphaDiversity, fill=CurrentStomaOrPouchType)) + geom_violin() + geom_boxplot(width = 0.1, fill = "white") + theme_minimal() + theme(legend.position="none") + theme(axis.text.x = element_text(hjust = 0.5,vjust = 1, size=8,color="black")) + scale_fill_manual(values=c('grey54','slateblue2','gold1','red2')) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test") + ylab ("Shannon Index") + xlab ("Bowel Phenotype") + scale_x_discrete(labels=c('HealthyControl' = 'general population\n N = 1178', 'No Resections' = 'IBD non-resected intestine\n N = 309', 'Resections' = 'IBD resected intestine\n N = 169', 'IleostomaPouch' = 'IBD small intestine\n N = 57'))


# --- Influence of confounders on alpha diversity ---

SIvsColon <- LLD.IBD_MetaSp15.asin[c("DiagnosisCurrent", "NumberOfResectionsAny", "IleocecalValveInSitu", "AgeAtFecalSampling", 
                                     "Sex", "PPI", "SmokeExposedToSmokeForMoreThanOneHour", "BowelMovementADayDef",
                                     "P_plant_en","SUMOFKCAL", "BMI", "P_animal_en", "PFReads", "Alc_en" ,"antibiotics_merged",
                                     "metformin", "beta_blockers", "FecalCalprotectinOver200yesno", "calcium", "vitamin_D", "folic_acid",
                                     "paracetamol","mesalazines", "vitamin_B12", "oral_steroid", "ferrum", "CurrentStomaOrPouchType",colnames(LLD.IBD_MetaSp15.asin)[165:298])] # See section 6 for "LLD.IBD_MetaSp15.asin" dataframe

#SIvsColon$CurrentStomaOrPouchType <- fct_collapse(SIvsColon$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))

### --- Impute missing values ---

SIvsColon$NumberOfResectionsAny <- ifelse(SIvsColon$DiagnosisCurrent == "generalpopulation" & is.na(SIvsColon[2])== TRUE, 0, SIvsColon[,2])
SIvsColon$NumberOfResectionsAny <- as.numeric(SIvsColon$NumberOfResectionsAny)

SIvsColon$IleocecalValveInSitu <- ifelse(SIvsColon$DiagnosisCurrent == "generalpopulation" & is.na(SIvsColon[3])== TRUE, "yes", as.character(SIvsColon[,3]))
SIvsColon$IleocecalValveInSitu <- as.factor(SIvsColon$IleocecalValveInSitu)


# Impute median for numerical variables & least presumptive for categorical
for (i in 1:ncol(SIvsColon)) {
  if(is.numeric(SIvsColon[,i])){
    #list_na <- colnames(SIvsColon)[apply(SIvsColon, 2, anyNA)]
    SIvsColon[,i] <- ifelse(is.na(SIvsColon[,i]),
                            median(SIvsColon[,i], na.rm = TRUE),
                            as.numeric(SIvsColon[,i]))
  }
  if(length(levels(SIvsColon[,i])) == 2){
    SIvsColon[,i] <- ifelse(is.na(SIvsColon[,i]),
                            "no",
                            as.character(SIvsColon[,i]))
    SIvsColon[,i] <- as.factor(SIvsColon[,i])
  }
}

SIvsColon$SmokeExposedToSmokeForMoreThanOneHour <- ifelse(is.na(SIvsColon$SmokeExposedToSmokeForMoreThanOneHour),
                          "RarelyOrNever",
                          as.character(SIvsColon[,7]))
SIvsColon$SmokeExposedToSmokeForMoreThanOneHour <- as.factor(SIvsColon$SmokeExposedToSmokeForMoreThanOneHour)

# --- Linear regression ---

Alpha1 <- as.data.frame(diversity(LLD_IBD_MetTax15[,c(747:880)], index="shannon"))
colnames(Alpha1)[1] <- 'AlphaDiversity' 

SIvsColon2 <- cbind(Alpha1,SIvsColon)

SIvsColon2$CurrentStomaOrPouchType <- factor(SIvsColon2$CurrentStomaOrPouchType, levels = c("IleostomaPouch","HealthyControl", "No Resections", "Resections"))

SIvsColon2$SmokeExposedToSmokeForMoreThanOneHour <- factor(SIvsColon2$SmokeExposedToSmokeForMoreThanOneHour, levels = c("RarelyOrNever","OncePerWeek", "MoreThanOncePerWeek", "Daily"))

x = SIvsColon2[,-c(2,3,4,24,29:162)]
my_Lm_sum <- Linear_model_2(x) 

write.table(my_Lm_sum, file = "AlphaDiversity_AfterCovariateCorrection.txt", sep = "\t", col.names = T)

***************
BETA DIVERSITY
***************

#With centroids
gd <- PCoA_values %>% 
  group_by(CurrentStomaOrPouchType) %>% 
  summarise(V1 = mean(V1),
            V2  = mean(V2))


ggplot(PCoA_values, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + 
  labs(y="PCo2", x="PCo1", color = "Bowel Phenotype") + geom_point(size = 2) + 
  scale_color_manual(labels = c('HealthyControl' = 'general population', 'No Resections' = 'IBD non-resected bowel', 'Resections' = 'IBD resected bowel', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + 
  theme_classic() + geom_point(data = gd, size = 5, shape=24, fill = c('grey44','slateblue3', 'gold2', 'red3'), color = 'black')


# Add variance explained percentages to the x and y axis

beta <- vegdist(LLD_IBD_TaxaS1, method="bray")
my_pcoa <- cmdscale(beta, eig = T)
Var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)


ggplot(PCoA_values, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + 
  labs(y=paste("PCo2 (", Var_expl[2],"%)", sep="") , x=paste("PCo1 (", Var_expl[1],"%)", sep=""), color = "Bowel Phenotype") + geom_point(size = 2) + 
  scale_color_manual(labels = c('HealthyControl' = 'General Population', 'No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1', 'red2')) + 
  theme_classic() + geom_point(data = gd, size = 5, shape=24, fill = c('grey44','slateblue3', 'gold2', 'red3'), color = 'black') 

## --- IBD Only ---

PCoA_valuesIBD <- read.table(file = "../RScripts/pcoa_tableIBD.txt", header = TRUE, sep = "\t",row.names = 1)
PCoA_valuesIBD$CurrentStomaOrPouchType <- factor(PCoA_valuesIBD$CurrentStomaOrPouchType, levels=c("No Resections", "Resections"))

StomaIBD <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = CurrentStomaOrPouchType)
StomaIBD$CurrentStomaOrPouchType <- as.numeric(StomaIBD$CurrentStomaOrPouchType)
IBD_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))

beta1 <- vegdist(IBD_taxa, method="bray")
my_pcoa1 <- cmdscale(beta1, eig = T)
Var_expl1 <- round(my_pcoa1$eig/sum(my_pcoa1$eig)*100,digits = 1)

gd1 <- PCoA_valuesIBD %>% 
  group_by(CurrentStomaOrPouchType) %>% 
  summarise(V1 = mean(V1),
            V2  = mean(V2))
            
#Plot PCoA graphs
#PCo1 vs PCo2

ggplot(PCoA_valuesIBD, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + 
  labs(y=paste("PCo2 (", Var_expl1[2],"%)", sep="") , x=paste("PCo1 (", Var_expl1[1],"%)", sep=""), color = "Bowel Phenotype") + geom_point(size = 2) + 
  scale_color_manual(labels = c('No Resections' = 'IBD non-resected intestine', 'Resections' = 'IBD resected intestine'), values=c('slateblue2', 'gold1')) + 
  theme_classic() + geom_point(data = gd1, size = 5, shape=24, fill = c('slateblue3', 'gold2'), color = 'black') 

            
## --- Small Intestine only ---

StomaSI <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = CurrentStomaOrPouchType)
StomaSI$CurrentStomaOrPouchType <- as.numeric(StomaSI$CurrentStomaOrPouchType)
SmallIntestine_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))

PCoA_valuesSI <- read.table(file = "../RScripts/pcoa_tableSI.txt", header = TRUE, sep = "\t",row.names = 1)

beta2 <- vegdist(SmallIntestine_taxa, method="bray")
my_pcoa2 <- cmdscale(beta2, eig = T)
Var_expl2 <- round(my_pcoa2$eig/sum(my_pcoa2$eig)*100,digits = 1)

gd2 <- PCoA_valuesSI %>% 
  group_by(CurrentStomaOrPouchType) %>% 
  summarise(V1 = mean(V1),
            V2  = mean(V2))

#Plot PCoA graphs
#PCoA1 vs PCoA2

ggplot(PCoA_valuesSI, aes(x=V1, y=V2)) + labs(y=paste("PCo2 (", Var_expl2[2],"%)", sep="") , x=paste("PCo1 (", Var_expl2[1],"%)", sep="")) + 
  geom_point(size = 2, colour = "red2") + theme_classic() + geom_point(data = gd2, size = 5, shape=24, fill = c('red2'), color = 'black') 


## --- IBD-Res stratified according to ileocecal valve status ---

#Import table generated from 'plot_pcoa' function
PCoA_values <- read.table(file = "../RScripts/pcoa_tableALL.txt", header = TRUE, sep = "\t", quote = "\\") 
#rownames(PCoA_values) <- PCoA_values$X
#PCoA_values$X <- NULL

Valve <- LLD_IBD_MetTax15[c(231)] 
merge1 <- merge(Valve,PCoA_values, by = "row.names")
merge1$IleocecalValveInSitu[is.na(merge1$IleocecalValveInSitu)] <- "yes"

PCoA_valve <- merge1 %>% 
  mutate(CurrentStomaOrPouchType = ifelse(IleocecalValveInSitu =="no", gsub('Resections', 'NotInSitu', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

PCoA_valve$CurrentStomaOrPouchType <- factor(PCoA_valve$CurrentStomaOrPouchType, levels=c("HealthyControl", "No Resections","Resections","NotInSitu", "IleostomaPouch"))
rownames(PCoA_valve) <- PCoA_valve$Row.names
PCoA_valve[,1:2] <- NULL

#With centroids
gd <- PCoA_valve %>% 
  group_by(CurrentStomaOrPouchType) %>% 
  summarise(V1 = mean(V1),
            V2  = mean(V2))

beta <- vegdist(LLD_IBD_TaxaS1, method="bray")
my_pcoa <- cmdscale(beta, eig = T)
Var_expl <- round(my_pcoa$eig/sum(my_pcoa$eig)*100,digits = 1)

ggplot(PCoA_valve, aes(x=V1, y=V2, color= CurrentStomaOrPouchType)) + 
  labs(y=paste("PCo2 (", Var_expl[2],"%)", sep="") , x=paste("PCo1 (", Var_expl[1],"%)", sep=""), color = "Bowel Phenotype") + geom_point(size = 2) + 
  scale_color_manual(labels = c('HealthyControl' = 'General population', 'No Resections' = 'IBD-NoRes', 'Resections' = 'IBD-Res_InSitu', 'NotInSitu' = 'IBD-Res_NotInSitu', 'IleostomaPouch' = 'IBD small intestine'), values=c('grey54','slateblue2', 'gold1','#FA8072', 'red2')) + 
  theme_classic() + geom_point(data = gd, size = 5, shape=24, fill = c('grey44','slateblue3', 'gold2','#FA8072', 'red3'), color = 'black')


## --- PERMANOVA with covariates --- 

SI_NonSI_taxa <- LLD.IBD_MetaSp15.asin[c(165:ncol(LLD.IBD_MetaSp15.asin))]
colnames(SIvsColon)[colSums(is.na(SIvsColon)) > 0]
table(is.na(Test$P_plant_en))

SIvsColon <- LLD.IBD_MetaSp15.asin[c("DiagnosisCurrent", "NumberOfResectionsAny", "IleocecalValveInSitu", "AgeAtFecalSampling", 
                                       "Sex", "PPI", "SmokeExposedToSmokeForMoreThanOneHour", "BowelMovementADayDef",
                                       "P_plant_en","SUMOFKCAL", "BMI", "P_animal_en", "PFReads", "Alc_en" ,"antibiotics_merged",
                                       "metformin", "beta_blockers", "FecalCalprotectinOver200yesno", "calcium", "vitamin_D", "folic_acid",
                                     "paracetamol","mesalazines", "vitamin_B12", "oral_steroid", "ferrum", "CurrentStomaOrPouchType",colnames(LLD.IBD_MetaSp15.asin)[165:298])]

#SIvsColon$CurrentStomaOrPouchType <- fct_collapse(SIvsColon$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))

#Impute missing values 

SIvsColon$NumberOfResectionsAny <- ifelse(SIvsColon$DiagnosisCurrent == "generalpopulation" & is.na(SIvsColon[2])== TRUE, 0, SIvsColon[,2])
SIvsColon$NumberOfResectionsAny <- as.numeric(SIvsColon$NumberOfResectionsAny)

SIvsColon$IleocecalValveInSitu <- ifelse(SIvsColon$DiagnosisCurrent == "generalpopulation" & is.na(SIvsColon[3])== TRUE, "yes", as.character(SIvsColon[,3]))
SIvsColon$IleocecalValveInSitu <- as.factor(SIvsColon$IleocecalValveInSitu)


table(is.na(SIvsColon$IleocecalValveInSitu ),SIvsColon$DiagnosisCurrent)
table(SIvsColon$DiagnosisCurrent,SIvsColon$IleocecalValveInSitu)
table(is.na(SIvsColon$SmokeExposedToSmokeForMoreThanOneHour ),SIvsColon$DiagnosisCurrent)
table(SIvsColon$DiagnosisCurrent,SIvsColon$SmokeExposedToSmokeForMoreThanOneHour)


# Impute median for numerical variables & least presumptive for categorical
for (i in 1:ncol(SIvsColon)) {
  if(is.numeric(SIvsColon[,i])){
    #list_na <- colnames(SIvsColon)[apply(SIvsColon, 2, anyNA)]
    SIvsColon[,i] <- ifelse(is.na(SIvsColon[,i]),
                            median(SIvsColon[,i], na.rm = TRUE),
                            as.numeric(SIvsColon[,i]))
  }
  if(length(levels(SIvsColon[,i])) == 2){
    SIvsColon[,i] <- ifelse(is.na(SIvsColon[,i]),
                            "no",
                            as.character(SIvsColon[,i]))
  }
}

SIvsColon$SmokeExposedToSmokeForMoreThanOneHour <- ifelse(is.na(SIvsColon$SmokeExposedToSmokeForMoreThanOneHour),
                                                          "RarelyOrNever",
                                                          as.character(SIvsColon[,7]))

# PERMANOVA without IBD specific variables
ad1 <- adonis(formula = SI_NonSI_taxa ~ AgeAtFecalSampling+ 
               Sex + PPI + SmokeExposedToSmokeForMoreThanOneHour+ BowelMovementADayDef+
               P_plant_en+SUMOFKCAL+BMI+P_animal_en+PFReads+Alc_en +antibiotics_merged+
               metformin+beta_blockers+FecalCalprotectinOver200yesno+calcium+vitamin_D+folic_acid+
               paracetamol+vitamin_B12+oral_steroid+ferrum+CurrentStomaOrPouchType, data = SIvsColon, permutations = 9999, method = "bray")
aov_table1 <- ad1$aov.tab
adonis_All_3 <- matrix(ncol = 4, nrow=23) 
adonis_All_3[,4]=nrow(SIvsColon)
adonis_All_3[,1]=aov_table1[1:23,1]
adonis_All_3[,2]=aov_table1[1:23,5]
adonis_All_3[,3]=aov_table1[1:23,6]

adonis_All_3 <- as.data.frame(adonis_All_3)

colnames(adonis_All_3) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_All_3)=rownames(aov_table1)[1:23]
adonis_All_3$FDR <- p.adjust(c(adonis_All_3[,3]), method = 'fdr')

write.table(adonis_All_3, file = "PERMANOVA_withCovariates.txt", sep = "\t", col.names = T)

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


*****************************************
Within SI: ...Ileostomy vs Ileoanal pouch
*****************************************

LLD.IBD_Final$NumberIndicationabcessOrfistula <- as.numeric(as.character(LLD.IBD_Final$NumberIndicationabcessOrfistula))

#Remove 'correlated/related' variables
LLD.IBD_Final_x <- LLD.IBD_Final[c(1,2,5:13,15,18,21:31,33,36:43,45,47,49,52:108,110,111,114:134,302,310:317,319,321,323:370,373,377,379,381:385,387,389,403:1054)]
LLD.IBD_Final_x <- LLD.IBD_Final_x[-c(74,76:78,80:93,96,98,113)]

# Remove bacteria from dataframe with a prevalence < 15%
LLD.IBD_MetaSp15_x <- Prevalence.filter(LLD.IBD_Final_x, RowSpecies = 165)

#Transform species abundance values
LLD.IBD_MetaSp15.asin_x <-transform_and_filter_taxa(LLD.IBD_MetaSp15_x[,c(165:298)])
LLD.IBD_MetaSp15.asin_x <-cbind(LLD.IBD_MetaSp15_x[c(1:164)],LLD.IBD_MetaSp15.asin_x)

Ileo_Pouch <- subset(LLD.IBD_MetaSp15.asin_x, CurrentStomaOrPouchType == 'ileostoma' | CurrentStomaOrPouchType == 'pouch', select = c(76,165:298))
Ileo_Pouch$CurrentStomaOrPouchType <- fct_drop(Ileo_Pouch$CurrentStomaOrPouchType)


#For loop to correlate Resections type with each bacterial species

my_matrix <-  matrix(nrow = ncol(Ileo_Pouch) - 1, ncol = 10)
colnames(my_matrix) <- c('Taxonomy', 'Group1', 'Group2', 'NGroup1', 'NGroup2', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'MeanGroup1', 'MeanGroup2')

combos <- combn(levels(Ileo_Pouch[,1]),2)
a=1
for (i in 2:ncol(Ileo_Pouch)) {
  wilcox <- (pairwise.wilcox.test(Ileo_Pouch[,i], Ileo_Pouch[,1], p.adjust.method = 'none'))$p.value
  for (z in 1:ncol(combos)) {
    my_matrix[a,1] <- colnames(Ileo_Pouch)[i]
    my_matrix[a,2] <- combos[1,z]
    my_matrix[a,3] <- combos[2,z]
    my_matrix[a,4] <- sum(Ileo_Pouch[,1]==combos[1,z], na.rm = T) #Number of samples in group1
    my_matrix[a,5] <- sum(Ileo_Pouch[,1]==combos[2,z], na.rm = T)
    my_matrix[a,6] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test coefficient values 
    my_matrix[a,9] <- mean(subset(Ileo_Pouch, Ileo_Pouch[c(1)] == combos[1,z], select = i)[,1]) #Mean taxonomy relative abundances group 3
    my_matrix[a,10] <- mean(subset(Ileo_Pouch, Ileo_Pouch[c(1)] == combos[2,z], select = i)[,1]) #Mean taxonomy relative abundances group 2
    a=a+1
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')

my_matrix$Taxonomy <- as.character(gsub("_", " ", my_matrix$Taxonomy))

write.table(my_matrix, file = 'Univariate_Ileostoma_vs_Pouch.txt', sep = "\t",col.names = TRUE)


**********************
Within SI: ...CD vs UC 
**********************

#Univariate Analysis
IleostomaPouch <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == 'IleostomaPouch', select = c(15,165:298))
IleostomaPouch[10,1] <- "UC"  #temporarily change IBDU to UC because by definition should not have small intestinal inflammation 
IleostomaPouch$DiagnosisCurrent <- fct_drop(IleostomaPouch$DiagnosisCurrent)

#For loop to correlate IBD type with each bacterial species

my_matrix <-  matrix(nrow = ncol(IleostomaPouch) - 1, ncol = 10)
colnames(my_matrix) <- c('Taxonomy', 'Group1', 'Group2', 'NGroup1', 'NGroup2', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'MeanGroup1', 'MeanGroup2')

combos <- combn(levels(IleostomaPouch[,1]),2)
a=1
for (i in 2:ncol(IleostomaPouch)) {
  wilcox <- (pairwise.wilcox.test(IleostomaPouch[,i], IleostomaPouch[,1], p.adjust.method = 'none'))$p.value
  for (z in 1:ncol(combos)) {
    my_matrix[a,1] <- colnames(IleostomaPouch)[i]
    my_matrix[a,2] <- combos[1,z]
    my_matrix[a,3] <- combos[2,z]
    my_matrix[a,4] <- sum(IleostomaPouch[,1]==combos[1,z], na.rm = T) #Number of samples in group1
    my_matrix[a,5] <- sum(IleostomaPouch[,1]==combos[2,z], na.rm = T)
    my_matrix[a,6] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test coefficient values 
    my_matrix[a,9] <- mean(subset(IleostomaPouch, IleostomaPouch[c(1)] == combos[1,z], select = i)[,1]) #Mean taxonomy relative abundances group 3
    my_matrix[a,10] <- mean(subset(IleostomaPouch, IleostomaPouch[c(1)] == combos[2,z], select = i)[,1]) #Mean taxonomy relative abundances group 2
    a=a+1
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')

my_matrix$Taxonomy <- as.character(gsub("_", " ", my_matrix$Taxonomy))

write.table(my_matrix, file = 'Univariate_CD_vs_UC.txt', sep = "\t",col.names = TRUE)


***********************************************
Within SI: ...Ileal vs Colonic disease location
***********************************************

IleostomaPouch1 <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == 'IleostomaPouch', select = c(18,165:298))
IleostomaPouch1[c(44,54),1] <- "both"  #temporarily change ileum to both  
IleostomaPouch1$DiseaseLocation <- fct_drop(IleostomaPouch1$DiseaseLocation)

#For loop to correlate disease location with each bacterial species

my_matrix <-  matrix(nrow = ncol(IleostomaPouch1) - 1, ncol = 10)
colnames(my_matrix) <- c('Taxonomy', 'Group1', 'Group2', 'NGroup1', 'NGroup2', 'P-value', 'BonferroniAdjust', 'FDRadjust', 'MeanGroup1', 'MeanGroup2')

combos <- combn(levels(IleostomaPouch1[,1]),2)
a=1
for (i in 2:ncol(IleostomaPouch1)) {
  wilcox <- (pairwise.wilcox.test(IleostomaPouch1[,i], IleostomaPouch1[,1], p.adjust.method = 'none'))$p.value
  for (z in 1:ncol(combos)) {
    my_matrix[a,1] <- colnames(IleostomaPouch1)[i]
    my_matrix[a,2] <- combos[1,z]
    my_matrix[a,3] <- combos[2,z]
    my_matrix[a,4] <- sum(IleostomaPouch1[,1]==combos[1,z], na.rm = T) #Number of samples in group1
    my_matrix[a,5] <- sum(IleostomaPouch1[,1]==combos[2,z], na.rm = T)
    my_matrix[a,6] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test coefficient values 
    my_matrix[a,9] <- mean(subset(IleostomaPouch1, IleostomaPouch1[c(1)] == combos[1,z], select = i)[,1]) #Mean taxonomy relative abundances group 3
    my_matrix[a,10] <- mean(subset(IleostomaPouch1, IleostomaPouch1[c(1)] == combos[2,z], select = i)[,1]) #Mean taxonomy relative abundances group 2
    a=a+1
  }
}
my_matrix <- as.data.frame(my_matrix)
my_matrix$`P-value` <- as.numeric(as.character(my_matrix$`P-value`))
my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')

my_matrix$Taxonomy <- as.character(gsub("_", " ", my_matrix$Taxonomy))

write.table(my_matrix, file = 'Univariate_DiseaseLocation_ileumcolon_vs_colon.txt', sep = "\t",col.names = TRUE)

```


### 7. Maaslin Multivariate analyses

```{r}

#Packages required:
#install.packages(c('agricolae','gam','gamlss','gbm','glmnet','inlinedocs','logging','MASS','nlme','optparse','outliers','penalized','pscl','robustbase'))
install.packages('agricolae')
install.packages('gam')
install.packages('gamlss')
install.packages('gbm')
install.packages('glmnet')
install.packages('inlinedocs')
install.packages('logging')
install.packages('MASS')
install.packages('nlme')
install.packages('optparse')
install.packages('outliers')
install.packages('penalized')
install.packages('pscl')
install.packages('robustbase')

install.packages("https://bitbucket.org/biobakery/maaslin/downloads/Maaslin_0.0.5.tar.gz", repo=NULL, type="source")

library(Maaslin)
library(forcats)
library(ggplot2)
library(plotly)
library(wesanderson)
library(RColorBrewer)
library(reshape2)
FactorVariableHeatmapMaaslin #see 'Functions' file, point 11.


#***IleostomaPouch vs Healthy Controls***

#Prepare input data for Maaslin 

# 1. Metadata & Taxa
# Dataframe needed: LLD.IBD_MetaSp15 (i.e. non-transformed) (see rscript 'UnivariateAnalysis260320' for dataframe)

Maaslin_MetaTaxa1 <- subset(LLD.IBD_MetaSp15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'HealthyControl')
Maaslin_MetaTaxa1 <- Maaslin_MetaTaxa1[c(1,2,3,9,14,74,76,96,98,99,102,113,115,119,120,121,128,133,135,138,147,148,163,165:298)]
Maaslin_MetaTaxa1$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaTaxa1$CurrentStomaOrPouchType)

#Set sample IDs as column 1
#Maaslin_MetaTaxa1$SampleID <- rownames(Maaslin_MetaTaxa1) 
#Maaslin_MetaTaxa1 <- Maaslin_MetaTaxa1[c(159,1:158)]

for (i in 24:157) {
  Maaslin_MetaTaxa1[,i] <- Maaslin_MetaTaxa1[,i]/100
}     # Relative abundance should be between 0 - 1


write.table(Maaslin_MetaTaxa1, file = "MaaslinMetaTaxa_HealthControl_vs_IleostomaPouch.tsv", sep = "\t", quote = F)
Maaslin(strInputTSV = "MaaslinMetaTaxa_HealthControl_vs_IleostomaPouch.tsv", strOutputDIR = "Maaslin_HC_vs_Ileostoma_Output", strInputConfig = "MaasLinConfig_HC_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0)

#MaaslinTest <- read.table(file = "../Maaslin/MaaslinMetaTaxa_HealthControl_vs_IleostomaPouch.tsv", header = TRUE, sep = "\t", quote = "\\") #Metadata



#***Within IBD __ IleostomaPouch vs No Resection***

# 1. Metadata & Taxa
Maaslin_MetaTaxa2 <- subset(LLD.IBD_MetaSp15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'No Resections')
Maaslin_MetaTaxa2 <- Maaslin_MetaTaxa2[c(1,2,3,9,14,15,74,76,96,98,99,102,113,115,119,120,121,127,128,133,135,138,147,148,163,165:298)]
Maaslin_MetaTaxa2$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaTaxa2$CurrentStomaOrPouchType)

Maaslin_MetaTaxa2$CurrentStomaOrPouchType <- relevel(Maaslin_MetaTaxa2$CurrentStomaOrPouchType, ref = 'No Resections') #Change Reference (Baseline) Category

Maaslin_MetaTaxa2$DiagnosisCurrent <- fct_drop(Maaslin_MetaTaxa2$DiagnosisCurrent)

for (i in 26:159) {
  Maaslin_MetaTaxa2[,i] <- Maaslin_MetaTaxa2[,i]/100
} 

write.table(Maaslin_MetaTaxa2, file = "MaaslinMetaTaxa_WithinIBD_NoResections_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinConfig_NoResection_vs_IleostomaPouch.read.config.R

Maaslin(strInputTSV = "MaaslinMetaTaxa_WithinIBD_NoResections_vs_IleostomaPouch.tsv", strOutputDIR = "Maaslin_NoResection_vs_Ileostoma_Output", strInputConfig = "MaaslinConfig_NoResection_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0)


#***Within IBD __ IleostomaPouch vs Resections***

# 1. Metadata & Taxa
Maaslin_MetaTaxa3 <- subset(LLD.IBD_MetaSp15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'Resections')
Maaslin_MetaTaxa3 <- Maaslin_MetaTaxa3[c(1,2,3,9,14,15,74,76,79,82,96,98,99,102,113,115,119,120,121,127,128,133,135,138,147,148,163,165:298)]
Maaslin_MetaTaxa3$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaTaxa3$CurrentStomaOrPouchType)

Maaslin_MetaTaxa3$CurrentStomaOrPouchType <- relevel(Maaslin_MetaTaxa3$CurrentStomaOrPouchType, ref = 'Resection NoStoma') #Change Reference (Baseline) Category

#Remove levels with 0 samples
Maaslin_MetaTaxa3$DiagnosisCurrent <- fct_drop(Maaslin_MetaTaxa3$DiagnosisCurrent)

for (i in 28:161) {
  Maaslin_MetaTaxa3[,i] <- Maaslin_MetaTaxa3[,i]/100
} 

write.table(Maaslin_MetaTaxa3, file = "MaaslinMetaTaxa_WithinIBD_Resections_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinConfig_ResectionNoStoma_IleostomaPouch.read.config.R

Maaslin(strInputTSV = "MaaslinMetaTaxa_WithinIBD_Resections_vs_IleostomaPouch.tsv", strOutputDIR = "Maaslin_Resections_vs_Ileostoma_Output", strInputConfig = "MaaslinConfig_ResectionNoStoma_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0)



#***Within IBD __ IleostomaPouch vs Resections IleoCecalValveInSitu***

# 1. Metadata & Taxa
Maaslin_MetaTaxa3a <- subset(Maaslin_MetaTaxa3, CurrentStomaOrPouchType == 'Resections' & IleocecalValveInSitu == 'yes')
Maaslin_MetaTaxa3b <- rbind(subset(Maaslin_MetaTaxa3, CurrentStomaOrPouchType == 'IleostomaPouch'), Maaslin_MetaTaxa3a)
Maaslin_MetaTaxa3b$IleocecalValveInSitu <- NULL

write.table(Maaslin_MetaTaxa3b, file = "MaaslinMetaTaxa_WithinIBD_ResectionsWithIleoCecalValveInSitu_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaasLinConfig_HC_vs_IleostomaPouch.read.config.R


Maaslin(strInputTSV = "MaaslinMetaTaxa_WithinIBD_ResectionsWithIleoCecalValveInSitu_vs_IleostomaPouch.tsv", strOutputDIR = "Maaslin_ResectionsCecalValveInSitu_vs_Ileostoma_Output", strInputConfig = "MaasLinConfig_HC_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0)



#Display results
#Heatmap

Healthy_IBD.SI.Maaslin <- read.table(file = "../Maaslin/Maaslin_HC_vs_Ileostoma_Output/MaaslinMetaTaxa_HealthControl_vs_IleostomaPouch-CurrentStomaOrPouchType.txt", header = TRUE, sep = "\t") 

IBD.NoResec_SI.Maaslin <- read.table(file = "../Maaslin/Maaslin_NoResection_vs_Ileostoma_Output/MaaslinMetaTaxa_WithinIBD_NoResections_vs_IleostomaPouch-CurrentStomaOrPouchType.txt", header = TRUE, sep = "\t")

IBD.Resec_SI.Maaslin <- read.table(file = "../Maaslin/Maaslin_Resections_vs_Ileostoma_Output/MaaslinMetaTaxa_WithinIBD_Resections_vs_IleostomaPouch-CurrentStomaOrPouchType.txt", header = TRUE, sep = "\t") 

IBD.IleoValve_SI.Maaslin <- read.table(file = "../Maaslin/Maaslin_ResectionsCecalValveInSitu_vs_Ileostoma_Output/MaaslinMetaTaxa_WithinIBD_ResectionsWithIleoCecalValveInSitu_vs_IleostomaPouch-CurrentStomaOrPouchType.txt", header = TRUE, sep = "\t") 

Healthy_IBD.SI.Maaslin$Value <- as.factor(gsub("CurrentStomaOrPouchTypeIleostomaPouch", "IleostomaPouch01", Healthy_IBD.SI.Maaslin$Value ))
IBD.NoResec_SI.Maaslin$Value <- as.factor(gsub("CurrentStomaOrPouchTypeNo Resections", "IleostomaPouch02", IBD.NoResec_SI.Maaslin$Value))
IBD.Resec_SI.Maaslin$Value <- as.factor(gsub("CurrentStomaOrPouchTypeResections", "IleostomaPouch03", IBD.Resec_SI.Maaslin$Value))
IBD.IleoValve_SI.Maaslin$Value <- as.factor(gsub("CurrentStomaOrPouchTypeResections", "IleostomaPouch04", IBD.IleoValve_SI.Maaslin$Value))

Healthy_IBD.SI.Maaslin <- FactorVariableHeatmapMaaslin(Healthy_IBD.SI.Maaslin)
IBD.NoResec_SI.Maaslin <- FactorVariableHeatmapMaaslin(IBD.NoResec_SI.Maaslin)
IBD.Resec_SI.Maaslin <- FactorVariableHeatmapMaaslin(IBD.Resec_SI.Maaslin)
IBD.IleoValve_SI.Maaslin <- FactorVariableHeatmapMaaslin(IBD.IleoValve_SI.Maaslin)

colnames(Healthy_IBD.SI.Maaslin)[9] <- 'Maaslin_01'
colnames(IBD.NoResec_SI.Maaslin)[9] <- 'Maaslin_02'
colnames(IBD.Resec_SI.Maaslin)[9] <- 'Maaslin_03'
colnames(IBD.IleoValve_SI.Maaslin)[9] <- 'Maaslin_04'

# the direction for the last 4 Maaslin test values refer to the non IleostomaPouch group - multiply logFDR column by -1 
IBD.NoResec_SI.Maaslin$Maaslin_02 <- IBD.NoResec_SI.Maaslin$Maaslin_02 * -1
IBD.Resec_SI.Maaslin$Maaslin_03 <- IBD.Resec_SI.Maaslin$Maaslin_03 * -1
IBD.IleoValve_SI.Maaslin$Maaslin_04 <- IBD.IleoValve_SI.Maaslin$Maaslin_04 * -1

Merge10 <- merge(Healthy_IBD.SI.Maaslin, IBD.NoResec_SI.Maaslin, by = 'Feature', all = T)
Merge10 <- Merge10[c(1,9,17)]

Merge11 <- merge(Merge10, IBD.Resec_SI.Maaslin, by = 'Feature', all = T)
Merge11<- Merge11[c(1,2,3,11)]

#Merge12 <- merge(Merge11, IBD.IleoValve_SI.Maaslin, by = 'Feature', all = T)
#Merge12<- Merge12[c(1,2,3,4,12)]

#plot heatmap
Merge11$Feature <- as.character(gsub("_", " ", Merge11$Feature))

Heatmap_Maaslin <- melt(Merge11)
Heatmap_Maaslin$value <- as.factor(Heatmap_Maaslin$value)
Heatmap_Maaslin[is.na(Heatmap_Maaslin)] <- 0

cols <- c("-3" = "#2171B5", "-2" = "#6BAED6", "-1" = "#C6DBEF", "0" = "white",  "1" = "#FCBBA1", "2" = "#FB6A4A", "3" = "#CB181D")

w <- ggplot(Heatmap_Maaslin, aes(x = variable, y = Feature, fill= value)) + geom_tile(color="gray80") + scale_fill_manual(values = cols) + theme_classic() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size=5))

ggplotly(w)


================================================================================================================================================
REVISED Manuscript - 17.12.20
================================================================================================================================================

#IBD Non-Resection vs IBD Resection
# 1. Metadata & Taxa
Maaslin_MetaTaxa4 <- subset(LLD.IBD_MetaSp15, CurrentStomaOrPouchType == 'No Resections' | CurrentStomaOrPouchType == 'Resections')
Maaslin_MetaTaxa4 <- Maaslin_MetaTaxa4[c(1,2,3,9,14,15,74,76,96,98,99,102,113,115,119,120,121,127,128,133,135,138,147,148,163,165:298)]
Maaslin_MetaTaxa4$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaTaxa4$CurrentStomaOrPouchType)

Maaslin_MetaTaxa4$DiagnosisCurrent <- fct_drop(Maaslin_MetaTaxa4$DiagnosisCurrent)

for (i in 26:159) {
  Maaslin_MetaTaxa4[,i] <- Maaslin_MetaTaxa4[,i]/100
} 
write.table(Maaslin_MetaTaxa4, file = "MaaslinMetaTaxa_WithinIBD_NoResections_vs_Resections.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinConfig_ResectionNoStoma_IleostomaPouch.read.config.R

Maaslin(strInputTSV = "MaaslinMetaTaxa_WithinIBD_NoResections_vs_Resections.tsv", strOutputDIR = "Maaslin_NoResections_vs_Resections_Output", strInputConfig = "MaaslinConfig_ResectionNoStoma_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0)


## -- Display results --
## LOLLIPOP PLOT - IBD-Res vs IBD-NoRes (Maaslin)

IBD.NoRes_Res.Maaslin <- read.table(file = "../Maaslin/Maaslin_NoResections_vs_Resections_Output/MaaslinMetaTaxa_WithinIBD_NoResections_vs_Resections-CurrentStomaOrPouchType.txt", header = TRUE, sep = "\t")

IBD.NoRes_Res.Maaslin$Feature <- as.character(gsub("_", " ", IBD.NoRes_Res.Maaslin$Feature))

IBD.NoRes_Res.Maaslin <- FactorVariableHeatmapMaaslin(IBD.NoRes_Res.Maaslin)

colnames(IBD.NoRes_Res.Maaslin)[9] <- 'Significance'
IBD.NoRes_Res.Maaslin$flag <- ifelse(IBD.NoRes_Res.Maaslin$Significance < 0, "Negative", ifelse(IBD.NoRes_Res.Maaslin$Significance > 0, "Positive", "Zero"))

ggplot(IBD.NoRes_Res.Maaslin, aes(x = Significance, y = Feature, color = flag)) +
        geom_segment(aes(x = 0, y = Feature, xend = Significance, yend = Feature)) +
        geom_point() +
        scale_colour_manual(values = c("blue", "red", "black")) + theme_minimal() + theme (legend.position = "none") + xlab("IBD-Res vs IBD-NoRes") + ylab("Species")

#Significant associations only

IBD.NoRes_Res.Maaslin_signif <- subset(IBD.NoRes_Res.Maaslin, Significance > 0 | Significance < 0)

ggplot(IBD.NoRes_Res.Maaslin_signif, aes(x = Significance, y = Feature, color = flag)) +
        geom_segment(aes(x = 0, y = Feature, xend = Significance, yend = Feature)) +
        geom_point() +
        scale_colour_manual(values = c("blue", "red")) + theme_minimal() + theme (legend.position = "none") + xlab("IBD-Res vs IBD-NoRes") + ylab("Species")


```

### 8. ADONIS analysis

```{r}

#Packages & functions required:
library(vegan)
library (dplyr)
RemoveColumnsNA #see 'Functions' file, point 10.

Final_phenotypes1 <- read.table(file = "../UnivariateAssociations_asin.txt",header = TRUE, sep = "\t", row.names = 1)
Final_phenotypes1 <- as.character(unique(Final_phenotypes1$Phenotype))

HealthyControls_meta <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "HealthyControl", select = c(Final_phenotypes1))
IBD_meta <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = c(Final_phenotypes1))
SmallIntestine_meta <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = c(Final_phenotypes1))

HealthyControls_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "HealthyControl", c(165:ncol(LLD.IBD_MetaSp15.asin)))
IBD_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))
SmallIntestine_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "IleostomaPouch" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))


write.table(HealthyControls_meta, file = "HealthyControls_meta.txt", sep = "\t", col.names = TRUE)
write.table(IBD_meta, file = "IBD_meta.txt", sep = "\t", col.names = TRUE)
write.table(SmallIntestine_meta, file = "SmallIntestine_meta.txt", sep = "\t", col.names = TRUE)

write.table(HealthyControls_taxa, file = "HealthyControls_taxa.txt", sep = "\t", col.names = TRUE)
write.table(IBD_taxa, file = "IBD_taxa.txt", sep = "\t", col.names = TRUE)
write.table(SmallIntestine_taxa, file = "SmallIntestine_taxa.txt", sep = "\t", col.names = TRUE)




#Run Adonis

#LLD
adonis_HC <- matrix(ncol = 4, nrow=ncol(HealthyControls_meta)) 
#HealthyControls_meta2 <- HealthyControls_meta[rownames(si_taxa),]
for ( i in 1:ncol(HealthyControls_meta)){
  HealthyControls_meta2=HealthyControls_meta[complete.cases(HealthyControls_meta[,i]),]
  adonis_HC[i,4]=nrow(HealthyControls_meta2)
  HealthyControls_taxa2=subset(HealthyControls_taxa, row.names(HealthyControls_taxa) %in% row.names(HealthyControls_meta2))
  if (length(unique(HealthyControls_meta2[,i]))!=1 & nrow(HealthyControls_meta2)>1){
    ad <- adonis(formula = HealthyControls_taxa2 ~ HealthyControls_meta2[,i] , data = HealthyControls_meta2, permutations = 1000, method = "bray")
    aov_table <- ad$aov.tab
    adonis_HC[i,1]=aov_table[1,1]
    adonis_HC[i,2]=aov_table[1,5]
    adonis_HC[i,3]=aov_table[1,6]
  } else{
    adonis_HC[i,1]="Only one level"
    adonis_HC[i,2]="Only one level"
    adonis_HC[i,3]="Only one level"
  }
}
colnames(adonis_HC) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_HC)=colnames(HealthyControls_meta)

adonis_HC <- as.data.frame(adonis_HC)
adonis_HC[,3] <- as.numeric(as.character(adonis_HC[,3]))

adonis_HC$BonferroniAdjust <- p.adjust(c(adonis_HC[,3]), method = 'bonferroni')
adonis_HC$FDRadjust <- p.adjust(c(adonis_HC[,3]), method = 'fdr')

write.table(adonis_HC, file = 'Adonis_HealthyControls1.txt', sep = "\t", col.names = T)




#IBD excl. ileostoma/ileoanal pouch
#IBD_meta1 <- IBD_meta #[c(1:335,337:483),] #exclude sample IBDFEC1024

adonis_IBD <- matrix(ncol = 4, nrow=ncol(IBD_meta)) 
#IBD_meta2 <- IBD_meta[rownames(si_taxa),]
for ( i in 1:ncol(IBD_meta)){
  IBD_meta2=IBD_meta[complete.cases(IBD_meta[,i]),]
  adonis_IBD[i,4]=nrow(IBD_meta2)
  IBD_taxa2=subset(IBD_taxa, row.names(IBD_taxa) %in% row.names(IBD_meta2))
  if (length(unique(IBD_meta2[,i]))!=1 & nrow(IBD_meta2)>1){
    ad <- adonis(formula = IBD_taxa2 ~ IBD_meta2[,i] , data = IBD_meta2, permutations = 1000, method = "bray")
    aov_table <- ad$aov.tab
    adonis_IBD[i,1]=aov_table[1,1]
    adonis_IBD[i,2]=aov_table[1,5]
    adonis_IBD[i,3]=aov_table[1,6]
  } else{
    adonis_IBD[i,1]="Only one level"
    adonis_IBD[i,2]="Only one level"
    adonis_IBD[i,3]="Only one level"
  }
}
colnames(adonis_IBD) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_IBD)=colnames(IBD_meta)

adonis_IBD <- as.data.frame(adonis_IBD)
adonis_IBD[,3] <- as.numeric(as.character(adonis_IBD[,3]))

adonis_IBD$BonferroniAdjust <- p.adjust(c(adonis_IBD[,3]), method = 'bonferroni')
adonis_IBD$FDRadjust <- p.adjust(c(adonis_IBD[,3]), method = 'fdr')

write.table(adonis_IBD, file = 'Adonis_IBD_excl.ileostoma&pouch1.txt', sep = "\t", col.names = T)



#__Small Intestine (leostoma/ileoanal pouch)

adonis_SI <- matrix(ncol = 4, nrow=ncol(SmallIntestine_meta)) 
#SmallIntestine_meta2 <- SmallIntestine_meta[rownames(si_taxa),]
for ( i in 1:ncol(SmallIntestine_meta)){
  SmallIntestine_meta2=SmallIntestine_meta[complete.cases(SmallIntestine_meta[,i]),]
  adonis_SI[i,4]=nrow(SmallIntestine_meta2)
  SmallIntestine_taxa2=subset(SmallIntestine_taxa, row.names(SmallIntestine_taxa) %in% row.names(SmallIntestine_meta2))
  if (length(unique(SmallIntestine_meta2[,i]))!=1 & nrow(SmallIntestine_meta2)>1){
    ad <- adonis(formula = SmallIntestine_taxa2 ~ SmallIntestine_meta2[,i] , data = SmallIntestine_meta2, permutations = 1000, method = "bray")
    aov_table <- ad$aov.tab
    adonis_SI[i,1]=aov_table[1,1]
    adonis_SI[i,2]=aov_table[1,5]
    adonis_SI[i,3]=aov_table[1,6]
  } else{
    adonis_SI[i,1]="Only one level"
    adonis_SI[i,2]="Only one level"
    adonis_SI[i,3]="Only one level"
  }
}
colnames(adonis_SI) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_SI)=colnames(SmallIntestine_meta)

adonis_SI <- as.data.frame(adonis_SI)
adonis_SI[,3] <- as.numeric(as.character(adonis_SI[,3]))

adonis_SI$BonferroniAdjust <- p.adjust(c(adonis_SI[,3]), method = 'bonferroni')
adonis_SI$FDRadjust <- p.adjust(c(adonis_SI[,3]), method = 'fdr')

write.table(adonis_SI, file = 'Adonis_SmallIntestine1.txt', sep = "\t", col.names = T)



#Display results

#variables not applicable to/not differing in SI removed
Outcome_function = RemoveColumnsNA(SmallIntestine_meta, Threshold = 0.92) 
Features_removed = Outcome_function[[2]]

adonis_HC_sub <- adonis_HC[-c(22,23,45,46,47,49,50,51,54,55,108:113),]
adonis_IBD_sub <- adonis_IBD[-c(22,23,45,46,47,49,50,51,54,55,108:113),]
adonis_SI_sub <- adonis_SI[-c(22,23,45,46,47,49,50,51,54,55,108:113),]

adonis_HC_sub$FDRadjust <- p.adjust(c(adonis_HC_sub[,3]), method = 'fdr')
adonis_IBD_sub$FDRadjust <- p.adjust(c(adonis_IBD_sub[,3]), method = 'fdr')
adonis_SI_sub$FDRadjust <- p.adjust(c(adonis_SI_sub[,3]), method = 'fdr')

adonis_HC_sub$BonferroniAdjust <- p.adjust(c(adonis_HC_sub[,3]), method = 'bonferroni')
adonis_IBD_sub$BonferroniAdjust <- p.adjust(c(adonis_IBD_sub[,3]), method = 'bonferroni')
adonis_SI_sub$BonferroniAdjust <- p.adjust(c(adonis_SI_sub[,3]), method = 'bonferroni')

adonis_HC_sig <- subset(adonis_HC_sub, `Pr(>F)` < 0.05)
adonis_IBD_sig <- subset(adonis_IBD_sub, `Pr(>F)` < 0.05)
adonis_SI_sig <- subset(adonis_SI_sub, `Pr(>F)` < 0.05)

R2_HC <- sort(adonis_HC_sig$R2, decreasing = T) #order
R2_IBD <- sort(adonis_IBD_sig$R2,decreasing = T)
R2_SI <- sort(adonis_SI_sig$R2,decreasing = T)

adonis_HC_Top_pval <- adonis_HC_sub[unique(c(TopNTaxa(R2_HC,10),TopNTaxa(R2_IBD,10),TopNTaxa(R2_SI,length(R2_SI)))),]
adonis_IBD_Top_pval <- adonis_IBD_sub[unique(c(TopNTaxa(R2_HC,10),TopNTaxa(R2_IBD,10),TopNTaxa(R2_SI,length(R2_SI)))),]
adonis_SI_Top_pval <- adonis_SI_sub[unique(c(TopNTaxa(R2_HC,10),TopNTaxa(R2_IBD,10),TopNTaxa(R2_SI,length(R2_SI)))),]

adonis_HC_Top_pval$R2 <- as.numeric(as.character(adonis_HC_Top_pval$R2))*100
adonis_IBD_Top_pval$R2 <- as.numeric(as.character(adonis_IBD_Top_pval$R2))*100
adonis_SI_Top_pval$R2 <- as.numeric(as.character(adonis_SI_Top_pval$R2))*100

adonis_HC_Top_pval$ID <- rownames(adonis_HC_Top_pval)
adonis_IBD_Top_pval$ID <- rownames(adonis_IBD_Top_pval)
adonis_SI_Top_pval$ID <- rownames(adonis_SI_Top_pval)

adonis_HC_Top_pval <- adonis_HC_Top_pval[complete.cases(adonis_HC_Top_pval),] 
adonis_IBD_Top_pval <- adonis_IBD_Top_pval[complete.cases(adonis_IBD_Top_pval),]
adonis_SI_Top_pval <- adonis_SI_Top_pval[complete.cases(adonis_SI_Top_pval),]


ggplot(adonis_HC_Top_pval, aes(x=reorder(ID, -R2), y=R2, fill = FDRadjust < 0.05)) + labs(x = "Phenotype", y= "Variance Explained") + geom_bar(stat="identity",color = "black")+scale_fill_manual(values=c("white","grey24")) + theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 55, hjust = 1, size=5))+ylim(0,8)
ggplot(adonis_IBD_Top_pval, aes(x=reorder(ID, -R2), y=R2, fill = FDRadjust < 0.05)) + labs(x = "Phenotype", y= "Variance Explained") + geom_bar(stat="identity",color = "black") +scale_fill_manual(values=c("white","grey24")) + theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 55, hjust = 1, size=5))+ylim(0,8)
ggplot(adonis_SI_Top_pval, aes(x=reorder(ID, -R2), y=R2, fill = FDRadjust < 0.05)) + labs(x = "Phenotype", y= "Variance Explained") + geom_bar(stat="identity",color = "black") +scale_fill_manual(values=c("white","grey24")) + theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 55, hjust = 1, size=5))+ylim(0,8)

ggplot(adonis_HC_Top_pval, aes(x=reorder(ID, -R2), y=R2,label = ifelse(FDRadjust < 0.05, "*", ""))) + labs(x = "Phenotype", y= "Variance Explained") + geom_bar(stat="identity", colour = "black", fill = "grey50") + geom_text(vjust = 0) + theme_classic()+theme(legend.position="none",axis.text.x = element_text(angle = 55, hjust = 1, size=5))+ylim(0,8)


write.table(adonis_HC_sub, file = 'Adonis_HC_1.txt', sep = "\t", col.names = T)
write.table(adonis_IBD_sub, file = 'Adonis_IBD_1.txt', sep = "\t", col.names = T)
write.table(adonis_SI_sub, file = 'Adonis_SI_1.txt', sep = "\t", col.names = T)

write.table(adonis_HC_Top_pval, file = 'Adonis_HC_TopPheno1.txt', sep = "\t", col.names = T)
write.table(adonis_IBD_Top_pval, file = 'Adonis_IBD_TopPheno1.txt', sep = "\t", col.names = T)
write.table(adonis_SI_Top_pval, file = 'Adonis_SI_TopPheno1.txt', sep = "\t", col.names = T)



#____ SI versus Rest (IBD & LLD faecal)

SIvsColon <- LLD.IBD_MetSp_Below15[c(76,1,3,165:ncol(LLD.IBD_MetSp_Below15))]
SIvsColon$CurrentStomaOrPouchType <- fct_collapse(SIvsColon$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))

SI_NonSI <- SIvsColon[c(1)]
SI_NonSI_taxa <- LLD.IBD_MetaSp15.asin[c(165:ncol(LLD.IBD_MetaSp15.asin))]

adonis_All <- matrix(ncol = 4, nrow=ncol(SI_NonSI)) 
  adonis_All[1,4]=nrow(SI_NonSI)
    ad <- adonis(formula = SI_NonSI_taxa ~ SI_NonSI[,1] , data = SI_NonSI, permutations = 1000, method = "bray")
    aov_table <- ad$aov.tab
    adonis_All[1,1]=aov_table[1,1]
    adonis_All[1,2]=aov_table[1,5]
    adonis_All[1,3]=aov_table[1,6]

colnames(adonis_All) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_All)=colnames(SI_NonSI)

write.table(adonis_All, file = 'Adonis_SI_NonSI.txt', sep = "\t", col.names = T)


#__ Resection vs Non-Resec
NoResec_Resec <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections", select = "CurrentStomaOrPouchType")
NoResec_Resec_taxa <- subset(LLD.IBD_MetaSp15.asin, CurrentStomaOrPouchType == "Resections" | CurrentStomaOrPouchType == "No Resections" , select = c(165:ncol(LLD.IBD_MetaSp15.asin)))

adonis_Resec <- matrix(ncol = 4, nrow=ncol(NoResec_Resec)) 
adonis_Resec[1,4]=nrow(NoResec_Resec)
ad <- adonis(formula = NoResec_Resec_taxa ~ NoResec_Resec[,1] , data = NoResec_Resec, permutations = 1000, method = "bray")
aov_table <- ad$aov.tab
adonis_Resec[1,1]=aov_table[1,1]
adonis_Resec[1,2]=aov_table[1,5]
adonis_Resec[1,3]=aov_table[1,6]

colnames(adonis_Resec) <- c("Df", "R2", "Pr(>F)", "Num.samples")
rownames(adonis_Resec)=colnames(NoResec_Resec)

write.table(adonis_Resec, file = 'Adonis_Resec_NoResections.txt', sep = "\t", col.names = T)

```

### 9. Microbial Pathways - Data Import

```{r}

#Packages required:
library(stringr)


IBD_pathways <- read.table(file = "../Data/SI_project/PathwayData/IBD_humann2_pathways_uniref90_082017.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) # Humann2 Pathways
LLD_pathways <- read.table(file = "../Data/SI_project/PathwayData/LLD_humann2_Uniref90_092017.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) # Humann2 Pathways
IBD_ID.Key <- read.table(file = "../Data/SI_project/PathwayData/rename_IBD.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) # Humann2 Pathways
LLD_ID.Key <- read.table(file = "../Data/SI_project/PathwayData/rename_LLD.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) # Humann2 Pathways

#Remove 'X' at the beginning of ID numbers
#colnames(IBD_pathways) <- gsub("^X", "", colnames(IBD_pathways))
colnames(LLD_pathways) <- gsub("^X", "", colnames(LLD_pathways))

#Remove '_kneaddata_merged_Abundance' at the end of ID numbers
colnames(IBD_pathways) <- gsub("_kneaddata_merged_Abundance", "", colnames(IBD_pathways))
colnames(LLD_pathways) <- gsub("_kneaddata_merged_Abundance", "", colnames(LLD_pathways))

#Remove duplicate pathways
LLD_pathways <- LLD_pathways[!grepl('\\|',rownames(LLD_pathways),ignore.case = T),]
IBD_pathways <- IBD_pathways[!grepl('\\|',rownames(IBD_pathways),ignore.case = T),]


#Translate IDs
#IBD
IBD_pathwaysT <- as.data.frame(t(IBD_pathways))
IBD_pathwaysIDmerge <- merge(IBD_ID.Key, IBD_pathwaysT, by= 'row.names', all = T)

IBD_pathwaysIDmerge$Row.names <- as.character(IBD_pathwaysIDmerge$Row.names)
IBD_pathwaysIDmerge$Classic <- as.character(IBD_pathwaysIDmerge$Classic)

IBD_pathwaysIDmerge$Classic[107:545] <- IBD_pathwaysIDmerge$Row.names[107:545]
rownames(IBD_pathwaysIDmerge) <- IBD_pathwaysIDmerge$Classic
IBD_pathwaysIDmerge$Classic <- NULL
IBD_pathwaysIDmerge$Row.names <- NULL
IBD_pathways_Final <- as.data.frame(t(IBD_pathwaysIDmerge))

#LLD
LLD_pathwaysT <- as.data.frame(t(LLD_pathways))
LLD_pathwaysIDmerge <- merge(LLD_ID.Key, LLD_pathwaysT, by= 'row.names')

LLD_pathwaysIDmerge$Row.names <- as.character(LLD_pathwaysIDmerge$Row.names)
LLD_pathwaysIDmerge$Classic <- as.character(LLD_pathwaysIDmerge$Classic)

LLD_pathwaysIDmerge$Classic[107:545] <- LLD_pathwaysIDmerge$Row.names[107:545]

rownames(LLD_pathwaysIDmerge) <- LLD_pathwaysIDmerge$ID.1
LLD_pathwaysIDmerge$ID.1 <- NULL
LLD_pathwaysIDmerge$Row.names <- NULL
LLD_pathways_Final <- as.data.frame(t(LLD_pathwaysIDmerge))

#Merge both cohorts
IBD_LLD_pathways <- merge(LLD_pathways_Final, IBD_pathways_Final, by = 'row.names', all = T)
rownames(IBD_LLD_pathways) <- IBD_LLD_pathways$Row.names
IBD_LLD_pathways$Row.names <- NULL
IBD_LLD_pathways <- as.data.frame(t(IBD_LLD_pathways))

IBD_LLD_pathways_RelAbun <- as.data.frame(t(apply(IBD_LLD_pathways, 1, function(x)
  {x/sum(x, na.rm = T)
})))

# log transformation

#IBD_LLD_pathways_log10 <- as.data.frame(t(apply(IBD_LLD_pathways_RelAbun, 1, function(x){
#  -log10(x)
#})))      - produces inf. value

IBD_LLD_pathways_log <- normalize(as.data.frame(t(IBD_LLD_pathways_RelAbun)), to.abundance = F)
IBD_LLD_pathways_log <- as.data.frame(t(IBD_LLD_pathways_log))

#Replace special characters in column name
colnames(IBD_LLD_pathways_log) <- str_replace_all(colnames(IBD_LLD_pathways_log), "[^[:alnum:]]", "_")


#_Metadata
#metadata (LLD and IBD faecal ID's combined in one column to allow for merge (see following steps))
LLD.IBD_Meta <- read.table(file = "../Data/SI_project/Metadata16032020.txt", header = TRUE, sep = "\t", quote = "\\") #Metadata

#To distinguish patients with no stoma and no resections from patients with resections and no stoma
LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="yes", gsub('none', 'Resection NoStoma', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta <- LLD.IBD_Meta %>%
  mutate(CurrentStomaOrPouchType = ifelse(ResectionAny=="no", gsub('none', 'No Resections', CurrentStomaOrPouchType), as.character(CurrentStomaOrPouchType)))

LLD.IBD_Meta$CurrentStomaOrPouchType <- as.factor(LLD.IBD_Meta$CurrentStomaOrPouchType)

#Add additional level to 'CurrentStomaOrPouchType' column to represent LLD individuals 
levels <- levels(LLD.IBD_Meta$CurrentStomaOrPouchType)
levels[length(levels) + 1] <- "HealthyControl"

# refactor CurrentStomaOrPouchType to include "HealthyControl" as a factor level
# and replace NA with "HealthyControl"
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = levels)
LLD.IBD_Meta$CurrentStomaOrPouchType[is.na(LLD.IBD_Meta$CurrentStomaOrPouchType)] <- "HealthyControl"

#Change order of the levels for variable CurrentStomaOrPouchType
LLD.IBD_Meta$CurrentStomaOrPouchType <- factor(LLD.IBD_Meta$CurrentStomaOrPouchType, levels = c('HealthyControl', 'No Resections', 'Resection NoStoma', 'colostoma', 'pouch', 'ileostoma'))

#Remove columns with 'Date' in the column names
LLD.IBD_Meta1 <- LLD.IBD_Meta[,!grepl('Date', colnames(LLD.IBD_Meta),ignore.case = T)]

#Remove rows with unknown ID
LLD.IBD_Meta2 <- subset(LLD.IBD_Meta1, is.na(CombinedIDs) == 'FALSE')

#Set IDs as rownames
rownames(LLD.IBD_Meta2) <- LLD.IBD_Meta2$CombinedIDs
LLD.IBD_Meta2$CombinedIDs <- NULL

#Taxonomy
LLD_Taxa <- read.table(file = "../Data/SI_project/LLD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #LLD Taxonomy 
rownames(LLD_Taxa) <- LLD_Taxa$ID
LLD_Taxa$ID <- NULL

LLD <- gsub("^X", "", colnames(LLD_Taxa)) #Remove 'X' at the beginning of ID numbers
colnames(LLD_Taxa) <- LLD

IBD_Taxa <- read.table(file = "../Data/SI_project/IBD_taxonomy_unstrat_clean.txt", header = TRUE, sep = "\t", quote = "\\") #IBD Taxonomy 
rownames(IBD_Taxa) <- IBD_Taxa$ID
IBD_Taxa$ID <- NULL

#Merge the cohort taxa's
LLD_IBD_Taxa <- merge(LLD_Taxa, IBD_Taxa, by = 'row.names', all = T)
rownames(LLD_IBD_Taxa) <- LLD_IBD_Taxa$Row.names
LLD_IBD_Taxa$Row.names <- NULL
LLD_IBD_Taxa <- as.data.frame(t(LLD_IBD_Taxa))

#Extract Species only
LLD_IBD_TaxaS <- LLD_IBD_Taxa[,grepl('s__',colnames(LLD_IBD_Taxa),ignore.case = T)]
LLD_IBD_TaxaS <- LLD_IBD_TaxaS[,!grepl('t__',colnames(LLD_IBD_TaxaS),ignore.case = T)]
LLD_IBD_TaxaS[is.na(LLD_IBD_TaxaS)] <- 0

#Merge Species Taxo with stoma/pouch type and ID's
LLD_IBD_MetTaxS <- merge(LLD.IBD_Meta2, LLD_IBD_TaxaS, by = 'row.names')

x <- gsub('.*\\|s__','',colnames(LLD_IBD_MetTaxS)) #Keep taxonomy level name only i.e. 'class name' 
colnames(LLD_IBD_MetTaxS) <- x

rownames(LLD_IBD_MetTaxS) <- LLD_IBD_MetTaxS$Row.names
LLD_IBD_MetTaxS$Row.names <- NULL


write.table(LLD_IBD_MetTaxS, file = 'LLD_IBD_MetaSpecies150420.txt', sep = "\t", col.names = TRUE)
LLD_IBD_MetTaxS <- read.table(file = "../RScripts/LLD_IBD_MetaSpecies150420.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) #contains updated medicine  


#Merge MetaData with pathways (same individuals as Taxonomy data)
IBD_LLD_MetaPathways <- merge(LLD_IBD_MetTaxS[,1:680], IBD_LLD_pathways_log, by = 'row.names')  
rownames(IBD_LLD_MetaPathways) <- IBD_LLD_MetaPathways$Row.names
IBD_LLD_MetaPathways$Row.names <- NULL


write.table(IBD_LLD_MetaPathways, file = "IBD_LLD_MetaPathways.txt", sep = "\t", col.names = TRUE)

```

### 10. Microbial Pathways - Univariate analyses

```{r}
#Packages & Functions required:
univariate.analysis2 #see 'Functions' file, point 7.

IBD_LLD_MetaPathways <- read.table(file = "IBD_LLD_MetaPathways.txt", header = TRUE, sep = "\t", quote = "\\",row.names = 1) #see PathwayMetaDataIBD&LLD script for details

colnames(IBD_LLD_MetaPathways) <- gsub("^X", "", colnames(IBD_LLD_MetaPathways))

#remove individual with no pathway abundances ...
which(rowSums(IBD_LLD_MetaPathways[681:1168]) == 0)
IBD_LLD_MetaPathways <- IBD_LLD_MetaPathways[-c(212,218,379,1554,1671),]

#Merge with the metadata included in the taxa analyses
IBD_LLD_MetaPathwaysx <- merge(LLD.IBD_MetaSp15[c(1:164)], IBD_LLD_MetaPathways[c(681:1168)], by = "row.names")
rownames(IBD_LLD_MetaPathwaysx) <- IBD_LLD_MetaPathwaysx$Row.names
IBD_LLD_MetaPathwaysx$Row.names <- NULL


# Remove pathways from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 165:ncol(IBD_LLD_MetaPathwaysx)) {
  cname <- colnames(IBD_LLD_MetaPathwaysx)[i]
  if(colSums(IBD_LLD_MetaPathwaysx[i] != 0) / nrow(IBD_LLD_MetaPathwaysx) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD.IBD_MetaPwy15 <- IBD_LLD_MetaPathwaysx %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

LLD.IBD_MetaPwy15$NumberIndicationabcessOrfistula <- as.numeric(as.character(LLD.IBD_MetaPwy15$NumberIndicationabcessOrfistula))

LLD.IBD.UniA_Pathways <- univariate.analysis2(LLD.IBD_MetaPwy15,164,165,341) 

LLD.IBD.UniA_Pathways$PrevalenceGroup1 <- as.numeric(as.character(LLD.IBD.UniA_Pathways$PrevalenceGroup1)) 
LLD.IBD.UniA_Pathways$PrevalenceGroup2 <- as.numeric(as.character(LLD.IBD.UniA_Pathways$PrevalenceGroup2))
LLD.IBD.UniA_Pathways$`NA` <- as.numeric(as.character(LLD.IBD.UniA_Pathways$`NA`))
LLD.IBD.UniA_Pathways$`NSamplesG1/Nnonzero` <- as.numeric(as.character(LLD.IBD.UniA_Pathways$`NSamplesG1/Nnonzero`))
LLD.IBD.UniA_Pathways$`NSampleG2/Nzero` <- as.numeric(as.character(LLD.IBD.UniA_Pathways$`NSampleG2/Nzero`))

LLD.IBD.UniA_Pathways1 <- subset(LLD.IBD.UniA_Pathways, `NA` < 1720*0.95) #remove rows with NA > 95%
LLD.IBD.UniA_Pathways2 <- LLD.IBD.UniA_Pathways1


#Filter out tests with a group sample number < 20
remove.rows <- function(x){ # to complete function replace LLD.IBD.UniA2 with 'x' 
  remove_rows <- vector()
  for (i in 1:nrow(LLD.IBD.UniA_Pathways2)) {
    names <- row.names(LLD.IBD.UniA_Pathways2)[i] 
    if(LLD.IBD.UniA_Pathways2[i,5] == 'Category' & LLD.IBD.UniA_Pathways2[i,7] < 20 | LLD.IBD.UniA_Pathways2[i,5] == 'Category' & LLD.IBD.UniA_Pathways2[i,8] < 20){
      remove_rows <- c(remove_rows, names)
    }
  }
  x_new <- subset(LLD.IBD.UniA_Pathways2, !(rownames(LLD.IBD.UniA_Pathways2) %in% remove_rows))
  print(paste(c("Function removed a total of ",length(remove_rows), "variables"), collapse= ""))
  return(x_new)
}

write.table(x_new, file = 'Pathways_UnivariateAnalysis.txt', sep = "\t", col.names = T)

### Top 25 phenotypes with most bacterial associations
Significant_pwy <- subset(x_new, BelowAdjustPval == 'yes') #subset associations with p-values below adjusted p-value

Significant_pwy <- transform(Significant_pwy, Combinations=paste(Phenotype, Group1, Group2, sep=":")) #merge Phenotype, Group1, Group2 variables into a new variable

colnames(Significant_pwy)[1] <- "Pathway"
 
sorted_pwy_TopPhe <- Significant_pwy %>% #create dataframe with 'combinations' variable in descending order of frequency
  count(Combinations) %>%
  arrange(-n) %>% 
  mutate(Combinations = factor(Combinations, levels = unique(Combinations)))

#plot top 25 most frequent variables 
ggplot(sorted_pwy_TopPhe[1:25,], aes(x = Combinations, y = n)) +
  geom_bar(stat="identity", fill = "steelblue") + xlab('Phenotype') + geom_text(aes(label=n), vjust=-0.3, size=2) + 
  ylab('Significant Associations (n)') + theme_classic() + theme(axis.text.x = element_text(angle = 55, hjust = 1, size=5)) +
  scale_x_discrete(labels=c('CurrentStomaOrPouchType:HealthyControl:IleostomaPouch' = 'Bowel Phenotype\n (LLD General Population vs IBD Small Intestine)',
                            'CurrentStomaOrPouchType:IleostomaPouch:No Resections' = 'Bowel Phenotype\n (IBD Small Intestine vs IBD Intact Colon)',
                            'EverHadStomaOrPouch:no:yes' = 'Ever Had A Stoma Or Pouch\n (No/Yes)',
                            'DiagnosisCurrent:CD:generalpopulation' = 'Current IBD Diagnosis\n (CD vs General Population)',
                            'NumberOfResectionColonic:NumberOfResectionColonic:NumberOfResectionColonic' = 'Number Of Colonic Resections',
                            'ResectionColonicAny:no:yes' = 'Colonic Resection Any\n (No/Yes)',
                            'CurrentStomaOrPouchType:HealthyControl:Resections' = 'Bowel Phenotype\n (LLD General Population vs IBD Partially Intact Colon)',
                            'IleocecalValveInSitu:no:yes' = 'Ileocecal Valve In Situ\n (No/Yes)',
                            'ResectionAny:no:yes' = 'Resection Any\n (No/Yes)',
                            'NumberOfResectionsAny:NumberOfResectionsAny:NumberOfResectionsAny' = 'Total Number Of Resections',
                            'vitamin_B12:no:yes' = 'Vitamin B12\n (No/Yes)', 
                            'PPI:no:yes' = 'PPI\n (No/Yes)',
                            'CurrentStomaOrPouchType:IleostomaPouch:Resections' ='Bowel Phenotype\n (IBD Small Intestine vs IBD Partially Intact Colon)',
                            'vitamin_D:no:yes' = 'Vitamin D\n (No/Yes)',
                            'FecalCalprotectinOver200yesno:no:yes' = 'Fecal Calprotectin Over 200\n (No/Yes)',
                            'Alc_en:Alc_en:Alc_en' = 'Alcohol\n (% Total Energy Intake)',
                            'calcium:no:yes' = 'Calcium\n (No/Yes)',
                            'CurrentStomaOrPouchType:No Resections:Resections' = 'Bowel Phenotype\n (IBD Intact Colon vs IBD Partially Intact Colon)',
                            'ResectionIlealCecalAny:no:yes' = 'IlealCecal Resection Any\n (No/Yes)',
                            'NumberOfResetionsIleoCecal:NumberOfResetionsIleoCecal:NumberOfResetionsIleoCecal' = 'Number Of IleoCecal Resections',
                            'DiagnosisCurrent:generalpopulation:UC' = 'Current IBD Diagnosis\n (General Population vs UC)',
                            'PFReads:PFReads:PFReads' = 'Reads Depth',
                            'FrequencyFormedStoolDay:FrequencyFormedStoolDay:FrequencyFormedStoolDay' = 'Frequency Stool Formed Per Day',
                            'oral_steroid:no:yes'= 'Oral Steroids\n (No/Yes)',
                            'CurrentStomaOrPouchType:HealthyControl:No Resections' = 'Bowel Phenotype\n (LLD General Population vs IBD Intact Colon)'))
                   


## Top 25 species according to number of significant phenotype associations  
sorted_TopPwy <- Significant_pwy %>% #create dataframe with 'Pathway' variable in descending order of frequency
  count(Pathway) %>%
  arrange(-n) %>% 
  mutate(Pathway = factor(Pathway, levels = unique(Pathway)))

sorted_TopPwy$Pathway <- as.character(gsub("_", " ", sorted_TopPwy$Pathway))

ggplot(sorted_TopPwy[1:25,], aes(x = reorder(Pathway, -n), y = n)) +
  geom_bar(stat="identity", fill = "steelblue") +  theme_classic() + geom_text(aes(label=n), vjust=-0.3, size=2) +
  theme(axis.text.x = element_text(angle = 55, hjust = 1, size=5)) + ylab('Significant Associations (n)') + xlab('Pathways') #plot top 25 most frequent variables


write.table(Significant_pwy, file = 'Significant_pwy15_associations3.txt', sep = "\t", col.names = T)

```

### 11. Microbial pathways - Multivariate analysis

```{r}

#Packages required:
library(Maaslin) 
library(forcats)


#***IleostomaPouch vs Healthy Controls***

#Prepare input data for Maaslin 

# 1. Metadata & Pathways
# Dataframe needed: IBD_LLD_MetaPathwaysx (i.e. non-transformed) (see rscript 'Univariate_Analysis_Pathways' for dataframe)

Maaslin_MetaPwy1 <- subset(LLD.IBD_MetaPwy15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'HealthyControl')
Maaslin_MetaPwy1 <- Maaslin_MetaPwy1[c(1,2,3,9,14,74,76,96,98,99,102,113,115,119,120,121,128,133,135,138,147,148,163,165:505)]
Maaslin_MetaPwy1$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaPwy1$CurrentStomaOrPouchType)


#for (i in 24:157) {
#  Maaslin_MetaPwy1[,i] <- Maaslin_MetaPwy1[,i]/100
#}     # Relative abundance should be between 0 - 1


write.table(Maaslin_MetaPwy1, file = "MaaslinMetaPwy_HealthControl_vs_IleostomaPouch.tsv", sep = "\t", quote = F)
Maaslin(strInputTSV = "MaaslinMetaPwy_HealthControl_vs_IleostomaPouch.tsv", strOutputDIR = "MaaslinPwy_HC_vs_Ileostoma_Output", strInputConfig = "MaaslinPwyConfig_HC_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0, strTransform = 'none')




#***Within IBD __ IleostomaPouch vs No Resection***

# 1. Metadata & Taxa
Maaslin_MetaPwy2 <- subset(LLD.IBD_MetaPwy15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'No Resections')
Maaslin_MetaPwy2 <- Maaslin_MetaPwy2[c(1,2,3,9,14,15,74,76,96,98,99,102,113,115,119,120,121,127,128,133,135,138,147,148,163,165:505)]
Maaslin_MetaPwy2$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaPwy2$CurrentStomaOrPouchType)

Maaslin_MetaPwy2$CurrentStomaOrPouchType <- relevel(Maaslin_MetaPwy2$CurrentStomaOrPouchType, ref = 'No Resections') #Change Reference (Baseline) Category

Maaslin_MetaPwy2$DiagnosisCurrent <- fct_drop(Maaslin_MetaPwy2$DiagnosisCurrent)

write.table(Maaslin_MetaPwy2, file = "MaaslinMetaPwy_WithinIBD_NoResections_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinPwyConfig_NoResection_vs_IleostomaPouch.read.config.R

Maaslin(strInputTSV = "MaaslinMetaPwy_WithinIBD_NoResections_vs_IleostomaPouch.tsv", strOutputDIR = "MaaslinPwy_NoResection_vs_Ileostoma_Output", strInputConfig = "MaaslinPwyConfig_NoResection_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0, strTransform = 'none')


#***Within IBD __ IleostomaPouch vs Resections ***

# 1. Metadata & Taxa
Maaslin_MetaPwy3 <- subset(LLD.IBD_MetaPwy15, CurrentStomaOrPouchType == 'IleostomaPouch' | CurrentStomaOrPouchType == 'Resections')
Maaslin_MetaPwy3 <- Maaslin_MetaPwy3[c(1,2,3,9,14,15,74,76,79,82,96,98,99,102,113,115,119,120,121,127,128,133,135,138,147,148,163,165:505)]
Maaslin_MetaPwy3$CurrentStomaOrPouchType <- fct_drop(Maaslin_MetaPwy3$CurrentStomaOrPouchType)


#Remove levels with 0 samples
Maaslin_MetaPwy3$DiagnosisCurrent <- fct_drop(Maaslin_MetaPwy3$DiagnosisCurrent)


write.table(Maaslin_MetaPwy3, file = "MaaslinMetaPwy_WithinIBD_Resections_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinPwyConfig_ResectionNoStoma_IleostomaPouch.read.config.R

Maaslin(strInputTSV = "MaaslinMetaPwy_WithinIBD_Resections_vs_IleostomaPouch.tsv", strOutputDIR = "MaaslinPwy_Resection_vs_Ileostoma_Output", strInputConfig = "MaaslinPwyConfig_ResectionNoStoma_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0, strTransform = 'none')



#***Within IBD __ IleostomaPouch vs ResectionsIleoCecalValveInSitu***

# 1. Metadata & Taxa
Maaslin_MetaPwy3a <- subset(Maaslin_MetaPwy3, CurrentStomaOrPouchType == 'Resections' & IleocecalValveInSitu == 'yes')
Maaslin_MetaPwy3b <- rbind(subset(Maaslin_MetaPwy3, CurrentStomaOrPouchType == 'IleostomaPouch'), Maaslin_MetaPwy3a)
Maaslin_MetaPwy3b$IleocecalValveInSitu <- NULL

write.table(Maaslin_MetaPwy3b, file = "MaaslinMetaPwy_WithinIBD_ResectionsWithIleoCecalValveInSitu_vs_IleostomaPouch.tsv", sep = "\t", quote = F)

#2. Config file
# Using R, name: MaaslinPwyConfig_ResectionNoStoma_IleostomaPouch.read.config.R


Maaslin(strInputTSV = "MaaslinMetaPwy_WithinIBD_ResectionsWithIleoCecalValveInSitu&NoStoma_vs_IleostomaPouch.tsv", strOutputDIR = "MaaslinPwy_ResectionsCecalValveInSitu_vs_Ileostoma_Output", strInputConfig = "MaaslinPwyConfig_HC_vs_IleostomaPouch.read.config.R", dMinAbd = 0, dMinSamp = 0, strTransform = 'none')


================================================================================================================================================
REVISED Manuscript - 17.12.20
================================================================================================================================================

#PACKAGES & FUNCTIONS

install.packages("pheatmap")
install.packages("dendextend")
install.packages("circlize")
install.packages("ape")

# load packages
library(pheatmap)
library(dendextend)
library(circlize)
library(ape)


#IMPORT & PROCESS DATA
IBD_LLD_MetaPathways <- read.table(file = "IBD_LLD_MetaPathways.txt", header = TRUE, sep = "\t", quote = "\\",row.names = 1) #see PathwayMetaDataIBD&LLD script for details

colnames(IBD_LLD_MetaPathways) <- gsub("^X", "", colnames(IBD_LLD_MetaPathways))

#remove individual with no pathway abundances ...
which(rowSums(IBD_LLD_MetaPathways[681:1168]) == 0)
IBD_LLD_MetaPathways <- IBD_LLD_MetaPathways[-c(212,218,379,1554,1671),]

#Merge with the metadata included in the taxa analyses
IBD_LLD_MetaPathwaysx <- merge(LLD.IBD_MetaSp15[c(1:164)], IBD_LLD_MetaPathways[c(681:1168)], by = "row.names")
rownames(IBD_LLD_MetaPathwaysx) <- IBD_LLD_MetaPathwaysx$Row.names
IBD_LLD_MetaPathwaysx$Row.names <- NULL


# Remove pathways from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 165:ncol(IBD_LLD_MetaPathwaysx)) {
  cname <- colnames(IBD_LLD_MetaPathwaysx)[i]
  if(colSums(IBD_LLD_MetaPathwaysx[i] != 0) / nrow(IBD_LLD_MetaPathwaysx) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}
LLD.IBD_MetaPwy15 <- IBD_LLD_MetaPathwaysx %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))

LLD.IBD_MetaPwy15$NumberIndicationabcessOrfistula <- as.numeric(as.character(LLD.IBD_MetaPwy15$NumberIndicationabcessOrfistula))

LLD.IBD_MetaPwy15$CurrentStomaOrPouchType <- factor(LLD.IBD_MetaPwy15$CurrentStomaOrPouchType, levels = c("HealthyControl", "No Resections", "Resections", "IleostomaPouch"))

**********************************
DENDROGRAM WITH HEATMAP (Figure 4)
**********************************

#Pathways only (Heatmap Rows)
Pathways <- LLD.IBD_MetaPwy15[,165:505]
colnames(Pathways) <- gsub("\\__.*", "",colnames(Pathways)) 
Pathways$UNINTEGRATED <- NULL
Pathways$UNMAPPED <- NULL
Pathways <- as.data.frame(t(Pathways))

#Column annotation (samples by bowel phenotype)                                            
my_sample_col <- LLD.IBD_MetaPwy15[,c(76:77)]
my_sample_col$EverHadStomaOrPouch <- NULL
colnames(my_sample_col)[1] <- "BowelPhenotype"
my_sample_col$a <- my_sample_col$BowelPhenotype
my_sample_col$a <- ifelse(my_sample_col$a != "IleostomaPouch", "Other", "IleostomaPouch")
my_sample_col[,c("b","c")] <- my_sample_col$BowelPhenotype
my_sample_col$b <- ifelse(my_sample_col$b != "Resections", "Other", "Resections")
my_sample_col$c <- ifelse(my_sample_col$c != "No Resections", "Other", "No Resections")

my_colour = list(
  a = c("IleostomaPouch" = 'red2', "Other" = "white"),
  b = c("Other" = "white","Resections" = 'gold1'),
  c = c("Other" = "white", "No Resections" = 'slateblue2'))

# Plot Heatmap
p <- pheatmap(Pathways,clustering_distance_rows = "euclidean", annotation_colors = my_colour, annotation_col = my_sample_col[,2:4], show_colnames = F, fontsize_row = 2, annotation_names_col = F)                                         

Pwy_Order_All <- as.data.frame(rownames(Pathways[p$tree_row[["order"]],]))
colnames(Pwy_Order_All)[1] <- "Pathway"

#Merge Pwy order with Maaslin results
Pwy_Order_All$Pathway <- as.character(Pwy_Order_All$Pathway)

HC_IBD.SI.Maaslin.pwy$Feature <- gsub("\\__.*", "",HC_IBD.SI.Maaslin.pwy$Feature)  #See Previous
merge1 <- merge(Pwy_Order_All, HC_IBD.SI.Maaslin.pwy[,c(2,9)], by.x = "Pathway", by.y = "Feature", all.x = T)

NoResec_SI.Maaslin.pwy$Feature <- gsub("\\__.*", "",NoResec_SI.Maaslin.pwy$Feature) 
merge2 <- merge(merge1, NoResec_SI.Maaslin.pwy[,c(2,9)], by.x = "Pathway", by.y = "Feature", all.x = T)

Resec_SI.Maaslin.pwy$Feature <- gsub("\\__.*", "",Resec_SI.Maaslin.pwy$Feature) 
merge3 <- merge(merge2, Resec_SI.Maaslin.pwy[,c(2,9)], by.x = "Pathway", by.y = "Feature", all.x = T)

target <- Pwy_Order_All$Pathway
Dendrogram_legends <- merge3[match(target, merge3$Pathway),]
Dendrogram_legends[,2:4][is.na(Dendrogram_legends[,2:4])] <- 0


## --- Plot dendrogram with the three comparison tests as a heatmap around the rim ---

col_fun = colorRamp2(breaks= c(-3,-2,-1,0,1,2,3), 
                     colors = c("#2171B5","#6BAED6","#C6DBEF","#FFFFFF","#FCBBA1", "#FB6A4A", "#CB181D")) 

circos.clear()
circos.par("start.degree" = 90,cell.padding = c(0, 0, 0, 0), gap.degree = 15) 
circos.initialize("a", xlim =c(0,339)) 
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.15, 
             panel.fun = function(x, y) {
               for(i in seq_len(nrow(Dendrogram_legends))) {
                 circos.text(i-0.5, 0, Dendrogram_legends$Pathway[i], adj = c(0, 0.5), 
                             facing = "clockwise", niceFacing = TRUE,
                             cex = 0.3)                
               }
             })

circos.track(ylim = c(0, 3), bg.border = NA,track.height = 0.1, panel.fun = function(x, y) {
  mm = Dendrogram_legends
  rownames(mm) <- mm$Pathway
  mm$Pathway <- NULL
  mm$Cluster <- NULL
  colnames(mm)[1:3] <- c("IBD-SI vs GP", "IBD-SI vs IBD-NoRes", "IBD-SI vs IBD-Res")
  mm2 <- as.matrix(t(mm))
  dend1 = as.dendrogram(p$tree_row)
  #changed variable for 1 heatmap setting
  col_mat = col_fun(mm2)
  nr = nrow(mm2)
  nc = ncol(mm2)
  for(i in 1:nr) {
    circos.rect(1:nc - 1, rep(nr - i, nc), 
                1:nc, rep(nr - i + 1, nc), 
                border = col_mat[i, ], col = col_mat[i, ])
  }
  #adding row label
  circos.text(rep(1, 3), 3:1, 
              rownames(mm2), 
              facing = "downward", adj = c(1.2, 1), cex = 0.3) 
  
})

#Dendrogram
dend1 = as.dendrogram(p$tree_row)
max_height = attr(dend1, "height") 
circos.track(ylim = c(0, max_height), bg.border = NA, track.height = 0.3, 
             panel.fun = function(x, y) {
               circos.dendrogram(dend1, max_height = max_height)
             })



```


### 12. Logistic Regression

```{r}

#Packages required:
library(forcats)
Remove_FactCols #see 'Functions' file, point 2.
LogisticRegression.function #see 'Functions' file, point 12.
LgRgr.DescriptiveStats #see 'Functions' file, point 13.


#See 'pathways' R script for origin 
LLD_IBD_MetTaxS <- read.table(file = "../RScripts/LLD_IBD_MetaSpecies150420.txt", header = TRUE, sep = "\t", quote = "\\", row.names = 1) # 
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

LLD.IBD_Final <- Remove_FactCols(LLD.IBD_MetTaxS.F) 
y <- lapply(LLD.IBD_Final[,1:428], class)

LLD.IBD_Final <- LLD.IBD_Final[c(1:114,142:1080)] #remove logical variables, not informative and cannot be applied to the 'univariate.analysis' function
y <- lapply(LLD.IBD_Final[,1:401], class)      


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

#Remove 'correlated/releted' variables
LLD.IBD_Finalxx <- LLD.IBD_Finalx[c(1,2,5:13,15,18,21:31,33,36:43,45,47,49,52:108,110,111,114:134,302,310:317,319,321,323:370,373,377,379,381:385,387,389,403:1054)]
LLD.IBD_Finalxx <- LLD.IBD_Finalxx[-c(74,76:78,80:93,96,98,113)]

# Remove bacteria from dataframe with a prevalence < 15%
remove_cols <- vector()
for (i in 165:ncol(LLD.IBD_Finalxx)) {
  cname <- colnames(LLD.IBD_Finalxx)[i]
  if(colSums(LLD.IBD_Finalxx[i] != 0) / nrow(LLD.IBD_Finalxx) *100 < 15){
    remove_cols <- c(remove_cols, cname) 
  }
}

LLD.IBD_Species_Below15 <- LLD.IBD_Finalxx %>% select(remove_cols)
LLD.IBD_Species_Below15 <- as.data.frame(ifelse(LLD.IBD_Species_Below15 > 0, 1, 0))

LLD.IBD_MetSp_Below15 <- cbind(LLD.IBD_Finalxx[,c(1:164)], LLD.IBD_Species_Below15)

#LLD.IBD_MetSp_Below15$CurrentStomaOrPouchType <- relevel(LLD.IBD_MetSp_Below15$CurrentStomaOrPouchType, ref = 'IleostomaPouch') #Change Reference (Baseline) Category



# --- Logistic regression with Ileostomy/Pouch vs Rest (i.e. SI vs Colon) ---
SIvsColon <- LLD.IBD_MetSp_Below15[c(76,1,3,165:ncol(LLD.IBD_MetSp_Below15))]
SIvsColon$CurrentStomaOrPouchType <- fct_collapse(SIvsColon$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))

my_preds=c("CurrentStomaOrPouchType","Sex","AgeAtFecalSampling")

LgR <- LogisticRegression.function(SIvsColon,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(SIvsColon, 3, 4)

LogisticRegrSIvsColon <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrSIvsColon, "LogisticRegression_SIvsColon1.txt", sep = "\t", col.names = T)


#Bowel Phenotype only
SIvsColon_only <- LLD.IBD_MetSp_Below15[c(76,165:ncol(LLD.IBD_MetSp_Below15))]
SIvsColon_only$CurrentStomaOrPouchType <- fct_collapse(SIvsColon$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))

my_preds=c("CurrentStomaOrPouchType")

LgR <- LogisticRegression.function(SIvsColon_only,1,2,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(SIvsColon_only, 1, 2)

LogisticRegrSIvsColon_only <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrSIvsColon_only, "LogisticRegression_SIvsColon2.txt", sep = "\t", col.names = T)

================================================================================================================================================
REVISED Manuscript - 17.12.20
================================================================================================================================================

#Packages/Functions required:

library(reshape2)
Prev  #see 'Functions' #see 'Functions' file, point 12.

***************************
Display Results (Figure 3b) 
***************************

LogRegr <- read.table("../RScripts/LogisticRegression_SIvsColon1.txt", header = TRUE, sep = "\t", row.names = 1)
LogRegr <- subset(LogRegr, Phenotype == "CurrentStomaOrPouchType")
LogRegr <- subset(LogRegr, FDRAdjust < 0.05)

SIvsColon_prev <- LLD.IBD_MetSp_Below15[c(76,165:ncol(LLD.IBD_MetSp_Below15))]
SIvsColon_prev$CurrentStomaOrPouchType <- fct_collapse(SIvsColon_prev$CurrentStomaOrPouchType, IleostomaPouch = c("IleostomaPouch"),Colon = c("Resections","HealthyControl","No Resections"))
SIvsColon_prev$CurrentStomaOrPouchType <- fct_drop(SIvsColon_prev$CurrentStomaOrPouchType)

a <- LogRegr$Taxonomy

SIvsColon_prev1 <- SIvsColon_prev[colnames(SIvsColon_prev) %in% a]
SIvsColon_prev1 <- cbind(SIvsColon_prev[,1], SIvsColon_prev1)
colnames(SIvsColon_prev1)[1] <- "CurrentStomaOrPouchType"

agg = aggregate(SIvsColon_prev1[2:111],
                by = list(SIvsColon_prev1$CurrentStomaOrPouchType),
                FUN = Prev)
                
melt <- melt(agg, 'Group.1')

melt1 <- melt %>% mutate(value = if_else(Group.1 == "Colon", -value, value))
melt1$variable <- gsub("_", " ", melt1$variable)

## find the order
temp_df <-
  melt1 %>% 
  filter(Group.1 == "IleostomaPouch") %>% 
  arrange(value)
the_order <- temp_df$variable

#Remove species with a prevalence below 10/15%
temp_df_15 <- subset(temp_df, value > 15)
the_order <- temp_df_15$variable

melt1 %>%  
  ggplot(aes(x = variable, y = value, group = Group.1, fill = Group.1)) +
  geom_bar(stat = "identity", width = 0.4) + #thinner bars
  coord_flip() + theme_classic()+ scale_x_discrete(limits = the_order) +
  scale_y_continuous(limits = c(-80, 80),breaks = seq(-100, 100, 10), 
                     labels = abs(seq(-100, 100, 10))) +
  labs(x = "Species", y = "Prevalence (%)") +  theme(legend.title = element_blank(),text = element_text(size=7)) + 
  scale_fill_manual(values=c("#808080", "#D72626"),labels=c("Colon", "Small Intestine"))


```

