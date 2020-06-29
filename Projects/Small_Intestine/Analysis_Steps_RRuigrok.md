Small Intestine Microbiota Project
==================================

Creator: Renate Ruigrok
Year: 2020

### 1. Cohort Descriptive Statistics

```{r}
# Functions Required:

library(forcats)
summary_statistics_metadata   #see 'Functions' folder, 1.
Remove_FactCols #see 'Functions' folder, 2.
kruskal.test_function. #see 'Functions' folder, 3.


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




