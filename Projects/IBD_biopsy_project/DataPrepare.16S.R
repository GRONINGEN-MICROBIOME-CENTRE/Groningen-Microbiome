# this script is for analysis
library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
library(RColorBrewer)
library(ggalluvial)
require(cowplot)
library(ggpubr)
library("dplyr")
library(edgeR)
library(lme4)
library(limma)
library(mixOmics)
library(robCompositions)
library(microbiome)
source("Microbiome.function.R")
library(lmerTest)

# ===========================================================================================
# prepare data
# ===========================================================================================

# import bacteria
genus=read.table("Decontamination_dada2/Genus.2000reads.txt",stringsAsFactors = F,header = T,sep = "\t",check.names = F)
family=read.table("Decontamination_dada2/Family.2000reads.txt",stringsAsFactors = F,header = T,sep = "\t",check.names = F)
phylum=read.table("Decontamination_dada2/Phylum.2000reads.txt",stringsAsFactors = F,header = T,sep = "\t",check.names = F)
orders=read.table("Decontamination_dada2/Order.2000reads.txt",stringsAsFactors = F,header = T,sep = "\t",check.names = F)
class=read.table("Decontamination_dada2/Class.2000reads.txt",stringsAsFactors = F,header = T,sep = "\t",check.names = F)
bacteria=rbind(phylum,orders,class,family,genus)

# import phenotype data
phenotype_bac=read.table("MetaData/Metadaata.16S.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# harmonize data, use old batch rnaseq data for overlapped samples (new batch has some failed)
# two new batch rnaseq samples have zero counts, so also removed
# 758 16S samples, 824 RNAseq samples
genus=genus[,colnames(genus) %in% phenotype_bac$`16SSampleID`,]
family=family[,colnames(family) %in% phenotype_bac$`16SSampleID`,]
phylum=phylum[,colnames(phylum) %in% phenotype_bac$`16SSampleID`,]
orders=orders[,colnames(orders) %in% phenotype_bac$`16SSampleID`,]
class=class[,colnames(class) %in% phenotype_bac$`16SSampleID`,]
bacteria=bacteria[,colnames(bacteria) %in% phenotype_bac$`16SSampleID`,]
phenotype_bac=phenotype_bac[phenotype_bac$`16SSampleID` %in% colnames(genus),]

genus=as.data.frame(t(genus))
family=as.data.frame(t(family))
phylum=as.data.frame(t(phylum))
orders=as.data.frame(t(orders))
class=as.data.frame(t(class))
bacteria=as.data.frame(t(bacteria))

genus=genus[order(rownames(genus)),]
family=family[order(rownames(family)),]
phylum=phylum[order(rownames(phylum)),]
orders=orders[order(rownames(orders)),]
class=class[order(rownames(class)),]
bacteria=bacteria[order(rownames(bacteria)),]
phenotype_bac=phenotype_bac[order(phenotype_bac$`16SSampleID`),]
rownames(genus)=phenotype_bac$biopsy_number
rownames(family)=phenotype_bac$biopsy_number
rownames(phylum)=phenotype_bac$biopsy_number
rownames(orders)=phenotype_bac$biopsy_number
rownames(class)=phenotype_bac$biopsy_number
rownames(bacteria)=phenotype_bac$biopsy_number

phenotype_bac$Batch=NA
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run1"]="batch1"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run2"]="batch2"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run3"]="batch3"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run4"]="batch4"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run5"]="batch5"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run6"]="batch6"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run7"]="batch7"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run8"]="batch8"
phenotype_bac$Batch[phenotype_bac$`16SSampleID` %like% "run9"]="batch9"

phenotype_bac$BMI=gsub(",",".",phenotype_bac$BMI)
phenotype_bac$BMI=as.numeric(phenotype_bac$BMI)

phenotype_bac$Aminosalicylates=as.numeric(phenotype_bac$Aminosalicylates)
phenotype_bac$Aminosalicylates[is.na(phenotype_bac$Aminosalicylates)]=0
phenotype_bac$Thiopurines=as.numeric(phenotype_bac$Thiopurines)
phenotype_bac$Thiopurines[is.na(phenotype_bac$Thiopurines)]=0
phenotype_bac$Methotrexaat=as.numeric(phenotype_bac$Methotrexaat)
phenotype_bac$Methotrexaat[is.na(phenotype_bac$Methotrexaat)]=0
phenotype_bac$Steroids=as.numeric(phenotype_bac$Steroids)
phenotype_bac$Steroids[is.na(phenotype_bac$Steroids)]=0
phenotype_bac$resec_part_colon=as.numeric(phenotype_bac$resec_part_colon)
phenotype_bac$resec_part_colon[is.na(phenotype_bac$resec_part_colon)]=0
phenotype_bac$resec_ileocec=as.numeric(phenotype_bac$resec_ileocec)
phenotype_bac$resec_ileocec[is.na(phenotype_bac$resec_ileocec)]=0
phenotype_bac$resec_part_small=as.numeric(phenotype_bac$resec_part_small)
phenotype_bac$resec_part_small[is.na(phenotype_bac$resec_part_small)]=0

covariate_assess=c("ResearchID","Diagnosis","Cohort","Inflammation","Location_rough","Batch","age_at_biopsy",
                   "sex","BMI","smoking_DB","MontrealA","MontrealL","MontrealB","MontrealE","MontrealS","Aminosalicylates",
                   "Thiopurines","Steroids","Methotrexaat","resec_part_colon","resec_ileocec","resec_part_small",
                   "SCCAI","HBI")
covariate_bac=phenotype_bac[,colnames(phenotype_bac) %in% covariate_assess]
rownames(covariate_bac)=phenotype_bac$biopsy_number

write.table(covariate_bac,file = "Covariate.bacteria.txt",sep = "\t",quote = F,row.names = T)
write.table(genus,file = "Decontamination_dada2/Input.genus.txt",row.names = T,quote = F,sep = "\t")
write.table(family,file = "Decontamination_dada2/Input.family.txt",row.names = T,quote = F,sep = "\t")
write.table(phylum,file = "Decontamination_dada2/Input.phylum.txt",row.names = T,quote = F,sep = "\t")
write.table(class,file = "Decontamination_dada2/Input.class.txt",row.names = T,quote = F,sep = "\t")
write.table(orders,file = "Decontamination_dada2/Input.order.txt",row.names = T,quote = F,sep = "\t")
write.table(bacteria,file = "Decontamination_dada2/Input.bacteria.txt",row.names = T,quote = F,sep = "\t")
