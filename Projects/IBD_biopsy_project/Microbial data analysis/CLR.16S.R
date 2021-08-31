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
# bacteria abundance transformation, clr
# ===========================================================================================

# fucntions are from Johannes, rows are samples and colums are taxa, and the input is abundance table
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
family=read.table("Decontamination_dada2/Input.family.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
orders=read.table("Decontamination_dada2/Input.order.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
phylum=read.table("Decontamination_dada2/Input.phylum.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
class=read.table("Decontamination_dada2/Input.class.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

genus=genus[,colSums(genus)>0]
genus=genus[rowSums(genus)>0,]
genus=apply(genus,1,function(x){
  x=x/sum(x)
})
genus=as.data.frame(t(genus))

family=family[,colSums(family)>0]
family=family[rowSums(family)>0,]
family=apply(family,1,function(x){
  x=x/sum(x)
})
family=as.data.frame(t(family))

phylum=phylum[,colSums(phylum)>0]
phylum=phylum[rowSums(phylum)>0,]
phylum=apply(phylum,1,function(x){
  x=x/sum(x)
})
phylum=as.data.frame(t(phylum))

class=class[,colSums(class)>0]
class=class[rowSums(class)>0,]
class=apply(class,1,function(x){
  x=x/sum(x)
})
class=as.data.frame(t(class))

orders=orders[,colSums(orders)>0]
orders=orders[rowSums(orders)>0,]
orders=apply(orders,1,function(x){
  x=x/sum(x)
})
orders=as.data.frame(t(orders))

genus_clr <- zCompositions::cmultRepl(genus, method="CZM", label=0)
genus_clr = compositions::clr(genus_clr)
genus_clr=as.data.frame(genus_clr)
genus_clr = genus_clr[,colSums(genus>0)>nrow(genus) * 0.1]

family_clr <- zCompositions::cmultRepl(family, method="CZM", label=0)
family_clr = compositions::clr(family_clr)
family_clr=as.data.frame(family_clr)
family_clr = family_clr[,colSums(family>0)>nrow(family) * 0.1]

order_clr <- zCompositions::cmultRepl(orders, method="CZM", label=0)
order_clr = compositions::clr(order_clr)
order_clr=as.data.frame(order_clr)
order_clr = order_clr[,colSums(orders>0)>nrow(orders) * 0.1]

phylum_clr <- zCompositions::cmultRepl(phylum, method="CZM", label=0)
phylum_clr = compositions::clr(phylum_clr)
phylum_clr=as.data.frame(phylum_clr)
phylum_clr = phylum_clr[,colSums(phylum>0)>nrow(phylum) * 0.1]

class_clr <- zCompositions::cmultRepl(class, method="CZM", label=0)
class_clr = compositions::clr(class_clr)
class_clr=as.data.frame(class_clr)
class_clr = class_clr[,colSums(class>0)>nrow(class) * 0.1]

write.table(genus_clr,"OutputTable/CLR.genus.txt",row.names = T,quote = F,sep = "\t")
write.table(family_clr,"OutputTable/CLR.family.txt",row.names = T,quote = F,sep = "\t")
write.table(order_clr,"OutputTable/CLR.order.txt",row.names = T,quote = F,sep = "\t")
write.table(phylum_clr,"OutputTable/CLR.phylum.txt",row.names = T,quote = F,sep = "\t")
write.table(class_clr,"OutputTable/CLR.class.txt",row.names = T,quote = F,sep = "\t")
write.table(cbind(genus_clr,family_clr,class_clr,order_clr,phylum_clr),"OutputTable/CLR.bacteria.txt",row.names = T,quote = F,sep = "\t")

covariate_bac=covariate_bac[!is.na(covariate_bac$Diagnosis),]
covariate_bac=covariate_bac[!is.na(covariate_bac$Inflammation),]
covariate_bac$BMI[is.na(covariate_bac$BMI)]=median(covariate_bac$BMI[!is.na(covariate_bac$BMI)])
covariate_bac$smoking_DB[is.na(covariate_bac$smoking_DB)]=median(covariate_bac$smoking_DB[!is.na(covariate_bac$smoking_DB)])

write.table(covariate_bac,"OutputTable/Covariate_bac.organized.txt",sep = "\t",row.names = T,quote = F)

