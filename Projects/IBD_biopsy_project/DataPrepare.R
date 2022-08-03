# this script is for cleaning 16S+RNAseq project data

library(edgeR)
library(limma)
library(mixOmics)
library(HTSFilter)
library(foreach)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
require(cowplot)
library(ggpubr)

# 16S data
bacteria=read.table("Decontamination_dada2/Genus.txt",sep = "\t",row.names = 1,stringsAsFactors = F,check.names = F)
bacteria=as.data.frame(t(bacteria))
phenotype_bac=read.table("MetaData/Metadaata.16S.txt",fill = T,stringsAsFactors = F,header = T,sep = "\t",check.names = F)

# RNA-seq data
gene1=read.table("RNAseq/NewRelease.GeneCount.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
gene2=read.table("RNAseq/OldRelease.ExpressionTable.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
rownames(gene1)=gsub("\\..*","",rownames(gene1))
overlap=intersect(rownames(gene1),rownames(gene2))
gene1=gene1[rownames(gene1) %in% overlap,]
gene2=gene2[rownames(gene2) %in% overlap,]
colnames(gene1)=gsub(".count","",colnames(gene1))
coupling=read.table("MetaData/couplingFile.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t")
coupling=coupling[coupling$`RNA-seqID` %in% colnames(gene1),]
gene1=gene1[,colnames(gene1) %in% coupling$`RNA-seqID`]

coupling=coupling[order(coupling$`RNA-seqID`),]
gene1=gene1[,order(colnames(gene1))]
colnames(gene1)=coupling$ChangedID
gene=cbind(gene2,gene1)

phenotype_rna=read.table("MetaData/Metadaata.RNAseq.txt",header = T,stringsAsFactors = F,fill = T,check.names = F,sep = "\t")

write.table(gene,file = "RNAseq/Merged.txt",quote = F,row.names = T,sep = "\t")

# normalize and transform
count1=as.data.frame(t(gene1))
count1=count1[rowSums(count1)!=0,colSums(count1)!=0]
dgeFull <- DGEList(count1, remove.zeros = TRUE) # generate project
dgeFull <- calcNormFactors(dgeFull, method="TMM") # TMM normalization
expression1=cpm(dgeFull,log = TRUE, prior.count = 1) 
expression1=as.data.frame(expression1)

count2=as.data.frame(t(gene2))
count2=count2[rowSums(count2)!=0,colSums(count2)!=0]
dgeFull <- DGEList(count2, remove.zeros = TRUE) # generate project
dgeFull <- calcNormFactors(dgeFull, method="TMM") # TMM normalization
expression2=cpm(dgeFull,log = TRUE, prior.count = 1) 
expression2=as.data.frame(expression2)

expression1=expression1[,colnames(expression1) %in% colnames(expression2)]
expression2=expression2[,colnames(expression2) %in% colnames(expression1)]

expression=rbind(expression2,expression1)
write.table(expression,file = "RNAseq/Merged.normalized.txt",quote = F,row.names = T,sep = "\t")

# ===============================================================================
# paired samples
# ===============================================================================

phenotype_bac=read.table("MetaData/Metadaata.16S.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
phenotype_rna=read.table("MetaData/Metadaata.RNAseq.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# colon vs. ileum
unique_individual=unique(phenotype_bac$ResearchID)

phenotype_bac.tmp=phenotype_bac[phenotype_bac$Inflammation=="No",]
multi_location=as.data.frame(matrix(ncol = ncol(phenotype_bac.tmp)))
colnames(multi_location)=colnames(phenotype_bac.tmp)
for(i in unique_individual){
  tmp.individual=i
  tmp.subset=phenotype_bac.tmp[phenotype_bac.tmp$ResearchID==i,]
  
  if(length(unique(tmp.subset$Location_rough))>=2){
    multi_location=rbind(multi_location,tmp.subset)
  }
  
}
multi_location=multi_location[-1,]

multi_location_CD=multi_location[multi_location$Diagnosis=="CD",]
multi_location_UC=multi_location[multi_location$Diagnosis=="UC",]
multi_location_HC=multi_location[multi_location$Diagnosis=="Control",]
length(unique(multi_location_HC$ResearchID))

tmp1=data.frame(ResearchID=multi_location_CD$ResearchID,BiopsyID=multi_location_CD$biopsy_number,Disease="CD",Inflammation="Yes")
tmp2=data.frame(ResearchID=multi_location_UC$ResearchID,BiopsyID=multi_location_UC$biopsy_number,Disease="UC",Inflammation="Yes")
tmp3=data.frame(ResearchID=multi_location_HC$ResearchID,BiopsyID=multi_location_HC$biopsy_number,Disease="HC",Inflammation="Yes")
tmp=rbind(tmp1,tmp2)
write.table(tmp,"OutputTable/Paired.multipleLcation.noninflamed.txt",sep="\t",quote = F,row.names = F)

# inflamed vs. non-inflamed
unique_individual=unique(phenotype_bac$ResearchID)

phenotype_bac.tmp=phenotype_bac[phenotype_bac$Location_rough=="ileum",]
multi_inflammation=as.data.frame(matrix(ncol = ncol(phenotype_bac)))
colnames(multi_inflammation)=colnames(phenotype_bac.tmp)
for(i in unique_individual){
  tmp.individual=i
  tmp.subset=phenotype_bac.tmp[phenotype_bac.tmp$ResearchID==i,]
  
  if(length(unique(tmp.subset$Inflammation))>=2){
    multi_inflammation=rbind(multi_inflammation,tmp.subset)
  }
  
}
multi_inflammation=multi_inflammation[-1,]

multi_inflammation_CD=multi_inflammation[multi_inflammation$Diagnosis=="CD",]
multi_inflammation_UC=multi_inflammation[multi_inflammation$Diagnosis=="UC",]
multi_inflammation_HC=multi_inflammation[multi_inflammation$Diagnosis=="Control",]
length(unique(multi_inflammation_CD$ResearchID))

tmp1=data.frame(ResearchID=multi_inflammation_CD$ResearchID,Disease="CD",BiopsyID=multi_inflammation_CD$biopsy_number,Location="ileum")
tmp2=data.frame(ResearchID=multi_inflammation_UC$ResearchID,Disease="UC",BiopsyID=multi_inflammation_UC$biopsy_number,Location="ileum")
tmp3=data.frame(ResearchID=multi_inflammation_HC$ResearchID,Disease="HC",BiopsyID=multi_inflammation_HC$biopsy_number,Location="ileum")
tmp=rbind(tmp1,tmp2)
write.table(tmp,"OutputTable/Paired.multipleInflammation.ileum.txt",sep="\t",quote = F,row.names = F)


# annotation new releease RNAssseq
gene1=read.table("RNAseq/NewRelease.GeneCount.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
rownames(gene1)=gsub("\\..*","",rownames(gene1))
expression=gene1

annotation=read.table("annotation.file.txt",stringsAsFactors = F,header = T)
annotation=annotation[!duplicated(annotation$id),]
expression=expression[rownames(expression) %in% annotation$Gene,]
annotation=annotation[annotation$Gene %in% rownames(expression),]

annotation=annotation[order(annotation$Gene),]
expression=expression[order(rownames(expression)),]

rownames(expression)[20000]
annotation$Gene[20000]

rownames(expression)=annotation$id
covariate_rna=covariate_rna[covariate_rna$ResearchID %like% "UMCGIBD",]
coupling=merge(coupling,covariate_rna,by.x="ChangedID",by.y="row.names",all=F)
expression=expression[,colnames(expression) %in% rownames(covariate_rna)]
write.table(expression,"RNAseq/New.release.annot.txt",sep = "\t",row.names = T,quote = F)





























