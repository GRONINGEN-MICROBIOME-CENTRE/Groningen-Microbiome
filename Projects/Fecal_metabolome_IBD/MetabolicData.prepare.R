# =============================================================================================
#  This script is for preparing metabolites data for GWAS
# =============================================================================================
library(ggplot2)
library(factoextra)
library(vegan)
library(ggsci)
source("Microbiome.function.R")

# ====================================== quantitative trait ==============================================
# prepare covariates
phenotype_quan_IBD=read.table("2.Input/IBD/ibd_phenos2.txt",sep = "\t",row.names = 1,check.names = F,header = T,stringsAsFactors = F)
phenotype_quan_HC=read.table("2.Input/Controls/controls_phenos2.txt",sep = "\t",row.names = 1,check.names = F,header = T,stringsAsFactors = F)

covariate_quan_IBD=phenotype_quan_IBD[,c("ibd_Diagnosis","host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer"),drop=F]
covariate_quan_CD=covariate_quan_IBD[covariate_quan_IBD$ibd_Diagnosis=="CD",]
covariate_quan_UC=covariate_quan_IBD[covariate_quan_IBD$ibd_Diagnosis=="UC",]
covariate_quan_HC=phenotype_quan_HC[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer"),drop=F]

coupling=read.table("2.Input/All.coupling.txt",header = F,sep = "\t")

pca_CD=read.table("2.Input/Genetics.PC//CD.pca.txt",header = T,sep = " ")
rownames(pca_CD)=pca_CD$FID
pca_CD=pca_CD[,c("PC1","PC2","PC3","PC4","PC5")]
pca_CD=merge(pca_CD,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_CD)=pca_CD$V2
pca_CD$Row.names=NULL
pca_CD$V2=NULL

pca_UC=read.table("2.Input/Genetics.PC/UC.pca.txt",header = T,sep = " ")
rownames(pca_UC)=pca_UC$FID
pca_UC=pca_UC[,c("PC1","PC2","PC3","PC4","PC5")]
pca_UC=merge(pca_UC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_UC)=pca_UC$V2
pca_UC$Row.names=NULL
pca_UC$V2=NULL

pca_HC=read.table("2.Input/Genetics.PC/HC.pca.txt",header = T,sep = " ")
rownames(pca_HC)=pca_HC$IID
pca_HC=pca_HC[,c("PC1","PC2","PC3","PC4","PC5")]
pca_HC=merge(pca_HC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_HC)=pca_HC$V2
pca_HC$Row.names=NULL
pca_HC$V2=NULL

pca=rbind(pca_CD,pca_UC,pca_HC)

quantitative_metabolic=read.table("2.Input/Case_control/cc_v4_in.txt",header = T,sep = "\t",check.names = F,row.names = 1)

# split by groups
covariate_quan_CD=merge(covariate_quan_CD,pca_CD,by="row.names",all=F)
covariate_quan_UC=merge(covariate_quan_UC,pca_UC,by="row.names",all=F)
covariate_quan_HC=merge(covariate_quan_HC,pca_HC,by="row.names",all=F)
covariate_quan_CD$ibd_Diagnosis=NULL
covariate_quan_UC$ibd_Diagnosis=NULL
covariate_quan_HC$clinical_BowelMovementADayDef[is.na(covariate_quan_HC$clinical_BowelMovementADayDef)]=median(!is.na(covariate_quan_HC$clinical_BowelMovementADayDef))

quantitative_CD=quantitative_metabolic[rownames(quantitative_metabolic) %in% covariate_quan_CD$Row.names,]
quantitative_UC=quantitative_metabolic[rownames(quantitative_metabolic) %in% covariate_quan_UC$Row.names,]
quantitative_HC=quantitative_metabolic[rownames(quantitative_metabolic) %in% covariate_quan_HC$Row.names,]

# re-order samples, remove first 14 columns which are not features
covariate_quan_CD=covariate_quan_CD[order(covariate_quan_CD$Row.names),]
covariate_quan_UC=covariate_quan_UC[order(covariate_quan_UC$Row.names),]
covariate_quan_HC=covariate_quan_HC[order(covariate_quan_HC$Row.names),]
quantitative_CD=quantitative_CD[order(rownames(quantitative_CD)),]
quantitative_UC=quantitative_UC[order(rownames(quantitative_UC)),]
quantitative_HC=quantitative_HC[order(rownames(quantitative_HC)),]
quantitative_CD=quantitative_CD[,c(15:ncol(quantitative_CD))]
quantitative_UC=quantitative_UC[,c(15:ncol(quantitative_UC))]
quantitative_HC=quantitative_HC[,c(15:ncol(quantitative_HC))]

# PCA for separate cohorts
quantitative_CD_permutated=as.data.frame(apply(quantitative_CD,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_CD=prcomp((quantitative_CD_permutated),scale = TRUE)
eigenvalue=get_eig(pca_CD)
ind <- get_pca_ind(pca_CD)
pca_CD_matrix=as.data.frame(ind$coord)
ggplot (pca_CD_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.CD.PCA.pdf",width = 10,height = 4)
summary_pc_CD=as.data.frame(t(data.frame(summary(pca_CD)$importance)))
summary_pc_CD$PCs=rownames(summary_pc_CD)
summary_pc_CD$PCs=factor(summary_pc_CD$PCs,levels = summary_pc_CD$PCs)
ggplot(summary_pc_CD, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.CD.PCA.accumulate.pdf",width = 10,height = 4)

quantitative_UC_permutated=as.data.frame(apply(quantitative_UC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_UC=prcomp((quantitative_UC_permutated),scale = TRUE)
eigenvalue=get_eig(pca_UC)
ind <- get_pca_ind(pca_UC)
pca_UC_matrix=as.data.frame(ind$coord)
ggplot (pca_UC_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.UC.PCA.pdf",width = 10,height = 4)
summary_pc_UC=as.data.frame(t(data.frame(summary(pca_UC)$importance)))
summary_pc_UC$PCs=rownames(summary_pc_UC)
summary_pc_UC$PCs=factor(summary_pc_UC$PCs,levels = summary_pc_UC$PCs)
ggplot(summary_pc_UC, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.UC.PCA.accumulate.pdf",width = 10,height = 4)

quantitative_HC_permutated=as.data.frame(apply(quantitative_HC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_HC=prcomp((quantitative_HC_permutated),scale = TRUE)
eigenvalue=get_eig(pca_HC)
ind <- get_pca_ind(pca_HC)
pca_HC_matrix=as.data.frame(ind$coord)
ggplot (pca_HC_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.HC.PCA.pdf",width = 10,height = 4)
summary_pc_HC=as.data.frame(t(data.frame(summary(pca_HC)$importance)))
summary_pc_HC$PCs=rownames(summary_pc_HC)
summary_pc_HC$PCs=factor(summary_pc_HC$PCs,levels = summary_pc_HC$PCs)
ggplot(summary_pc_HC, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.HC.PCA.accumulate.pdf",width = 10,height = 4)

# linear correction for quantitative traits
covariate_quan_CD$host_Sex[covariate_quan_CD$host_Sex=="female"]=1
covariate_quan_CD$host_Sex[covariate_quan_CD$host_Sex=="male"]=2
covariate_quan_UC$host_Sex[covariate_quan_UC$host_Sex=="female"]=1
covariate_quan_UC$host_Sex[covariate_quan_UC$host_Sex=="male"]=2
covariate_quan_HC$host_Sex[covariate_quan_HC$host_Sex=="female"]=1
covariate_quan_HC$host_Sex[covariate_quan_HC$host_Sex=="male"]=2
covariate_quan_CD$host_Sex=as.numeric(covariate_quan_CD$host_Sex)
covariate_quan_UC$host_Sex=as.numeric(covariate_quan_UC$host_Sex)
covariate_quan_HC$host_Sex=as.numeric(covariate_quan_HC$host_Sex)

rownames(covariate_quan_CD)=covariate_quan_CD$Row.names
covariate_quan_CD$Row.names=NULL
covariate_quan_CD=as.data.frame(apply(covariate_quan_CD,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
quantitative_CD_correct = apply(quantitative_CD,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_quan_CD[!is.na(x),,drop = FALSE]

  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
quantitative_CD_correct=as.data.frame(quantitative_CD_correct)

rownames(covariate_quan_UC)=covariate_quan_UC$Row.names
covariate_quan_UC$Row.names=NULL
covariate_quan_UC=as.data.frame(apply(covariate_quan_UC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
quantitative_UC_correct = apply(quantitative_UC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_quan_UC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
quantitative_UC_correct=as.data.frame(quantitative_UC_correct)

rownames(covariate_quan_HC)=covariate_quan_HC$Row.names
covariate_quan_HC$Row.names=NULL
covariate_quan_HC=as.data.frame(apply(covariate_quan_HC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
quantitative_HC_correct = apply(quantitative_HC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_quan_HC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
quantitative_HC_correct=as.data.frame(quantitative_HC_correct)

# write out quantitative phenotype files
quantitative_CD_correct = as.data.frame(t(quantitative_CD_correct))
quantitative_CD_correct = cbind(rownames(quantitative_CD_correct),quantitative_CD_correct)
colnames(quantitative_CD_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(quantitative_CD_correct),
                   Gene = rownames(quantitative_CD_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(quantitative_CD_correct, file = "OutputWithCorrection/CD_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/CD_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
write.table(colnames(quantitative_CD_correct),file = "OutputWithCorrection/CD.list.txt",sep="\t",row.names=F,quote = F)

quantitative_UC_correct = as.data.frame(t(quantitative_UC_correct))
quantitative_UC_correct = cbind(rownames(quantitative_UC_correct),quantitative_UC_correct)
colnames(quantitative_UC_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(quantitative_UC_correct),
                   Gene = rownames(quantitative_UC_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(quantitative_UC_correct, file = "OutputWithCorrection/UC_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/UC_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
write.table(colnames(quantitative_UC_correct),file = "OutputWithCorrection/UC.list.txt",sep="\t",row.names=F,quote = F)

quantitative_HC_correct = as.data.frame(t(quantitative_HC_correct))
quantitative_HC_correct = cbind(rownames(quantitative_HC_correct),quantitative_HC_correct)
colnames(quantitative_HC_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(quantitative_HC_correct),
                   Gene = rownames(quantitative_HC_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(quantitative_HC_correct, file = "OutputWithCorrection/HC_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/HC_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
write.table(colnames(quantitative_HC_correct),file = "OutputWithCorrection/HC.list.txt",sep="\t",row.names=F,quote = F)



# ====================================== binary trait ==============================================
# prepare covariates
phenotype_binary_IBD=read.table("2.Input/IBD/ibd_phenos_prev.txt",sep = "\t",row.names = 1,check.names = F,header = T,stringsAsFactors = F)
phenotype_binary_HC=read.table("2.Input/Controls/controls_phenos_prev2.txt",sep = "\t",row.names = 1,check.names = F,header = T,stringsAsFactors = F)

covariate_binary_IBD=phenotype_binary_IBD[,c("ibd_Diagnosis","host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",
                                             "run_day_cat"),drop=F]
covariate_binary_CD=covariate_binary_IBD[covariate_binary_IBD$ibd_Diagnosis=="CD",]
covariate_binary_UC=covariate_binary_IBD[covariate_binary_IBD$ibd_Diagnosis=="UC",]
covariate_binary_HC=phenotype_binary_HC[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",
                                           "run_day_cat"),drop=F]

coupling=read.table("2.Input/All.coupling.txt",header = F,sep = "\t")

pca_CD=read.table("2.Input/Genetics.PC//CD.pca.txt",header = T,sep = " ")
rownames(pca_CD)=pca_CD$FID
pca_CD=pca_CD[,c("PC1","PC2","PC3","PC4","PC5")]
pca_CD=merge(pca_CD,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_CD)=pca_CD$V2
pca_CD$Row.names=NULL
pca_CD$V2=NULL

pca_UC=read.table("2.Input/Genetics.PC/UC.pca.txt",header = T,sep = " ")
rownames(pca_UC)=pca_UC$FID
pca_UC=pca_UC[,c("PC1","PC2","PC3","PC4","PC5")]
pca_UC=merge(pca_UC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_UC)=pca_UC$V2
pca_UC$Row.names=NULL
pca_UC$V2=NULL

pca_HC=read.table("2.Input/Genetics.PC/HC.pca.txt",header = T,sep = " ")
rownames(pca_HC)=pca_HC$IID
pca_HC=pca_HC[,c("PC1","PC2","PC3","PC4","PC5")]
pca_HC=merge(pca_HC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_HC)=pca_HC$V2
pca_HC$Row.names=NULL
pca_HC$V2=NULL

binary_metabolic=read.table("2.Input/Case_control/cc_v4_prev2.txt",header = T,sep = "\t",check.names = F,row.names = 1)

# split by groups
covariate_binary_CD=merge(covariate_binary_CD,pca_CD,by="row.names",all=F)
covariate_binary_UC=merge(covariate_binary_UC,pca_UC,by="row.names",all=F)
covariate_binary_HC=merge(covariate_binary_HC,pca_HC,by="row.names",all=F)
covariate_binary_CD$ibd_Diagnosis=NULL
covariate_binary_UC$ibd_Diagnosis=NULL
covariate_binary_HC$clinical_BowelMovementADayDef[is.na(covariate_binary_HC$clinical_BowelMovementADayDef)]=median(!is.na(covariate_binary_HC$clinical_BowelMovementADayDef))

binary_CD=binary_metabolic[rownames(binary_metabolic) %in% covariate_binary_CD$Row.names,]
binary_UC=binary_metabolic[rownames(binary_metabolic) %in% covariate_binary_UC$Row.names,]
binary_HC=binary_metabolic[rownames(binary_metabolic) %in% covariate_binary_HC$Row.names,]

# re-order samples, remove first 14 columns which are not features
covariate_binary_CD=covariate_binary_CD[order(covariate_binary_CD$Row.names),]
covariate_binary_UC=covariate_binary_UC[order(covariate_binary_UC$Row.names),]
covariate_binary_HC=covariate_binary_HC[order(covariate_binary_HC$Row.names),]
binary_CD=binary_CD[order(rownames(binary_CD)),]
binary_UC=binary_UC[order(rownames(binary_UC)),]
binary_HC=binary_HC[order(rownames(binary_HC)),]
binary_CD=binary_CD[,c(325:ncol(binary_CD))]
binary_UC=binary_UC[,c(325:ncol(binary_UC))]
binary_HC=binary_HC[,c(325:ncol(binary_HC))]

binary_CD[!is.na(binary_CD)]=1
binary_CD[is.na(binary_CD)]=0
binary_UC[!is.na(binary_UC)]=1
binary_UC[is.na(binary_UC)]=0
binary_HC[!is.na(binary_HC)]=1
binary_HC[is.na(binary_HC)]=0
# PCA for separate cohorts
binary_CD_permutated=as.data.frame(apply(binary_CD,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_CD=prcomp((binary_CD_permutated),scale = TRUE)
eigenvalue=get_eig(pca_CD)
ind <- get_pca_ind(pca_CD)
pca_CD_matrix=as.data.frame(ind$coord)
ggplot (pca_CD_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.CD.PCA.pdf",width = 10,height = 4)
summary_pc_CD=as.data.frame(t(data.frame(summary(pca_CD)$importance)))
summary_pc_CD$PCs=rownames(summary_pc_CD)
summary_pc_CD$PCs=factor(summary_pc_CD$PCs,levels = summary_pc_CD$PCs)
ggplot(summary_pc_CD, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.CD.PCA.accumulate.pdf",width = 10,height = 4)

binary_UC_permutated=as.data.frame(apply(binary_UC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_UC=prcomp((binary_UC_permutated),scale = TRUE)
eigenvalue=get_eig(pca_UC)
ind <- get_pca_ind(pca_UC)
pca_UC_matrix=as.data.frame(ind$coord)
ggplot (pca_UC_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.UC.PCA.pdf",width = 10,height = 4)
summary_pc_UC=as.data.frame(t(data.frame(summary(pca_UC)$importance)))
summary_pc_UC$PCs=rownames(summary_pc_UC)
summary_pc_UC$PCs=factor(summary_pc_UC$PCs,levels = summary_pc_UC$PCs)
ggplot(summary_pc_UC, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.UC.PCA.accumulate.pdf",width = 10,height = 4)

binary_HC_permutated=as.data.frame(apply(binary_HC,2,function(x){
  x[which(is.na(x))]=median(!is.na(x))
  return(x)
}))
pca_HC=prcomp((binary_HC_permutated),scale = TRUE)
eigenvalue=get_eig(pca_HC)
ind <- get_pca_ind(pca_HC)
pca_HC_matrix=as.data.frame(ind$coord)
ggplot (pca_HC_matrix, aes(Dim.1,Dim.2)) + geom_point(size=2) + theme_bw()
ggsave("Plot/Metabolites.HC.PCA.pdf",width = 10,height = 4)
summary_pc_HC=as.data.frame(t(data.frame(summary(pca_HC)$importance)))
summary_pc_HC$PCs=rownames(summary_pc_HC)
summary_pc_HC$PCs=factor(summary_pc_HC$PCs,levels = summary_pc_HC$PCs)
ggplot(summary_pc_HC, aes(x = PCs,y = `Cumulative Proportion`)) +
  geom_hline(yintercept=0.8,color="red")+
  geom_bar(stat="identity",fill="white",color="black") +theme_bw()+guides(fill=F)
ggsave("Plot/Metabolites.HC.PCA.accumulate.pdf",width = 10,height = 4)

# glm correction for binary traits
covariate_binary_CD$host_Sex[covariate_binary_CD$host_Sex=="female"]=1
covariate_binary_CD$host_Sex[covariate_binary_CD$host_Sex=="male"]=2
covariate_binary_UC$host_Sex[covariate_binary_UC$host_Sex=="female"]=1
covariate_binary_UC$host_Sex[covariate_binary_UC$host_Sex=="male"]=2
covariate_binary_HC$host_Sex[covariate_binary_HC$host_Sex=="female"]=1
covariate_binary_HC$host_Sex[covariate_binary_HC$host_Sex=="male"]=2
covariate_binary_CD$host_Sex=as.numeric(covariate_binary_CD$host_Sex)
covariate_binary_UC$host_Sex=as.numeric(covariate_binary_UC$host_Sex)
covariate_binary_HC$host_Sex=as.numeric(covariate_binary_HC$host_Sex)

rownames(covariate_binary_CD)=covariate_binary_CD$Row.names
covariate_binary_CD$Row.names=NULL
covariate_binary_CD$host_BMI[is.na(covariate_binary_CD$host_BMI)]=median(!is.na(covariate_binary_CD$host_BMI))
binary_CD_correct = apply(binary_CD,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_binary_CD[!is.na(x),,drop = FALSE]
  
  x.resid = resid(glm(x.subset ~ .,data = covariate.subset,family=binomial))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
binary_CD_correct=as.data.frame(binary_CD_correct)

rownames(covariate_binary_UC)=covariate_binary_UC$Row.names
covariate_binary_UC$Row.names=NULL
covariate_binary_UC$host_BMI[is.na(covariate_binary_UC$host_BMI)]=median(!is.na(covariate_binary_UC$host_BMI))
binary_UC_correct = apply(binary_UC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_binary_UC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(glm(x.subset ~ .,data = covariate.subset,family=binomial))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
binary_UC_correct=as.data.frame(binary_UC_correct)

rownames(covariate_binary_HC)=covariate_binary_HC$Row.names
covariate_binary_HC$Row.names=NULL
binary_HC_correct = apply(binary_HC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = covariate_binary_HC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(glm(x.subset ~ .,data = covariate.subset,family=binomial))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
binary_HC_correct=as.data.frame(binary_HC_correct)

# write out binary phenotype files
binary_CD_correct = as.data.frame(t(binary_CD_correct))
binary_CD_correct = cbind(rownames(binary_CD_correct),binary_CD_correct)
colnames(binary_CD_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(binary_CD_correct),
                   Gene = rownames(binary_CD_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(binary_CD_correct, file = "OutputWithCorrection/CD_binary.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/CD_binary.metabolic.txt.annot",sep="\t",row.names=F,quote = F)

binary_UC_correct = as.data.frame(t(binary_UC_correct))
binary_UC_correct = cbind(rownames(binary_UC_correct),binary_UC_correct)
colnames(binary_UC_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(binary_UC_correct),
                   Gene = rownames(binary_UC_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(binary_UC_correct, file = "OutputWithCorrection/UC_binary.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/UC_binary.metabolic.txt.annot",sep="\t",row.names=F,quote = F)

binary_HC_correct = as.data.frame(t(binary_HC_correct))
binary_HC_correct = cbind(rownames(binary_HC_correct),binary_HC_correct)
colnames(binary_HC_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(binary_HC_correct),
                   Gene = rownames(binary_HC_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(binary_HC_correct, file = "OutputWithCorrection/HC_binary.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputWithCorrection/HC_binary.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
