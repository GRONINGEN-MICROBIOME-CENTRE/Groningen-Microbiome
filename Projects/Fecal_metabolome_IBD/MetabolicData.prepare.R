# =============================================================================================
#  This script is for preparing metabolites data for GWAS
# =============================================================================================

# prepare covariates
covariate=read.table("Input/Covariate.txt",header = T,sep = "\t",stringsAsFactors = F,check.names = F)
rownames(covariate)=covariate$ID
covariate$ID=NULL
covariate=covariate[,c("LC.COLUMN","run_day_cat","Amount_sample_gram","metabolon_Month_in_freezer",
                       "host_Age_sampling","host_Sex","ibd_Diagnosis")]
coupling=read.table("Input/All.coupling.txt",header = F,sep = "\t")

pca_CD=read.table("Input/CD.pca.txt",header = T,sep = " ")
rownames(pca_CD)=pca_CD$FID
pca_CD=pca_CD[,c("PC1","PC2","PC3")]
pca_CD=merge(pca_CD,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_CD)=pca_CD$V2
pca_CD$Row.names=NULL
pca_CD$V2=NULL

pca_UC=read.table("Input/UC.pca.txt",header = T,sep = " ")
rownames(pca_UC)=pca_UC$FID
pca_UC=pca_UC[,c("PC1","PC2","PC3")]
pca_UC=merge(pca_UC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_UC)=pca_UC$V2
pca_UC$Row.names=NULL
pca_UC$V2=NULL

pca_HC=read.table("Input/HC.pca.txt",header = T,sep = " ")
rownames(pca_HC)=pca_HC$IID
pca_HC=pca_HC[,c("PC1","PC2","PC3")]
pca_HC=merge(pca_HC,coupling,by.x="row.names",by.y="V1",all=F)
rownames(pca_HC)=pca_HC$V2
pca_HC$Row.names=NULL
pca_HC$V2=NULL

pca=rbind(pca_CD,pca_UC,pca_HC)

quantitative=read.table("Input/Quantitative.metabolites.txt",header = T,sep = "\t",check.names = F,row.names = 1)
binary=read.table("Input/Binary.metabolites.txt",header = T,sep = "\t",check.names = F,row.names = 1)

covariate=merge(covariate,pca,by="row.names",all = F)
rownames(covariate)=covariate$Row.names
covariate$Row.names=NULL
covariate$run_day_cat=as.factor(covariate$run_day_cat)
covariate$host_Sex=as.factor(covariate$host_Sex)

# split by groups
Covariate_CD=covariate[covariate$ibd_Diagnosis=="CD",]
Covariate_UC=covariate[covariate$ibd_Diagnosis=="UC",]
Covariate_HC=covariate[covariate$ibd_Diagnosis=="A_Control",]
Covariate_CD$ibd_Diagnosis=NULL
Covariate_UC$ibd_Diagnosis=NULL
Covariate_HC$ibd_Diagnosis=NULL

quantitative_CD=quantitative[rownames(quantitative) %in% rownames(Covariate_CD),]
quantitative_UC=quantitative[rownames(quantitative) %in% rownames(Covariate_UC),]
quantitative_HC=quantitative[rownames(quantitative) %in% rownames(Covariate_HC),]

binary_CD=binary[rownames(binary) %in% rownames(Covariate_CD),]
binary_UC=binary[rownames(binary) %in% rownames(Covariate_UC),]
binary_HC=binary[rownames(binary) %in% rownames(Covariate_HC),]

# re-order samples
Covariate_CD=Covariate_CD[order(rownames(Covariate_CD)),]
Covariate_UC=Covariate_UC[order(rownames(Covariate_UC)),]
Covariate_HC=Covariate_HC[order(rownames(Covariate_HC)),]
quantitative_CD=quantitative_CD[order(rownames(quantitative_CD)),]
quantitative_UC=quantitative_UC[order(rownames(quantitative_UC)),]
quantitative_HC=quantitative_HC[order(rownames(quantitative_HC)),]
binary_CD=binary_CD[order(rownames(binary_CD)),]
binary_UC=binary_UC[order(rownames(binary_UC)),]
binary_HC=binary_HC[order(rownames(binary_HC)),]

# linear correction for quantitative traits
quantitative_CD_correct = apply(quantitative_CD,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_CD[!is.na(x),,drop = FALSE]

  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
quantitative_CD_correct=as.data.frame(quantitative_CD_correct)

quantitative_UC_correct = apply(quantitative_UC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_UC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
quantitative_UC_correct=as.data.frame(quantitative_UC_correct)

quantitative_HC_correct = apply(quantitative_HC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_HC[!is.na(x),,drop = FALSE]
  
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

# glm correction for binary traits
binary_CD_correct = apply(binary_CD,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_CD[!is.na(x),,drop = FALSE]
  
  x.resid = resid(glm(x.subset ~ .,data = covariate.subset,family=binomial))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
binary_CD_correct=as.data.frame(binary_CD_correct)

binary_UC_correct = apply(binary_UC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_UC[!is.na(x),,drop = FALSE]
  
  x.resid = resid(glm(x.subset ~ .,data = covariate.subset,family=binomial))
  x[!is.na(x)] = x.resid
  
  return(x)
  
})
binary_UC_correct=as.data.frame(binary_UC_correct)

binary_HC_correct = apply(binary_HC,2,function(x){
  
  x.subset=x[!is.na(x)]
  covariate.subset = Covariate_HC[!is.na(x),,drop = FALSE]
  
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

write.table(Covariate_CD,file = "Covariate.CD.txt",row.names = T,quote = F,sep = "\t")
write.table(Covariate_UC,file = "Covariate.UC.txt",row.names = T,quote = F,sep = "\t")
write.table(Covariate_HC,file = "Covariate.CT.txt",row.names = T,quote = F,sep = "\t")


