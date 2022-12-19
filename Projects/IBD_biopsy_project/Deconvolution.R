# ====================================================================================================
# cell type deconvolution
# ====================================================================================================

library(xCell)
library(foreach)
library(lme4)
library(nlme)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(RColorBrewer)
library(randomcoloR)
library(crayon)
library(ggforce)
library(edgeR)
library(matrixStats)
library(ComplexHeatmap)
library(circlize)
library(vegan)
library(pheatmap)
library(reshape2)

RNAseq=read.table("RNAseq/Merged.normalized.txt",row.names = 1,header = T,stringsAsFactors = F,check.names = F,sep = "\t")
gene=as.data.frame((RNAseq))
deconvolution=xCellAnalysis(gene)

# ====================================================================================================
# start analysis
# ====================================================================================================
covariate_rna=read.table("Covariate.rna.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
deconvolution=read.table("OutputTable/Deconvolution.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)

# PCA
beta_diversity=vegdist((deconvolution_sub),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,covariate_rna,all = F,by="row.names")

# cell type basic correction, + Th1/Th2 cell ratio
deconvolution_sub = apply(deconvolution_sub,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="cell"
  
  fit=glmer(cell~age_at_biopsy+sex+BMI+Aminosalicylates+Thiopurines+Steroids+Batch+Biological_use +(1|ResearchID),data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
deconvolution_sub=as.data.frame(deconvolution_sub)
deconvolution_sub$ThCellRatio=deconvolution_sub$`Th1 cells`/deconvolution_sub$`Th2 cells`
  
# compare CDi vs. UCi
covariate1=covariate_rna[covariate_rna$Diagnosis=="CD" | covariate_rna$Diagnosis=="UC",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
nn=0
result0 = foreach(i=1:ncol(deconvolution_sub1),.combine = rbind) %do%  {
  tmp.genus=colnames(deconvolution_sub1)[i]
  x=deconvolution_sub1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Diagnosis,data=tmp.data) 
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== counting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,diagnosis.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="DiagnosisUC"],
                           diagnosis.beta=coefs.fit$Estimate[rownames(coefs.fit)=="DiagnosisUC"],
                           diagnosis.se=coefs.fit$Std..Error[rownames(coefs.fit)=="DiagnosisUC"])
}
result0$FDR=p.adjust(result0$diagnosis.P,method = "BH")

compare=matrix(nrow = ncol(deconvolution_sub),ncol = 2)
compare=as.data.frame(compare)
colnames(compare)=c("CellType","Pvalue")
for(i in 1:ncol(deconvolution_sub)){
  cell=colnames(deconvolution_sub)[i]
  inflamed=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$Inflammation=="Yes"],]
  noninf=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$Inflammation=="No"],]
  mm=wilcox.test(inflamed[,cell],noninf[,cell])
  Pvalue=mm$p.value
  compare$CellType[i]=cell
  compare$Pvalue[i]=Pvalue
}
compare$FDR=p.adjust(compare$Pvalue)
compare=compare[order(compare$FDR),]

inflamed=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$Inflammation=="Yes"],]
noninf=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$Inflammation=="No"],]
compare_all=data.frame(Score=NA,Group=NA,Cell=NA,Sample=NA)
for(i in 1:nrow(compare)){
  cell=compare$CellType[i]
  sub.in=inflamed[,cell,drop=F]
  sub.non=noninf[,cell,drop=F]
  sub.in$Group="Inflammation"
  sub.non$Group="Non_inflammation"
  sub.in$Cell=cell
  sub.non$Cell=cell
  colnames(sub.in)[1]="Score"
  colnames(sub.non)[1]="Score"
  sub.data=rbind(sub.in,sub.non)
  sub.data$Sample=rownames(sub.data)
  rownames(sub.data)=NULL
  compare_all=rbind(compare_all,sub.data)
}
compare_all=na.omit(compare_all)
compare_all$Cell=factor(compare_all$Cell,levels = rev(unique(compare$CellType)))

# ====================================================================================================
# cell-host, use relative abundance instead of 0/1
# ====================================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
result1=matrix(nrow = 1,ncol = 4)
result1=as.data.frame(result1)
colnames(result1)=c("Bacteria","Cells","Pvalue","Estimate")
for (i in 1:ncol(bacteria)){
  
  tmp.genus=colnames(bacteria)[i]
  tmp.data=(bacteria[,i,drop=F])
  tmp.data=na.omit(tmp.data)
  tmp.cell=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(tmp.data),]
  tmp.cell[tmp.cell==0]=NA
  tmp.data=tmp.data[rownames(tmp.data) %in% rownames(tmp.cell),,drop=F]
  tmp.cell=tmp.cell[order(rownames(tmp.cell)),]
  tmp.data=tmp.data[order(rownames(tmp.data)),,drop=F]
  
  tmp.result=matrix(nrow = ncol(tmp.cell),ncol = 4)
  tmp.result=as.data.frame(tmp.result)
  colnames(tmp.result)=c("Bacteria","Cells","Pvalue","Estimate")
  tmp.result$Bacteria=tmp.genus
  for(n in 1:ncol(tmp.cell)){
    tmp.result$Cells[n]=colnames(tmp.cell)[n]
    tmp=merge(tmp.cell[,n,drop=F],tmp.data,by="row.names",all=F)
    rownames(tmp)=tmp$Row.names
    tmp$Row.names=NULL
    tmp=merge(tmp,basic_factors,by="row.names",all=F)
    tmp$Row.names=NULL
    tmp=na.omit(tmp)
    
    fit=lm(((tmp[,1])) ~ tmp[,2]+Inflammation+Location_rough,data=tmp)
    mm=as.data.frame(summary(fit)$coef)
    tmp.result$Pvalue[n]=mm$`Pr(>|t|)`[2]
    tmp.result$Estimate[n]=mm$Estimate[2]
    
    cat(i,tmp.genus,"-------",colnames(tmp.cell)[n],"is done","\n")
  }
  
  result1=rbind(result1,tmp.result)
  cat(tmp.genus,"DONE=======","\n")
}
result1=result1[-1,]
result1$FDR=p.adjust(result1$Pvalue,method = "BH")

# ====================================== ########## ============================================
# ====================================== predictors ============================================
# ====================================== ########## ============================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
lasso_fit=function(feature,outcome){
  if(length(colnames(feature))==0){
    result=data.frame(Protein=colnames(outcome),Feature="No",Lasso.beta=0)
  }else if(length(colnames(feature))==1){
    model=lm(outcome[,1]~ feature[,1])
    beta = summary(model)$coefficients[2]
    result=data.frame(Protein=colnames(outcome),Feature=colnames(feature),Lasso.beta=beta)
  }else{
    cv=cv.glmnet(as.matrix(feature),as.matrix(outcome), alpha = 1, nfolds = 10, type.measure="mse",standardize=T)
    beta <- (coef(cv, s = "lambda.min"))
    beta=as.data.frame(beta[-1,1])
    beta$variable=rownames(beta)
    colnames(beta)[1]="beta"
    result=data.frame(Protein=colnames(outcome),Feature=beta$variable,Lasso.beta=beta$beta)
  }
}

# ====================================================================================
# models (residules after correcting for repearted measurements)
# model 1: cell ~ factor_group1: basic factors (age, sex, BMI, batch)
# model 2: cell ~ factor_group2: medication (Aminosalicylates, Thiopurines, Steroids, Biological_use)
# model 3: cell ~ factor_group3: Inflammation 
# model 4: cell ~ factor_group4: location
# model 5: cell ~ factor_group5: bacteria
# model 6: cell ~ factor_group6: full model with all factors
# ====================================================================================

factor_group1=basic_factors[,c("age_at_biopsy","sex","BMI","Batch")]
factor_group2=basic_factors[,c("Aminosalicylates","Thiopurines","Steroids","Biological_use")]
factor_group3=basic_factors[,c("Inflammation"),drop=F]
factor_group4=basic_factors[,c("Location_rough"),drop=F]
factor_group5=bacteria
factor_group6=merge(basic_factors,bacteria,by="row.names",all=F)
rownames(factor_group6)=factor_group6$Row.names
factor_group6$Row.names=NULL

# model(s): cell ~ factor_group1/2/3/4/5/6
set.seed(1)
groups=split(deconvolution_sub.correct, sample(1:10, nrow(deconvolution_sub.correct), replace=T))
factors=factor_group5
deconvolution_sub.correct=deconvolution_sub.correct[rownames(deconvolution_sub.correct) %in% rownames(factors),]

variation=list()

for(mm in 1:10){
  
  deconvolution_sub_test <- groups[[mm]]
  deconvolution_sub_train=deconvolution_sub.correct[!rownames(deconvolution_sub.correct) %in% rownames(deconvolution_sub_test),]
  
  trainset= foreach(i=1:ncol(deconvolution_sub_train),.combine = rbind) %do%  {
    
    tmp.protein=colnames(deconvolution_sub_train)[i]
    tmp.feature=as.character(colnames(factors))
    
    cat(yellow(tmp.protein,"lasso running","\n"))
    
    feature=na.omit(factors[,tmp.feature,drop=F])
    outcome=na.omit(deconvolution_sub_train[,tmp.protein,drop=F])
    
    feature=feature[rownames(feature) %in% rownames(outcome),,drop=F]
    outcome=outcome[rownames(outcome) %in% rownames(feature),,drop=F]
    
    feature=feature[order(rownames(feature)),,drop=F]
    outcome=outcome[order(rownames(outcome)),,drop=F]
    
    result=lasso_fit(feature,outcome)
    return.string = result
  }
  trainset=trainset[trainset$Lasso.beta!=0,]
  
  testset.variation=foreach(i=1:length(unique(trainset$Protein)),.combine = rbind) %do%  {
    pro=as.character(unique(trainset$Protein)[i])
    cov=as.character(trainset$Feature[trainset$Protein==pro])
    
    tmp.pro=deconvolution_sub_test[,pro,drop=F]
    tmp.cov=factors[,colnames(factors) %in% cov,drop=F]
    tmp.pro=na.omit(tmp.pro)
    tmp.cov=na.omit(tmp.cov)
    tmp.pro=(tmp.pro[rownames(tmp.pro) %in% rownames(tmp.cov),,drop=F])
    tmp.cov=(tmp.cov[rownames(tmp.cov) %in% rownames(tmp.pro),,drop=F])
    tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
    tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
    
    tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
    tmp.coef=as.data.frame(summary(tmp.glm)$coef)
    tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
    tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
    tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
    tmp.av=anova(tmp.glm)
    tmp.av$Explain=NA
    for(j in 2:nrow(tmp.av)){
      tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
    }
    tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
    
    if(nrow(tmp.av)==0){
      return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA)
    }else{
      return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                               MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,Estimate=tmp.av$Estimate)
    }
  }
  
  variation[[mm]]=testset.variation
  cat(green(mm,"fold is done","=====================","\n"))
  
}

mm1=variation[[1]]
mm2=variation[[2]]
mm3=variation[[3]]
mm4=variation[[4]]
mm5=variation[[5]]
mm6=variation[[6]]
mm7=variation[[7]]
mm8=variation[[8]]
mm9=variation[[9]]
mm10=variation[[10]]
mm1$Fold=1
mm2$Fold=2
mm3$Fold=3
mm4$Fold=4
mm5$Fold=5
mm6$Fold=6
mm7$Fold=7
mm8$Fold=8
mm9$Fold=9
mm10$Fold=10
mm=rbind(mm1,mm2,mm3,mm4,mm5,mm6,mm7,mm8,mm9,mm10)

mm_select=foreach(i=1:length(unique(mm$Protein)), .combine = rbind) %dopar% { # select stable features
  tmp.cell=as.character(unique(mm$Protein)[i])
  tmp.data=mm[mm$Protein==tmp.cell,]
  
  freq=data.frame(table(tmp.data$MultiVariate))
  freq=freq[freq$Freq>8,]
  tmp.data=tmp.data[tmp.data$MultiVariate %in% freq$Var1,]
  
  return.string=data.frame(tmp.data)
}
