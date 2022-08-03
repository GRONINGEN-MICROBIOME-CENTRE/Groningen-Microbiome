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

annotation=read.table("annotation.file.txt",sep = "\t",stringsAsFactors = F,header = T)
gene=gene[,colnames(gene) %in% annotation$Gene]
annotation=annotation[annotation$Gene %in% colnames(gene),]
gene=gene[,order(colnames(gene))]
annotation=annotation[order(annotation$Gene),]
colnames(gene)=annotation$id
gene=as.data.frame(t(gene))

deconvolution=xCellAnalysis(gene)
deconvolution=as.data.frame(deconvolution)
deconvolution=as.data.frame(t(deconvolution))
write.table(deconvolution,file = "OutputTable/Deconvolution.txt",sep = "\t",row.names = T,quote = F)


# ====================================================================================================
# start analysis
# ====================================================================================================
covariate_rna=read.table("Covariate.rna.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
deconvolution=read.table("OutputTable/Deconvolution.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
rownames(deconvolution)=gsub("X","",rownames(deconvolution))
deconvolution_sub=deconvolution[,c("cDC",	"Macrophages M1",	"NK cells",	"pDC",	"Macrophages M2",
                                   "CD4+ naive T-cells",	"CD4+ Tcm",	"CD8+ naive T-cells",	"CD8+ Tcm",	
                                   "Tgd cells",	"Th2 cells",	"Tregs",	"Th1 cells",
                                   "NKT",	"CD8+ Tem",	"CD4+ Tem",	"Class-switched memory B-cells",
                                   "Plasma cells",	"naive B-cells",	"Memory B-cells",	"Basophils",
                                   "Mast cells",	"Neutrophils",	"Eosinophils",	"Endothelial cells",
                                   "Epithelial cells",	"Fibroblasts"),drop=F]
deconvolution_sub=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna),]
covariate_rna=covariate_rna[rownames(covariate_rna) %in% rownames(deconvolution_sub),]
covariate_rna$Inflammation[covariate_rna$Inflammation=="Light"]="No"

covariate_rna$Inflammation[covariate_rna$Inflammation=="Yes"]=3
covariate_rna$Inflammation[covariate_rna$Inflammation=="No" & covariate_rna$Cohort=="Control"]=1
covariate_rna$Inflammation[covariate_rna$Inflammation=="No"]=2
covariate_rna$BMI[is.na(covariate_rna$BMI)]=median(covariate_rna$BMI[!is.na(covariate_rna$BMI)])

# PCA
beta_diversity=vegdist((deconvolution_sub),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,covariate_rna,all = F,by="row.names")

pca_analysis[pca_analysis==""]=NA
ggplot (pca_analysis, aes(V1,V2,fill=Diagnosis,alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  stat_ellipse(type = "norm")+
  scale_fill_nejm()+
  theme(legend.position = 'top')+guides(alpha=F)
ggsave("OutputPlot/Inflammation.diagnosis.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=Inflammation,alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  scale_fill_aaas()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')+guides(alpha=F)
ggsave("OutputPlot/Inflammation.cell.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=Location_rough,alpha=0.5)) + 
  geom_point(shape = 21, size = 2)+ theme_bw() +
  guides(size=F)+
  scale_fill_simpsons()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')+guides(alpha=F)
ggsave("OutputPlot/Location.cell.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=Batch,alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  scale_fill_nejm()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')
ggsave("OutputPlot/Location.batch.PCA.pdf",width = 3,height = 3)

# cluster
annotation=covariate_rna[,c("Diagnosis","Inflammation","Location_rough","resec_ileocec")]
annotation$resec_ileocec=as.character(annotation$resec_ileocec)
annotation$Inflammation[annotation$Inflammation=="Light"]="No"
heatmap.data=as.data.frame((deconvolution))
heatmap.data=heatmap.data[,!colnames(heatmap.data) %in% "MPP"]
heatmap.data=heatmap.data[rownames(heatmap.data) %in% rownames(annotation),]
my_colour = list(Diagnosis = c(Control = "#4EB265", CD = "#E8601C", UC="grey"),Inflammation=c(Yes="#BB5566",No="#004488"),
  Location_rough=c(ileum="#AAAA00",colon="lightblue"),resec_ileocec=c("1"="darkred","0"="#DDDDDD"))

pdf("OutputPlot/heatmap.cell.pdf",width = 20,height = 15)
pheatmap((heatmap.data),cluster_cols = T, cluster_rows = T,scale = "column",
         show_rownames=F, show_colnames=T,
         annotation_row  =annotation,annotation_colors = my_colour,
         cellheight = 0.5,cellwidth = 2,fontsize_number=5,fontsize_row=10,fontsize_col = 3,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

# cell type basic correction, + Th1/Th2 cell ratio
basic_factors=covariate_rna[,c("age_at_biopsy","sex","BMI","Aminosalicylates","Thiopurines","Steroids","Batch","Biological_use")]
str(basic_factors)
which(is.na(basic_factors))
which(is.na(deconvolution_sub))

deconvolution_sub=deconvolution_sub[order(rownames(deconvolution_sub)),]
basic_factors=basic_factors[order(rownames(basic_factors)),]
stopifnot(all(rownames(deconvolution_sub) == rownames(basic_factors)))

deconvolution_sub = apply(deconvolution_sub,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="cell"
  
  fit=lm(cell~age_at_biopsy+sex+BMI+Aminosalicylates+Thiopurines+Steroids+Batch+Biological_use,data=subset.data)
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
covariate1=covariate1[covariate1$Inflammation=="3",]
covariate1=covariate1[,c("Diagnosis","Inflammation","Location_rough")]
which(is.na(covariate1))
deconvolution_sub1=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate1),]
deconvolution_sub1=deconvolution_sub1[order(rownames(deconvolution_sub1)),]
covariate1=covariate1[order(rownames(covariate1)),]

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
write.table(result0,file = "OutputTable/Cell.CDi.UC.txt",sep = "\t",quote = F,row.names = F)

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

ggplot (compare_all,aes(x= Cell, y=Score, fill=Group))+
  geom_boxplot(alpha=0.8,outlier.shape = NA)+
  geom_point(aes(color = Group),alpha=0.1,position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5),size=1)+
  theme_bw()+
  scale_color_manual(values=c("#BB5566","#004488"))+
  scale_fill_manual(values=c("#BB5566","#004488"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab("CellTypes")+ylab("EnrichmentScore")+
  coord_flip()+guides(fill=F)+guides(color=F)
ggsave("OutputPlot/Cell.compare.pdf",width = 8,height = 8)

# ====================================================================================================
# cell-host, use relative abundance instead of 0/1
# ====================================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
basic_factors=covariate_rna[,c("Inflammation","Location_rough")]
basic_factors$Inflammation=as.numeric(basic_factors$Inflammation)

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

bac="Subdoligranulum"
cell="Eosinophils"
tmp.data=merge(deconvolution_sub[,cell,drop=F],bacteria[,bac,drop=F],all=F,by="row.names")

ggplot(tmp.data, aes(tmp.data[,bac], (tmp.data[,cell]))) +
  geom_point(shape = 21, 
             color = "black", size = 2)+
  geom_smooth(method = lm)+
  scale_color_jama()+
  theme_bw()+ theme(legend.position="bottom")+guides(fill=F)

# ====================================== ########## ============================================
# ====================================== predictors ============================================
# ====================================== ########## ============================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)

library (hdi)
library(tidyverse)
library (dplyr)
library(glmnet)
library(glmnetUtils)
library(crayon)
library(caret)
library(pROC)
library(plyr)
library(readr)
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
# models
# model 1: cell ~ factor_group1: basic factors (age, sex, BMI, batch)
# model 2: cell ~ factor_group2: medication (Aminosalicylates, Thiopurines, Steroids, Biological_use)
# model 3: cell ~ factor_group3: Inflammation 
# model 4: cell ~ factor_group4: location
# model 5: cell ~ factor_group5: bacteria
# model 6: cell ~ factor_group6: full model with all factors
# ====================================================================================

deconvolution_sub.correct=deconvolution[,c("cDC",	"Macrophages M1",	"NK cells",	"pDC",	"Macrophages M2",
                                           "CD4+ naive T-cells",	"CD4+ Tcm",	"CD8+ naive T-cells",	"CD8+ Tcm",	
                                           "Tgd cells",	"Th2 cells",	"Tregs",	"Th1 cells",
                                           "NKT",	"CD8+ Tem",	"CD4+ Tem",	"Class-switched memory B-cells",
                                           "Plasma cells",	"naive B-cells",	"Memory B-cells",	"Basophils",
                                           "Mast cells",	"Neutrophils",	"Eosinophils",	"Endothelial cells",
                                           "Epithelial cells",	"Fibroblasts"),drop=F]
deconvolution_sub.correct=deconvolution_sub.correct[rownames(deconvolution_sub.correct) %in% rownames(covariate_rna),]
basic_factors=covariate_rna[,c("age_at_biopsy","sex","BMI","Aminosalicylates","Thiopurines","Steroids","Batch","Biological_use","Inflammation","Location_rough")]
basic_factors$Inflammation=as.numeric(basic_factors$Inflammation)
str(basic_factors)
which(is.na(basic_factors))
basic_factors$Batch[basic_factors$Batch=="batch1"]=1
basic_factors$Batch[basic_factors$Batch=="batch2"]=2
basic_factors$Location_rough[basic_factors$Location_rough=="colon"]=1
basic_factors$Location_rough[basic_factors$Location_rough=="ileum"]=0
basic_factors$Location_rough=as.numeric(basic_factors$Location_rough)
basic_factors$Batch=as.numeric(basic_factors$Batch)

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

model6_select= foreach(i=1:length(unique(mm_select$Protein)),.combine = rbind) %do%  {
  
  tmp.pro=as.character(unique(mm_select$Protein)[i])
  tmp.var=mm_select[mm_select$Protein==tmp.pro,]
  
  tmp.summary=data.frame(Fold=c(1:10),Bacteria=tmp.pro,Variation=NA)
  for(j in 1:10){
    tmp.var.fold=tmp.var[tmp.var$Fold==j,]
    
    if(nrow(tmp.var.fold)==0){
      tmp.summary$Variation[j]=0
    }else{
      tmp.summary$Variation[j]=sum(tmp.var.fold$MultiVariate.explain)
    }
    
  }
  cat(red(tmp.pro,"done","\n"))
  
  return.string = data.frame(Fold=tmp.summary$Fold,Cell=tmp.pro,Variation=tmp.summary$Variation)
}

model6_summary= foreach(i=1:length(unique(model6_select$Cell)),.combine = rbind) %do%  {
  tmp.cell=as.character(unique(model6_select$Cell)[i])
  tmp.data=model6_select[model6_select$Cell==tmp.cell,]
  tmp.variation=mean(tmp.data$Variation)
  
  return.string=data.frame(Models="model6",Cell=tmp.cell,Factors="full",Variation=tmp.variation)
}

model_summary=rbind(model1_summary,model2_summary,model3_summary,model4_summary,model5_summary)
model_summary=as.data.frame(acast(model_summary, Cell~Factors, value.var="Variation"))
model_summary=as.data.frame(t(model_summary))
model_summary[is.na(model_summary)]=0

model6_summary=model6_summary[order(model6_summary$Variation,decreasing = T),]
orders=c(as.character(model6_summary$Cell))
model_summary=model_summary[,orders]

paletteLength <- 20
myColor <- colorRampPalette(c("white","#993404"))(paletteLength)

pdf("OutputPlot//Cell.variation.heatmap.pdf",width = 8,height = 3)
pheatmap(model_summary,cluster_cols = F, cluster_rows = F,show_rownames=T, show_colnames=T,
         cellheight = 15,cellwidth = 15,fontsize_number=12,fontsize_row=8,fontsize_col = 8,
         color = myColor)
dev.off()

model6_select$Cell=factor(model6_select$Cell,levels = orders)
ggplot(model6_select, aes(x=Cell, y=Variation,fill=Cell,na.rm = F)) +
  geom_boxplot(width=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  guides(fill=F)
ggsave("OutputPlot/Cell.variation.pdf",width = 8,height = 4)

model5_summary=model5_summary[order(model5_summary$Variation,decreasing = T),]
model5_select$Cell=factor(model5_select$Cell,levels = model5_summary$Cell)
ggplot(model5_select, aes(x=Cell, y=Variation,fill=Cell,na.rm = F)) +
  geom_boxplot(width=0.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  guides(fill=F)
ggsave("OutputPlot/Cell.variation.bacteria.pdf",width = 7,height = 3)

tmp.cell="NK cells"
tmp.data=mm_select[mm_select$Protein==tmp.cell,]
tmp.data$MultiVariate=as.factor(tmp.data$MultiVariate)
ggplot(tmp.data, aes(x=reorder(MultiVariate, MultiVariate.explain, FUN = median), y=MultiVariate.explain,fill=MultiVariate,na.rm = F)) +
  geom_boxplot(width=0.3)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  guides(fill=F)+ coord_flip()+ylab("Variation")+xlab("")
ggsave("OutputPlot/NK cells.bacteria.pdf",width = 3,height = 5)





















var_factor=mm_select[mm_select$MultiVariate %in% colnames(basic_factors),]
var_bacteria=mm_select[!mm_select$MultiVariate %in% tmp.data1$MultiVariate,]
tmp <- foreach(i=1:length(unique(var_bacteria$Protein)), .combine = rbind) %dopar% {

  tmp.cell=as.character(unique(var_bacteria$Protein))[i]
  tmp.data=var_bacteria[var_bacteria$Protein==tmp.cell,]
  freq=data.frame(table(tmp.data$MultiVariate))
  tmp.data=tmp.data[tmp.data$MultiVariate %in% freq$Var1[freq$Freq>8],]
  
  if(nrow(tmp.data)==0){
    cat(blue(tmp.cell,"is empty","\n"))
    return.string=data.frame(cell=tmp.cell,fold=0,var=0)
  }else{
    tmp.count=data.frame(cell=tmp.cell,fold=unique(tmp.data$Fold),var=0)
    for(n in unique(tmp.data$Fold)){
      tmp.var=tmp.data[tmp.data$Fold==n,]
      count_bacteria=nrow(tmp.var)
      var_sum=sum(tmp.var$MultiVariate.explain)
      
      tmp.count$var[tmp.count$fold==n]=var_sum
    }
    cat(yellow(tmp.cell,"is done","\n")) 
    
    return.string=as.data.frame(tmp.count)
    }
}
var_bacteria=tmp
ggplot(var_bacteria, aes(x=cell, y=var,fill=cell,na.rm = F)) +
  geom_boxplot(width=0.1, fill="white")+
  theme_bw()+
  scale_fill_jama()+
  theme(axis.text.x = element_text(angle = 90))+
  guides(fill=F)

var_bacteria$category="Bacteria"
var_factor=var_factor[,c("Protein","Fold","MultiVariate.explain","MultiVariate")]
colnames(var_factor)=c("cell","fold","var","category")
nn=nn[,c("Cell","Fold","Variation")]
nn$category="Total"
colnames(nn)=c("cell","fold","var","category")
orders=data.frame(cell=unique(nn$cell),mean=0)
for(i in unique(nn$cell)){
  tmp.data=nn[nn$cell==i,]
  mean=mean(tmp.data$var)
  orders$mean[orders$cell==i]=mean
}
orders=orders[order(orders$mean,decreasing = T),]

var_all=rbind(var_bacteria,var_factor,nn)
var_all$cell=factor(var_all$cell,levels = orders$cell)
var_all=var_all[var_all$category %in% c("Bacteria","Inflammation","Location_rough","Total"),]
ggplot(var_all, aes(x=cell, y=var,fill=cell,na.rm = F)) +
  geom_boxplot(width=0.1)+
  theme_bw()+
  facet_grid(category ~ .) +
  theme(axis.text.x = element_text(angle = 45,vjust = 0.5,size = 10))+
  theme(strip.text.y.right = element_text(angle = 360))+
  guides(fill=F)
ggsave("OutputPlot/cell.variation.pdf",width = 15,height = 8)
write.table(var_all,file = "OutputTable/Cell.variation.txt",sep = "\t",quote = F,row.names = F)
























































