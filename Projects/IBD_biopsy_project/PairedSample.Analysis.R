# ====================================================================================
# paired sample analysis
# ====================================================================================
library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
library(RColorBrewer)
library(randomcoloR)
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
library("ape")
library(dendextend)

bacteria=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)

pair_tissue=read.table("OutputTable/Paired.multipleLcation.noninflamed.txt",sep = "\t",header = T,stringsAsFactors = F)
pair_inflammation_colon=read.table("OutputTable/Paired.multipleInflammation.colon.txt",sep = "\t",header = T,stringsAsFactors = F)
pair_inflammation_ileum=read.table("OutputTable/Paired.multipleInflammation.Ileum.txt",sep = "\t",header = T,stringsAsFactors = F)

# individual difference, vs. tissue difference
bacteria_tissue=bacteria[rownames(bacteria) %in% pair_tissue$BiopsyID,,]
beta_diversity=vegdist((bacteria_tissue),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=2))
pca_analysis=merge(pca_analysis,covariate_bac,by="row.names",all=F)

dd <- dist(scale(bacteria_tissue), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
dend <- as.dendrogram(hc)
groupCodes=pca_analysis$ResearchID
colorCodes <- distinctColorPalette(368)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

distmat <- (as.matrix(beta_diversity))
distmat=data.frame(as.table(distmat))[lower.tri(distmat, diag = F), ]
distmat=merge(distmat,pair_tissue[,c("ResearchID","BiopsyID")],by.x="Var1",by.y="BiopsyID")
distmat=as.data.frame(distmat)
colnames(distmat)[4]="Var1ID"
distmat=merge(distmat,pair_tissue[,c("ResearchID","BiopsyID")],by.x="Var2",by.y="BiopsyID")
colnames(distmat)[5]="Var2ID"
distmat=merge(distmat,covariate_bac[,"Location_rough",drop=F],by.x="Var1",by.y="row.names",all=F)
colnames(distmat)[6]="Var1Group"
distmat=merge(distmat,covariate_bac[,"Location_rough",drop=F],by.x="Var2",by.y="row.names",all=F)
colnames(distmat)[7]="Var2Group"

distmat$group=NA
distmat$group[distmat$Var1ID==distmat$Var2ID]="WithinIndividual"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="colon" & distmat$Var2Group=="colon"]="BetweenIndividual.colon"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="ileum" & distmat$Var2Group=="ileum"]="BetweenIndividual.ileum"
distmat=na.omit(distmat)

wilcox.test(distmat$Freq[distmat$group=="WithinIndividual"],distmat$Freq[distmat$group=="BetweenIndividual.colon"])
wilcox.test(distmat$Freq[distmat$group=="WithinIndividual"],distmat$Freq[distmat$group=="BetweenIndividual.ileum"])
wilcox.test(distmat$Freq[distmat$group=="BetweenIndividual.colon"],distmat$Freq[distmat$group=="BetweenIndividual.ileum"])

# individual difference, vs. inflammation difference (colon)
bacteria_colon=bacteria[rownames(bacteria) %in% pair_inflammation_colon$BiopsyID,,]
beta_diversity=vegdist((bacteria_colon),method = "euclidean")

distmat <- (as.matrix(beta_diversity))
distmat=data.frame(as.table(distmat))[lower.tri(distmat, diag = F), ]
distmat=merge(distmat,pair_inflammation_colon[,c("ResearchID","BiopsyID")],by.x="Var1",by.y="BiopsyID")
distmat=as.data.frame(distmat)
colnames(distmat)[4]="Var1ID"
distmat=merge(distmat,pair_inflammation_colon[,c("ResearchID","BiopsyID")],by.x="Var2",by.y="BiopsyID")
colnames(distmat)[5]="Var2ID"
distmat=merge(distmat,covariate_bac[,"Inflammation",drop=F],by.x="Var1",by.y="row.names",all=F)
colnames(distmat)[6]="Var1Group"
distmat=merge(distmat,covariate_bac[,"Inflammation",drop=F],by.x="Var2",by.y="row.names",all=F)
colnames(distmat)[7]="Var2Group"

distmat$group=NA
distmat$group[distmat$Var1ID==distmat$Var2ID]="WithinIndividual"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="No" & distmat$Var2Group=="No"]="BetweenIndividual.nonInf"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="Yes" & distmat$Var2Group=="Yes"]="BetweenIndividual.Inf"
distmat=na.omit(distmat)

# individual difference, vs. inflammation difference (ileum)
bacteria_ileum=bacteria[rownames(bacteria) %in% pair_inflammation_ileum$BiopsyID,,]
beta_diversity=vegdist((bacteria_ileum),method = "euclidean")

distmat <- (as.matrix(beta_diversity))
distmat=data.frame(as.table(distmat))[lower.tri(distmat, diag = F), ]
distmat=merge(distmat,pair_inflammation_ileum[,c("ResearchID","BiopsyID")],by.x="Var1",by.y="BiopsyID")
distmat=as.data.frame(distmat)
colnames(distmat)[4]="Var1ID"
distmat=merge(distmat,pair_inflammation_ileum[,c("ResearchID","BiopsyID")],by.x="Var2",by.y="BiopsyID")
colnames(distmat)[5]="Var2ID"
distmat=merge(distmat,covariate_bac[,"Inflammation",drop=F],by.x="Var1",by.y="row.names",all=F)
colnames(distmat)[6]="Var1Group"
distmat=merge(distmat,covariate_bac[,"Inflammation",drop=F],by.x="Var2",by.y="row.names",all=F)
colnames(distmat)[7]="Var2Group"

distmat$group=NA
distmat$group[distmat$Var1ID==distmat$Var2ID]="WithinIndividual"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="No" & distmat$Var2Group=="No"]="BetweenIndividual.nonInf"
distmat$group[distmat$Var1ID!=distmat$Var2ID & distmat$Var1Group=="Yes" & distmat$Var2Group=="Yes"]="BetweenIndividual.Inf"
distmat=na.omit(distmat)

# tissue location, colon vs ileum, in noninflamed tissue
tmp1=covariate_bac[rownames(covariate_bac) %in% pair_tissue$BiopsyID,]
tmp1$ID=paste(tmp1$Location_rough,tmp1$Inflammation,tmp1$ResearchID)
tmp1=tmp1[!duplicated(tmp1$ID),]
tmp1$ResearchID=as.character(tmp1$ResearchID)
count=as.data.frame(table(tmp1$ResearchID))
tmp1=tmp1[tmp1$ResearchID %in% count$Var1[count$Freq==2],]
shannon_tmp=shannon[shannon$Row.names %in% rownames(tmp1),]
ggplot(shannon_tmp, aes(x=Location_rough, y=shannon, fill=Location_rough)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_npg()+
  theme(legend.position = 'top')+
  xlab("")+geom_line(aes(group = ResearchID), alpha = 0.1, colour = "black", data = shannon_tmp)+theme_bw()
ggsave("OutputPlot/Paired.shannon.location.pdf",width = 3,height = 2)
id=unique(as.character(shannon_tmp$ResearchID))
tmp.data = foreach(i=1:length(id),.combine = rbind) %do%  {
  tmp.id=id[i]
  tmp.sub=shannon_tmp[shannon_tmp$ResearchID==tmp.id,]
  mm=tmp.sub$shannon[tmp.sub$Location_rough=="colon"]
  nn=tmp.sub$shannon[tmp.sub$Location_rough=="ileum"]
  
  return.string=data.frame(ID=tmp.id,colon=mm,ileum=nn)

  }
wilcox.test(tmp.data$colon, tmp.data$ileum, paired = TRUE)

covariate1=covariate_bac[rownames(covariate_bac) %in% rownames(tmp1),]
bacteria1=bacteria_clr[rownames(bacteria_clr) %in% rownames(tmp1),]
bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result1 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(covariate1[,c("ResearchID","Location_rough"),drop=F],x,by="row.names",all=F)
  tmp.data$Row.names=NULL
  tmp.data=dcast(tmp.data,ResearchID~Location_rough)


  mm=wilcox.test(tmp.data$colon, tmp.data$ileum, paired = TRUE)
  return.string=data.frame(Bacteria=tmp.genus,PairedPvalue=mm$p.value,MeanColon=mean(tmp.data$colon),MeanIleum=mean(tmp.data$ileum))
  
}
result1$FDR=p.adjust(result1$PairedPvalue,method = "BH")

# inflamed vs non-inflamedd, in ileum
tmp1=covariate_bac[rownames(covariate_bac) %in% pair_inflammation_ileum$BiopsyID,]
tmp1$ID=paste(tmp1$Location_rough,tmp1$Inflammation,tmp1$ResearchID)
tmp1=tmp1[!duplicated(tmp1$ID),]
tmp1$ResearchID=as.character(tmp1$ResearchID)
count=as.data.frame(table(tmp1$ResearchID))
tmp1=tmp1[tmp1$ResearchID %in% count$Var1[count$Freq==2],]
shannon_tmp=shannon[shannon$Row.names %in% rownames(tmp1),]
shannon_tmp$Inflammation=as.factor(shannon_tmp$Inflammation)
ggplot(shannon_tmp, aes(x=Inflammation, y=shannon, fill=Inflammation)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_npg()+
  theme(legend.position = 'top')+
  xlab("")+geom_line(aes(group = ResearchID), alpha = 0.1, colour = "black", data = shannon_tmp)+theme_bw()
ggsave("OutputPlot/Paired.shannon.inflammation.iluem.pdf",width = 3,height = 2)
id=unique(as.character(shannon_tmp$ResearchID))
tmp.data = foreach(i=1:length(id),.combine = rbind) %do%  {
  tmp.id=id[i]
  tmp.sub=shannon_tmp[shannon_tmp$ResearchID==tmp.id,]
  mm=tmp.sub$shannon[tmp.sub$Inflammation=="2"]
  nn=tmp.sub$shannon[tmp.sub$Inflammation=="3"]
  
  return.string=data.frame(ID=tmp.id,noninf=mm,inf=nn)
  
}
wilcox.test(tmp.data$noninf, tmp.data$inf, paired = TRUE)

covariate1=covariate_bac[rownames(covariate_bac) %in% rownames(tmp1),]
bacteria1=bacteria_clr[rownames(bacteria_clr) %in% rownames(tmp1),]
bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result2 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(covariate1[,c("ResearchID","Inflammation"),drop=F],x,by="row.names",all=F)
  tmp.data$Row.names=NULL
  tmp.data=dcast(tmp.data,ResearchID~Inflammation)
  
  
  mm=wilcox.test(tmp.data$`2`, tmp.data$`3`, paired = TRUE)
  return.string=data.frame(Bacteria=tmp.genus,PairedPvalue=mm$p.value,MeanNonInf=mean(tmp.data$`2`),MeanInf=mean(tmp.data$`3`))
  
}
result2$FDR=p.adjust(result2$PairedPvalue,method = "BH")

# inflamed vs non-inflamedd, in colon
tmp1=covariate_bac[rownames(covariate_bac) %in% pair_inflammation_colon$BiopsyID,]
tmp1$ID=paste(tmp1$Location_rough,tmp1$Inflammation,tmp1$ResearchID)
tmp1=tmp1[!duplicated(tmp1$ID),]
tmp1$ResearchID=as.character(tmp1$ResearchID)
count=as.data.frame(table(tmp1$ResearchID))
tmp1=tmp1[tmp1$ResearchID %in% count$Var1[count$Freq==2],]
shannon_tmp=shannon[shannon$Row.names %in% rownames(tmp1),]
shannon_tmp$Inflammation=as.factor(shannon_tmp$Inflammation)
ggplot(shannon_tmp, aes(x=Inflammation, y=shannon, fill=Inflammation)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_npg()+
  theme(legend.position = 'top')+
  xlab("")+geom_line(aes(group = ResearchID), alpha = 0.1, colour = "black", data = shannon_tmp)+theme_bw()
ggsave("OutputPlot/Paired.shannon.inflammation.location.pdf",width = 3,height = 2)
id=unique(as.character(shannon_tmp$ResearchID))
tmp.data = foreach(i=1:length(id),.combine = rbind) %do%  {
  tmp.id=id[i]
  tmp.sub=shannon_tmp[shannon_tmp$ResearchID==tmp.id,]
  mm=tmp.sub$shannon[tmp.sub$Inflammation=="2"]
  nn=tmp.sub$shannon[tmp.sub$Inflammation=="3"]
  
  return.string=data.frame(ID=tmp.id,noninf=mm,inf=nn)
  
}
wilcox.test(tmp.data$noninf, tmp.data$inf, paired = TRUE)

covariate1=covariate_bac[rownames(covariate_bac) %in% rownames(tmp1),]
bacteria1=bacteria_clr[rownames(bacteria_clr) %in% rownames(tmp1),]
bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result3 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(covariate1[,c("ResearchID","Inflammation"),drop=F],x,by="row.names",all=F)
  tmp.data$Row.names=NULL
  tmp.data=dcast(tmp.data,ResearchID~Inflammation)
  
  
  mm=wilcox.test(tmp.data$`2`, tmp.data$`3`, paired = TRUE)
  return.string=data.frame(Bacteria=tmp.genus,PairedPvalue=mm$p.value,MeanNonInf=mean(tmp.data$`2`),MeanInf=mean(tmp.data$`3`))
  
}
result3$FDR=p.adjust(result3$PairedPvalue,method = "BH")

