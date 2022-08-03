# dysbiosis analysis
library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(gridExtra)
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
library(jcolors)
library(crayon)

covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]
genus=genus[!rownames(genus) %in% remove,]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(genus),]
genus=genus[rownames(genus) %in% rownames(covariate_bac),]
bacteria=read.table("OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]

# PCA plot
beta_diversity=vegdist((genus),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,covariate_bac,all = F,by="row.names")
pca_analysis[pca_analysis==""]=NA
pca_analysis=pca_analysis[!is.na(pca_analysis$Diagnosis),]
pca_analysis$resec_ileocec=as.factor(pca_analysis$resec_ileocec)
ggplot (pca_analysis, aes(V1,V2,color=Diagnosis)) + 
  geom_point() + theme_bw() +
  guides(size=F)+
  scale_color_npg()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")
ggsave("OutputPlot//PCA.bac.diagnosis.pdf",width = 8,height = 5)

# calculate median distance score as dysbiosis
beta_diversity=vegdist(genus,method = "euclidean")
distmat <- (as.matrix(beta_diversity))
distmat=data.frame(as.table(distmat))
distmat_IBD=distmat[distmat$Var1 %in% rownames(covariate_bac)[covariate_bac$Cohort=="1000IBD"],]
distmat_IBD=distmat_IBD[distmat_IBD$Var2 %in% rownames(covariate_bac)[covariate_bac$Cohort=="Control"],]
distmat_control=distmat[distmat$Var1 %in% rownames(covariate_bac)[covariate_bac$Cohort=="Control"],]
distmat_control=distmat_control[distmat_control$Var2 %in% rownames(covariate_bac)[covariate_bac$Cohort=="Control"],]
distmat_control=distmat_control[distmat_control$Freq!=0,]

sample_IBD=as.character(unique(distmat_IBD$Var1))
dysbiosis_IBD = foreach(i=1:length(sample_IBD),.combine = rbind) %do%  {
  tmp.sample=sample_IBD[i]
  tmp.cal=distmat_IBD[distmat_IBD$Var1==tmp.sample,]
  tmp.median=median(tmp.cal$Freq)
  
  return.string=data.frame(ID=tmp.sample,Dysbiosis=tmp.median,Group="IBD")
  }
sample_control=as.character(unique(distmat_control$Var1))
dysbiosis_control = foreach(i=1:length(sample_control),.combine = rbind) %do%  {
  tmp.sample=sample_control[i]
  tmp.cal=distmat_control[distmat_control$Var1==tmp.sample,]
  tmp.median=median(tmp.cal$Freq)
  
  return.string=data.frame(ID=tmp.sample,Dysbiosis=tmp.median,Group="Control")
}
dysbiosis_all=rbind(dysbiosis_control,dysbiosis_IBD)

pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,dysbiosis_all,all = F,by.x="row.names",by.y="ID")

ggplot (pca_analysis, aes(V1,V2,color=Dysbiosis)) + 
  geom_point() + theme_bw() +
  guides(size=F)+scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+
  #theme(legend.position = 'top')+
  xlab("PCA1")+ylab("PCA2")+theme_classic()
ggsave("OutputPlot/Dysbiosis.PCA.pdf",width = 3,height = 2)
write.table(dysbiosis_all,file = "OutputTable/Dysbiosis.score.all.txt",sep = "\t",row.names = F,quote = F)

# dysbiosis vs. disease group
dysbiosis_all_group=dysbiosis_all
dysbiosis_all_group$Group=as.character(dysbiosis_all_group$Group)
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"]]="CD"
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"]]="UC"
dysbiosis_all_group=dysbiosis_all_group[dysbiosis_all_group$Group!="IBD",]
ggplot(dysbiosis_all_group, aes(x=Dysbiosis)) + geom_density(aes(group=Group, fill=Group), alpha=0.5) +
  scale_fill_manual(values=c("#EE8866","#44BB99","#AA4499"))+theme_classic()
ggsave("OutputPlot/Dysbiosis.density.pdf",width = 3,height = 2)
wilcox.test(dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="Control"],dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="CD"])
wilcox.test(dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="Control"],dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="UC"])

# check tech parameters
reads_count=read.table("DadaOutput/Reads.count.txt",sep = "\t",header = T,stringsAsFactors = F)
covariate_all=read.table("MetaData/Metadaata.16S.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
coupling=covariate_all[,c(2,4)]
reads_count=merge(reads_count,coupling,by.x="SampleID",by.y="16SSampleID",all=F)
pca_analysis=merge(pca_analysis,reads_count,by.x="Row.names",by.y="biopsy_number",all=F)
pca_analysis=merge(pca_analysis,covariate_bac,by.x="Row.names",by.y="row.names",all=F)
m=ggplot (pca_analysis, aes(V1,V2,color=log(AlignedReads))) + 
  geom_point() + theme_bw() +
  guides(size=F)+scale_color_jcolors_contin("pal2")+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")
ggMarginal(m)
ggsave("OutputPlot/Dysbiosis.PCA.reads.pdf",width = 5,height = 5)

# covariate
pca_analysis$Inflammation[pca_analysis$Inflammation=="Yes"]=1
pca_analysis$Inflammation[pca_analysis$Inflammation=="No"]=0
pca_analysis$Location_rough[pca_analysis$Location_rough=="ileum"]=0
pca_analysis$Location_rough[pca_analysis$Location_rough=="colon"]=1

tmp.covariate=c("AlignedReads","Inflammation","Location_rough","age_at_biopsy","sex","BMI","resec_part_colon","resec_ileocec","resec_part_small","SCCAI","HBI")
dys_covariate = foreach(i=1:length(tmp.covariate),.combine = rbind) %do%  {
  tmp.factor=tmp.covariate[i]
  tmp.data=cbind(pca_analysis[,tmp.factor,drop=F],pca_analysis[,"Dysbiosis",drop=F])
  tmp.data=na.omit(tmp.data)
  tmp.data[,1]=as.numeric(tmp.data[,1])
  
  mm=cor.test(tmp.data[,1],tmp.data[,2],method = "spearman")
  
  return.string=data.frame(factor=tmp.factor,cor=mm$estimate,p=mm$p.value)

}
pca_analysis$resec_ileocec=as.factor(pca_analysis$resec_ileocec)
pca_analysis$resec_part_small=as.factor(pca_analysis$resec_part_small)
pca_analysis$resec_part_colon=as.factor(pca_analysis$resec_part_colon)
pca_analysis$sex=as.factor(pca_analysis$sex)
pca_analysis$Inflammation=as.factor(pca_analysis$Inflammation)
pca_analysis$Inflammation[pca_analysis$Inflammation==""]=0
pca_analysis$Location_rough=as.factor(pca_analysis$Location_rough)
p=ggplot(pca_analysis) +
  geom_point(aes(V1,V2,fill=Location_rough),shape = 21, color = "black", size = 2)+
  theme(legend.position="bottom")+
  scale_fill_manual(values=c("#EECC66", "#0077BB"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
y_box = axis_canvas(p, axis = "y", coord_flip = TRUE) +
  geom_density(data = pca_analysis, aes(x=V1,fill=Location_rough),alpha=0.5) +
  scale_fill_manual(values = c("#EECC66", "#0077BB"))+
  coord_flip()
combined_plot <- insert_yaxis_grob(p, y_box, position = "right")
x_box = axis_canvas(p, axis = "x") +
  scale_fill_manual(values = c("#EECC66", "#0077BB"))+
  geom_density(data = pca_analysis, aes(x=V2,fill=Location_rough),alpha=0.5)
combined_plot <- insert_xaxis_grob(combined_plot, x_box, position = "top")
ggdraw(combined_plot)
pdf("OutputPlot/dys_tmp.pdf",height = 4,width = 5)
ggdraw(combined_plot)
dev.off()

# dysbiosis vs. bacteria
dys_genus = foreach(i=1:ncol(bacteria),.combine = rbind) %do%  {
  tmp.genus=bacteria[,i,drop=F]
  tmp.data=merge(dysbiosis_all,tmp.genus,by.x="ID",by.y="row.names",all=F)
  tmp.test=cor.test(tmp.data$Dysbiosis,tmp.data[,4],method = "spearman")
  
  return.string=data.frame(Bacteria=colnames(bacteria)[i],Cor=tmp.test$estimate,Pvalue=tmp.test$p.value)
}
dys_genus$FDR=p.adjust(dys_genus$Pvalue,method = "BH")
write.table(dys_genus,file = "OutputTable/Dys.bacteria.txt",sep = "\t",row.names = F,quote = F)

tmp.genus="Veillonella"
tmp.data=merge(dysbiosis_all,genus[,tmp.genus,drop=F],by.x="ID",by.y="row.names",all=F)
ggplot(tmp.data, aes(`Veillonella`,Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2,fill="#DDAA33")+
  geom_smooth(method = lm,color="black")+theme_bw()
ggsave("OutputPlot/Dysbiosis.bacteria.pdf",width = 3,height = 3)

# dysbiosis vs. shannon
tmp.genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F)
tmp.genus=tmp.genus[rownames(tmp.genus) %in% rownames(genus),]
#tmp.genus=tmp.genus[,colSums(tmp.genus)>0.05*nrow(tmp.genus)]
tmp.genus=tmp.genus[,colnames(tmp.genus) %in% colnames(genus)]
shannon <-microbiome::alpha(t(tmp.genus), index = "shannon")
colnames(shannon)="shannon"

tmp.data=merge(dysbiosis_all,shannon,by.x="ID",by.y="row.names",all=F)
ggplot(tmp.data, aes(shannon,Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2)+
  geom_smooth(method = lm,color="black")+theme_bw()
cor.test(tmp.data$Dysbiosis,tmp.data$shannon)

# dysbiosis vs. montreal classification
table(covariate_bac$MontrealB)
covariate_bac$MontrealBB[covariate_bac$MontrealB==0 | covariate_bac$MontrealB==3]=0
covariate_bac$MontrealBB[covariate_bac$MontrealB==1 | covariate_bac$MontrealB==4]=1
covariate_bac$MontrealBB[covariate_bac$MontrealB==2 | covariate_bac$MontrealB==5]=2
covariate_bac$perianal[covariate_bac$MontrealB==0 | covariate_bac$MontrealB==1 | covariate_bac$MontrealB==2]=0
covariate_bac$perianal[covariate_bac$MontrealB==3 | covariate_bac$MontrealB==4 | covariate_bac$MontrealB==5]=1

montreal=merge(dysbiosis_all_group,covariate_bac[,c("MontrealA","perianal","MontrealBB","MontrealL","MontrealE","MontrealS")],by.x="ID",by.y="row.names",all=F)
cor.test(montreal$Dysbiosis,montreal$MontrealA)
cor.test(montreal$Dysbiosis,montreal$perianal)
cor.test(montreal$Dysbiosis,montreal$MontrealBB)
kruskal.test(Dysbiosis~MontrealL,montreal)
cor.test(montreal$Dysbiosis,montreal$MontrealE)
cor.test(montreal$Dysbiosis,montreal$MontrealS)

tmp=montreal[,c("Dysbiosis","MontrealL")]
tmp=na.omit(tmp)
tmp$MontrealL=as.factor(tmp$MontrealL)
ggplot(tmp, aes(perianal,Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2)+
  geom_smooth(method = lm,color="black",linetype="dashed")+
  theme(legend.position="bottom")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
ggsave("OutputPlot/dys.MontrealL.pdf",width = 5,height = 5)
ggplot(tmp, aes(x=MontrealL, y=Dysbiosis,fill=MontrealL)) + 
  geom_boxplot(color="black")+scale_fill_npg()+theme_bw()+guides(fill=F)

# dysbiosis vs. gene expression
gene=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
samples_count=intersect(rownames(gene),intersect(rownames(genus),intersect(dysbiosis_all$ID,rownames(covariate_rna))))

gene=gene[rownames(gene) %in% samples_count,]
genus=genus[rownames(genus) %in% samples_count,]
dysbiosis_all=dysbiosis_all[(dysbiosis_all$ID) %in% samples_count,]
covariate_rna=covariate_rna[rownames(covariate_rna) %in% samples_count,]
write.table(covariate_rna,file = "OutputTable/Final.used.samples.16S+RNAseq.txt",sep = "\t",row.names = T,quote = F)

# RNA data correction (age, gender, BMI, bacth, inflammation, location (no diagnosis, location == diagnosis), and medication)
covariate_rna=covariate_rna[,c("Inflammation","Location_rough","age_at_biopsy","sex","BMI","Batch","Aminosalicylates","Thiopurines","Steroids")]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]
gene=gene[order(rownames(gene)),]
genes_correct = apply(gene,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_rna[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="tmp.gene"
  
  fit=lm(tmp.gene~age_at_biopsy+sex+BMI+Batch+Inflammation+Location_rough+Aminosalicylates+Thiopurines+Steroids,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
genes_correct=as.data.frame(genes_correct)

# select inflammed genes
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation_genes=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
genes_correct=genes_correct[,colnames(genes_correct) %in% inflammation_genes]

# select overlapped with HMP2
#hmp2.genus=read.table("OutputTable/1000IBD,HMP2.consistent.genus.txt",sep="\t",header = T)
#genus=genus[,colnames(genus) %in% hmp2.genus$x]

# model: gene ~ shannon 
tmp.genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F)
tmp.genus=tmp.genus[,colSums(tmp.genus)>0]
tmp.genus=tmp.genus[rowSums(tmp.genus)>0,]

tmp.genus=tmp.genus[rownames(tmp.genus) %in% rownames(genus),]
tmp.genus=tmp.genus[,colSums(tmp.genus>0)>0.05*nrow(tmp.genus)]
tmp.genus=tmp.genus[,colnames(tmp.genus) %in% colnames(genus)]
shannon <- as.data.frame(diversity(t(tmp.genus), index="shannon"))
colnames(shannon)="shannon"
modelfit_shannon = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],shannon,by="row.names",all=F)
  
  colnames(tmp.data)=c("ID","Gene","Shannon")
  tmp.model=lm(Gene ~ Shannon,data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
  cat(green(i,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,Bac="Shannon",Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}

# model: gene ~ dys
quantile(dysbiosis_all$Dysbiosis)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=14.92912]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>14.92912]=1
modelfit00 = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=(genes_correct[,tmp.gene,drop=F])
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="row.names",by.x="ID",all=F)
  colnames(tmp.data)=c("ID","Dys","Group","Gene")

  tmp.model=lm(Gene ~ Dys ,data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coefficients)

  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=data.frame(Gene=tmp.gene,Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}

# model: gene ~ bac 
modelfit_taxa = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],bacteria,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta","SE","Zscore","P")
  
  for(n in 1:ncol(bacteria)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=lm(Gene ~ Bac ,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$SE[n]=tmp.model$`Std. Error`[rownames(tmp.model)=="Bac"]
    tmp.result$Zscore[n]=tmp.model$`t value`[rownames(tmp.model)=="Bac"]
    tmp.result$P[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    
  }
  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=as.data.frame(tmp.result)
}

# model: gene ~ bac + dys
quantile(dysbiosis_all$Dysbiosis)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=14.92912]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>14.92912]=1
modelfit0 = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],genus,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(genus),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_dys","P_taxa","P_dys")
  
  for(n in 1:ncol(genus)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=lm(Gene ~ Bac + Dys ,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$Beta_dys[n]=tmp.model$Estimate[rownames(tmp.model)=="Dys"]
    tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    tmp.result$P_dys[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Dys"]
    
  }
  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=as.data.frame(tmp.result)
}

# model: gene ~ bac + dysbiosis + bac * dysbiosis
quantile(dysbiosis_all$Dysbiosis,0.95)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=17.66]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>17.66]=1
modelfit = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],bacteria,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  #tmp.data$Dysbiosis=sample(tmp.data$Dysbiosis)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 8))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_dys","Beta_interact","P_taxa","P_dys","P_interact")
  
  for(n in 1:ncol(bacteria)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=lm(Gene ~ Bac + Dys + Bac * Dys,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$Beta_dys[n]=tmp.model$Estimate[rownames(tmp.model)=="Dys"]
    tmp.result$Beta_interact[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac:Dys"]
    tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    tmp.result$P_dys[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Dys"]
    tmp.result$P_interact[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac:Dys"]
    
  }
  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=as.data.frame(tmp.result)
}

modelfit$FDR_interact=p.adjust(modelfit$P_interact,method = "BH")
#modelfit0$FDR_taxa=p.adjust(modelfit0$P_taxa,method = "BH")
modelfit_shannon$FDR=p.adjust(modelfit_shannon$Pvalue,method = 'BH')
modelfit00$FDR=p.adjust(modelfit00$Pvalue,method = "BH")
modelfit_taxa$FDR=p.adjust(modelfit_taxa$P,method = "BH")

write.table(modelfit_taxa,file = "OutputTable/gene.bacteria.txt",quote = F,sep = "\t",row.names = F)
write.table(modelfit,file = "OutputTable/Dysbiosis.interaction.result.txt",quote = F,sep = "\t",row.names = F)
write.table(modelfit0,file = "OutputTable/Dysbiosis.Nointeraction.result.txt",quote = F,sep = "\t",row.names = F)
write.table(modelfit_shannon,file = "OutputTable/Dysbiosis.shannon.result.txt",quote = F,sep = "\t",row.names = F)

library("ggExtra")
tmp.genus="Lachnospiraceae"
tmp.gene="S100A8"
tmp.data=merge(dysbiosis_quantitle,bacteria[,tmp.genus,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data=merge(tmp.data,genes_correct[,tmp.gene,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data=merge(tmp.data,montreal,by="ID",all=F)
colnames(tmp.data)[2]="Dysbiosis"
tmp.data$Dysbiosis[tmp.data$Dysbiosis==0]="Dysbiosis score [0~95%]"
tmp.data$Dysbiosis[tmp.data$Dysbiosis==1]="Dysbiosis score [95~100%]"
tmp.data$perianal=as.factor(tmp.data$perianal)

ggplot(tmp.data, aes(Lachnospiraceae,S100A8,fill=Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2)+
  geom_smooth(method = lm,color="black",linetype="dashed")+
  theme_classic()+
  facet_wrap(~ Dysbiosis)+
  theme(strip.background = element_blank(), strip.text = element_blank())+scale_fill_manual(values=c("#BB5566", "#004488"))+
  guides(fill=F)
ggsave("OutputPlot/test.pdf",width = 3,height = 2)
ggplot(tmp.data, aes(x=Dysbiosis, y=PLAUR, fill=Dysbiosis)) +theme_bw()+
  geom_boxplot()+theme(axis.title.x=element_blank(),
                          axis.text.x=element_blank(),
                          axis.ticks.x=element_blank())+
  scale_fill_manual(values=c("#BB5566", "#004488"))+guides(fill=F)+xlab("")
ggsave("OutputPlot/test.boxplot.pdf",width = 2,height = 4)
wilcox.test(tmp.data$S100A8[tmp.data$Dysbiosis=="Dysbiosis score [0~75%]"],
            tmp.data$S100A8[tmp.data$Dysbiosis=="Dysbiosis score [75~100%]"])


p=ggplot(tmp.data) +
  geom_point(aes(`Lachnospiraceae`,LCN15,fill=Dysbiosis),shape = 21, color = "black", size = 2)+
  geom_smooth(aes(`Lachnospiraceae`,LCN15),method = lm,color="#DDAA33",linetype="dashed") +
  theme(legend.position="bottom")+scale_fill_manual(values=c("#BB5566", "#004488"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
y_box = axis_canvas(p, axis = "y", coord_flip = TRUE) +
  geom_density(data = tmp.data, aes(x=LCN15,fill=Dysbiosis),alpha=0.5) +
  scale_fill_manual(values = c("#BB5566", "#004488"))+
  coord_flip()
combined_plot <- insert_yaxis_grob(p, y_box, position = "right")
x_box = axis_canvas(p, axis = "x") +
  geom_density(data = tmp.data, aes(x=`Lachnospiraceae`,fill=Dysbiosis),alpha=0.5) +
  scale_fill_manual(values = c("#BB5566", "#004488"))
combined_plot <- insert_xaxis_grob(combined_plot, x_box, position = "top")
ggdraw(combined_plot)
ggsave("OutputPlot/Test.pdf",width = 2,height = 3)

library(visNetwork)
library(geomnet)
library(igraph)
# network of non-interaction results
edge=modelfit0[modelfit0$FDR_taxa<0.05,c(1:3)]
colnames(edge)=c("from","to","width")
edge$color=NA
edge$color[edge$width<0]="#CCEEFF"
edge$color[edge$width>0]="#EEEEBB"
edge$width=100*abs(edge$width)

node=c(unique(edge$from),unique(edge$to))
node=data.frame(id=node,lable=node)
node$group=NA
node$group[1:43]="Gene"
node$group[44:56]="Bacteria"
node$color=NA
node$color[node$group=="Gene"]="#C2A5CF"
node$color[node$group=="Bacteria"]="#5AAE61"
node$font.size=50
node$shape=NA
node$shape[node$group=="Gene"]="dot"
node$shape[node$group=="Bacteria"]="square"

network=visNetwork(node, edge, width = 3000,height = 1500) %>%
  visIgraphLayout() %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges(smooth = list(roundness = 0.3)) %>%
  visLayout(randomSeed = 1111)
network
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.no.interaction.html")
tmp.genus="Faecalibacterium"
tmp.gene="CREB3L3"
tmp.data=merge(tmp.data,genes_correct[,tmp.gene,drop=F],by.x="ID",by.y="row.names",all=F)
ggplot(tmp.data, aes(Faecalibacterium,CREB3L3)) +
  geom_point(shape = 21, color = "#6699CC", size = 2)+
  geom_smooth(method = lm,color="#6699CC",linetype="dashed")+
  theme(legend.position="bottom")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
ggsave("OutputPlot/tmpplot.pdf",width = 5,height = 4)

library(pheatmap)
library(reshape2)
tmp.gene=modelfit0$Gene[modelfit0$FDR_taxa<0.05]
tmp.bac=modelfit0$Taxa[modelfit0$FDR_taxa<0.05]
associations_cor=acast(modelfit0[modelfit0$Gene %in% tmp.gene & modelfit0$Taxa %in% tmp.bac,], Gene~Taxa, value.var="Beta_taxa")
associations_sig=acast(modelfit0[modelfit0$Gene %in% tmp.gene & modelfit0$Taxa %in% tmp.bac,], Gene~Taxa, value.var="FDR_taxa")
heatmap.data=(associations_cor)
heatmap.data.lable=(associations_sig)
heatmap.data.lable[heatmap.data.lable<0.05 & heatmap.data.lable>0]="*"
heatmap.data.lable[heatmap.data.lable<0.1 & heatmap.data.lable>0]="*"
heatmap.data.lable[heatmap.data.lable>0.05]=""

paletteLength <- 50
myColor <- colorRampPalette(c("#4393C3", "white", "#FEDA8B"))(paletteLength)
myBreaks <- c(seq(min(heatmap.data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap.data)/paletteLength, max(heatmap.data), length.out=floor(paletteLength/2)))


pdf("OutputPlot//host-microbe.association.pdf",width = 5,height = 10)
pheatmap(heatmap.data,cluster_cols = T, cluster_rows = T,
         show_rownames=T, show_colnames=T,
         cellheight = 8,cellwidth = 10,fontsize_number=12,fontsize_row=6,fontsize_col = 8,
         display_numbers = heatmap.data.lable,color = myColor,breaks=myBreaks)
dev.off()

# network of interaction results
edge=modelfit[modelfit$FDR_interact<0.05,c(1:3)]
colnames(edge)=c("from","to","width")
edge$color="#DDDDDD"
edge$width=8

node=c(unique(edge$from),unique(edge$to))
node=data.frame(id=node,lable=node)
node$group=NA
node$group[1:101]="Gene"
node$group[102:129]="Bacteria"
node$color=NA
node$color[node$group=="Gene"]="#C2A5CF"
node$color[node$group=="Bacteria"]="#5AAE61"
node$font.size=50
node$shape=NA
node$shape[node$group=="Gene"]="dot"
node$shape[node$group=="Bacteria"]="square"

network=visNetwork(node, edge, width = 3000,height = 1500) %>%
  visIgraphLayout() %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges(smooth = list(roundness = 0.3)) %>%
  visLayout(randomSeed = 1111)
network
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.interaction.html")

# network of gene ~ taxa results
edge=modelfit_taxa[modelfit_taxa$FDR<0.05,c(1:3)]
colnames(edge)=c("from","to","width")
edge$color="#DDDDDD"
edge$width=8

node=c(unique(edge$from),unique(edge$to))
node=data.frame(id=node,lable=node)
node$group=NA
node$group[1:120]="Gene"
node$group[121:162]="Bacteria"
node$color=NA
node$color[node$group=="Gene"]="#C2A5CF"
node$color[node$group=="Bacteria"]="#5AAE61"
node$font.size=50
node$shape=NA
node$shape[node$group=="Gene"]="dot"
node$shape[node$group=="Bacteria"]="square"

network=visNetwork(node, edge, width = 3000,height = 1500) %>%
  visIgraphLayout() %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges(smooth = list(roundness = 0.3)) %>%
  visLayout(randomSeed = 1111)
network
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.gene-taxa.html")

# =================================================================================
# context-specific analysis
##### As Montreal E, Montreal B and biological_usage have significant effects on
##### both bacteria taxa and dysbiosis, therefore to identify the how the bacteria
##### modulates host inflammation-related gene expressions, we focus on conext-sepcific 
##### host-microbe interactions
# =================================================================================

write.table(genes_correct,file = "OutputTable/Context.specific.genes.correct.data.txt",sep = "\t",row.names = T,quote = F)
write.table(dysbiosis_all,file = "OutputTable/Context.specific.dysbiosis.data.txt",sep = "\t",row.names = F,quote = F)
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
covariate_context=merge(covariate_rna[,c("ResearchID","biological_use")],montreal,by.x="row.names",by.y="ID",all=F)
write.table(covariate_context,file = "OutputTable/Context.specific.context.data.txt",sep = "\t",row.names = F,quote = F)

# import files
genes_correct=read.table("OutputTable/Context.specific.genes.correct.data.txt",sep = "\t",header = T,row.names = 1,stringsAsFactors = F,check.names = F)
bacteria=read.table("OutputTable/CLR.bacteria.txt",header = T,row.names = 1,stringsAsFactors = F,check.names = F,sep = "\t")
covariate_context=read.table("OutputTable/Context.specific.context.data.txt",header = T,row.names = 1,stringsAsFactors = F,check.names = F,sep = "\t")
bacteria=bacteria[rownames(bacteria) %in% rownames(genes_correct),]
covariate_context=covariate_context[rownames(covariate_context) %in% rownames(bacteria),]
dysbiosis_all=read.table("OutputTable/Context.specific.dysbiosis.data.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t")

# individual interaction
# model: gene ~ bac + Montreal + bac * Montreal
quantile(dysbiosis_all$Dysbiosis)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=14.92912]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>14.92912]=1
modelfit_context = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],bacteria,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 26))
  colnames(tmp.result)=c("Gene","Taxa",
                         "Beta_TaxaBiological","Beta_Biological","Beta_interactBiological","P_taxaBiological","P_Biological","P_interactBiological",
                         "Beta_TaxaE","Beta_MontrealE","Beta_interactE","P_taxaE","P_MontrealE","P_interactE",
                         "Beta_TaxaB1B2","Beta_MontrealB1B2","Beta_interactB1B2","P_taxaB1B2","P_MontrealB1B2","P_interactB1B2",
                         "Beta_TaxaB1B3","Beta_MontrealB1B3","Beta_interactB1B3","P_taxaB1B3","P_MontrealB1B3","P_interactB1B3")
  
  for(n in 1:ncol(bacteria)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    
    # Biological_usage
    cat(yellow(i,"+++",tmp.gene,"bacteria","+++",n,"Biological_usage","\n"))
    
    tmp.biological_use=merge(tmp.test,covariate_context[,c("biological_use"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.biological_use=na.omit(tmp.biological_use)
    tmp.model.biological_use=lm(Gene ~ Bac + biological_use + Bac * biological_use,data = tmp.biological_use)
    tmp.model.biological_use=as.data.frame(summary(tmp.model.biological_use)$coefficients)
    
    # montreal E
    cat(green(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealE","\n"))
    
    tmp.MontrealE=merge(tmp.test,covariate_context[,c("MontrealE"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.MontrealE=na.omit(tmp.MontrealE)
    tmp.model.MontrealE=lm(Gene ~ Bac + MontrealE + Bac * MontrealE,data = tmp.MontrealE)
    tmp.model.MontrealE=as.data.frame(summary(tmp.model.MontrealE)$coefficients)
    
    # montreal B1B2
    cat(red(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealB1B2","\n"))
    
    tmp.B1B2=merge(tmp.test,covariate_context[,c("MontrealBB"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.B1B2=na.omit(tmp.B1B2)
    tmp.B1B2=tmp.B1B2[tmp.B1B2$MontrealBB==0 | tmp.B1B2$MontrealBB==1,]
    colnames(tmp.B1B2)[6]="B1B2"
    tmp.model.B1B2=lm(Gene ~ Bac + B1B2 + Bac * B1B2,data = tmp.B1B2)
    tmp.model.B1B2=as.data.frame(summary(tmp.model.B1B2)$coefficients)
    
    # montreal B1B3
    cat(cyan(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealB1B3","\n"))
    
    tmp.B1B3=merge(tmp.test,covariate_context[,c("MontrealBB"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.B1B3=na.omit(tmp.B1B3)
    tmp.B1B3=tmp.B1B3[tmp.B1B3$MontrealBB==0 | tmp.B1B3$MontrealBB==2,]
    colnames(tmp.B1B3)[6]="B1B3"
    tmp.model.B1B3=lm(Gene ~ Bac + B1B3 + Bac * B1B3,data = tmp.B1B3)
    tmp.model.B1B3=as.data.frame(summary(tmp.model.B1B3)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    
    tmp.result$Beta_TaxaBiological[n]=tmp.model.biological_use$Estimate[rownames(tmp.model.biological_use)=="Bac"]
    tmp.result$Beta_Biological[n]=tmp.model.biological_use$Estimate[rownames(tmp.model.biological_use)=="biological_use"]
    tmp.result$Beta_interactBiological[n]=tmp.model.biological_use$Estimate[rownames(tmp.model.biological_use)=="Bac:biological_use"]
    tmp.result$P_taxaBiological[n]=tmp.model.biological_use$`Pr(>|t|)`[rownames(tmp.model.biological_use)=="Bac"]
    tmp.result$P_Biological[n]=tmp.model.biological_use$`Pr(>|t|)`[rownames(tmp.model.biological_use)=="biological_use"]
    tmp.result$P_interactBiological[n]=tmp.model.biological_use$`Pr(>|t|)`[rownames(tmp.model.biological_use)=="Bac:biological_use"]
    
    tmp.result$Beta_TaxaE[n]=tmp.model.MontrealE$Estimate[rownames(tmp.model.MontrealE)=="Bac"]
    tmp.result$Beta_MontrealE[n]=tmp.model.MontrealE$Estimate[rownames(tmp.model.MontrealE)=="MontrealE"]
    tmp.result$Beta_interactE[n]=tmp.model.MontrealE$Estimate[rownames(tmp.model.MontrealE)=="Bac:MontrealE"]
    tmp.result$P_taxaE[n]=tmp.model.MontrealE$`Pr(>|t|)`[rownames(tmp.model.MontrealE)=="Bac"]
    tmp.result$P_MontrealE[n]=tmp.model.MontrealE$`Pr(>|t|)`[rownames(tmp.model.MontrealE)=="MontrealE"]
    tmp.result$P_interactE[n]=tmp.model.MontrealE$`Pr(>|t|)`[rownames(tmp.model.MontrealE)=="Bac:MontrealE"]
    
    tmp.result$Beta_TaxaB1B2[n]=tmp.model.B1B2$Estimate[rownames(tmp.model.B1B2)=="Bac"]
    tmp.result$Beta_MontrealB1B2[n]=tmp.model.B1B2$Estimate[rownames(tmp.model.B1B2)=="B1B2"]
    tmp.result$Beta_interactB1B2[n]=tmp.model.B1B2$Estimate[rownames(tmp.model.B1B2)=="Bac:B1B2"]
    tmp.result$P_taxaB1B2[n]=tmp.model.B1B2$`Pr(>|t|)`[rownames(tmp.model.B1B2)=="Bac"]
    tmp.result$P_MontrealB1B2[n]=tmp.model.B1B2$`Pr(>|t|)`[rownames(tmp.model.B1B2)=="B1B2"]
    tmp.result$P_interactB1B2[n]=tmp.model.B1B2$`Pr(>|t|)`[rownames(tmp.model.B1B2)=="Bac:B1B2"]
    
    tmp.result$Beta_TaxaB1B3[n]=tmp.model.B1B3$Estimate[rownames(tmp.model.B1B3)=="Bac"]
    tmp.result$Beta_MontrealB1B3[n]=tmp.model.B1B3$Estimate[rownames(tmp.model.B1B3)=="B1B3"]
    tmp.result$Beta_interactB1B3[n]=tmp.model.B1B3$Estimate[rownames(tmp.model.B1B3)=="Bac:B1B3"]
    tmp.result$P_taxaB1B3[n]=tmp.model.B1B3$`Pr(>|t|)`[rownames(tmp.model.B1B3)=="Bac"]
    tmp.result$P_MontrealB1B3[n]=tmp.model.B1B3$`Pr(>|t|)`[rownames(tmp.model.B1B3)=="B1B3"]
    tmp.result$P_interactB1B3[n]=tmp.model.B1B3$`Pr(>|t|)`[rownames(tmp.model.B1B3)=="Bac:B1B3"]
    
  }
  cat(black$bgWhite$bold("==========================================","\n"))
  cat(black$bgWhite$bold(i,"+++",tmp.gene,"\n"))
  cat(black$bgWhite$bold("==========================================","\n"))
  
  return.string=as.data.frame(tmp.result)
}

modelfit_context$FDR_interactBiological_use=p.adjust(modelfit_context$P_interactBiological,method = "BH")
modelfit_context$FDR_interactE=p.adjust(modelfit_context$P_interactE,method = "BH")
modelfit_context$FDR_interactB1B2=p.adjust(modelfit_context$P_interactB1B2,method = "BH")
modelfit_context$FDR_interactB1B3=p.adjust(modelfit_context$P_interactB1B3,method = "BH")
write.table(modelfit_context,"OutputTable/Context.interaction.txt",row.names = F,quote = F)

covariate_context$biological_use=as.factor(covariate_context$biological_use)
ggplot(covariate_context, aes(x=biological_use, y=Dysbiosis, fill=biological_use)) +theme_bw()+
  geom_boxplot()+
  scale_fill_manual(values=c("#BB5566", "#004488"))+guides(fill=F)
ggsave("OutputPlot/dys_biological.pdf")
wilcox.test(covariate_context$Dysbiosis[covariate_context$biological_use==0],covariate_context$Dysbiosis[covariate_context$biological_use==1])



















