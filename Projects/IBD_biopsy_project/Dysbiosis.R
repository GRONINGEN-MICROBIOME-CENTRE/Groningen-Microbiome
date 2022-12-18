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

# import data
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=read.table("OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

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


# calculate median distance as dysbiosis scores
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

# dysbiosis vs. disease group
dysbiosis_all_group=dysbiosis_all
dysbiosis_all_group$Group=as.character(dysbiosis_all_group$Group)
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"]]="CD"
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"]]="UC"
dysbiosis_all_group=dysbiosis_all_group[dysbiosis_all_group$Group!="IBD",]
wilcox.test(dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="Control"],dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="CD"])
wilcox.test(dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="Control"],dysbiosis_all_group$Dysbiosis[dysbiosis_all_group$Group=="UC"])

# check if any tech parameter confounded dysbiosis, e.g. sequencing reads
covariate_all=read.table("MetaData/Metadaata.16S.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
pca_analysis=merge(pca_analysis,covariate_all,by.x="Row.names",by.y="row.names",all=F)
m=ggplot (pca_analysis, aes(V1,V2,color=log(AlignedReads))) + 
  geom_point() + theme_bw() +
  guides(size=F)+scale_color_jcolors_contin("pal2")+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")

# check if any clinical covariate confounded dysbiosis, e.g. inflammation
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

# dysbiosis vs. bacteria, generalized linear model
dys_genus = foreach(i=1:ncol(bacteria),.combine = rbind) %do%  {
  tmp.genus=bacteria[,i,drop=F]
  tmp.data=merge(dysbiosis_all,tmp.genus,by.x="ID",by.y="row.names",all=F)
  tmp.test( Dysbiosis ~ tmp.data[,4] + (1|ResearchID),  data=tmp.data)
  
  return.string=data.frame(Bacteria=colnames(bacteria)[i],Cor=tmp.test$estimate,Pvalue=tmp.test$p.value)
}
dys_genus$FDR=p.adjust(dys_genus$Pvalue,method = "BH")

# dysbiosis vs. gene expression after correcting for potential confounders
gene=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
covariate_rna=covariate_rna[,c("Inflammation","Location_rough","age_at_biopsy","sex","BMI","Batch","Aminosalicylates","Thiopurines","Steroids")]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]
gene=gene[order(rownames(gene)),]
genes_correct = apply(gene,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_rna[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="tmp.gene"
  
  fit=glmer(tmp.gene~age_at_biopsy+sex+BMI+Batch+Inflammation+Location_rough+Aminosalicylates+Thiopurines+Steroids + (1|ResearchID),data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
genes_correct=as.data.frame(genes_correct)

# model: gene ~ dys, generalized linear model
quantile(dysbiosis_all$Dysbiosis,.90)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=14.92912]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>14.92912]=1
modelfit00 = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=(genes_correct[,tmp.gene,drop=F])
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="row.names",by.x="ID",all=F)
  colnames(tmp.data)=c("ID","Dys","Group","Gene")

  tmp.model=glmer(Gene ~ Dys +  (1|ResearchID) ,data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coefficients)

  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=data.frame(Gene=tmp.gene,Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}

# model: gene ~ bac, generalized linear model
modelfit_taxa = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],bacteria,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta","SE","Zscore","P")
  
  for(n in 1:ncol(bacteria)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=glm(Gene ~ Bac + (1|ResearchID),data = tmp.test)
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

# model: gene ~ bac + dys, generalized linear model
modelfit0 = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],genus,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(genus),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_dys","P_taxa","P_dys")
  
  for(n in 1:ncol(genus)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=glmer(Gene ~ Bac + Dys + (1|ResearchID),data = tmp.test)
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

# model: gene ~ bac + dysbiosis + bac * dysbiosis, generalized linear model
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
    tmp.model=glmer(Gene ~ Bac + Dys + Bac * Dys + (1|ResearchID),data = tmp.test)
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

# =================================================================================
# context-specific analysis
##### As Montreal E, Montreal B and biological_usage have significant effects on
##### both bacteria taxa and dysbiosis, therefore to identify the how the bacteria
##### modulates host inflammation-related gene expressions, we focus on conext-sepcific 
##### host-microbe interactions
# =================================================================================
# import files
genes_correct=read.table("OutputTable/Context.specific.genes.correct.data.txt",sep = "\t",header = T,row.names = 1,stringsAsFactors = F,check.names = F)
bacteria=read.table("OutputTable/CLR.bacteria.txt",header = T,row.names = 1,stringsAsFactors = F,check.names = F,sep = "\t")
covariate_context=read.table("OutputTable/Context.specific.context.data.txt",header = T,row.names = 1,stringsAsFactors = F,check.names = F,sep = "\t")
dysbiosis_all=read.table("OutputTable/Context.specific.dysbiosis.data.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t")

# individual interaction
# model: gene ~ bac + Montreal + bac * Montreal, generalized linear model
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
    tmp.model.biological_use=glmer(Gene ~ Bac + biological_use + Bac * biological_use + (1|ResearchID),data = tmp.biological_use)
    tmp.model.biological_use=as.data.frame(summary(tmp.model.biological_use)$coefficients)
    
    # montreal E
    cat(green(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealE","\n"))
    
    tmp.MontrealE=merge(tmp.test,covariate_context[,c("MontrealE"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.MontrealE=na.omit(tmp.MontrealE)
    tmp.model.MontrealE=glmer(Gene ~ Bac + MontrealE + Bac * MontrealE + (1|ResearchID),data = tmp.MontrealE)
    tmp.model.MontrealE=as.data.frame(summary(tmp.model.MontrealE)$coefficients)
    
    # montreal B1B2
    cat(red(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealB1B2","\n"))
    
    tmp.B1B2=merge(tmp.test,covariate_context[,c("MontrealBB"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.B1B2=na.omit(tmp.B1B2)
    tmp.B1B2=tmp.B1B2[tmp.B1B2$MontrealBB==0 | tmp.B1B2$MontrealBB==1,]
    colnames(tmp.B1B2)[6]="B1B2"
    tmp.model.B1B2=glmer(Gene ~ Bac + B1B2 + Bac * B1B2 + (1|ResearchID),data = tmp.B1B2)
    tmp.model.B1B2=as.data.frame(summary(tmp.model.B1B2)$coefficients)
    
    # montreal B1B3
    cat(cyan(i,"+++",tmp.gene,"bacteria","+++",n,"MontrealB1B3","\n"))
    
    tmp.B1B3=merge(tmp.test,covariate_context[,c("MontrealBB"),drop=F],by.x="ID",by.y="row.names",all=F)
    tmp.B1B3=na.omit(tmp.B1B3)
    tmp.B1B3=tmp.B1B3[tmp.B1B3$MontrealBB==0 | tmp.B1B3$MontrealBB==2,]
    colnames(tmp.B1B3)[6]="B1B3"
    tmp.model.B1B3=glmer(Gene ~ Bac + B1B3 + Bac * B1B3 + (1|ResearchID),data = tmp.B1B3)
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


