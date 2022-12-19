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
library(RRHO)  
require(cowplot)
library(ggpubr)
library("biomaRt")
library("dplyr")
library(edgeR)
library(crayon)
library(DOSE)
library(ReactomePA)
library(clusterProfiler)
library(limma)
library(mixOmics)
library(robCompositions)
library(ComplexHeatmap)
library(circlize)
library(biomaRt)
source("Microbiome.function.R")
library(stringr)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# import data

covariate_rna=read.table("Covariate.rna.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
expression=read.table("RNAseq/Input.RNAseq.CLR.annot.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
basic_factors=c("ResearchID","age_at_biopsy","sex","BMI","Aminosalicylates","Thiopurines","Steroids","Batch","Biological_use")

pca=prcomp(genes,scale = TRUE)
fviz_eig(pca)
eigenvalue=get_eig(pca)
ind <- get_pca_ind(pca)
pca_matrix=as.data.frame(ind$coord)
pca_matrix=pca_matrix[,1:6]

pca_matrix=merge(pca_matrix,covariate_rna,by="row.names",all=F)
pca_matrix$Inflammation=as.factor(pca_matrix$Inflammation)

ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Location_rough)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+guides(color=F)+
  theme(legend.position = 'top')
ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Batch)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+
  theme(legend.position = 'top')

which(is.na(basic_factors))
which(is.na(genes))
# basic factors correction
genes_correct = apply(genes,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="gene"
  
  fit=glmer(gene~age_at_biopsy+sex+BMI+Aminosalicylates+Thiopurines+Steroids+Batch + (1|ResearchID),data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
genes_correct=as.data.frame(genes_correct)

# ===========================================================================================
# gene expression characterization, inflammation effects
# ===========================================================================================

# CD vs. HC at ileum
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
nn=0
result1 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== counting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                             inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                             inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])

}

# CD vs. HC at colon
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
nn=0
result2 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== couting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])
  
}

# UC vs. HC at colon
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="UC" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
nn=0
result3 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== couting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])
  
}

result1$group=1
result2$group=2
result3$group=3
result1$FDR=p.adjust(result1$inflammation.P)
result2$FDR=p.adjust(result2$inflammation.P)
result3$FDR=p.adjust(result3$inflammation.P)
result=rbind(result1[result1$FDR<0.05,],result2[result2$FDR<0.05,],result3[result3$FDR<0.05,])

# enrichment analysis
result1=result[result$group==1,]
result2=result[result$group==2,]
result3=result[result$group==3,]
tmp_genes=as.character(result1$Gene[result1$FDR<0.05])
tmp_genes = bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich= summary(enrichPathway(tmp_genes$ENTREZID,readable=T,pvalueCutoff = 0.05))
ratio=as.matrix(as.data.frame(str_split_fixed(enrich$GeneRatio,"/",2)))
enrich$ratio=as.numeric(ratio[,1])/as.numeric(ratio[,2])
enrich=enrich[order(enrich$qvalue,decreasing = F),]
enrich$Description=factor(enrich$Description,levels = rev(enrich$Description))
enrich=enrich[1:10,]
ggplot(enrich, aes(x=Description, y=ratio, fill=-log10(qvalue)))+
  geom_bar(stat="identity", color="black")+scale_fill_gradient(low="#4DBBD5B2", high="#0077BB")+
  theme_classic()+ coord_flip()+xlab("")+theme(axis.text.y=element_text(angle=-10,vjust=0.8,hjust = 0.8))

tmp_genes=as.character(result2$Gene[result2$FDR<0.05])
tmp_genes = bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich= summary(enrichPathway(tmp_genes$ENTREZID,readable=T,pvalueCutoff = 0.05))
ratio=as.matrix(as.data.frame(str_split_fixed(enrich$GeneRatio,"/",2)))
enrich$ratio=as.numeric(ratio[,1])/as.numeric(ratio[,2])
enrich=enrich[order(enrich$qvalue,decreasing = F),]
enrich$Description=factor(enrich$Description,levels = rev(enrich$Description))
enrich=enrich[1:10,]
ggplot(enrich, aes(x=Description, y=ratio, fill=-log10(qvalue)))+
  geom_bar(stat="identity", color="black")+scale_fill_gradient(low="#00A087B2", high="#009988")+
  theme_classic()+ coord_flip()+xlab("")+theme(axis.text.y=element_text(angle=-10,vjust=0.8,hjust = 0.8))

tmp_genes=as.character(result3$Gene[result3$FDR<0.05])
tmp_genes = bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich= summary(enrichPathway(tmp_genes$ENTREZID,readable=T,pvalueCutoff = 0.05))
ratio=as.matrix(as.data.frame(str_split_fixed(enrich$GeneRatio,"/",2)))
enrich$ratio=as.numeric(ratio[,1])/as.numeric(ratio[,2])
enrich=enrich[order(enrich$qvalue,decreasing = F),]
enrich$Description=factor(enrich$Description,levels = rev(enrich$Description))
enrich=enrich[1:20,]
ggplot(enrich, aes(x=Description, y=ratio, fill=-log10(qvalue)))+
  geom_bar(stat="identity", color="black")+scale_fill_gradient(low="#F39B7FB2", high="#EE7733")+
  theme_bw()+ coord_flip()

gene_selection=unique(union(result2$Gene[result3$FDR<0.05],result1$Gene[result1$FDR<0.05]))
list2=result2[result2$Gene %in% gene_selection,c("Gene","inflammation.beta")]
list3=result3[result3$Gene %in% gene_selection,c("Gene","inflammation.beta")]
list1=result1[result1$Gene %in% gene_selection,c("Gene","inflammation.beta")]

list2=list2[order(list2$inflammation.beta,decreasing = T),]
list3=list3[order(list3$inflammation.beta,decreasing = T),]
list1=list1[order(list1$inflammation.beta,decreasing = T),]
list2$rank2=1:nrow(list2)
list3$rank3=1:nrow(list3)
list1$rank1=1:nrow(list1)

list2=list2[order(rownames(list2)),]
list3=list3[order(rownames(list3)),]
list1=list1[order(rownames(list1)),]

dataset=cbind(list3[,"rank3",drop=F],list1[,"rank1",drop=F])
ggplot(dataset, aes(rank3, rank1,color="#00A087B2")) +
  geom_point(shape = 21, 
             color = "black", size = 1,fill="darkgreen")+
  geom_smooth(method = lm)+
  theme_bw()+ theme(legend.position="bottom")
cor.test(dataset$rank3,dataset$rank1)

x=list(A=result$Gene[result$group==1],B=result$Gene[result$group==2],C=result$Gene[result$group==3])
pdf("OutputPlot/Venn.inflammation.gene.pdf",width = 6,height = 6)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 8
)
dev.off()

# ===========================================================================================
# gene expression characterization, CDi colon vs. UCi colon
# ===========================================================================================

# CDi colon vs. UCi colon
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="UC",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
covariate1=covariate1[covariate1$Inflammation=="3",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]
genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]

nn=0
result0 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
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

# ===========================================================================================
# gene expression characterization, systematic inflammation: IBD non-inf vs. HC
# ===========================================================================================

# CD vs. HC at ileum
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
nn=0
result1.1 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== counting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])
  
}

# CD vs. HC at colon
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
nn=0
result2.1 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== couting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])
  
}

# UC vs. HC at colon
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="UC" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
nn=0
result3.1 = foreach(i=1:ncol(genes_correct_protein1),.combine = rbind) %do%  {
  tmp.genus=colnames(genes_correct_protein1)[i]
  x=genes_correct_protein1[,i,drop=F]
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  fit<-lm(tmp.data[,tmp.genus]~Inflammation,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  
  nn=nn+1
  cat(tmp.genus,"is done","\n")
  cat(yellow(nn, "========== couting", "\n"))
  
  return.string=data.frame(Gene=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"])
  
}

result1.1$group=1
result2.1$group=2
result3.1$group=3
result1.1$FDR=p.adjust(result1.1$inflammation.P)
result2.1$FDR=p.adjust(result2.1$inflammation.P)
result3.1$FDR=p.adjust(result3.1$inflammation.P)
result.1=rbind(result1.1[result1.1$FDR<0.05,],result2.1[result2.1$FDR<0.05,],result3.1[result3.1$FDR<0.05,])
