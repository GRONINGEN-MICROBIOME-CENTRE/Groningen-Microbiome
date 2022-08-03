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

genes=as.data.frame((expression))
genes=genes[order(rownames(genes)),]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]
covariate_rna$Inflammation[covariate_rna$Inflammation=="Light"]="No"

covariate_rna$Inflammation[covariate_rna$Inflammation=="Yes"]=3
covariate_rna$Inflammation[covariate_rna$Inflammation=="No" & covariate_rna$Diagnosis=="CD"]=2
covariate_rna$Inflammation[covariate_rna$Inflammation=="No" & covariate_rna$Diagnosis=="UC"]=2
covariate_rna$Inflammation[covariate_rna$Inflammation=="No" & covariate_rna$Diagnosis=="Control"]=1
covariate_rna$Inflammation=as.numeric(covariate_rna$Inflammation)

covariate_rna$BMI[is.na(covariate_rna$BMI)]=median(covariate_rna$BMI[!is.na(covariate_rna$BMI)])
covariate_rna$Inflammation[is.na(covariate_rna$Inflammation)]=median(covariate_rna$Inflammation[!is.na(covariate_rna$Inflammation)])
covariate_rna$smoking_DB[is.na(covariate_rna$smoking_DB)]=median(covariate_rna$smoking_DB[!is.na(covariate_rna$smoking_DB)])
write.table(covariate_rna,file = "Covariate.rna.organized.txt",row.names = T,quote = F,sep = "\t")

basic_factors=c("ResearchID","age_at_biopsy","sex","BMI","Aminosalicylates","Thiopurines","Steroids","Batch","Biological_use")
basic_factors=covariate_rna[,colnames(covariate_rna) %in% basic_factors]
basic_factors$ResearchID=as.factor(basic_factors$ResearchID)
genes=genes[rownames(genes) %in% rownames(basic_factors),]
basic_factors=basic_factors[rownames(basic_factors) %in% rownames(genes),]
basic_factors=basic_factors[order(rownames(basic_factors)),]
genes=genes[order(rownames(genes)),]

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
  
  fit=lm(gene~age_at_biopsy+sex+BMI+Aminosalicylates+Thiopurines+Steroids+Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
genes_correct=as.data.frame(genes_correct)
write.table(genes_correct,file = "OutputTable/Genes.basic.correction.txt",sep = "\t",row.names = T,quote = F)

# ===========================================================================================
# protein coding gene annotation
# ===========================================================================================

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
geneID <-  colnames(genes)
protein_coding <- getBM(attributes= c("hgnc_symbol","gene_biotype"),filters = c("hgnc_symbol","biotype"),
                        values=list(geneID,"protein_coding"), mart= mart)
genes_correct_protein=genes_correct[,colnames(genes_correct) %in% protein_coding$hgnc_symbol]
genes=genes[,colnames(genes) %in% protein_coding$hgnc_symbol]
write.table(genes_correct_protein,file = "OutputTable/Genes.basic.correction.protein.coding.txt",sep = "\t",row.names = T,quote = F)
write.table(genes,file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = T,quote = F)

genes=genes[,colSums(genes>0)>(0.1*nrow(genes))]

gene_dist <- dist(genes_correct_protein)
gene_hclust <- hclust(gene_dist, method = "complete")
plot(gene_hclust, labels = FALSE)

sample_order=rownames(genes_correct_protein)[gene_hclust$order]

col.annot=tmp_covariate
col.annot=col.annot[match(sample_order, rownames(col.annot)),]

column_ha = HeatmapAnnotation(Diagnosis=col.annot$Diagnosis,
                              Inflammation = col.annot$Inflammation, 
                              Location = col.annot$Location_rough,
                              surgery=col.annot$surgery)
abundance=as.matrix(t(genes_correct_protein))
set.seed(20)
pdf("heatmap.pdf",width = 10,height = 10)
Heatmap(abundance,column_order = sample_order, name = "mat", top_annotation = column_ha)
dev.off()

# ===========================================================================================
# import data, start analysis
# ===========================================================================================
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
genes=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
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
  theme(legend.position = 'top')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
ggplot (pca_matrix, aes(Dim.1,Dim.2,color=Inflammation)) + geom_point(size=2) + 
  theme_bw()+scale_color_npg()+
  theme(legend.position = 'top')+
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

ggscatter(pca_matrix, x = "Dim.1", y = "Dim.2",
          color = "Location_rough", palette = "jco",
          ellipse = TRUE, 
          mean.point = TRUE,
          star.plot = TRUE,alpha=0.5,xlab = "PC1",ylab = "PC2")
ggsave("OutputPlot/PCA.CLR.RNAseq.location.pdf",width = 3,height = 3)
mm=wilcox.test(pca_matrix$Dim.2[pca_matrix$Location_rough=="colon"],pca_matrix$Dim.2[pca_matrix$Location_rough=="ileum"])

ggscatter(pca_matrix, x = "Dim.1", y = "Dim.2",
          color = "Inflammation", palette = c("#117733", "#88CCEE", "#BB5566"),
          ellipse = TRUE, 
          mean.point = TRUE,
          star.plot = TRUE,alpha=0.5,xlab = "PC1",ylab = "PC2")
ggsave("OutputPlot/PCA.CLR.RNAseq.inflammation.pdf",width = 3,height = 3)
mm=kruskal.test(Dim.1 ~ Inflammation, data = pca_matrix)
kruskal.test(Dim.2 ~ Inflammation, data = pca_matrix)

ggscatter(pca_matrix, x = "Dim.1", y = "Dim.2",
          color = "Diagnosis", palette = c("#EE8866","#44BB99","#AA4499"),
          ellipse = TRUE, 
          mean.point = TRUE,
          star.plot = TRUE,alpha=0.5,xlab = "PC1",ylab = "PC2")
ggsave("OutputPlot/PCA.CLR.RNAseq.IBD.pdf",width = 3,height = 3)
kruskal.test(Dim.1 ~ Diagnosis, data = pca_matrix)
mm=kruskal.test(Dim.2 ~ Diagnosis, data = pca_matrix)

# ===========================================================================================
# gene expression characterization, inflammation effects
# ===========================================================================================
genes_correct_protein=read.table("OutputTable/Genes.basic.correction.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
genes_correct_protein=genes_correct_protein[,colSums(genes_correct_protein>0)>0.1*nrow(genes_correct_protein)]
tmp_covariate=covariate_rna[,c("Diagnosis","Cohort","Inflammation","Location_rough","resec_part_colon","resec_ileocec","resec_part_small")]
tmp_covariate=tmp_covariate[rownames(tmp_covariate) %in% rownames(genes_correct_protein),]
tmp_covariate=tmp_covariate[order(rownames(tmp_covariate)),]

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
write.table(result,file="OutputTable/RNAseq.inflammation.compare.txt",row.names = F,sep = "\t",quote = F)

result=read.table("OutputTable//RNAseq.inflammation.compare.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)

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
ggsave("OutputPlot/RNAseq.inflammation.result1.pdf",width = 6,height = 6)

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
ggsave("OutputPlot/RNAseq.inflammation.result2.pdf",width = 6,height = 6)

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
ggsave("OutputPlot/RNAseq.inflammation.result3.pdf",width = 8,height = 8)

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
ggsave("OutputPlot/list3list1.compare.pdf",width = 5,height = 5)
cor.test(dataset$rank3,dataset$rank1)

result=read.table("OutputTable//RNAseq.inflammation.compare.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
library(ggvenn)
x=list(A=result$Gene[result$group==1],B=result$Gene[result$group==2],C=result$Gene[result$group==3])
pdf("OutputPlot/Venn.inflammation.gene.pdf",width = 6,height = 6)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#CD534CFF"),
  stroke_size = 1, set_name_size = 8
)
dev.off()

tmp.gene="IL17A"
tmp.values=genes[,tmp.gene,drop=F]
tmp.values=merge(tmp.values,covariate_rna,by="row.names",all=F)
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
covariate2=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate2=covariate2[covariate2$Location_rough=="colon",]
covariate3=tmp_covariate[tmp_covariate$Diagnosis=="UC" | tmp_covariate$Diagnosis=="Control",]
covariate3=covariate3[covariate3$Location_rough=="colon",]
tmp.values$Group=NA
tmp.values$Group[tmp.values$Row.names %in% rownames(covariate1)]="Group1"
tmp.values$Group[tmp.values$Row.names %in% rownames(covariate2)]="Group2"
tmp.values$Group[tmp.values$Row.names %in% rownames(covariate3)]="Group3"
tmp.data=tmp.values[tmp.values$Group=="Group3",]
tmp.data=tmp.data[tmp.data$Inflammation==1,]
tmp.data$Group="Group2"
tmp.data=tmp.data[!is.na(tmp.data$Row.names),]
tmp.values=rbind(tmp.values,tmp.data)
tmp.values$Inflammation=as.factor(tmp.values$Inflammation)
tmp.values=tmp.values[!is.na(tmp.values$Group),]
ggplot(tmp.values, aes(x=Inflammation, y=IL17A,fill=Inflammation)) + 
  geom_boxplot(aes(fill=Inflammation),position=position_dodge(2),size=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("#117733", "#88CCEE", "#BB5566"))+
  guides(color=FALSE)+guides(fill=FALSE)+
  facet_grid(.~Group)+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())+xlab("")
ggsave("OutputPlot/IL17A.pdf",width = 4,height = 2)

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

tmp_genes=as.character(result0$Gene[result0$FDR<0.05 & result0$diagnosis.beta>0])
tmp_genes = bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich= summary(enrichPathway(tmp_genes$ENTREZID,readable=T,pvalueCutoff = 0.05))
ratio=as.matrix(as.data.frame(str_split_fixed(enrich$GeneRatio,"/",2)))
enrich$ratio=as.numeric(ratio[,1])/as.numeric(ratio[,2])
enrich=enrich[order(enrich$qvalue,decreasing = F),]
enrich$Description=factor(enrich$Description,levels = rev(enrich$Description))
enrich=enrich[1:10,]
ggplot(enrich, aes(x=Description, y=ratio, fill=-log10(qvalue)))+
  geom_bar(stat="identity", color="black")+scale_fill_gradient(low="#D1BBD7", high="#882E72")+
  theme_bw()+ coord_flip()
ggsave("OutputPlot/RNAseq.diagnosis.result0.UCup.pdf",width = 8,height = 3)

tmp_genes=as.character(result0$Gene[result0$FDR<0.05 & result0$diagnosis.beta<0])
tmp_genes = bitr(tmp_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
enrich= summary(enrichPathway(tmp_genes$ENTREZID,readable=T,pvalueCutoff = 0.05))
ratio=as.matrix(as.data.frame(str_split_fixed(enrich$GeneRatio,"/",2)))
enrich$ratio=as.numeric(ratio[,1])/as.numeric(ratio[,2])
enrich=enrich[order(enrich$qvalue,decreasing = F),]
enrich$Description=factor(enrich$Description,levels = rev(enrich$Description))
enrich=enrich[1:10,]
ggplot(enrich, aes(x=Description, y=ratio, fill=-log10(qvalue)))+
  geom_bar(stat="identity", color="black")+scale_fill_gradient(low="#EC7014", high="#993404")+
  theme_bw()+ coord_flip()
ggsave("OutputPlot/RNAseq.diagnosis.result0.UCdown.pdf",width = 8,height = 3)

tmp.gene="APOB"
tmp.values=genes_correct_protein[,tmp.gene,drop=F]
tmp.values=merge(tmp.values,covariate_rna,by="row.names",all=F)
tmp.values=tmp.values[tmp.values$Inflammation=="3",]
tmp.values=tmp.values[tmp.values$Diagnosis=="CD" | tmp.values$Diagnosis=="UC",]
tmp.values=tmp.values[tmp.values$Location_rough=="colon",]

ggplot(tmp.values, aes(x=Diagnosis, y=APOB,fill=Diagnosis)) + 
  geom_boxplot(aes(fill=Diagnosis),position=position_dodge(2),size=0.6)+
  theme_classic()+
  scale_fill_manual(values = c("#FDB366", "#9970AB"))+
  guides(color=FALSE)+guides(fill=FALSE)+
  theme(strip.background =element_rect(fill="white"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("")
ggsave("OutputPlot/APOB.pdf",width = 1.8,height = 1.8)

library(pheatmap)
library(reshape2)

heatmap.data=as.data.frame(t(genes_correct_protein1))
heatmap.data=heatmap.data[rownames(heatmap.data) %in% result0$Gene[result0$FDR<0.05],]
paletteLength <- 20
myColor <- colorRampPalette(c("#993404", "white", "#2166AC"))(paletteLength)
myBreaks <- c(seq(min(heatmap.data), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(heatmap.data)/paletteLength, max(heatmap.data), length.out=floor(paletteLength/2)))
annotation <- covariate1[,c("Diagnosis"),drop=F]
Diagnosis        <- c("#FDB366", "#9970AB")
names(Diagnosis) <- c("CD", "UC")
anno_colors <- list(Diagnosis = Diagnosis)

pdf("OutputPlot//CDiUCi.pdf",width = 10,height = 10)
pheatmap(heatmap.data,cluster_cols = T, cluster_rows = T,scale = "row",
         show_rownames=F, show_colnames=F,
         cellheight = 0.2,cellwidth = 1,fontsize_number=12,fontsize_row=6,fontsize_col = 8,
         color = myColor,breaks=myBreaks,annotation = annotation, annotation_colors = anno_colors)
dev.off()
write.table(result0,"OutputTable/CDi.vs.UCi.txt",sep = "\t",quote = F,row.names = F)
result0=read.table("OutputTable/CDi.vs.UCi.txt",sep = "\t",header = T)

# ===========================================================================================
# gene expression characterization, systematic inflammation: IBD non-inf vs. HC
# ===========================================================================================

# CD vs. HC at ileum
covariate1=tmp_covariate[tmp_covariate$Diagnosis=="CD" | tmp_covariate$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
covariate1=covariate1[covariate1$Inflammation!="3",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
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
covariate1=covariate1[covariate1$Inflammation!="3",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
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
covariate1=covariate1[covariate1$Inflammation!="3",]
which(is.na(covariate1))
genes_correct_protein1=genes_correct_protein[rownames(genes_correct_protein) %in% rownames(covariate1),]

genes_correct_protein1=genes_correct_protein1[order(rownames(genes_correct_protein1)),]
covariate1=covariate1[order(rownames(covariate1)),]
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
write.table(result.1,file="OutputTable/RNAseq.inflammation.trend.txt",row.names = F,sep = "\t",quote = F)

trend_genes=intersect(result.1$Gene[result.1$group==3],intersect(result.1$Gene[result.1$group==1],result.1$Gene[result.1$group==2]))

# ===========================================================================================
# inflammation genes enrichment analysis
# ===========================================================================================
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation_genes=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))

IBD_genes=read.table("IBDloci/IBD.locus.genes.NG2017.txt",header = F,stringsAsFactors = F)
IBD_genes=unique(IBD_genes$V1)
all_genes=colnames(genes)

IBD_genes=IBD_genes[IBD_genes %in% all_genes]
length(IBD_genes)
length(inflammation_genes)
length(all_genes)
length(inflammation_genes[inflammation_genes %in% IBD_genes])

data <-
  matrix(c(75, 595, 1212, 15934),
         nrow = 2,
         dimnames = list(Guess = c("Milk", "Tea"),
                         Truth = c("Milk", "Tea")))
fisher.test(data, alternative = "greater")























