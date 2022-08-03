
.libPaths("/groups/umcg-weersma/tmp01/Shixian/Rpackage")
library(data.table)
library(foreach)

# load bacteria data and covariate
covariate_bac=read.table("data/Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]
bacteria=read.table("data/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]

# load gene data and covariate
gene=read.table(file = "data/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
metadata_rna=read.table("data/Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
samples_count=intersect(rownames(gene),intersect(rownames(bacteria),rownames(metadata_rna)))

gene=gene[rownames(gene) %in% samples_count,]
bacteria=bacteria[rownames(bacteria) %in% samples_count,]
metadata_rna=metadata_rna[rownames(metadata_rna) %in% samples_count,]

# RNA data correction (age, gender, BMI, bacth, inflammation, location (no diagnosis, location == diagnosis), and medication)
covariate_rna=metadata_rna[,c("Inflammation","Location_rough","age_at_biopsy","sex","BMI","Batch","Aminosalicylates","Thiopurines","Steroids")]
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

# ==============================================================
# Network construction
# ==============================================================
library("metagMisc")
library(NetCoMi,verbose=F)
library("Rcpp")
library(igraph)
library("Hmisc")
library(visNetwork)
library(geomnet)
library(reshape2)

MontrealB_bacteria=read.table("data/montrealB_bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
MontrealB_gene=read.table("data/montrealB_gene.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# Montreal B Network construction in B1 samples
matrix_B1B2=merge(bacteria[,colnames(bacteria) %in% MontrealB_bacteria$bacteria[MontrealB_bacteria$FDR_B1B2<0.05]],
                  gene[,colnames(gene) %in% MontrealB_gene$Gene[MontrealB_gene$FDR_B1B2<0.05]],by="row.names",all=F)
rownames(matrix_B1B2)=matrix_B1B2$Row.names
matrix_B1B2$Row.names=NULL
matrix_B1B2=matrix_B1B2[rownames(matrix_B1B2) %in% rownames(metadata_rna)[metadata_rna$MontrealB==0 | metadata_rna$MontrealB==3],]
matrix_B1B2_bacteria=bacteria[rownames(bacteria) %in% rownames(matrix_B1B2),colnames(bacteria) %in% colnames(matrix_B1B2)]
matrix_B1B2__gene=gene[rownames(gene) %in% rownames(matrix_B1B2),colnames(gene) %in% colnames(matrix_B1B2)]

net_B1B2_bacteria <- netConstruct(matrix_B1B2_bacteria,  
                                  measure = "spearman",
                                  normMethod = "none", 
                                  zeroMethod = "none",
                                  sparsMethod = "t-test", 
                                  adjust="BH",
                                  alpha=0.05,
                                  verbose = 3,seed = 100)
props_B1B2_bacteria <- netAnalyze(net_B1B2_bacteria, centrLCC = TRUE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector",connectivity=F,normNatConnect=F,
                                  weightDeg = FALSE, normDeg = FALSE)
sink(file = "props_B1B2_bacteria.txt")
summary(props_B1B2_bacteria, numbNodes = 5L)
sink()

net_B1B2_gene <- netConstruct(matrix_B1B2__gene,  
                              measure = "spearman",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "t-test", 
                              adjust="BH",
                              alpha=0.05,
                              verbose = 3,seed = 100)
props_B1B2_gene <- netAnalyze(net_B1B2_gene, centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector",
                              weightDeg = FALSE, normDeg = FALSE)
sink(file = "props_B1B2_gene.txt")
summary(props_B1B2_gene, numbNodes = 5L)
sink()











