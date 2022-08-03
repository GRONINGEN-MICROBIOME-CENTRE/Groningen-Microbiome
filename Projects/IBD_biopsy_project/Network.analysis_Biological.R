
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
library("Hmisc")
library(geomnet)
library(reshape2)

Biological_gene=read.table("data/biological_gene.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
Biological_bacteria=read.table("data/biological_bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# Biological Network construction in non-user samples
matrix_biological=merge(bacteria[,colnames(bacteria) %in% Biological_bacteria$Gene[Biological_bacteria$FDR<0.05]],
                  gene[,colnames(gene) %in% Biological_gene$Gene[Biological_gene$FDR<0.05]],by="row.names",all=F)
rownames(matrix_biological)=matrix_biological$Row.names
matrix_biological$Row.names=NULL
matrix_biological=matrix_biological[rownames(matrix_biological) %in% rownames(metadata_rna)[metadata_rna$biological_use==0],]
matrix_biological_bacteria=bacteria[rownames(bacteria) %in% rownames(matrix_biological),colnames(bacteria) %in% colnames(matrix_biological)]
matrix_biological_gene=gene[rownames(gene) %in% rownames(matrix_biological),colnames(gene) %in% colnames(matrix_biological)]

net_Biological_bacteria <- netConstruct(matrix_biological_bacteria,  
                                  measure = "spearman",
                                  normMethod = "none", 
                                  zeroMethod = "none",
                                  sparsMethod = "t-test", 
                                  adjust="BH",
                                  alpha=0.05,
                                  verbose = 3,seed = 101)
props_Biological_bacteria <- netAnalyze(net_Biological_bacteria, centrLCC = TRUE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector",connectivity=F,normNatConnect=F,
                                  weightDeg = FALSE, normDeg = FALSE)
sink(file = "props_Biological_bacteria.txt")
summary(props_Biological_bacteria, numbNodes = 5L)
sink()

net_Biological_gene <- netConstruct(matrix_biological_gene,  
                              measure = "spearman",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "t-test", 
                              adjust="BH",
                              alpha=0.05,
                              verbose = 3,seed = 101)
props_Biological_gene <- netAnalyze(net_Biological_gene, centrLCC = TRUE,
                              clustMethod = "cluster_fast_greedy",
                              hubPar = "eigenvector",
                              weightDeg = FALSE, normDeg = FALSE)
sink(file = "props_Biological_gene.txt")
summary(props_Biological_gene, numbNodes = 5L)
sink()













































