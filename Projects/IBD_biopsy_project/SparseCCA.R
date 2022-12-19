# SparseCCA

library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
library(ggalluvial)
library(ggpubr)
library("dplyr")
library(PMA) ##for sparse CCA
library(data.table)

###################### Functions ################
load_gene_expr <- function(filename){
  genes <- data.frame(fread(filename,sep="\t",head=T), row.names = 1, check.names = F)
  genes <- as.matrix(genes)
  
}

load_microbiome_abnd <- function(filename){
  microbes <- read.table(filename,sep="\t",head=T, row.names = 1, check.names =F)
  microbes <- as.matrix(microbes)
  
}

filter_genes <- function(genes, qt){
  genes.t <- t(genes)
  genes.sd <- transform(as.data.frame(genes.t), SD=apply(as.data.frame(genes.t),1, sd, na.rm = TRUE))
  ## select top genes with high SD (~ variability) across samples
  SD_quantile <- quantile(genes.sd$SD) ## identical to summary(genes.sd$SD)
  SD_cutoff <- SD_quantile[qt] ## 2nd quantile -- 25th quantile.
  genes.sd <- genes.sd[order(genes.sd$SD, decreasing = T),]
  top.variable.genes <- rownames(genes.sd[genes.sd$SD > SD_cutoff,])
  ## subset these genes from gene table
  select <- which(colnames(genes) %in% top.variable.genes)
  genes <- genes[,select]
}

get_tuning_params <- function(X, Z, outputFile = NULL, pdfFile = NULL, num_perm=25){
  perm.out <- CCA.permute(X,Z,typex="standard",typez="standard", nperms = num_perm)
  ## can tweak num. of permutations if needed (default = 25)
  
  ## get p-value for first left and right canonical covariates resulting from the selected
  ## tuning parameter value
  perm.out.pval <- perm.out$pvals[which(perm.out$penaltyxs == perm.out$bestpenaltyx)]
  
  if(!is.null(outputFile)){
    sink(outputFile)
    print(perm.out)
    print(paste0("P-value for selected tuning params:",perm.out.pval)) #0
    sink()
  }
  
  if(!is.null(pdfFile)){
    pdf(pdfFile)
    plot(perm.out)
    dev.off()
  }
  
  return(perm.out)
}

run_sparseCCA <- function(X, Z, CCA.K, penaltyX, penaltyZ, vInit=NULL, outputFile=NULL){
  CCA.out <-  CCA(X,Z,typex="standard",typez="standard",K=CCA.K,
                  penaltyx=penaltyX,penaltyz=penaltyZ,
                  v=vInit) ## standardize=T by default
  if(!is.null(outputFile)){
    sink(outputFile)
    print(CCA.out)
    sink()
  }
  
  ## add rownames to output factors
  rownames(CCA.out$u) <- colnames(X)
  rownames(CCA.out$v) <- colnames(Z)
  ## compute contribution of selected features to each of the samples.
  CCA_var_genes <- X %*% CCA.out$u ## canonical variance for genes 
  CCA_var_microbes <- Z %*% CCA.out$v ## canonical variance for microbes
  
  return(list(CCA.out, CCA_var_genes, CCA_var_microbes))
  
}

get_avg_features <- function(cca_cov, CCA.K){
  num_features <- 0
  for(k in 1:CCA.K){
    num_features <- num_features + length(which(cca_cov[,k]!=0))
  }
  avg_features <- num_features/CCA.K
}

save_CCA_components <- function(CCA.out, CCA.K, dirname){
  ## Print canonical covariates in files 
  for(i in CCA.K){
    #i <- 2 ##debug
    print(paste0("Writing significant component = ", i))
    selected_X <- which(CCA.out$u[,i]!=0) 
    selected_X <- rownames(CCA.out$u)[selected_X]
    coeff_X <- unname(CCA.out$u[selected_X,i])
    selected_Z <- which(CCA.out$v[,i]!=0)
    selected_Z <- rownames(CCA.out$v)[selected_Z]
    coeff_Z <- unname(CCA.out$v[selected_Z,i])
    ## Make all vectors of same length to avoid repetition of elements from shorter vectors.
    n <- max(length(selected_X), length(selected_Z))
    length(selected_X) <- n                      
    length(selected_Z) <- n
    length(coeff_X) <- n
    length(coeff_Z) <- n
    selected_XZ <- as.data.frame(cbind(gene = selected_X, gene_coeff = coeff_X,
                                       taxa = selected_Z, taxa_coeff = coeff_Z))
    write.table(selected_XZ, file=paste0(dirname,"gene_taxa_component_",i,".txt"), sep = "\t", col.names = NA)
  }
  
}

# ======================================================================================================================================
# select protein-coding, inflammation-related genes
# ======================================================================================================================================

# import inflammation-related
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
# import clr gene table
gene=read.table("OutputTable/Genes.factor_corrected.txt",sep = "\t",header = T,row.names = 1)
covariate_rna=read.table("Covariate.rna.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)

# ======================================================================================================================================
# all the data are used after factors correction
# ======================================================================================================================================

bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)
inflammation=read.table("sparseCCA/Inflammation.genes.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)

dim(bacteria)
dim(inflammation)

keep=as.data.frame(matrix(ncol = ncol(bacteria),nrow = 1))
for(i in 1:ncol(bacteria)){
  keep[,i]=length(unique(bacteria[,i]))
}
keep=which(keep>2)
bacteria=bacteria[,keep]

keep=as.data.frame(matrix(ncol = ncol(inflammation),nrow = 1))
for(i in 1:ncol(inflammation)){
  keep[,i]=length(unique(inflammation[,i]))
}
keep=which(keep>2)
inflammation=inflammation[,keep]
stopifnot(all(rownames(inflammation) == rownames(bacteria)))

## select tuning parameters using grid-search
X <- as.matrix(inflammation)
Y <- as.matrix(bacteria)
scoreXcv <- c()
scoreYcv <- c()
penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)
corr_all <- matrix(nrow = length(penaltyX), ncol =  length(penaltyY))
num_samples <- nrow(inflammation)
start_time <- Sys.time()
for( i in 1:length(penaltyX)){
  for(j in 1:length(penaltyY)){

    for(k in 1:num_samples){

      print(paste0("Index: i = ",i,", j =", j," k = ",k)); flush.console()

      res <- CCA(X[-k,],Y[-k,], penaltyx = penaltyX[i], penaltyz = penaltyY[j], K=1, niter = 5, trace = F)

      scoreXcv[k] <- X[k,]%*%res$u 
      scoreYcv[k] <- Y[k,]%*%res$v
    }
    ## correlation between scores for X and Y for all held out samples.
    corr_all[i,j] = cor(scoreXcv,scoreYcv) 
  }
}
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste0("Time elapsed for all = ", time_elapsed))

row.names(corr_all) <- as.character(penaltyX)
colnames(corr_all) <- as.character(penaltyY)

corr_all_df <- as.data.frame(corr_all)
rownames(corr_all_df)
colnames(corr_all_df)

## save to file
write.table(corr_all_df, file = ("sparseCCA//corr_all.txt"), sep="\t", row.names = T, col.names = NA )
save(corr_all, file = paste0("sparseCCA/corr_all.RData"))

penaltyX <- seq(0.1,0.4,length=10)
penaltyY <- seq(0.15,0.4,length=10)

# find index with max absolute corr
bestpenalty <- which(abs(corr_all) == max(abs(corr_all)), arr.ind = TRUE)
bestpenalty

bestpenaltyX <- penaltyX[bestpenalty[1]]
bestpenaltyX 
bestpenaltyY <- penaltyY[bestpenalty[2]]
bestpenaltyY

## order abs corr to get top 5 corr
index <- order(abs(corr_all), decreasing = T)
abs(corr_all)[index][1:5] ## top 5 absolute corr

cca.k = 10 ## number of desired components

## Run sparse CCA using selected tuning param using permutation search
cca <- run_sparseCCA(as.matrix(inflammation), as.matrix(bacteria), cca.k, bestpenaltyX, bestpenaltyY,
                     outputFile=paste0("sparseCCA//CCA.output.",bestpenaltyX,"_",bestpenaltyY,".txt"))

## canonical correlation for each component:
cca[[1]]$cors

## average number of genes and microbes in resulting components
avg_genes <- get_avg_features(cca[[1]]$u, cca.k)
avg_genes

avg.microbes <- get_avg_features(cca[[1]]$v, cca.k)
avg.microbes


## Test significance of correlation using LOOCV
X <- as.matrix(inflammation)
Y <- as.matrix(bacteria)
cca.k = 10
scoresXcv <- matrix(nrow = nrow(X), ncol = cca.k)
scoresYcv <-  matrix(nrow = nrow(Y), ncol = cca.k)
corr_pval <- c()
corr_r <- c()
for(i in 1:nrow(inflammation)){ 
  res <- CCA(X[-i,],Y[-i,], penaltyx=bestpenaltyX, penaltyz=bestpenaltyY, K=cca.k, trace = F) 
  for(j in 1:cca.k){
    print(paste0("i = ", i," K = ", j)); flush.console()
    scoresXcv[i,j] <- X[i,]%*%res$u[,j]
    scoresYcv[i,j] <- Y[i,]%*%res$v[,j]
  }
}
for(j in 1:cca.k){
  # plot(scoresXcv,scoresYcv)
  corr <- cor.test(scoresXcv[,j],scoresYcv[,j])
  corr_pval[j] <- corr$p.value
  corr_r[j] <- corr$estimate
}
corr_pval

length(which(corr_pval < 0.1)) 
which(corr_pval < 0.1)

length(which(corr_pval < 0.05))
which(corr_pval < 0.05)

corr_padj <- p.adjust(corr_pval, method = "BH")
corr_padj
which(corr_padj < 0.1)

length(which(corr_padj < 0.1))

## LOOCV corr
corr_r

outputFile <- paste0("sparseCCA/crc_sparseCCA_summary_",bestpenaltyX,"_",bestpenaltyY,".txt")
sink(outputFile)
cat(paste0(" bestpenaltyX = ", bestpenaltyX, ", bestpenaltyY = ", bestpenaltyY))
cat(paste0("\n cor(Xu,Yv): \n"))
cat(paste0(signif(cca[[1]]$cors, digits = 4)))
cat(paste0("\n Avg. no. of genes across components = ",avg_genes))
cat(paste0("\n Avg. no. of microbes across components= ", avg.microbes))
cat(paste0("\n P-value for components (LOOCV): \n"))
cat(paste0(signif(corr_pval, digits = 4)))
cat(paste0("\n LOOCV corr: \n"))
cat(paste0(signif(corr_r, digits = 4)))
cat(paste0("\n No. of components with p-value < 0.1 = ", length(which(corr_pval < 0.1))))
cat(paste0("\n No. of components with p-value < 0.05 = ", length(which(corr_pval < 0.05))))
cat(paste0("\n No. of components with FDR < 0.1 = ", length(which(corr_padj < 0.1))))
cat(paste0("\n Significant components: \n" ))
cat(paste0(which(corr_padj < 0.1)))
sink()

write.table(cca[[2]], file = paste0("sparseCCA/CCA_var_genes.txt"), sep="\t", row.names = T, col.names = NA )
write.table(cca[[3]], file = paste0("sparseCCA/CCA_var_microbes.txt"), sep="\t", row.names = T, col.names = NA )

## only spit out significant components
sig <- which(corr_padj < 0.1)
dirname <- paste0("sparseCCA/sig_gene_taxa_components_",bestpenaltyX,"_", bestpenaltyY,"_padj/")
ifelse(!dir.exists(dirname), dir.create(dirname), FALSE)
save_CCA_components(cca[[1]],sig,dirname)

