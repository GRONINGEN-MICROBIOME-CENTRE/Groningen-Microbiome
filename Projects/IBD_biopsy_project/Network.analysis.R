#
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
library(Rcpp)
library(msigdbr)
library(data.table)
msigdb_enrichment <- function(pathway_DB, pathways, background_genes, genes_of_interest){
  
  enrichment_list <- list()
  for(i in 1:length(pathways)){
    
    pathway <- pathways[i]
    
    ## genes in this pathway
    pathway_gene_set <- pathway_DB[pathway_DB$gs_name == pathway,]$human_gene_symbol
    length(pathway_gene_set) 
    ## If the criteria for min and max #genes in a given pathway is not satified, 
    ## skip testing the current pathway
    if(length(pathway_gene_set) < min_genes || length(pathway_gene_set) > max_genes) next
    
    ## The contingency table
    ## Inspired by: https://www.pathwaycommons.org/guide/primers/statistics/fishers_exact_test/
    ##                    genes_of_interest genes_NOT_of_interest  Total
    ## In_pathway               x                 m - x              m
    ## Not_in_pathway         k - x             n - (k - x)          n
    ## Total                    k                 m + n - k          m + n
    
    ## m, #overlapping genes in this pathway and background genes for CRC
    m <- length(intersect(background_genes,pathway_gene_set))
    m ## 11
    
    ## n, #genes in background but not in pathway 
    n <- length(setdiff(background_genes,pathway_gene_set))
    n #12502
    
    ## x, #genes of interest in pathway 
    x <- length(intersect(pathway_gene_set,genes_of_interest))
    x ## 1
    ## If the overlap between genes of interest and the pathway is less than the overlap cut-off, 
    ## then skip testing the current pathway
    if(x < overlap_genes) next 
    
    ## Extract list of genes in the genes of interest that are included in the pathway. 
    gene_names_in_pathway = noquote(paste(intersect(pathway_gene_set,genes_of_interest), collapse = ","))
    
    ## k, total #genes in genes list of interest
    k <- length(genes_of_interest)
    k
    
    ## Build the contigency table
    contingency_table <- matrix(c(x, k-x, m-x, n-(k-x)), #matrix is filled along the column by default
                                nrow = 2, ncol = 2, 
                                dimnames = list(c("In_pathway","Not_in_pathway"),
                                                c("Genes_of_interest","Genes_NOT_of_interest"))
    )
    
    contingency_table
    
    fisher_result <- fisher.test( contingency_table, alternative = "greater")
    
    ## save details in a dataframe
    enrichment_result_df <- data.frame( pathway = pathway,
                                        BG_genes = length(background_genes),
                                        genes_in_pathway = length(pathway_gene_set),
                                        genes_in_path_and_BG = m,
                                        genes_of_interest = k,
                                        genes_of_interest_in_pathway = x,
                                        gene_names_in_pathway = gene_names_in_pathway,
                                        ## fill in contingency table entries: x, k-x, m-x, n-(k-x) for z-score computation.
                                        cont_n1 = x,
                                        cont_n2 = k-x,
                                        cont_n3 = m-x,
                                        cont_n4 = n-(k-x),
                                        CI_95 = paste0(signif(fisher_result$conf.int[1],5),"-",signif(fisher_result$conf.int[2],5)),
                                        odds_ratio = unname(fisher_result$estimate),
                                        p_val = fisher_result$p.value
    )
    
    enrichment_list[[i]] <- enrichment_result_df
    
  }
  
  return(enrichment_list)
  
}

## Do enrichment per component
perform_enrichment <- function(input_dir, filenames, msigdb_collection_name, pathway_DB, pathways, background_genes, output_dir){
  enriched_list <- list()
  count <- 0
  
  ## perform enrichment for each component 
  for(i in filenames){
    
    ## debug
    # i <- "gene_taxa_component_1.txt"
    
    count <- count + 1
    
    print(paste0("Enrichment for component: ", i));flush.console()
    ## load genes list for a given component
    sparseCCA_genes <- read.table(paste0(input_dir,"/",i),sep="\t",head=T, row.names = 1, check.names = F, stringsAsFactors = F)
    sparseCCA_genes <- sparseCCA_genes$gene
    genes_of_interest <- sparseCCA_genes
    
    ## perform enrichment using pathways in msigdb collection
    msigdb_pathways <- msigdb_enrichment(pathway_DB, pathways, background_genes, genes_of_interest) 
    ## convert the list to dataframe
    msigdb_pathways_df <- do.call(rbind.data.frame, msigdb_pathways)
    ## Add current component as a column to enrichment result
    msigdb_pathways_df$Component <- rep(paste0(i),nrow(msigdb_pathways_df))
    ## Add collection name as a column to the enrichemnt result
    msigdb_pathways_df$Collection <- msigdb_collection_name
    
    ## Make last column(i.e., component) as first column
    msigdb_pathways_df <- msigdb_pathways_df[,c(ncol(msigdb_pathways_df)-1, ncol(msigdb_pathways_df),1:(ncol(msigdb_pathways_df)-2))]
    
    ## update column name to reflect component, except pathway colname
    # colnames(msigdb_pathways_df[,-3]) <- paste(colnames(msigdb_pathways_df[,-3]), count, sep = "_")
    
    ## sort pathways by pval
    msigdb_pathways_df <- msigdb_pathways_df[order(msigdb_pathways_df$p_val, decreasing = F),]
    length(msigdb_pathways_df$pathway)
    
    ## save dataframe of pathways for a component to file
    if(!is.null(output_dir)){
      filename <- strsplit(i,"\\.")[[1]][1]
      write.table(msigdb_pathways_df, file = paste0(output_dir,"/msigdb_",filename,".txt") , sep="\t", row.names = F)
    }
    ## add df to list
    enriched_list[[count]] <- msigdb_pathways_df
  }
  
  ## This is the list of dataframes, where each dataframe holds enrichment result for a component.  
  return(enriched_list)
}


# ======================================================================================================================
# single network construction and network compare
# B2 -> B1; B3 -> B1, use diff genes and bacteria to construct two networks, and compare the two
# similary to montreal E and biological usage
# ======================================================================================================================

# load bacteria data and covariate
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]
bacteria=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]
covariate_bac$surgery=NA
covariate_bac$surgery=covariate_bac$resec_ileocec+covariate_bac$resec_part_colon+covariate_bac$resec_part_small
covariate_bac$surgery[covariate_bac$surgery!=0]=1

# load gene data and covariate
gene=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
metadata_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
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
#write.table(genes_correct, file = "OutputTable/Genes.basic.correction.protein.coding.txt",row.names = T,quote = F,sep = "\t")

# ==============================================================
# Montreal B analysis
# ==============================================================

# gene
MontrealB_gene = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],metadata_rna[,"MontrealB",drop=F],by="row.names",all=F)
  
  tmp.data$MontrealB_class=NA
  tmp.data$MontrealB_class[tmp.data$MontrealB==0]="B1"
  tmp.data$MontrealB_class[tmp.data$MontrealB==3]="B1"
  tmp.data$MontrealB_class[tmp.data$MontrealB==1]="B2"
  tmp.data$MontrealB_class[tmp.data$MontrealB==4]="B2"
  tmp.data$MontrealB_class[tmp.data$MontrealB==2]="B3"
  tmp.data$MontrealB_class[tmp.data$MontrealB==5]="B3"
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="gene"
  
  # B1 as reference
  tmp.model.B1B2=lm(gene ~ MontrealB_class, data = subset(tmp.data,MontrealB_class!="B3"))
  tmp.model.B1B2=as.data.frame(summary(tmp.model.B1B2)$coef)
  
  tmp.model.B1B3=lm(gene ~ MontrealB_class, data = subset(tmp.data,MontrealB_class!="B2"))
  tmp.model.B1B3=as.data.frame(summary(tmp.model.B1B3)$coef)
  
  cat(green(i,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,
                           Beta_B1B2=tmp.model.B1B2$Estimate[2],Pvalue_B1B2=tmp.model.B1B2$`Pr(>|t|)`[2],
                           Beta_B1B3=tmp.model.B1B3$Estimate[2],Pvalue_B1B3=tmp.model.B1B3$`Pr(>|t|)`[2])
}
MontrealB_gene$FDR_B1B2=p.adjust(MontrealB_gene$Pvalue_B1B2,method = "BH")
MontrealB_gene$FDR_B1B3=p.adjust(MontrealB_gene$Pvalue_B1B3,method = "BH")

# bacteria
MontrealB_bacteria = foreach(i=1:ncol(bacteria),.combine = rbind) %do%  {
  
  tmp.bacteria=colnames(bacteria)[i]
  tmp.data=merge(bacteria[,tmp.bacteria,drop=F],covariate_bac[,c("age_at_biopsy","sex","Inflammation","Location_rough","surgery","MontrealB"),drop=F],by="row.names",all=F)
  
  tmp.data$MontrealB_class=NA
  tmp.data$MontrealB_class[tmp.data$MontrealB==0]="B1"
  tmp.data$MontrealB_class[tmp.data$MontrealB==3]="B1"
  tmp.data$MontrealB_class[tmp.data$MontrealB==1]="B2"
  tmp.data$MontrealB_class[tmp.data$MontrealB==4]="B2"
  tmp.data$MontrealB_class[tmp.data$MontrealB==2]="B3"
  tmp.data$MontrealB_class[tmp.data$MontrealB==5]="B3"
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="bacteria"
  
  # B1 as reference
  tmp.model.B1B2=lm(bacteria ~ MontrealB_class + age_at_biopsy + sex + Inflammation + Location_rough + surgery, data = subset(tmp.data,MontrealB_class!="B3"))
  tmp.model.B1B2=as.data.frame(summary(tmp.model.B1B2)$coef)
  
  tmp.model.B1B3=lm(bacteria ~ MontrealB_class + age_at_biopsy + sex + Inflammation + Location_rough + surgery, data = subset(tmp.data,MontrealB_class!="B2"))
  tmp.model.B1B3=as.data.frame(summary(tmp.model.B1B3)$coef)
  
  cat(green(i,"+++",tmp.bacteria,"\n"))
  
  return.string=data.frame(bacteria=tmp.bacteria,
                           Beta_B1B2=tmp.model.B1B2$Estimate[7],Pvalue_B1B2=tmp.model.B1B2$`Pr(>|t|)`[7],
                           Beta_B1B3=tmp.model.B1B3$Estimate[7],Pvalue_B1B3=tmp.model.B1B3$`Pr(>|t|)`[7])
}
MontrealB_bacteria$FDR_B1B2=p.adjust(MontrealB_bacteria$Pvalue_B1B2,method = "BH")
MontrealB_bacteria$FDR_B1B3=p.adjust(MontrealB_bacteria$Pvalue_B1B3,method = "BH")

# ==============================================================
# Montreal E analysis
# ==============================================================

# gene
MontrealE_gene = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],metadata_rna[,"MontrealE",drop=F],by="row.names",all=F)
  
  tmp.data$MontrealE_class=NA
  tmp.data$MontrealE_class[tmp.data$MontrealE==0]="E1"
  tmp.data$MontrealE_class[tmp.data$MontrealE==1]="E2"
  tmp.data$MontrealE_class[tmp.data$MontrealE==2]="E3"
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="gene"
  
  # E1 as reference
  tmp.model.E1E2=lm(gene ~ MontrealE_class, data = subset(tmp.data,MontrealE_class!="E3"))
  tmp.model.E1E2=as.data.frame(summary(tmp.model.E1E2)$coef)
  
  tmp.model.E1E3=lm(gene ~ MontrealE_class, data = subset(tmp.data,MontrealE_class!="E2"))
  tmp.model.E1E3=as.data.frame(summary(tmp.model.E1E3)$coef)
  
  cat(green(i,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,
                           Eeta_E1E2=tmp.model.E1E2$Estimate[2],Pvalue_E1E2=tmp.model.E1E2$`Pr(>|t|)`[2],
                           Eeta_E1E3=tmp.model.E1E3$Estimate[2],Pvalue_E1E3=tmp.model.E1E3$`Pr(>|t|)`[2])
}
MontrealE_gene$FDR_E1E2=p.adjust(MontrealE_gene$Pvalue_E1E2,method = "BH")
MontrealE_gene$FDR_E1E3=p.adjust(MontrealE_gene$Pvalue_E1E3,method = "BH")

# bacteria
MontrealE_bacteria = foreach(i=1:ncol(bacteria),.combine = rbind) %do%  {
  
  tmp.bacteria=colnames(bacteria)[i]
  tmp.data=merge(bacteria[,tmp.bacteria,drop=F],covariate_bac[,c("age_at_biopsy","sex","Inflammation","Location_rough","MontrealE","surgery"),drop=F],by="row.names",all=F)
  
  tmp.data$MontrealE_class=NA
  tmp.data$MontrealE_class[tmp.data$MontrealE==0]="E1"
  tmp.data$MontrealE_class[tmp.data$MontrealE==1]="E2"
  tmp.data$MontrealE_class[tmp.data$MontrealE==2]="E3"
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="bacteria"
  
  # E1 as reference
  tmp.model.E1E2=lm(bacteria ~ MontrealE_class + age_at_biopsy + sex + Inflammation + Location_rough +surgery, data = subset(tmp.data,MontrealE_class!="E3"))
  tmp.model.E1E2=as.data.frame(summary(tmp.model.E1E2)$coef)
  
  tmp.model.E1E3=lm(bacteria ~ MontrealE_class + age_at_biopsy + sex + Inflammation + Location_rough +surgery, data = subset(tmp.data,MontrealE_class!="E2"))
  tmp.model.E1E3=as.data.frame(summary(tmp.model.E1E3)$coef)
  
  cat(green(i,"+++",tmp.bacteria,"\n"))
  
  return.string=data.frame(bacteria=tmp.bacteria,
                           Eeta_E1E2=tmp.model.E1E2$Estimate[2],Pvalue_E1E2=tmp.model.E1E2$`Pr(>|t|)`[2],
                           Eeta_E1E3=tmp.model.E1E3$Estimate[2],Pvalue_E1E3=tmp.model.E1E3$`Pr(>|t|)`[2])
}
MontrealE_bacteria$FDR_E1E2=p.adjust(MontrealE_bacteria$Pvalue_E1E2,method = "BH")
MontrealE_bacteria$FDR_E1E3=p.adjust(MontrealE_bacteria$Pvalue_E1E3,method = "BH")

# ==============================================================
# Biological usage analysis
# ==============================================================

# gene
Biological_gene = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  
  tmp.gene=colnames(genes_correct)[i]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],metadata_rna[,"biological_use",drop=F],by="row.names",all=F)
  
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="gene"
  
  # non-user as reference
  tmp.model=lm(gene ~ biological_use, data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coef)
  
  cat(yellow(i,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,
                           Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}
Biological_gene$FDR=p.adjust(Biological_gene$Pvalue,method = "BH")

# bacteria
Biological_bacteria = foreach(i=1:ncol(bacteria),.combine = rbind) %do%  {
  
  tmp.bacteria=colnames(bacteria)[i]
  tmp.data=merge(bacteria[,tmp.bacteria,drop=F],covariate_bac[,c("age_at_biopsy","sex","Inflammation","Location_rough","surgery"),drop=F],by="row.names",all=F)
  tmp.data=merge(tmp.data,metadata_rna[,"biological_use",drop=F],by.x="Row.names",by.y="row.names",all = F)
  
  tmp.data=na.omit(tmp.data)
  colnames(tmp.data)[2]="bacteria"
  
  # non-user as reference
  tmp.model=lm(bacteria ~ biological_use + surgery + age_at_biopsy + sex + Inflammation + Location_rough, data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coef)
  
  cat(yellow(i,"+++",tmp.bacteria,"\n"))
  
  return.string=data.frame(Gene=tmp.bacteria,
                           Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}
Biological_bacteria$FDR=p.adjust(Biological_bacteria$Pvalue,method = "BH")

write.table(Biological_gene,file = "OutputTable/biological_gene.txt",sep = "\t",quote = F,row.names = F)
write.table(Biological_bacteria,file = "OutputTable/biological_bacteria.txt",sep = "\t",quote = F,row.names = F)
write.table(MontrealB_bacteria,file = "OutputTable/montrealB_bacteria.txt",sep = "\t",quote = F,row.names = F)
write.table(MontrealB_gene,file = "OutputTable/montrealB_gene.txt",sep = "\t",quote = F,row.names = F)
write.table(MontrealE_bacteria,file = "OutputTable/montrealE_bacteria.txt",sep = "\t",quote = F,row.names = F)
write.table(MontrealE_gene,file = "OutputTable/montrealE_gene.txt",sep = "\t",quote = F,row.names = F)

# ==============================================================
# Montreal Network construction
# ==============================================================
library("metagMisc")
library(NetCoMi,verbose=F)
library("Rcpp")
library(igraph)
library("Hmisc")
library(visNetwork)
library(geomnet)
library(reshape2)

MontrealB_bacteria=read.table("OutputTable/montrealB_bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
MontrealB_gene=read.table("OutputTable/montrealB_gene.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# Montreal B Network construction in B1 samples and define core genes/bacteria
matrix_B1B2=merge(bacteria[,colnames(bacteria) %in% MontrealB_bacteria$bacteria[MontrealB_bacteria$FDR_B1B2<0.05]],
                  gene[,colnames(gene) %in% MontrealB_gene$Gene[MontrealB_gene$FDR_B1B2<0.05]],by="row.names",all=F)
rownames(matrix_B1B2)=matrix_B1B2$Row.names
matrix_B1B2$Row.names=NULL
matrix_B1=matrix_B1B2[rownames(matrix_B1B2) %in% rownames(metadata_rna)[metadata_rna$MontrealB==0 | metadata_rna$MontrealB==3],]
matrix_B2=matrix_B1B2[rownames(matrix_B1B2) %in% rownames(metadata_rna)[metadata_rna$MontrealB==1 | metadata_rna$MontrealB==4],]
matrix_B1_bacteria=bacteria[rownames(bacteria) %in% rownames(matrix_B1),colnames(bacteria) %in% colnames(matrix_B1)]
matrix_B1_gene=gene[rownames(gene) %in% rownames(matrix_B1),colnames(gene) %in% colnames(matrix_B1)]

net_B1_bacteria <- netConstruct(matrix_B1_bacteria,  
                           measure = "spearman",
                           normMethod = "none", 
                           zeroMethod = "none",
                           sparsMethod = "t-test", 
                           adjust="BH",
                           alpha=0.05,
                           verbose = 3,seed = 100)
props_B1_bacteria <- netAnalyze(net_B1_bacteria, centrLCC = TRUE,
                         clustMethod = "cluster_fast_greedy",
                         hubPar = "eigenvector",connectivity=F,normNatConnect=F,
                         weightDeg = FALSE, normDeg = FALSE)

sink("OutputTable/props_B1_bacteria.txt")
summary(props_B1_bacteria, numbNodes = 5L)
sink()
plot(props_B1_bacteria)

net_B1_gene <- netConstruct(matrix_B1_gene,  
                              measure = "spearman",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "threshold", 
                              thresh = 0.3,
                              verbose = 3,seed = 100)
#props_B1_gene <- netAnalyze(net_B1_gene, centrLCC = TRUE,
#                                  clustMethod = "cluster_fast_greedy",
#                                  hubPar = "eigenvector",
#                                  weightDeg = FALSE, normDeg = FALSE)
sink("OutputTable/props_B1_gene.txt")
summary(props_B1_gene, numbNodes = 5L)
sink()
props_B1_gene=read.table("OutputTable/props_montrealB.genes.extract.txt",header = T,stringsAsFactors = F)

tmp.gene="AK3"
tmp.bacteria="Faecalibacterium"
tmp.B1=matrix_B1[,colnames(matrix_B1) %in% c(tmp.gene,tmp.bacteria)]
tmp.B2=matrix_B2[,colnames(matrix_B2) %in% c(tmp.gene,tmp.bacteria)]
tmp.B1$Group="B1"
tmp.B2$Group="B2"
tmp.data=rbind(tmp.B1,tmp.B2)

ggplot(tmp.data, aes(Group,Faecalibacterium,fill=Group)) +
  geom_boxplot()+
  theme(legend.position="bottom")+scale_fill_npg()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
ggsave("OutputPlot/test.pdf",width = 4,height = 4)

# view networks in B1
cor_calculation=rcorr(as.matrix(matrix_B1),type="spearman")
cor_B1=cor_calculation$r
sig_B1=cor_calculation$P

cor_B1[upper.tri(cor_B1)] = NA
sig_B1[upper.tri(sig_B1)] = NA

sig_B1=setNames(melt(sig_B1, na.rm=T), c('from', 'to', 'sig'))
edge_B1=setNames(melt(cor_B1, na.rm=T), c('from', 'to', 'width'))
edge_B1=edge_B1[edge_B1$width!=1,]

edge_B1=cbind(edge_B1,sig_B1[,"sig",drop=F])
edge_B1$FDR=p.adjust(edge_B1$sig,method = "BH")
edge_B1=edge_B1[edge_B1$FDR<0.05,]

edge_B1$color=NA
edge_B1$color[edge_B1$width<0]="#99DDFF"
edge_B1$color[edge_B1$width>0]="#EEDD88"
edge_B1$width=abs(edge_B1$width)

tmp.bacteria=edge_B1[edge_B1$to %in% colnames(matrix_B1_bacteria),]
tmp.bacteria=tmp.bacteria[tmp.bacteria$from %in% colnames(matrix_B1_bacteria),]
tmp.gene=edge_B1[edge_B1$to %in% colnames(matrix_B1_bacteria),]
tmp.gene=tmp.gene[!tmp.gene$from %in% colnames(matrix_B1_bacteria),] # this is the host-microbe interaction network
write.table(tmp.gene,"OutputTable/Network.B1.host_microbe.txt",sep = "\t",row.names = F,quote = F)

edge=rbind(tmp.bacteria,tmp.gene)
edge$width=edge$width*10

node=unique(c(as.character(unique(edge$from)),as.character(unique(edge$to))))
node=data.frame(id=node,lable=node)
node$group=NA
node$group[node$lable %in% colnames(matrix_B1_bacteria)]="Bacteria"
node$group[node$lable %in% colnames(matrix_B1_gene)]="Gene"
node$group[node$lable %in% props_B1_gene$HubGene]="GeneHubs"

node$color=NA
node$color[node$group=="Gene"]="#BBBBBB"
node$color[node$group=="Bacteria"]="#994455"
node$color[node$group=="GeneHubs"]="#0077BB"

node$shape=NA
node$shape[node$group=="Gene"]="dot"
node$shape[node$group=="Bacteria"]="dot"
node$shape[node$group=="GeneHubs"]="dot"

node$font.size[node$group=="Bacteria"]=80
node$font.size[node$group=="GeneHubs"]=40
node$font.size[node$group=="Gene"]=0
node$size[node$group=="Bacteria"]=20
node$size[node$group=="Gene"]=5
node$size[node$group=="GeneHubs"]=20

network=visNetwork(node, edge) %>%
  visIgraphLayout(randomSeed = 1111,layout = "layout_in_circle",type = "full") %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges(smooth = list(roundness = 0.3))
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.B1.html")

#enrichment
msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
class(msigdb_C2) #[1] "tbl_df"     "tbl"        "data.frame"
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)

## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)
msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)

path_DB <- c("CP:REACTOME")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 

pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways

min_genes <- 10
max_genes <- 10000
overlap_genes <- 5
msigdb_collection_name <- "C2_CP" 

background_genes=colnames(gene)
candidate_gene=unique(as.character(MontrealB_gene$Gene[MontrealB_gene$FDR_B1B2<0.05]))

aa=(msigdb_enrichment(pathway_DB, pathways, background_genes,candidate_gene))
msigdb_pathways_df <- do.call(rbind.data.frame, aa)
msigdb_pathways_df=msigdb_pathways_df[order(msigdb_pathways_df$p_val),]

# network compare
cor_calculation=rcorr(as.matrix(matrix_B2),type="spearman")
cor_B2=cor_calculation$r
sig_B2=cor_calculation$P

cor_B2[upper.tri(cor_B2)] = NA
sig_B2[upper.tri(sig_B2)] = NA

sig_B2=setNames(melt(sig_B2, na.rm=T), c('from', 'to', 'sig'))
edge_B2=setNames(melt(cor_B2, na.rm=T), c('from', 'to', 'cor'))
edge_B2=edge_B2[edge_B2$cor!=1,]
edge_B2=cbind(edge_B2,sig_B2[,"sig",drop=F])

cor_calculation=rcorr(as.matrix(matrix_B1),type="spearman")
cor_B1=cor_calculation$r
sig_B1=cor_calculation$P

cor_B1[upper.tri(cor_B1)] = NA
sig_B1[upper.tri(sig_B1)] = NA

sig_B1=setNames(melt(sig_B1, na.rm=T), c('from', 'to', 'sig'))
edge_B1=setNames(melt(cor_B1, na.rm=T), c('from', 'to', 'cor'))
edge_B1=edge_B1[edge_B1$cor!=1,]
edge_B1=cbind(edge_B1,sig_B1[,"sig",drop=F])

write.table(edge_B1[edge_B1$FDR<0.05 & edge_B1$to %in% colnames(matrix_B1_bacteria),],"OutputTable/B1.microbiota-gene.clusters.txt",sep = "\t",quote = F,row.names = F)
write.table(edge_B2[edge_B2$FDR<0.05 & edge_B2$to %in% colnames(matrix_B1_bacteria),],"OutputTable/B2.microbiota-gene.clusters.txt",sep = "\t",quote = F,row.names = F)

library(psych)
edge_B1$Zscore=fisherz(edge_B1$cor)
edge_B2$Zscore=fisherz(edge_B2$cor)
edge_B1$FDR=p.adjust(edge_B1$sig,method = "BH")
edge_B2$FDR=p.adjust(edge_B2$sig,method = "BH")

node_compare=matrix(nrow = ncol(matrix_B1_bacteria),ncol = 2)
node_compare=as.data.frame(node_compare)
colnames(node_compare)=c("Bacteria","Pvalue")
for(i in 1:ncol(matrix_B1_bacteria)){
  tmp.bacteria=colnames(matrix_B1_bacteria)[i]
  tmp.B1=edge_B1[edge_B1$to==tmp.bacteria,]
  tmp.B1=tmp.B1[tmp.B1$from %in% colnames(matrix_B1_gene),]
  tmp.B1$pair=paste(tmp.B1$from,tmp.B1$to)
  
  tmp.B2=edge_B2[edge_B2$to==tmp.bacteria,]
  tmp.B2=tmp.B2[tmp.B2$from %in% colnames(matrix_B1_gene),]
  tmp.B2$pair=paste(tmp.B2$from,tmp.B2$to)
  
  tmp.B1=tmp.B1[order(tmp.B1$Zscore),]
  tmp.B2=tmp.B2[order(tmp.B2$Zscore),]
  
  tmp.B1$rank=1:nrow(tmp.B1)
  tmp.B2$rank=1:nrow(tmp.B2)
  
  tmp.gene=union(tmp.B1$from[tmp.B1$FDR<0.05],tmp.B2$from[tmp.B2$FDR<0.05])
  tmp.B1=tmp.B1[tmp.B1$from %in% tmp.gene,]
  tmp.B2=tmp.B2[tmp.B2$from %in% tmp.gene,]
  
  tmp.B1=tmp.B1[order(tmp.B1$from),]
  tmp.B2=tmp.B2[order(tmp.B2$from),]
  
  mm=(wilcox.test(tmp.B1$Zscore,tmp.B2$Zscore,paired=T))$p.value
  
  #tmp.B1.pre=bacteria[rownames(bacteria) %in% rownames(matrix_B1),tmp.bacteria]
  #tmp.B1.pre=length(tmp.B1.pre[tmp.B1.pre>0])/length(tmp.B1.pre)
  #tmp.B2.pre=bacteria[rownames(bacteria) %in% rownames(matrix_B2),tmp.bacteria]
  #tmp.B2.pre=length(tmp.B2.pre[tmp.B2.pre>0])/length(tmp.B2.pre)
  
  node_compare$Bacteria[i]=tmp.bacteria
  #node_compare$B1_pre[i]=tmp.B1.pre
  #node_compare$B2_pre[i]=tmp.B2.pre
  node_compare$Pvalue[i]=mm
}
node_compare$FDR=p.adjust(node_compare$Pvalue,method = "BH")

node_enrich = foreach(i=1:length(node_compare$Bacteria[node_compare$FDR<0.05]),.combine = rbind) %do%  {
  tmp.bacteria=node_compare$Bacteria[node_compare$FDR<0.05][i]
  tmp_B1=edge_B1[edge_B1$to==tmp.bacteria,]
  tmp_B1=tmp_B1[tmp_B1$FDR<0.05,]
  tmp_B1=tmp_B1[tmp_B1$from %in% colnames(matrix_B1_gene),]
  tmp_B2=edge_B2[edge_B2$to==tmp.bacteria,]
  tmp_B2=tmp_B2[tmp_B2$FDR<0.05,]
  tmp_B2=tmp_B2[tmp_B2$from %in% colnames(matrix_B1_gene),]
  
  background_genes=colnames(gene)
  tmp.candidate1=unique(as.character(tmp_B1$from))
  tmp.candidate2=unique(as.character(tmp_B2$from))
  
  tmp.msigdb1 <- do.call(rbind.data.frame, (msigdb_enrichment(pathway_DB, pathways, background_genes,tmp.candidate1)))
  tmp.msigdb2 <- do.call(rbind.data.frame, (msigdb_enrichment(pathway_DB, pathways, background_genes,tmp.candidate2)))
  
  if(dim(tmp.msigdb2)[1]==0){
    tmp.msigdb2=as.data.frame(matrix(ncol=ncol(tmp.msigdb1),nrow=1))
    colnames(tmp.msigdb2)=colnames(tmp.msigdb1)
  }else if(dim(tmp.msigdb1)[1]==0){
    tmp.msigdb1=as.data.frame(matrix(ncol=ncol(tmp.msigdb2),nrow=1))
    colnames(tmp.msigdb1)=colnames(tmp.msigdb2)
  }
  tmp.msigdb1$group="B1"
  tmp.msigdb2$group="B2"
  tmp.msigdb=rbind(tmp.msigdb1,tmp.msigdb2)
  cat(yellow(tmp.bacteria,"\n"))
  
  return.string=data.frame(Bacteria=tmp.bacteria,number=tmp.msigdb$cont_n1,
                           pathway=tmp.msigdb$pathway,genes=tmp.msigdb$gene_names_in_pathway,odds=tmp.msigdb$odds_ratio,Pvalue=tmp.msigdb$p_val,group=tmp.msigdb$group)
  
}
node_enrich$FDR=p.adjust(node_enrich$Pvalue,method = "BH")
node_enrich=node_enrich[node_enrich$FDR<0.05,]
write.table(node_enrich,file = "OutputTable/Node_enrich.B1B2.txt",quote = F,row.names = F) # need manual extraction, top3

tmp.node=node_enrich[node_enrich$Bacteria=="Lachnoclostridium",]
tmp.node=tmp.node[!duplicated(tmp.node$pathway),]
tmp.pathway="REACTOME_OPIOID_SIGNALLING"
tmp.gene=unlist(strsplit((tmp.node$genes[tmp.node$pathway==tmp.pathway]), ","))
tmp_B1=edge_B1[edge_B1$from %in% tmp.gene & edge_B1$to %in% tmp.gene,]
tmp_B1=rbind(tmp_B1,edge_B1[edge_B1$from %in% tmp.gene & edge_B1$to=="Lachnoclostridium",])
tmp_B1=tmp_B1[tmp_B1$FDR<0.05,]
tmp_B1=tmp_B1[tmp_B1$from %in% colnames(matrix_B1_gene),]
tmp_B2=edge_B2[edge_B2$from %in% tmp.gene & edge_B2$to %in% tmp.gene,]
tmp_B2=rbind(tmp_B2,edge_B2[edge_B2$from %in% tmp.gene & edge_B2$to=="Lachnoclostridium",])
tmp_B2=tmp_B2[tmp_B2$FDR<0.05,]
tmp_B2=tmp_B2[tmp_B2$from %in% colnames(matrix_B1_gene),]
tmp_B1$Group="B1"
tmp_B2$Group="B2"
#write.table(rbind(tmp_B1,tmp_B2),file = "OutputTable/B1B2.Christensenellaceae_R-7_group.geneCluster.txt",row.names = F,quote = F)

tmp_B1$color[tmp_B1$cor>0]="#CCBB44"
tmp_B1$color[tmp_B1$cor<0]="#66CCEE"
tmp_B2$color[tmp_B2$cor>0]="#CCBB44"
tmp_B2$color[tmp_B2$cor<0]="#66CCEE"

tmp_node=unique(c(as.character(unique(tmp_B2$from)),as.character(unique(tmp_B2$to))))
tmp_node=data.frame(id=tmp_node,lable=tmp_node)
tmp_node$group=NA
tmp_node$group[tmp_node$lable %in% colnames(matrix_B1_bacteria)]="Bacteria"
tmp_node$group[tmp_node$lable %in% colnames(matrix_B1_gene)]="Gene"
tmp_node$group[tmp_node$lable %in% props_B1_gene$HubGene]="GeneHubs"

tmp_node$color=NA
tmp_node$color[tmp_node$group=="Gene"]="#BBBBBB"
tmp_node$color[tmp_node$group=="Bacteria"]="#994455"
tmp_node$color[tmp_node$group=="GeneHubs"]="#0077BB"

tmp_node$shape=NA
tmp_node$shape[tmp_node$group=="Gene"]="dot"
tmp_node$shape[tmp_node$group=="Bacteria"]="dot"
tmp_node$shape[tmp_node$group=="GeneHubs"]="dot"

tmp_node$font.size[tmp_node$group=="Bacteria"]=0
tmp_node$font.size[tmp_node$group=="GeneHubs"]=80
tmp_node$font.size[tmp_node$group=="Gene"]=0
tmp_node$size[tmp_node$group=="Bacteria"]=40
tmp_node$size[tmp_node$group=="Gene"]=0
tmp_node$size[tmp_node$group=="GeneHubs"]=40

network=visNetwork(tmp_node, tmp_B2) %>%
  visIgraphLayout(randomSeed = 1111,type = "full") %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges()
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.tmp.html")

tmp.bacteria="Lachnoclostridium"
tmp.gene="CD8A"
tmp.bacteria=bacteria[,tmp.bacteria,drop=F]
tmp.gene=genes_correct[,tmp.gene,drop=F]
tmp.data=merge(tmp.bacteria,tmp.gene,by="row.names",all=F)
tmp.data=merge(tmp.data,metadata_rna[,"MontrealB",drop=F],by.x="Row.names",by.y = "row.names",all=F)
tmp.data=na.omit(tmp.data)
tmp.data$MontrealB_class=NA
tmp.data$MontrealB_class[tmp.data$MontrealB==0]="B1"
tmp.data$MontrealB_class[tmp.data$MontrealB==3]="B1"
tmp.data$MontrealB_class[tmp.data$MontrealB==1]="B2"
tmp.data$MontrealB_class[tmp.data$MontrealB==4]="B2"
tmp.data$MontrealB_class[tmp.data$MontrealB==2]="B3"
tmp.data$MontrealB_class[tmp.data$MontrealB==5]="B3"
tmp.data=tmp.data[tmp.data$MontrealB_class!="B3",]

ggplot(tmp.data, aes(Lachnoclostridium,CD8A,fill=MontrealB_class)) +
  geom_point(shape = 21, color = "black", size = 1)+
  geom_smooth(method = lm,color="black",linetype="dashed")+
  facet_wrap(~ MontrealB_class,scales="free")+
  theme(legend.position="bottom")+theme_classic()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+scale_fill_manual(values = c("#44BB99","#AA4499"))+
  theme(strip.background = element_blank(), strip.text = element_blank())+guides(fill=F)
ggsave("OutputPlot/test.pdf",width = 2.6,height = 2)

#edge_B1$FDR=p.adjust(edge_B1$sig,method = "BH")
#edge_B1=edge_B1[edge_B1$FDR<0.05,]

#edge_B1=edge_B1[edge_B1$to %in% colnames(matrix_B1_bacteria),]
#edge_B1=edge_B1[!edge_B1$from %in% colnames(matrix_B1_bacteria),]
#edge_B1$pair=paste(edge_B1$from,edge_B1$to)

#edge_B2$pair=paste(edge_B2$from,edge_B2$to)
#edge_B2=edge_B2[edge_B2$pair %in% edge_B1$pair,]

#hetero_test=as.data.frame(matrix(nrow = nrow(edge_B1),ncol = 5))
#colnames(hetero_test)=c("from","to","B1.cor","B2.cor","Probability")
#for(i in 1:nrow(edge_B1)){
#  tmp.B1=edge_B1$cor[i]
#  tmp.B2=edge_B2$cor[i]
  
#  hetero_test$B1.cor[i]=edge_B1$cor[i]
#  hetero_test$B2.cor[i]=edge_B2$cor[i]
#  hetero_test$from[i]=as.character(edge_B1$from[i])
#  hetero_test$to[i]=as.character(edge_B1$to[i])
#  hetero_test$Probability[i]=r.test(n=204,n2=103,r12=tmp.B1,r34=tmp.B2)$p
#  cat(i,"\n") 
#}
#aa=(msigdb_enrichment(pathway_DB, pathways, hetero_test$from,hetero_test$from[hetero_test$Probability<0.05]))
#msigdb_pathways_df <- do.call(rbind.data.frame, aa)

#tmp.gene="PRKACA"
#tmp.bacteria="Lachnospira"
#tmp.B1=matrix_B1[,colnames(matrix_B1) %in% c(tmp.gene,tmp.bacteria)]
#tmp.B2=matrix_B2[,colnames(matrix_B2) %in% c(tmp.gene,tmp.bacteria)]
#tmp.B1$Group="B1"
#tmp.B2$Group="B2"
#tmp.data=rbind(tmp.B1,tmp.B2)

#ggplot(tmp.data, aes(PRKACA,Lachnospira)) +
#  geom_point(shape = 21, color = "#6699CC", size = 3)+
#  geom_smooth(method = lm,color="#6699CC",linetype="dashed")+
#  theme(legend.position="bottom")+facet_wrap(~ Group,scales = "free_x")+
#  theme(panel.background = element_rect(fill = "white", colour = "black"))

# distupted whole network using z-score
library(dcanr)
cor_calculation=rcorr(as.matrix(matrix_B1B2),type="spearman")
condition_B1B1=metadata_rna[,"MontrealB",drop=F]
condition_B1B1$MontrealB[condition_B1B1$MontrealB==0]="B1"
condition_B1B1$MontrealB[condition_B1B1$MontrealB==3]="B1"
condition_B1B1$MontrealB[condition_B1B1$MontrealB==1]="B2"
condition_B1B1$MontrealB[condition_B1B1$MontrealB==4]="B2"
condition_B1B1=condition_B1B1[condition_B1B1$MontrealB=="B1" | condition_B1B1$MontrealB=="B2",,drop=F]
condition_B1B1=na.omit(condition_B1B1)
condition_B1B1$MontrealB[condition_B1B1=="B1"]=1
condition_B1B1$MontrealB[condition_B1B1=="B2"]=2
condition_B1B1=condition_B1B1[order(rownames(condition_B1B1)),,drop=F]
matrix_B1B2=matrix_B1B2[rownames(matrix_B1B2) %in% rownames(condition_B1B1),]
matrix_B1B2=matrix_B1B2[order(rownames(matrix_B1B2)),]
condition_B1B1=as.factor(condition_B1B1$MontrealB)
matrix_B1B2=as.data.frame(t(matrix_B1B2))

z_scores <- dcScore(matrix_B1B2, condition_B1B1, dc.method = 'zscore', cor.method = 'spearman')
print(z_scores[1:5, 1:5])

z_scores=z_scores[rownames(z_scores) %in% colnames(matrix_B1_gene),colnames(z_scores) %in% colnames(matrix_B1_bacteria)]
sig_B1B2=setNames(melt(z_scores, na.rm=T), c('from', 'to', 'Zscore'))
sig_B1B2$sig=2*pnorm(-abs(sig_B1B2$Zscore))
sig_B1B2=sig_B1B2[sig_B1B2$sig<0.05,]

# prediction
library(caret)
library(pROC)
library(glmnet)
outcome=metadata_rna[,"MontrealB",drop=F]
outcome$MontrealB[outcome$MontrealB==0]="B1"
outcome$MontrealB[outcome$MontrealB==1]="B2"
outcome$MontrealB[outcome$MontrealB==3]="B1"
outcome$MontrealB[outcome$MontrealB==4]="B2"
outcome=outcome[outcome$MontrealB=="B1" | outcome$MontrealB=="B2",,drop=F]
outcome=na.omit(outcome)
predictors=matrix_B1B2[rownames(matrix_B1B2) %in% rownames(outcome),]
covariate=metadata_rna[rownames(metadata_rna) %in% rownames(outcome),c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough"),drop=F]
covariate$Location_rough[covariate$Location_rough=="ileum"]=0
covariate$Location_rough[covariate$Location_rough=="colon"]=1
covariate$Location_rough=as.numeric(covariate$Location_rough)
predictors=merge(predictors,covariate,by="row.names",all=F)
rownames(predictors)=predictors$Row.names
predictors$Row.names=NULL

predictors=predictors[order(rownames(predictors)),]
outcome=outcome[order(rownames(outcome)),,drop=F]

outcome[outcome=="B1"]=0
outcome[outcome=="B2"]=1
outcome$MontrealB=as.numeric(outcome$MontrealB)

set.seed(19910930)
inTrain <- createDataPartition(y=outcome$MontrealB,p=0.7,list=F)
trainSet_predict <- predictors[inTrain,]
trainSet_outcome <- outcome[inTrain,,drop=F]
trainSet=cbind(trainSet_outcome,trainSet_predict)

testSet_predict <- predictors[-inTrain,]
testSet_outcome <- outcome[-inTrain,,drop=F]
testSet=cbind(testSet_outcome,testSet_predict)

library(doParallel )
registerDoParallel(cores = 4)

# basic factors
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$MontrealB)
feature=as.matrix(trainSet_predict[,c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough")])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$MontrealB)
new_feature=as.matrix(testSet_predict[,c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough")])
modelroc1_train=roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot1_train=data.frame(FPR=1-modelroc1_train$specificities,TPR=modelroc1_train$sensitivities)
modelroc1_test=roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot1_test=data.frame(FPR=1-modelroc1_test$specificities,TPR=modelroc1_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# gene features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$MontrealB)
feature=as.matrix(trainSet_predict[,!colnames(trainSet_predict) %in% colnames(matrix_B1_bacteria)])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$MontrealB)
new_feature=as.matrix(testSet_predict[,!colnames(testSet_predict) %in% colnames(matrix_B1_bacteria)])
modelroc2_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot2_train=data.frame(FPR=1-modelroc2_train$specificities,TPR=modelroc2_train$sensitivities)
modelroc2_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot2_test=data.frame(FPR=1-modelroc2_test$specificities,TPR=modelroc2_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# bacteria features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$MontrealB)
feature=as.matrix(trainSet_predict[,!colnames(trainSet_predict) %in% colnames(matrix_B1_gene)])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$MontrealB)
new_feature=as.matrix(testSet_predict[,!colnames(testSet_predict) %in% colnames(matrix_B1_gene)])
modelroc3_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot3_train=data.frame(FPR=1-modelroc3_train$specificities,TPR=modelroc3_train$sensitivities)
modelroc3_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot3_test=data.frame(FPR=1-modelroc3_test$specificities,TPR=modelroc3_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# all features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$MontrealB)
feature=as.matrix(trainSet_predict)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$MontrealB)
new_feature=as.matrix(testSet_predict)
modelroc4_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot4_train=data.frame(FPR=1-modelroc4_train$specificities,TPR=modelroc4_train$sensitivities)
modelroc4_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot4_test=data.frame(FPR=1-modelroc4_test$specificities,TPR=modelroc4_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

modelroc_plot1_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation"
modelroc_plot2_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Gene"
modelroc_plot3_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Bacteria"
modelroc_plot4_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Bacteria+Gene"
modelroc_plot_test=rbind(modelroc_plot1_test,modelroc_plot2_test,modelroc_plot3_test,modelroc_plot4_test)
modelroc_plot_test=modelroc_plot_test[order(modelroc_plot_test$TPR),]
ggplot(modelroc_plot_test,aes(FPR,TPR,color=Set))+geom_line(size=1)+geom_segment(aes(x = 1, y = 1, xend = 0,yend = 0,linetype="dotted"))+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)")+scale_color_aaas()+
  theme(legend.position="bottom")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
ggsave("OutputPlot/prediction.B1B2.pdf",width = 3,height = 3)

# ==============================================================
# Biological use Network construction
# ==============================================================
library("metagMisc")
library(NetCoMi,verbose=F)
library("Rcpp")
library(igraph)
library("Hmisc")
library(visNetwork)
library(geomnet)
library(reshape2)

Biological_gene=read.table("OutputTable/biological_gene.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
Biological_bacteria=read.table("OutputTable/biological_bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

# Network construction in biological non-users
matrix_biological=merge(bacteria[,colnames(bacteria) %in% Biological_bacteria$Gene[Biological_bacteria$FDR<0.05]],
                  gene[,colnames(gene) %in% Biological_gene$Gene[Biological_gene$FDR<0.05]],by="row.names",all=F)
rownames(matrix_biological)=matrix_biological$Row.names
matrix_biological$Row.names=NULL
matrix_user=matrix_biological[rownames(matrix_biological) %in% rownames(metadata_rna)[metadata_rna$biological_use==1],]
matrix_nonuser=matrix_biological[rownames(matrix_biological) %in% rownames(metadata_rna)[metadata_rna$biological_use==0],]
matrix_nonuser_bacteria=bacteria[rownames(bacteria) %in% rownames(matrix_nonuser),colnames(bacteria) %in% colnames(matrix_nonuser)]
matrix_nonuser_gene=gene[rownames(gene) %in% rownames(matrix_nonuser),colnames(gene) %in% colnames(matrix_nonuser)]

net_nonuser_bacteria <- netConstruct(matrix_nonuser_bacteria,  
                                  measure = "spearman",
                                  normMethod = "none", 
                                  zeroMethod = "none",
                                  sparsMethod = "t-test", 
                                  adjust="BH",
                                  alpha=0.05,
                                  verbose = 3,seed = 100)
props_nonuser_bacteria <- netAnalyze(net_nonuser_bacteria, centrLCC = TRUE,
                                  clustMethod = "cluster_fast_greedy",
                                  hubPar = "eigenvector",connectivity=F,normNatConnect=F,
                                  weightDeg = FALSE, normDeg = FALSE)

sink("OutputTable/props_nonuser_bacteria.txt")
summary(props_nonuser_bacteria, numbNodes = 5L)
sink()
plot(props_nonuser_bacteria)

net_nonuser_gene <- netConstruct(matrix_nonuser_gene,  
                              measure = "spearman",
                              normMethod = "none", 
                              zeroMethod = "none",
                              sparsMethod = "threshold", 
                              thresh = 0.3,
                              verbose = 3,seed = 100)
#props_B1B2_gene <- netAnalyze(net_B1B2_gene, centrLCC = TRUE,
#                              clustMethod = "cluster_fast_greedy",
#                              hubPar = "eigenvector",
#                              weightDeg = FALSE, normDeg = FALSE)
#sink("OutputTable/props_biological_gene.txt")
#summary(props_biological_gene, numbNodes = 5L)
#sink()
props_nonuser_gene=read.table("OutputTable/props_biological.genes.extract.txt",header = T,stringsAsFactors = F)

tmp.gene="CXCL13"
tmp.bacteria="Faecalibacterium"
tmp.nonuser=matrix_nonuser[,colnames(matrix_nonuser) %in% c(tmp.gene,tmp.bacteria)]
tmp.user=matrix_user[,colnames(matrix_user) %in% c(tmp.gene,tmp.bacteria)]
tmp.nonuser$Group="Anti-TNF non-user"
tmp.user$Group="Anti-TNF user"
tmp.data=rbind(tmp.nonuser,tmp.user)

ggplot(tmp.data, aes(Group,`CXCL13`,fill=Group)) +
  geom_boxplot()+
  theme(legend.position="bottom")+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  scale_fill_manual(values=c( "#E69F00", "#56B4E9"))
ggsave("OutputPlot/test.pdf",width = 4,height = 4)

# view networks in nonusers
cor_calculation=rcorr(as.matrix(matrix_nonuser),type = "spearman")
cor_nonuser=cor_calculation$r
sig_nonuser=cor_calculation$P

cor_nonuser[upper.tri(cor_nonuser)] = NA
sig_nonuser[upper.tri(sig_nonuser)] = NA

sig_nonuser=setNames(melt(sig_nonuser, na.rm=T), c('from', 'to', 'sig'))
edge_nonuser=setNames(melt(cor_nonuser, na.rm=T), c('from', 'to', 'width'))
edge_nonuser=edge_nonuser[edge_nonuser$width!=1,]

edge_nonuser=cbind(edge_nonuser,sig_nonuser[,"sig",drop=F])
edge_nonuser$FDR=p.adjust(edge_nonuser$sig,method = "BH")
edge_nonuser=edge_nonuser[edge_nonuser$FDR<0.05,]

edge_nonuser$color=NA
edge_nonuser$color[edge_nonuser$width<0]="#99DDFF"
edge_nonuser$color[edge_nonuser$width>0]="#EEDD88"
edge_nonuser$width=abs(edge_nonuser$width)

tmp.bacteria=edge_nonuser[edge_nonuser$to %in% colnames(matrix_nonuser_bacteria),]
tmp.bacteria=tmp.bacteria[tmp.bacteria$from %in% colnames(matrix_nonuser_bacteria),]
tmp.gene=edge_nonuser[edge_nonuser$to %in% colnames(matrix_nonuser_bacteria),]
tmp.gene=tmp.gene[!tmp.gene$from %in% colnames(matrix_nonuser_bacteria),]
write.table(tmp.gene,file = "OutputTable/Network.biological_nonuse.host_microbe.txt",sep = "\t",row.names = F,quote = F)

edge=rbind(tmp.bacteria,tmp.gene)
edge$width=edge$width*10

node=unique(c(as.character(unique(edge$from)),as.character(unique(edge$to))))
node=data.frame(id=node,lable=node)
node$group=NA
node$group[node$lable %in% colnames(matrix_nonuser_bacteria)]="Bacteria"
node$group[node$lable %in% colnames(matrix_nonuser_gene)]="Gene"
node$group[node$lable %in% props_nonuser_gene$HubGene]="GeneHubs"

node$color=NA
node$color[node$group=="Gene"]="#BBBBBB"
node$color[node$group=="Bacteria"]="#994455"
node$color[node$group=="GeneHubs"]="#0077BB"

node$shape=NA
node$shape[node$group=="Gene"]="dot"
node$shape[node$group=="Bacteria"]="dot"
node$shape[node$group=="GeneHubs"]="dot"

node$font.size[node$group=="Bacteria"]=80
node$font.size[node$group=="GeneHubs"]=40
node$font.size[node$group=="Gene"]=0
node$size[node$group=="Bacteria"]=20
node$size[node$group=="Gene"]=5
node$size[node$group=="GeneHubs"]=20

network=visNetwork(node, edge) %>%
  visIgraphLayout(randomSeed = 1111,layout = "layout_in_circle",type = "full") %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges(smooth = list(roundness = 0.3))
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.biological.html")

#enrichment
background_genes=colnames(gene)

msigdb_C2 <- msigdbr(species = "Homo sapiens", category = "C2")
class(msigdb_C2) #[1] "tbl_df"     "tbl"        "data.frame"
msigdb_C2 <- as.data.frame(msigdb_C2)
table(msigdb_C2$gs_subcat)

## Only keep canonical pathways, or in other words, remove CGP (Chemical and genetic perturbations)
msigdb_C2_CP <- msigdb_C2[msigdb_C2$gs_subcat!= "CGP",]
table(msigdb_C2_CP$gs_subcat)

path_DB <- c("CP:REACTOME")
msigdb_C2_CP <- msigdb_C2_CP[(msigdb_C2_CP$gs_subcat %in% path_DB), ]
table(msigdb_C2_CP$gs_subcat)

msigdb_C2_CP_unique_pathways <- unique(msigdb_C2_CP$gs_name)
length(msigdb_C2_CP_unique_pathways) 

pathway_DB <- msigdb_C2_CP
pathways <- msigdb_C2_CP_unique_pathways

min_genes <- 10
max_genes <- 10000
overlap_genes <- 5
msigdb_collection_name <- "C2_CP" 

candidate_gene=unique(as.character(tmp.gene$from))

aa=(msigdb_enrichment(pathway_DB, pathways, background_genes,candidate_gene))
msigdb_pathways_df <- do.call(rbind.data.frame, aa)
msigdb_pathways_df=msigdb_pathways_df[order(msigdb_pathways_df$p_val),]
write.table(msigdb_pathways_df,file = "OutputTable/enrichment.B1B2.gene.microbe.txt",row.names = F,quote = F)

# network compare
library(psych)
cor_calculation=rcorr(as.matrix(matrix_user),type="spearman")
cor_user=cor_calculation$r
sig_user=cor_calculation$P

cor_user[upper.tri(cor_user)] = NA
sig_user[upper.tri(sig_user)] = NA

sig_user=setNames(melt(sig_user, na.rm=T), c('from', 'to', 'sig'))
edge_user=setNames(melt(cor_user, na.rm=T), c('from', 'to', 'cor'))
edge_user=edge_user[edge_user$cor!=1,]
edge_user=cbind(edge_user,sig_user[,"sig",drop=F])

cor_calculation=rcorr(as.matrix(matrix_nonuser),type="spearman")
cor_nonuser=cor_calculation$r
sig_nonuser=cor_calculation$P

cor_nonuser[upper.tri(cor_nonuser)] = NA
sig_nonuser[upper.tri(sig_nonuser)] = NA

sig_nonuser=setNames(melt(sig_nonuser, na.rm=T), c('from', 'to', 'sig'))
edge_nonuser=setNames(melt(cor_nonuser, na.rm=T), c('from', 'to', 'cor'))
edge_nonuser=edge_nonuser[edge_nonuser$cor!=1,]
edge_nonuser=cbind(edge_nonuser,sig_nonuser[,"sig",drop=F])

edge_user$Zscore=fisherz(edge_user$cor)
edge_nonuser$Zscore=fisherz(edge_nonuser$cor)
edge_user$FDR=p.adjust(edge_user$sig,method = "BH")
edge_nonuser$FDR=p.adjust(edge_nonuser$sig,method = "BH")

write.table(edge_user[edge_user$FDR<0.05 & edge_user$to %in% colnames(matrix_nonuser_bacteria),],"OutputTable/Anti-TNF-user.microbiota-gene.clusters.txt",sep = "\t",quote = F,row.names = F)
write.table(edge_nonuser[edge_nonuser$FDR<0.05 & edge_nonuser$to %in% colnames(matrix_nonuser_bacteria),],"OutputTable/Anti-TNF-nonuser.microbiota-gene.clusters.txt",sep = "\t",quote = F,row.names = F)

node_compare=matrix(nrow = ncol(matrix_nonuser_bacteria),ncol = 4)
node_compare=as.data.frame(node_compare)
colnames(node_compare)=c("Bacteria","user.pre","nonuser.pre","Pvalue")
for(i in 1:ncol(matrix_nonuser_bacteria)){
  tmp.bacteria=colnames(matrix_nonuser_bacteria)[i]
  tmp.nonuser=edge_nonuser[edge_nonuser$to==tmp.bacteria,]
  tmp.nonuser=tmp.nonuser[tmp.nonuser$from %in% colnames(matrix_nonuser_gene),]
  tmp.nonuser$pair=paste(tmp.nonuser$from,tmp.nonuser$to)
  
  tmp.user=edge_user[edge_user$to==tmp.bacteria,]
  tmp.user=tmp.user[tmp.user$from %in% colnames(matrix_nonuser_gene),]
  tmp.user$pair=paste(tmp.user$from,tmp.user$to)
  
  tmp.nonuser=tmp.nonuser[order(tmp.nonuser$Zscore),]
  tmp.user=tmp.user[order(tmp.user$Zscore),]
  
  tmp.nonuser$rank=1:nrow(tmp.nonuser)
  tmp.user$rank=1:nrow(tmp.user)
  
  tmp.gene=union(tmp.nonuser$from[tmp.nonuser$FDR<0.05],tmp.user$from[tmp.user$FDR<0.05])
  tmp.nonuser=tmp.nonuser[tmp.nonuser$from %in% tmp.gene,]
  tmp.user=tmp.user[tmp.user$from %in% tmp.gene,]
  
  tmp.nonuser=tmp.nonuser[order(tmp.nonuser$from),]
  tmp.user=tmp.user[order(tmp.user$from),]
  
  mm=(wilcox.test(tmp.nonuser$Zscore,tmp.user$Zscore,paired=T))$p.value
  
  tmp.user.pre=bacteria[rownames(bacteria) %in% rownames(matrix_user),tmp.bacteria]
  tmp.user.pre=length(tmp.user.pre[tmp.user.pre>0])/length(tmp.user.pre)
  tmp.nonuser.pre=bacteria[rownames(bacteria) %in% rownames(matrix_nonuser),tmp.bacteria]
  tmp.nonuser.pre=length(tmp.nonuser.pre[tmp.nonuser.pre>0])/length(tmp.nonuser.pre)
  
  node_compare$Bacteria[i]=tmp.bacteria
  node_compare$user.pre[i]=tmp.user.pre
  node_compare$nonuser.pre[i]=tmp.nonuser.pre
  node_compare$Pvalue[i]=mm
  
}
node_compare$FDR=p.adjust(node_compare$Pvalue,method = "BH")

node_enrich = foreach(i=1:length(node_compare$Bacteria[node_compare$FDR<0.05]),.combine = rbind) %do%  {
  tmp.bacteria=node_compare$Bacteria[node_compare$FDR<0.05][i]
  tmp_nonuser=edge_nonuser[edge_nonuser$to==tmp.bacteria,]
  tmp_nonuser=tmp_nonuser[tmp_nonuser$FDR<0.05,]
  tmp_nonuser=tmp_nonuser[tmp_nonuser$from %in% colnames(matrix_nonuser_gene),]
  tmp_user=edge_user[edge_user$to==tmp.bacteria,]
  tmp_user=tmp_user[tmp_user$FDR<0.05,]
  tmp_user=tmp_user[tmp_user$from %in% colnames(matrix_nonuser_gene),]
  
  background_genes=colnames(gene)
  tmp.candidate1=unique(as.character(tmp_nonuser$from))
  tmp.candidate2=unique(as.character(tmp_user$from))
  
  tmp.msigdnonuser <- do.call(rbind.data.frame, (msigdb_enrichment(pathway_DB, pathways, background_genes,tmp.candidate1)))
  tmp.msigduser <- do.call(rbind.data.frame, (msigdb_enrichment(pathway_DB, pathways, background_genes,tmp.candidate2)))
  
  if(dim(tmp.msigduser)[1]==0){
    tmp.msigduser=as.data.frame(matrix(ncol=ncol(tmp.msigdnonuser),nrow=1))
    colnames(tmp.msigduser)=colnames(tmp.msigdnonuser)
  }else if(dim(tmp.msigdnonuser)[1]==0){
    tmp.msigdnonuser=as.data.frame(matrix(ncol=ncol(tmp.msigduser),nrow=1))
    colnames(tmp.msigdnonuser)=colnames(tmp.msigduser)
  }
  tmp.msigdnonuser$group="nonuser"
  tmp.msigduser$group="user"
  tmp.msigdb=rbind(tmp.msigdnonuser,tmp.msigduser)
  cat(yellow(tmp.bacteria,"\n"))
  
  return.string=data.frame(Bacteria=tmp.bacteria,
                           pathway=tmp.msigdb$pathway,genes=tmp.msigdb$gene_names_in_pathway,odds=tmp.msigdb$odds_ratio,Pvalue=tmp.msigdb$p_val,group=tmp.msigdb$group)
  
}
node_enrich$FDR=p.adjust(node_enrich$Pvalue,method = "BH")
node_enrich=node_enrich[node_enrich$FDR<0.05,]
write.table(node_enrich,file = "OutputTable/Node_enrich.biological.txt",quote = F,row.names = F) # need manual extraction, top3

node_enrich=read.table("OutputTable/Node_enrich.biological.txt",sep = "\t",header = T,stringsAsFactors = F)
node_enrich$path=factor(node_enrich$path,levels = c("path3","path2","path1"))
node_enrich$group=factor(node_enrich$group,levels = c("nonuser","user"))
tmp.data=node_enrich[node_enrich$Bacteria %in% c("Veillonella","Dorea","Subdoligranulum"),]
tmp.data$odds[tmp.data$group=="nonuser"]=-tmp.data$odds[tmp.data$group=="nonuser"]
ggplot(data=tmp.data, aes(x=Bacteria, y=odds,group=path,fill=Bacteria)) +
  geom_bar(stat = "identity", position=position_dodge(),color="black")+
  coord_flip()+
  theme_bw()+
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  theme(strip.background = element_blank(), strip.text = element_blank())+
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.background = element_blank())+
  scale_fill_brewer(palette="Oranges")
ggsave("OutputPlot/Node_enrich.biological.pdf",width = 8,height = 4)

tmp.node=node_enrich[ node_enrich$Bacteria=="Ruminococcaceae_UCG-002",]
tmp.node=tmp.node[!duplicated(tmp.node$pathway),]
tmp.pathway="REACTOME_GLYCEROPHOSPHOLIPID_BIOSYNTHESIS"
tmp.gene=unlist(strsplit((tmp.node$genes[tmp.node$pathway==tmp.pathway]), ","))
tmp_nonuser=edge_nonuser[edge_nonuser$from %in% tmp.gene & edge_nonuser$to %in% tmp.gene,]
tmp_nonuser=rbind(tmp_nonuser,edge_nonuser[edge_nonuser$from %in% tmp.gene & edge_nonuser$to=="Ruminococcaceae_UCG-002",])
tmp_nonuser=tmp_nonuser[tmp_nonuser$FDR<0.05,]
tmp_nonuser=tmp_nonuser[tmp_nonuser$from %in% colnames(matrix_nonuser_gene),]
tmp_user=edge_user[edge_user$from %in% tmp.gene & edge_user$to %in% tmp.gene,]
tmp_user=rbind(tmp_user,edge_user[edge_user$from %in% tmp.gene & edge_user$to=="Ruminococcaceae_UCG-002",])
tmp_user=tmp_user[tmp_user$FDR<0.05,]
tmp_user=tmp_user[tmp_user$from %in% colnames(matrix_nonuser_gene),]
tmp_nonuser$Group="nonuser"
tmp_user$Group="user"
#write.table(rbind(tmp_nonuser,tmp_user),file = "OutputTable/nonuseruser.Veillonella.geneCluster.txt",row.names = F,quote = F)

tmp_nonuser$color[tmp_nonuser$cor>0]="#CCBB44"
tmp_nonuser$color[tmp_nonuser$cor<0]="#66CCEE"
tmp_user$color[tmp_user$cor>0]="#CCBB44"
tmp_user$color[tmp_user$cor<0]="#66CCEE"

tmp_node=unique(c(as.character(unique(tmp_user$from)),as.character(unique(tmp_user$to))))
tmp_node=data.frame(id=tmp_node,lable=tmp_node)
tmp_node$group=NA
tmp_node$group[tmp_node$lable %in% colnames(matrix_nonuser_bacteria)]="Bacteria"
tmp_node$group[tmp_node$lable %in% colnames(matrix_nonuser_gene)]="Gene"
tmp_node$group[tmp_node$lable %in% props_nonuser_gene$HubGene]="GeneHubs"

tmp_node$color=NA
tmp_node$color[tmp_node$group=="Gene"]="#BBBBBB"
tmp_node$color[tmp_node$group=="Bacteria"]="#994455"
tmp_node$color[tmp_node$group=="GeneHubs"]="#0077BB"

tmp_node$shape=NA
tmp_node$shape[tmp_node$group=="Gene"]="dot"
tmp_node$shape[tmp_node$group=="Bacteria"]="dot"
tmp_node$shape[tmp_node$group=="GeneHubs"]="dot"

tmp_node$font.size[tmp_node$group=="Bacteria"]=0
tmp_node$font.size[tmp_node$group=="GeneHubs"]=80
tmp_node$font.size[tmp_node$group=="Gene"]=0
tmp_node$size[tmp_node$group=="Bacteria"]=40
tmp_node$size[tmp_node$group=="Gene"]=0
tmp_node$size[tmp_node$group=="GeneHubs"]=40

network=visNetwork(tmp_node, tmp_user) %>%
  visIgraphLayout(randomSeed = 1111,type = "full") %>%
  visLayout(improvedLayout = T) %>%
  visOptions(highlightNearest = list(enabled = T, degree = 1, hover = T)) %>% 
  visEdges()
visSave(network, file = "/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/OutputPlot/network.tmp.html")

tmp.bacteria="Ruminococcaceae_UCG-002"
tmp.gene="ACAA1"
tmp.bacteria=bacteria[,tmp.bacteria,drop=F]
tmp.gene=genes_correct[,tmp.gene,drop=F]
tmp.data=merge(tmp.bacteria,tmp.gene,by="row.names",all=F)
tmp.data=merge(tmp.data,metadata_rna[,"biological_use",drop=F],by.x="Row.names",by.y = "row.names",all=F)
tmp.data=na.omit(tmp.data)
tmp.data$biological_use[tmp.data$biological_use==0]="Nonuser"
tmp.data$biological_use[tmp.data$biological_use==1]="User"
ggplot(tmp.data, aes(`Ruminococcaceae_UCG-002`,ACAA1,fill=biological_use)) +
  geom_point(shape = 21, color = "black", size = 1)+
  geom_smooth(method = lm,color="black",linetype="dashed")+
  facet_wrap(~ biological_use,scales="free")+
  theme(legend.position="bottom")+theme_classic()+
  theme(panel.background = element_rect(fill = "white", colour = "black"))+scale_fill_manual(values = c("#44BB99","#AA4499"))+
  theme(strip.background = element_blank(), strip.text = element_blank())+guides(fill=F)
ggsave("OutputPlot/test.pdf",width = 2.6,height = 2)

# prediction
library(caret)
library(pROC)
library(glmnet)
outcome=metadata_rna[,"biological_use",drop=F]
outcome=na.omit(outcome)
predictors=matrix_biological[rownames(matrix_biological) %in% rownames(outcome),]
covariate=metadata_rna[rownames(metadata_rna) %in% rownames(outcome),c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough"),drop=F]
covariate$Location_rough[covariate$Location_rough=="ileum"]=0
covariate$Location_rough[covariate$Location_rough=="colon"]=1
covariate$Location_rough=as.numeric(covariate$Location_rough)
predictors=merge(predictors,covariate,by="row.names",all=F)
rownames(predictors)=predictors$Row.names
predictors$Row.names=NULL

predictors=predictors[order(rownames(predictors)),]
outcome=outcome[order(rownames(outcome)),,drop=F]

outcome$biological_use=as.numeric(outcome$biological_use)

set.seed(19910930)
inTrain <- createDataPartition(y=outcome$biological_use,p=0.7,list=F)
trainSet_predict <- predictors[inTrain,]
trainSet_outcome <- outcome[inTrain,,drop=F]
trainSet=cbind(trainSet_outcome,trainSet_predict)

testSet_predict <- predictors[-inTrain,]
testSet_outcome <- outcome[-inTrain,,drop=F]
testSet=cbind(testSet_outcome,testSet_predict)

library(doParallel )
registerDoParallel(cores = 4)

# basic factors
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$biological_use)
feature=as.matrix(trainSet_predict[,c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough")])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$biological_use)
new_feature=as.matrix(testSet_predict[,c("age_at_biopsy","sex","BMI","smoking_DB","Inflammation","Location_rough")])
modelroc1_train=roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot1_train=data.frame(FPR=1-modelroc1_train$specificities,TPR=modelroc1_train$sensitivities)
modelroc1_test=roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot1_test=data.frame(FPR=1-modelroc1_test$specificities,TPR=modelroc1_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# gene features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$biological_use)
feature=as.matrix(trainSet_predict[,!colnames(trainSet_predict) %in% colnames(matrix_nonuser_bacteria)])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$biological_use)
new_feature=as.matrix(testSet_predict[,!colnames(testSet_predict) %in% colnames(matrix_nonuser_bacteria)])
modelroc2_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot2_train=data.frame(FPR=1-modelroc2_train$specificities,TPR=modelroc2_train$sensitivities)
modelroc2_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot2_test=data.frame(FPR=1-modelroc2_test$specificities,TPR=modelroc2_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# bacteria features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$biological_use)
feature=as.matrix(trainSet_predict[,!colnames(trainSet_predict) %in% colnames(matrix_nonuser_gene)])
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$biological_use)
new_feature=as.matrix(testSet_predict[,!colnames(testSet_predict) %in% colnames(matrix_nonuser_gene)])
modelroc3_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot3_train=data.frame(FPR=1-modelroc3_train$specificities,TPR=modelroc3_train$sensitivities)
modelroc3_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot3_test=data.frame(FPR=1-modelroc3_test$specificities,TPR=modelroc3_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

# all features
a <- seq(0.1, 0.9, 0.05)
response=as.matrix(trainSet_outcome$biological_use)
feature=as.matrix(trainSet_predict)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(feature, response, family = "binomial", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(feature, response, family = "binomial", lambda = cv3$lambda.1se, alpha = cv3$alpha)
selection=as.matrix(coef(md3))
selection=selection[selection[,1]!=0,,drop=F]
print(selection[-1,])
new_response=as.matrix(testSet_outcome$biological_use)
new_feature=as.matrix(testSet_predict)
modelroc4_train=roc(as.numeric(response), as.numeric(predict(md3, feature, type = "response")))
modelroc_plot4_train=data.frame(FPR=1-modelroc4_train$specificities,TPR=modelroc4_train$sensitivities)
modelroc4_test=roc(as.numeric(new_response), as.numeric(predict(md3, new_feature, type = "response")))
modelroc_plot4_test=data.frame(FPR=1-modelroc4_test$specificities,TPR=modelroc4_test$sensitivities)
roc(as.vector(response), as.numeric(predict(md3, feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )
roc(as.vector(new_response), as.numeric(predict(md3, new_feature, type = "response")), percent=F, boot.n=1000, ci.alpha=0.9, plot=TRUE, grid=TRUE, show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
    print.auc = TRUE, print.thres.col = "blue", ci=TRUE )

modelroc_plot1_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation"
modelroc_plot2_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Gene"
modelroc_plot3_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Bacteria"
modelroc_plot4_test$Set="Age+Gender+BMI+Smoking+Inflammation+TissueLocation+Bacteria+Gene"
modelroc_plot_test=rbind(modelroc_plot1_test,modelroc_plot2_test,modelroc_plot3_test,modelroc_plot4_test)
modelroc_plot_test=modelroc_plot_test[order(modelroc_plot_test$TPR),]
ggplot(modelroc_plot_test,aes(FPR,TPR,color=Set))+geom_line(size=2)+geom_segment(aes(x = 1, y = 1, xend = 0,yend = 0,linetype="dotted"))+
  labs(title= "ROC curve", 
       x = "False Positive Rate (1-Specificity)", 
       y = "True Positive Rate (Sensitivity)")+theme_minimal()+scale_color_aaas()+
  theme(legend.position="bottom")
ggsave("OutputPlot/prediction.biological.pdf")











































