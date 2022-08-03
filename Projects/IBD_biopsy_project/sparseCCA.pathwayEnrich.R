# sparseCCA genes pathways enrichemnt
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


# ======================================================================================================================================
# enrichment
# ======================================================================================================================================

# backgroud genes
genes=read.table("OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",header = T,row.names = 1,check.names = F)
genes=unique(colnames(genes))

inflammation=read.table("sparseCCA/Inflammation.genes.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)
inflammation=colnames(inflammation)

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
background_genes <- inflammation # define backgroud genes

min_genes <- 10
max_genes <- 300
overlap_genes <- 5
msigdb_collection_name <- "C2_CP" 

aa=read.table("/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/sparseCCAoutput//sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_1.txt",sep = "\t",
              header = T)
msigdb_enrichment(pathway_DB, pathways, background_genes, aa$gene)

####
input_dir="/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/"
filenames <- list.files(input_dir)
output_dir="/Users/hushixian/Desktop/PhD/700 RNA+16s project/Analysis/sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/"

enrichment_list <- perform_enrichment(input_dir, filenames, 
                                      msigdb_collection_name, pathway_DB, pathways, background_genes, 
                                      output_dir)

# ======================================================================================================================================
# pathway compare
# ======================================================================================================================================
# 
tmp1=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_1.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp2=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_3.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp3=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_5.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp4=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_7.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp5=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_8.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp6=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/msigdb_gene_taxa_component_9.txt",sep = "\t",header = T,stringsAsFactors = F)

tmp1$FDR=p.adjust(tmp1$p_val,method = "BH")
tmp2$FDR=p.adjust(tmp2$p_val,method = "BH")
tmp3$FDR=p.adjust(tmp3$p_val,method = "BH")
tmp4$FDR=p.adjust(tmp4$p_val,method = "BH")
tmp5$FDR=p.adjust(tmp5$p_val,method = "BH")
tmp6$FDR=p.adjust(tmp6$p_val,method = "BH")
pathway=rbind(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6)

# ======================================================================================================================================
# sparseCCA plots
# ======================================================================================================================================
library(ggplot2)
library(pheatmap)
library(grid)  
draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = -0.5, hjust = 1, rot = 135, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

# re-load files
bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(bacteria),]
inflammation=read.table("sparseCCA/Inflammation.genes.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)
inflammation=inflammation[rownames(inflammation) %in% rownames(bacteria),]

bacteria=bacteria[rownames(bacteria) %in% rownames(inflammation),]
inflammation=inflammation[rownames(inflammation) %in% rownames(bacteria),]

dim(bacteria)
dim(inflammation)

bacteria=bacteria[order(rownames(bacteria)),]
inflammation=inflammation[order(rownames(inflammation)),]
stopifnot(all(rownames(inflammation) == rownames(bacteria)))

component_gene=read.table("sparseCCAoutput/CCA_var_genes.txt",sep = "\t",header = T,row.names = 1)
component_bacteria=read.table("sparseCCAoutput/CCA_var_microbes.txt",sep = "\t",header = T,row.names = 1)

# gene_taxa_component_5.txt
tmp_plot=cbind(component_gene[,"V9",drop=F],component_bacteria[,"V9",drop=F])
colnames(tmp_plot)=c("Gene_Component","Microbe_Component")
ggplot(tmp_plot, aes(Microbe_Component, Gene_Component)) +
  geom_point(shape = 21, color = "black", size = 2,fill="#F4A582")+
  geom_smooth(method = lm,color="#F4A582")+
  theme_bw()+ theme(legend.position="bottom")+guides(color=guide_legend(nrow = 1))+theme_classic()
ggsave("OutputPlot/sparseCCA.plot8.pdf",width = 3,height = 5)
cor.test(tmp_plot$Gene_Component,tmp_plot$Microbe_Component)

cor.coef=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_8.txt",header = T)
cor.coef.gene=cor.coef[,c("gene","gene_coeff"),drop=F]
rownames(cor.coef.gene)=cor.coef.gene$gene
cor.coef.gene$gene=NULL
cor.coef.gene=cor.coef.gene[order(cor.coef.gene$gene_coeff,decreasing = T),,drop=F]
pdf("OutputPlot/sparseCCA.plot5.1.pdf",width = 2,height = 10)
pheatmap(cor.coef.gene,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = F,
         cellheight = 2,cellwidth = 70,color = colorRampPalette(c("#F6C141", "white", "#4EB265"))(100))
dev.off()

cor.coef.gene$name=rownames(cor.coef.gene)
cor.coef.gene=cor.coef.gene[order(cor.coef.gene$gene_coeff),]
cor.coef.gene$name=factor(cor.coef.gene$name,levels = cor.coef.gene$name)
ggplot(cor.coef.gene, aes(x=name, y=gene_coeff)) +
  geom_bar(stat="identity", color="#CCBB44", fill="white",width=0.05)+xlab("")+ylab("Correlation Coefficient")+coord_flip()+theme_classic()+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
ggsave("OutputPlot/sparseCCA.plot8.1.pdf",width = 3,height = 5)

cor.coef.taxa=cor.coef[,c("taxa","taxa_coeff"),drop=F]
cor.coef.taxa=na.omit(cor.coef.taxa)
rownames(cor.coef.taxa)=cor.coef.taxa$taxa
cor.coef.taxa$taxa=NULL
cor.coef.taxa=cor.coef.taxa[order(cor.coef.taxa$taxa_coeff,decreasing = T),,drop=F]
cor.coef.taxa=as.data.frame(t(cor.coef.taxa))
pdf("OutputPlot/sparseCCA.plot5.2.pdf",width = 8,height = 3)
pheatmap(cor.coef.taxa,cluster_rows = F,cluster_cols = F,show_rownames = F,show_colnames = T,
         cellheight = 20,cellwidth = 20,color = colorRampPalette(c("#2166AC", "#D1E5F0"))(100),legend = F)
dev.off()

cor.coef.taxa=as.data.frame(t(cor.coef.taxa))
cor.coef.taxa$name=rownames(cor.coef.taxa)
cor.coef.taxa=cor.coef.taxa[cor.coef.taxa$name!="Actinobacteria.1",]
cor.coef.taxa=cor.coef.taxa[order(cor.coef.taxa$taxa_coeff),]
cor.coef.taxa$name=factor(cor.coef.taxa$name,levels = cor.coef.taxa$name)
ggplot(cor.coef.taxa, aes(x=name, y=taxa_coeff)) +
  geom_bar(stat="identity", color="grey", fill="#0077BB",width=0.8)+xlab("")+ylab("Correlation Coefficient")+scale_y_reverse()+theme_classic()
ggsave("OutputPlot/sparseCCA.plot8.2.pdf",width = 3,height = 2)

tmp_data=cbind(bacteria[,"Erysipelotrichaceae_UCG.003",drop=F],inflammation[,"COL1A2",drop=F])
ggplot(tmp_data, aes(Erysipelotrichaceae_UCG.003, COL1A2)) +
  geom_point(shape = 21, color = "black", size = 1,fill="#9970AB")+
  geom_smooth(method = lm,color="#9970AB")+
  theme_classic()
ggsave("OutputPlot/test.pdf",width = 1.6,height = 2)
cor.test(tmp_data[,1],tmp_data[,2])

# ======================================================================================================================================
# pathway， individual compare
# ======================================================================================================================================
# 
tmp1=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_1.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp2=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_3.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp3=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_5.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp4=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_7.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp5=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_8.txt",sep = "\t",header = T,stringsAsFactors = F)
tmp6=read.table("sparseCCAoutput/sig_gene_taxa_components_0.4_0.205555555555556_padj/gene_taxa_component_9.txt",sep = "\t",header = T,stringsAsFactors = F)

pathway_gene=c(tmp1$gene,tmp2$gene,tmp3$gene,tmp4$gene,tmp5$gene,tmp6$gene)
pathway_gene=unique(pathway_gene)

individual_gene=read.table("OutputTable/gene.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F)
individual_gene=individual_gene$Gene[individual_gene$FDR<0.05]
individual_gene=unique(individual_gene)

length(intersect(individual_gene,pathway_gene))/length(individual_gene)
write.table((intersect(individual_gene,pathway_gene)),file = "OutputTable/test.txt",quote = F,row.names = F)






























































