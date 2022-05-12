#!/usr/bin/env Rscript

################################## Load R packages ################################
print('Loading R packages...')
shhh(library(devtools))
shhh(library(MAST))
shhh(library(Seurat))
shhh(library(data.table))
shhh(library(lme4))
shhh(library(plyr))
shhh(library(dplyr))
shhh(library(reshape2))
shhh(library(stringi))
shhh(library(stringr))
shhh(library(RCurl))

#################################### Functions ####################################
# 1. Pre-processing function
preprocess.func <- function(so, freq_expressed){
  # Cleaning the metadata
  print('Cleaning the metadata slot from the seurat object...')
  rownames(so@meta.data) <- colnames(so)
  so@meta.data$bare_barcode_lane <- rownames(so@meta.data)
  so@meta.data <- so@meta.data %>% mutate_if(is.character,as.factor)
  
  # Pre-processing
  ## Convert Seurat object to SingleCellExperiment object; and then from SingleCellExperiment object to SingleCellAssay object
  print('Converting seurat object to SingleCellAssay...')
  DefaultAssay(so) <- "RNA"
  sca <- as(as.SingleCellExperiment(so), 'SingleCellAssay')
  scaRaw <- sca
  
  ## Filtering (based on https://www.bioconductor.org/help/course-materials/2016/BioC2016/ConcurrentWorkshops2/McDavid/MAITAnalysis.html)
  ### Filter out lowly variable genes
  print('Filtering out lowly variable genes...')
  sca <- sca[freq(sca)>0,] #returns the frequency of expression, i.e., the proportion of non-zero values
  
  ### The number of genes detected in a sample (cngeneson), also called the cellular detection rate (cdr), is often the first principal component. 
  print('Calculating the cellular detection rate (cdr)...')
  colData(sca)$wellKey <- colnames(sca)
  rowData(sca)$primerid <- rownames(rowData(sca))
  colData(sca) <- droplevels(colData(sca))
  cdr2 <- colSums(assay(sca)>0) #assay(sca) is working on assays(sca)$logcounts
  colData(sca)$cngeneson <- scale(cdr2)
  
  ### Filter out lowly expressed genes
  print(paste0('Filtering out lowly expressed genes... Minimum expression threshold: ', freq_expressed))
  expressed_genes <- freq(sca) > freq_expressed
  sca <- sca[expressed_genes,]
  print(paste0('# of cells: ', ncol(sca)))
  print(paste0('# of genes tested: ', nrow(sca)))
  return(sca)
}

# 2. single-cell differential expression analysis function (scDEA): scDEA with MAST glmer (zlm) and testing the significance of the interesting variable with LRT
de_glmer.func <- function(sca_object, contrast, fixed_effects, random_effects, nagq0){
  print(paste0('Testing: ', contrast))
  random_effects.fmla <- paste(paste0('(1|',random_effects,')'),collapse='+') #try other nomenclature (=option 1, default)
  contrast_fixed.fmla <- paste(c(contrast,fixed_effects),collapse='+')
  zlm_vars <- paste0('~',paste(c(contrast_fixed.fmla,random_effects.fmla), collapse='+'))
  zlm_formula <- as.formula(zlm_vars)
  
  print(paste0('Fitting glmer: ',zlm_vars))
  if(nagq0){
    print(paste0('Trying to mitigate model failing to converge for some genes by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))'))
    zlmCond <- zlm(zlm_formula, 
                   sca_object,
                   method='glmer', ebayes=FALSE, fitArgsD=list(nAGQ=0))
    summaryCond <- summary(zlmCond, doLRT=contrast, fitArgsD=list(nAGQ=0))
  }else{
    zlmCond <- zlm(zlm_formula, 
                   sca_object,
                   method='glmer', ebayes=FALSE)
    summaryCond <- summary(zlmCond, doLRT=contrast)
  }
  summaryDt <- summaryCond$datatable
  return(summaryDt)
}

# 3. Function to extract the Hurdle component scDEA summary statistics
de_stats.func <- function(summaryDt, contrast_in){
  fcHurdle <- merge(summaryDt[contrast==contrast_in & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast==contrast_in & component=='logFC', .(primerid, coef, ci.hi, ci.lo, z)], #logFC coefficients
                    by = 'primerid')
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  fcHurdleSig <- merge(fcHurdle[fdr<.05],
                       as.data.table(mcols(sca)),
                       by='primerid')
  setorder(fcHurdleSig, fdr)
  setorder(fcHurdle, fdr)
  fcHurdle_df <- as.data.frame(fcHurdle)
  fcHurdle_df <- fcHurdle_df[,c(1,3:6,2,7)]
  return(fcHurdle_df)
}

# 4. Function to read and reformat by cell
reformat_by_cell.func <- function(cell_type, in_dir, fn_pattern){
  print(cell_type)
  in_sdir <- paste0(in_dir, cell_type, '/')
  in_fn <- paste0(in_sdir, fn_pattern)
  df <- readRDS(in_fn)
  df$Cell_type <- cell_type
  df_out <- df[,c(1,6,2,7)]
  colnames(df_out) <- c('gene_symbol','p_val','logFC','p_val_adj')
  df_list <- list(raw=df,
                  formatted=df_out)
  return(df_list)
}

# 5. Function to read and reformat
reformat.func <- function(in_dir, nagq0){
  fn_pattern <- ifelse(nagq0, 'de_glmer_stats.nagq0.rds', 'de_glmer_stats.rds')
  cell_types <- list.dirs(in_dir, full.names = F, recursive = F)
  out_list <- sapply(cell_types, function(ct) reformat_by_cell.func(cell_type = ct, 
                                                                    in_dir = in_dir,
                                                                    fn_pattern = fn_pattern), simplify = F)
  return(out_list)
}
