#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("-i","--in_dir"), action="store", default="https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/", type='character',
              help="Input directory (path/to/in_dir) with the files structured as in the example from https://downloads.molgeniscloud.org/downloads/combio_andreu_2022/."),
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="Cell level of interest. Example: approach_I (high resolution level) or approach II (low resolution level)."),
  make_option(c("--cell_type"), action="store", default=NA, type='character',
              help="Cell type of interest. Example: CD8T_memory (approach I) or T_cells (approach II)."),
  make_option(c("--covs"), action="store", default=NA, type='character',
              help="Covariates file. Example: covariates.approach_I.tab (approach I) or covariates.approach_II.tab (approach I)."),
  make_option(c("--genes"), action="store", default=NULL, type='character',
              help="Subset the original seurat object to perform the scDEA only on those genes. Example: top10genes.approach_I.tab (approach I) or covariates.approach_II.tab (approach I)."),
  make_option(c("--nagq0"), action="store", default=TRUE, type='logical',
              help="Mitigate model fails to converge by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))."),
  make_option(c("--freq"), action="store", default=0.1, type='double',
              help="Minimum expression threshold to filter out lowly expressed genes."),
  make_option(c("-o","--out_dir"), action="store", default=NA, type='character',
              help="Output directory. Example: /path/to/out_dir"))
opt = parse_args(OptionParser(option_list=option_list))

################################ Set Variables #################################
# Main directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions_fn <- 'scripts/accessory_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Input files
## Seurat object
cell_level <- opt$cell_level
cell_type <- opt$cell_type
in_fn <- paste0(opt$in_dir, '/', cell_level, '/', cell_type,'.rds')

## Covariates
covs <- opt$covs
covs_fn <- paste0(opt$in_dir,'/',covs)

## Read files
if(url.exists(opt$in_dir)){
  print('From URL...')
  print(paste0('Reading PBMC seurat object file in: ',in_fn))
  system.time(pbmc <- readRDS(url(in_fn)))
  print(paste0('Reading covariates file in: ',covs_fn))
  covs.df <- read.table(url(covs_fn), header = T)
}else if(dir.exists(opt$in_dir)){
  print('From local input directory...')
  print(paste0('Reading PBMC seurat object file in: ',in_fn))
  system.time(pbmc <- readRDS(in_fn))
  print(paste0('Reading covariates file in: ',covs_fn))
  covs.df <- read.table(covs_fn, header = T)
}else{
  err_in <- 'Neither URL nor local input directory exists. Please provide a valid URL or local input directory.'
  stop(err_in) 
}

# Output directory
out.dir <- paste0(opt$out_dir,'/')
print(paste0('Main output directory in: ',out.dir))
out.dir <- paste0(out.dir, cell_level, '/')
out.dir <- paste0(out.dir, cell_type, '/')
print(paste0('Creating output subdirectory in: ',out.dir))
if(!dir.exists(out.dir)){dir.create(out.dir, recursive = T)}

###################################### Pre-processing ###################################
print('##################### Pre-processing #####################')
sca <- preprocess.func(so = pbmc, freq_expressed=opt$freq)

#################### single-cell Differential Expression Analysis (scDEA) ################
# Subset the original seurat object with a specific set of genes to perform the scDEA only on those genes (optional)
if(!is.null(opt$genes)){
  # Reading the subset of genes
  genes_fn <- paste0(opt$in_dir,'/',opt$genes, '.', cell_level, '.tab')
  print(paste0('Testing with a subset of genes in: ',genes_fn))
  if(url.exists(genes_fn)){
    print('From URL...')
    genes <- read.delim(url(genes_fn), header = F)$V1
  }else if(file.exists(genes_fn)){
    print('From local input directory...')
    genes <- read.delim(genes_fn, header = F)$V1
  }else{
    err_in <- 'Neither URL nor local file exists. Please provide a valid URL or local input file.'
    stop(err_in) 
  }
  # Subsetting the SCA object
  sca <- sca[which(rownames(sca)%in%genes),]
}

# Saving the final SCA object
sca.fn <- paste0(out.dir,'pbmc_sca.rds') #save the filtered sca object
print(paste0('Saving sca object in: ',sca.fn))
saveRDS(sca, sca.fn)

# Running scDEA
## parameters
contrast_var <- unique(covs.df[covs.df$type=='contrast',]$covariate)
fixed_effects_var <- unique(covs.df[covs.df$type=='fixed',]$covariate)
random_effects_var <- unique(covs.df[covs.df$type=='random',]$covariate)

## apply function
print('##################### Running scDEA #####################')
de_glmer.res <- de_glmer.func(sca_object = sca,
                              contrast = contrast_var,
                              fixed_effects = fixed_effects_var,
                              random_effects = random_effects_var,
                              nagq0 = opt$nagq0)
# Saving scDEA results
if(opt$nagq0){
  de_glmer.fn <- paste0(out.dir, 'de_glmer.nagq0.rds')
}else{
  de_glmer.fn <- paste0(out.dir, 'de_glmer.rds')
}
print(paste0('Saving DEA results with MAST glmer in: ', de_glmer.fn))
saveRDS(de_glmer.res, de_glmer.fn)

######################## Extract scDEA summary statistics ########################
print('##################### Extracting scDEA summary statistics #####################')
# apply 'de_stats.func' function
de_stats.res <- de_stats.func(summaryDt = de_glmer.res,
                              contrast_in = contrast_var)

# Saving scDEA summary statistics
if(opt$nagq0){
  de_glmer_stats.fn <- paste0(out.dir, 'de_glmer_stats.nagq0.rds')
}else{
  de_glmer_stats.fn <- paste0(out.dir, 'de_glmer_stats.rds')
}
print(paste0('Saving DEA summary stats with MAST glmer in: ', de_glmer_stats.fn))
saveRDS(de_stats.res, de_glmer_stats.fn)