#!/usr/bin/env Rscript

# options parser
shhh <- suppressPackageStartupMessages
shhh(library(optparse))
option_list = list(
  make_option(c("--in_dir"), action="store", default=NA, type='character',
              help="Output directory from dea_MAST_glmer.R"),
  make_option(c("--cell_level"), action="store", default=NA, type='character',
              help="Cell level of interest. Example: approach_I or approach_II."),
  make_option(c("-m","--nagq0"), action="store", default=TRUE, type='logical',
              help="Mitigate model fails to converge by passing nAGQ=0 to the fitting function zlm (fitArgsD=list(nAGQ=0))."))
opt = parse_args(OptionParser(option_list=option_list))

################################ Set Variables #################################
# Main directory
cwd <- getwd()
setwd(cwd)

# Loading functions
functions_fn <- 'scripts/accessory_functions.R'
print(paste0('Loading functions from: ',functions_fn))
source(functions_fn)

# Input directory
cell_level <- opt$cell_level
in.dir <- paste0(opt$in_dir, '/', cell_level, '/')

############################### Reformat files #################################
# Read scDEA summary statistics across cell types for one specific approach into a list
print('################ Reformatting scDEA summary statistics #####################')
reformat.list <- reformat.func(in_dir = in.dir,
                               nagq0 = opt$nagq0)

# Saving outputs
## reformatted list (input for Downstream_SC_DE.R)
reformat_formatted <- lapply(reformat.list, function(x) x$formatted)
fn <- paste0(in.dir, 'dea.list_formatted')
fn <- ifelse(opt$nagq0, paste0(fn,'.nagq0.rds'), paste0(fn,'.rds'))
print(paste0('Saving reformatted scDEA summary statistics in: ', fn))
saveRDS(reformat_formatted, fn)

## reformatted table in Table_S7.1 (approach I) and Table_S7.2 (approach II)
reformat_raw.list <- lapply(reformat.list, function(x) x$raw)
reformat_raw <- do.call("rbind",reformat_raw.list)
rownames(reformat_raw) <- NULL
fn <- paste0(in.dir, 'Table_S7.')
suff <- ifelse(cell_level=='approach_I', '1', '2')
fn <- paste0(fn, suff, '.csv')
print(paste0('Saving reformatted scDEA supplementary table in: ', fn))
write.csv(reformat_raw, fn, row.names = F, quote = F)
