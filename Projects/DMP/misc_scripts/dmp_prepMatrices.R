# By: R.Gacesa, UMCG (2020)
#
# DMP raw data processor
# =================================================================
#
# calculates Bray-Curtis dissimilarity matrices and
# performs principal coordinate analysis
#

# ==== load libraries ====
library(vegan)

# =================================================================
# =================================================================
# MAIN 
# =================================================================
# =================================================================
setwd('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_statistics')
source('codes/myLibs_v3/R_Microbiome_scripts.R')
inFolderMB <- '/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome'
inMicrobiome <- '/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/raw/DAG3_metaphlan_merged.txt'
inPhenos <- '/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v27.csv'

# ===================
# LOAD PHENOS
# note: we merge them with microbiome and use only matched samples
# ===================
inPhDf <- read.table(inPhenos,header=T,sep=',')

# ======================================================================================================================================
# > CALCULATE B/C matrices
# ======================================================================================================================================
# ======================================================================================================================================
DIVlvls <- c("S","G","F","O","C","P")
print('============================================================')
print(' >>>>         PREPARING B/C MATRICES                   <<<< ')
print('============================================================')

# filter raw microbiome to obtain "good" samples
print('   > loading microbiome ... ')
inMB <- read.table(inMicrobiome,header=T,sep='\t')
rownames(inMB) <- inMB$ID
inMB$ID <- NULL
inMB <- as.data.frame(t(inMB))
inMB$ID <- rownames(inMB)

inMBs <- subsetMicrobiomeDF(inMB,getPWYs=F,getTaxa=T,getPhenos=F,getVFs=F,getCARDs=F,getDivs=F,getID=T,idToRowName=T)
rownames(inMBs) <- gsub('_metaphlan','',rownames(inMBs))
inMBss <- inMBs[rownames(inMBs) %in% inPhDf$DAG3_sampleID,]
inMBss <- inMBss[inMBss$k__Bacteria > 80,]
inMBss$ID <- rownames(inMBss)
write.table(inMBss,'/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/microbiome_8208samples_rdy.csv',sep=',',row.names=F)
inMBss$ID <- NULL

# >> load filtered microbiome, make bray-curtis matrix per taxonomy level
for (div in DIVlvls[1:length(DIVlvls)]) {
  print('    >> subsetting taxa (no filtering')
  inMBsss <- filterMetaGenomeDF(inDF=inMBss,presPerc=-1,minMRelAb=-1,minMedRelAb=-1,keepDomains=c("Bacteria","Archaea"),keepLevels=div)
  inMBsss$ID <- rownames(inMBsss)
  write.table(inMBsss,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/microbiome_8208samples_',div,'_unfiltered.csv',sep=',',row.names=F))
  inMBsss$ID <- NULL
  print('     > calculating B/C distance matrix')  
  dv <- vegdist(inMBsss,method = "bray")
  print('     > saving matrix')
  saveRDS(dv,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/BC_matrices/BC_matrix_unfiltered_',div,'.RData'))
  print('    >> subsetting taxa (filtered)')
  inMBsss <- filterMetaGenomeDF(inDF=inMBss,presPerc=0.05,minMRelAb=0.00001,minMedRelAb=-1,keepDomains=c("Bacteria","Archaea"),keepLevels=div)
  inMBsss$ID <- rownames(inMBsss)
  write.table(inMBsss,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/microbiome_8208samples_',div,'_filtered_prab.csv',sep=',',row.names=F))
  inMBsss$ID <- NULL
  print('     > calculating B/C distance matrix')  
  dv <- vegdist(inMBsss,method = "bray")
  print('     > saving matrix')
  saveRDS(dv,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/BC_matrices/BC_matrix_filtered_prab_',div,'.RData'))
  print('    > done <')
}

# ======================================================================================================================================
# > CALCULATE PCoA
# ======================================================================================================================================
# ======================================================================================================================================
# >>> UNFILTERED <<<<
# iterate over stuff to calculate (taxa per taxon lvl, pathways, VFs, CARDs)
DIVlvls <- c("S","G","F","O","C","P")
print('============================================================')
print(' >>>>         PREPARING PCoA MATRICES      <<<< ')
print('============================================================')
for (div in DIVlvls[1:length(DIVlvls)]) {
  print(paste0('  >> UNFILTERED DATA, ',div))
  # I) load data
  print('   > loading data ... ')
  distMatrix <- readRDS(paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/BC_matrices/BC_matrix_unfiltered_',div,'.RData'))
  distMatrix[is.na(distMatrix)] <- 1.0
  # II) prep ordination (PcoA)
  print('   > calculating PCoA ... ')
  pCoa <- cmdscale( distMatrix, eig = T,k = 20 )
  # calculate variance explained
  varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:20]
  xVar <- as.integer(varExp[1]*100)
  yVar <- as.integer(varExp[2]*100)
  # get axes
  pCoaVecs <- as.data.frame(pCoa$points)
  pCoaVecsDF <- cbind.data.frame(data.frame(PCoA=c(1:20)),data.frame(VarExplained=varExp))
  colnames(pCoaVecs) <- paste0("PCo",c(1:20))
  pCoaVecs$ID <- row.names(pCoaVecs)
  print('   > saving PCoA results ... ')
  write.table(pCoaVecs,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/PCoA_matrices/PCoA_unfiltered_',div,'_PCs.csv'),sep = ',',row.names = F)
  write.table(pCoaVecsDF,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/PCoA_matrices/PCoA_unfiltered_',div,'_PCs_varExp.csv'),sep = ',',row.names = F) 
}

# >>> FILTERED (5% prevalence, 0.01% abundance <<<<
DIVlvls <- c("S","G","F","O","C","P")
for (div in DIVlvls[1:length(DIVlvls)]) {
  print(paste0('  >> UNFILTERED DATA, ',div))
  # I) load data
  print('   > loading data ... ')
  distMatrix <- readRDS(paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/BC_matrices/BC_matrix_filtered_prab_',div,'.RData'))
  distMatrix[is.na(distMatrix)] <- 1.0
  # II) prep ordination (PcoA)
  print('   > calculating PCoA ... ')
  pCoa <- cmdscale( distMatrix, eig = T,k = 20 )
  # calculate variance explained
  varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:20]
  xVar <- as.integer(varExp[1]*100)
  yVar <- as.integer(varExp[2]*100)
  # get axes
  pCoaVecs <- as.data.frame(pCoa$points)
  pCoaVecsDF <- cbind.data.frame(data.frame(PCoA=c(1:20)),data.frame(VarExplained=varExp))
  colnames(pCoaVecs) <- paste0("PCo",c(1:20))
  pCoaVecs$ID <- row.names(pCoaVecs)
  print('   > saving PCoA results ... ')
  write.table(pCoaVecs,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/PCoA_matrices/PCoA_filtered_prab_',div,'_PCs.csv'),sep = ',',row.names = F)
  write.table(pCoaVecsDF,paste0('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/PCoA_matrices/PCoA_filtered_prab_',div,'_PCs_varExp.csv'),sep = ',',row.names = F) 
}
print ('========================================')
print ('>>>>>>>>>>>>>> PCOA DONE <<<<<<<<<<<<<<<')
print ('========================================')

