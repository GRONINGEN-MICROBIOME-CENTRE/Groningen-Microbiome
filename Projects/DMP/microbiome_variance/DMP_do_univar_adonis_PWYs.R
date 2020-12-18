# =============================================================================
# By: Weersma Group, UMCG (2020)
#
# DMP Variance explained analysis, Species-level adonis
#
# script does univariate ADONIS analysis of microbiome function (MetaCyc pathways) 
# for one phenotype for DMP data
# it is intended to be run once per phenotype in parallel
#
# Note: must be invoked with 2 CL params:
# (1) = number of phenotype (column in phenotype file)
# (2) = number of permutations
# ============================================================================

# load libraries
library(vegan)
library(parallel)

# SETTINGS
# ==========================================
outFolder <- "output_adonis_PWYS"
# ==========================================

# ==========================================
# >> MAIN <<
# ==========================================

# === COLLECT COMMAND LINE PARAMS ===
args = commandArgs(trailingOnly=TRUE)
# (1) = number of phenotype
# (2) = number of permutations
phenoNR <- as.numeric(args[1])
permNR <- as.numeric(args[2])
# ===================================

# set WD
setwd('adonis_pwys')

# load helper scripts (see https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/tree/master/Projects/DMP/r_scripts_library)
source('../myLibs_v3/R_Microbiome_scripts.R')

# load phenotypes
print(' >> LOADING PHENOTYPES')
inPhenos <- read.table('../phenotypes/DAG3_metadata_merged_ready.csv',sep=',',header=T,quote='"')
rownames(inPhenos) <- inPhenos$DAG3_sampleID

# load data (microbiome pwys)
print(' >> LOADING MICROBIOME')
inMB <- read.table('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome/processed/DAG3_humann_Metacyc_nofiltering.txt',sep='\t',header=T)

# cleanup of folders
# ==========================================
if (!dir.exists(outFolder)) {
  dir.create(outFolder)
}

# load phenotypes
# ==========================================
# > remove IDs and other (irrelevant) phenotypes
adonisVarsTouse <- colnames(inPhenos)
toExclude = c("DAG3_sampleID","ID")
adonisVarsTouse <- adonisVarsTouse[!adonisVarsTouse %in% toExclude]

# fix names of columns
# ================================================================
if ( ("ID" %in% colnames(inPhenos)) & !("DAG3_sampleID" %in% colnames(inPhenos)) ) {
    colnames(inPhenos)[colnames(inPhenos)=="ID"] <- "DAG3_sampleID"
}
if ( ("ID" %in% colnames(inMB)) & !("DAG3_sampleID" %in% colnames(inMB)) ) {
    colnames(inMB)[colnames(inMB)=="ID"] <- "DAG3_sampleID"
}
if (!("DAG3_sampleID" %in% colnames(inMB))) {
   stop('ERROR: DAG3_sampleID not in columns of microbiome data')
}
if (!("DAG3_sampleID" %in% colnames(inPhenos))) {
   stop('ERROR: DAG3_sampleID not in columns of phenotype data')
}

# =================================================================
# do adonis
# =================================================================
print (' >> STARTING ADONIS CALCULATIONS')
print (paste0('   >> phenotype NR ', phenoNR))
print (paste0('   >> NR of premut ', permNR))
print (' ============================================ ')
print(timestamp())
adonisResults <- NULL
for (i in adonisVarsTouse[phenoNR]) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  print(timestamp())
  print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"DAG3_sampleID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="DAG3_sampleID",y=inMB,by.y="DAG3_sampleID")
  rownames(allDF) <- allDF$DAG3_sampleID
  allDF$DAG3_sampleID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  print ('  >> calculating B/C distance')
  inBC <- vegdist(allDF,method = "bray",parallel=4)
  print(timestamp())
  print ('  >> doing adonis')
  nrRows <- length(av)
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    #print(paste0(' NR NAs: ',sum(is.na(av))))
    ad <- adonis(inBC ~ av,permutations=permNR)
    aov_table <- ad$aov.tab
    # accumulate results
    oneRow <- data.frame(Var=i,
                         NR_nonNA=nrRows,
                         DF=aov_table[1,1],
                         SumsOfSqs=aov_table[1,2],
                         MeanSqs=aov_table[1,3],
                         FModel=aov_table[1,4],
                         R2=aov_table[1,5],
                         pval=aov_table[1,6],
                         FDR.BH=NA,
                         Significant=NA)
    print(oneRow)
    adonisResults <- rbind.data.frame(adonisResults,oneRow)
    write.table(adonisResults,paste0(outFolder,"/adonis_results_table_",phenoNR,".csv"),sep=",",row.names=F)
    print (paste0('--- ',i,' DONE! ---'))
  }
}
print (' >> DONE WITH CALCULATIONS')
#print(timestamp())
print ('  >> saving output')
rownames(adonisResults) = adonisResults$Var
adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant="No"
adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"
adonisResults <- adonisResults[order(adonisResults$pval),]

write.table(adonisResults,paste0(outFolder,"/adonis_PWYs_results_table_",phenoNR,".csv"),sep=",",row.names=F)
