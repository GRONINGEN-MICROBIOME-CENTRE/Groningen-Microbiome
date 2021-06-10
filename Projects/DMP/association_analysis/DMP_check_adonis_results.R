# ======================================
# By: Weersma Group, UMCG (2020)
#
# DMP Variance explained analysis,
#
# helper script for Adonis analysis,
# checks if adonis results are accounted for
# ===================================

# === COLLECT COMMAND LINE PARAMS ===
args = commandArgs(trailingOnly=TRUE)
# (1) = adonis results
adonisResults <- args[1]
# ===================================

print (adonisResults)
# load adonis results
inRes <- read.table(adonisResults,header=T,sep=',')
# load phenotypes
print(' >> LOADING PHENOTYPES')
inPhenos <- read.table('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v26.csv',sep=',',header=T,quote='"')
inPhenos$ID <- NULL
inPhenos$DAG3_sampleID <- NULL
# iterate over them and check which phenotypes are not accounted for
cn = 0
for (cc in colnames(inPhenos)) {
   cn = cn + 1
   if (!cc %in% inRes$Var) {
      print(paste0('pheno ',cn,' = ',cc,' missing!'))
   }
}
