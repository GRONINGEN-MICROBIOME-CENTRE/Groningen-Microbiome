# =============================================================================
# By: Weersma Group, UMCG (2020)
#
# DMP Variance explained analysis, multivariate adonis
#
# script does multivariate ADONIS analysis of microbiome
# function (MetaCyc pathways) or composition (species-level)
# for multiple phenotypes from DMP data
# it is intended to be run once per group of phenotypes per data layer
#
# Note: must be invoked with 4 CL params:
# (1) = phenotypes
#  > must be in format "phenotype1,phenotype2,...,phenotypeN" (with quotations included), ex: "AGE,SEX,BMI"
# (2) = number of permutations
# (3) = what to use (taxa or pwys)
# (4) = type of analysis: (sum,seq or marg): 
#  > sum : summary : total variance explained for all phenotypes put together (without individual contribution)
#  > seq : sequential : variables are added one by one, as entered, result is total variance explained + per phenotype
#  > marg : marginal : variabes are all included into the model, result is total variance explained + per phenotype
#  > NOTE: type of analysis does not matter for total variance explained for group of phenotypes, 
#          but does affect results for individual phenotypes
#          in DMP, this script was primarily used for summary analysis of group of phenotype (example: all diseases or all medications)
# (5) = output file
# ============================================================================

# load libraries
library(vegan)
# load helper scripts (see https://github.com/GRONINGEN-MICROBIOME-CENTRE/Groningen-Microbiome/tree/master/Projects/DMP/r_scripts_library)
source('../myLibs_v3/R_Microbiome_scripts.R')

# maxFreqClass function 
# ================
# returns most common class
# in vector or matrix of class factor
# ex: maxFreqClass(c("A","A","B","C","C","C")) = "C" 
#     maxFreqClass(matrix(c("A","A","B","A","B","C"),2)) = "A"
maxFreqClass <- function(InVec, mult = FALSE) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  if (isTRUE(mult)) {
    levels(InVec)[A == max(A)]
  }
  else levels(InVec)[which.max(A)]
}

# SETTINGS
# ==========================================
outFolder <- "output_adonis_multivariate"
dataTaxa <- 'data/DAG3_metaphlan_bacteria_archaea_nofiltering.txt'
dataPWYs <- 'data/DAG3_humann_Metacyc_nofiltering.txt'
# ==========================================

# === COLLECT COMMAND LINE PARAMS ===
args = commandArgs(trailingOnly=TRUE)
# (1) = phenotypes
# (2) = number of permutations
# (3) = what to use (taxa or pwys)
# (4) = type of analysis
# (5) = output file
if (length(args) != 5)  {
   stop('CL params: <list of phenotypes, comma-separeted in one line [ex: "AGE,SEX,BMI"]> <nr of permutations> <PWYS or TAXA> <analysis type [sum,seq or marg]> <output file>')
}

# parse phenotypes
phenosAll <- as.character(args[1])
phenosList <- unlist(strsplit(phenosAll,split=','))
#debug:
# print(phenosList)

permNR <- as.numeric(args[2])
dataLayer <- as.character(args[3])
if (!(dataLayer %in% c("PWYS","TAXA"))) {
   stop(' input data type (CL 3) must be PWYS or TAXA')
}
analysisType <- as.character(args[4])
if (!(analysisType %in% c("sum","seq","marg"))) {
   stop(' input data type (CL 4) must be sum, seq or marg')
}
outFile <- as.character(args[5])

# ===================================

# load phenotypes
print(' >> LOADING PHENOTYPES')
inPhenos <- read.table('/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v26.csv',sep=',',header=T,quote='"')
rownames(inPhenos) <- inPhenos$DAG3_sampleID
# imputation & stuff (special cases)
# > diseases
inPhenos$MED.DISEASES.Cancer.Prostate[is.na(inPhenos$MED.DISEASES.Cancer.Prostate)] <- "N"
inPhenos$MED.DISEASES.Cancer.Breast[is.na(inPhenos$MED.DISEASES.Cancer.Breast)] <- "N"
inPhenos$MED.DISEASES.Cancer.Cervix[is.na(inPhenos$MED.DISEASES.Cancer.Cervix)] <- "N"
inPhenos$MED.DISEASES.Neurological.Headaches.Duration[is.na(inPhenos$MED.DISEASES.Neurological.Headaches.Duration)] <- "less than 6 months"
inPhenos$MED.DISEASES.Neurological.Headaches.HowOften[is.na(inPhenos$MED.DISEASES.Neurological.Headaches.HowOften)] <- "not at all"
# > early life
inPhenos$EXP.EARLYLIFE.Birth.lt[is.na(inPhenos$EXP.EARLYLIFE.Birth.lt)] <- median(inPhenos$EXP.EARLYLIFE.Birth.lt,na.rm=T)
inPhenos$EXP.EARLYLIFE.Birth.wt.g[is.na(inPhenos$EXP.EARLYLIFE.Birth.wt.g)] <- median(inPhenos$EXP.EARLYLIFE.Birth.wt.g,na.rm=T)
inPhenos$EXP.EARLYLIFE.Preg.Mother.Smoking[is.na(inPhenos$EXP.EARLYLIFE.Preg.Mother.Smoking)] <- "No"
inPhenos$EXP.EARLYLIFE.PretermBorn[is.na(inPhenos$EXP.EARLYLIFE.PretermBorn)] <- "N"
# > other
inPhenos$META.Antibiotics_3m[is.na(inPhenos$META.Antibiotics_3m)] <- "N"

# impute everything else
for (cn in colnames(inPhenos)) {
   if (cn != "DAG3_sampleID") {
      if (class(inPhenos[[cn]]) == "numeric" | class(inPhenos[[cn]])=="integer") {
         inPhenos[[cn]][is.na(inPhenos[[cn]])] <- median(inPhenos[[cn]],na.rm=T)
      } else if (class(inPhenos[[cn]]) == "factor") {
         inPhenos[[cn]][is.na(inPhenos[[cn]])] <- maxFreqClass(inPhenos[[cn]][!is.na(inPhenos[[cn]])])
      } else {
        print(paste0(cn,' = ',class(inPhenos[[cn]])))
      }
   }
}

# load data (microbiome pwys)
if (dataLayer == "TAXA") {
   print(' >> LOADING MICROBIOME')
   inDF <- read.table(dataTaxa,sep='\t',header=T)
   # grab species only
   inMB <- filterMetaGenomeDF(inDF,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = T,verbose = T,keepLevels = c("S"))
} else if (dataLayer=="PWYS") {
   print(' >> LOADING PATHWAYS')
   inMB <- read.table(dataPWYs,sep='\t',header=T)
}

# cleanup of folders
# ==========================================
if (!dir.exists(outFolder)) {
  dir.create(outFolder)
}

# =================================================================
# do adonis
# =================================================================
print (' >> STARTING ADONIS CALCULATIONS')
#print (paste0('   >> phenotypes: '))
#print (phenosList)
print (paste0('   >> NR of premut ', permNR))
print (' ============================================ ')
print(timestamp())
adonisResults <- NULL

inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(phenosList,"DAG3_sampleID")]
colnames(inMB)[colnames(inMB)=="ID"] <- "DAG3_sampleID"
allDF <- merge(x=inPhenosOneVarID,by.x="DAG3_sampleID",y=inMB,by.y="DAG3_sampleID")
rownames(allDF) <- allDF$DAG3_sampleID
allDF$DAG3_sampleID <- NULL

# === VARIABLE CHECK ===
# drop vars with < X cases
xCases <- 0
toDrop <- c()
for (cc in phenosList) {
   nrNA <- sum(is.na(allDF[[cc]]))
   nrTot <- length(allDF[[cc]])
   nrNonNA <- nrTot - nrNA
   if (nrNonNA < xCases) {
      print(paste0(' > WARNING: ',cc,' is ',nrNA,' NA; ',nrNonNA,' non-NA; ',nrTot,' Total samples, dropping it!'))
      toDrop <- c(toDrop,cc)
   }    
}
phenosList <- phenosList[!(phenosList %in% toDrop)]
allDF <- allDF[,!(colnames(allDF) %in% toDrop)]

# use complete cases
allDF <- allDF[complete.cases(allDF),]
#debug: make smaller for testing
#allDF <- allDF[1:1000,]

# =============================
av <- allDF[,colnames(allDF) %in% phenosList]
allDF <- allDF[,!(colnames(allDF) %in% phenosList)]

print ('  >> calculating B/C distance')
inBC <- vegdist(allDF,method = "bray")

print ('  >> doing adonis')
if (length(av) < 3 | length(unique(av)) < 2) {
  print(paste0(' >> WARNING: ',i,' has no useful data!!'))
} else {
    frm <- reformulate(termlabels=colnames(av),response='inBC')
    if (analysisType == 'seq') {
       print ('  >> running sequential adonis')
    ad <- adonis(formula=frm,data=av,permutations=permNR)
       print(ad)
    } else if (analysisType == 'sum') {
       print ('  >> running summary adonis')
       ad <- adonis2(formula=frm,data=av,permutations=permNR,by=NULL)
       print(ad)
    } else if (analysisType == 'marg') {
       print ('  >> running marginal adonis')
       ad <- adonis2(formula=frm,data=av,permutations=permNR,by="margin")
       print(ad)
    }
    print (paste0('--- DONE! ---'))
    saveRDS(ad,paste0(outFolder,"/",paste0("ad_",outFile)))
}
