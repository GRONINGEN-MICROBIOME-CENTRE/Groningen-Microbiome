# ==================================================
# By R.Gacesa, UMCG (2020)

#  Do analysis for IBD vs IBS prediction
# ==================================================
# > load R packages
library(caret)
library(vegan)
library(doSNOW)
# > prep multi-threading
cl <- makeCluster(12)
registerDoSNOW(cl)
#registerDoSEQ()

# > load my functions for Machine learning and microbiome 
setwd("D:/UMCG/ML_IBD_v5")
source('../myLibs_v3/R_ML_scripts_v4.R')
source('../myLibs_v3/R_Microbiome_scripts.R')
source('../myLibs_v3/R_Misc.R')

# color blind palette
#cbPalette <- c("#0072B2", "#D55E00", "#CC79A7")

# output folder
#outFolderName = "plots_final"

# Construction of models based on clinical phenotypes (Null models)
# =====================================
# > load data and prepare models for Age, Sex and BMI
allDataR <- read.table('data_step1_v5/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv',stringsAsFactors = F,sep=',',header = T)
allDataPhenos <- subsetMicrobiomeDF(inDF = allDataR,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,getDivs = F,verbose = T)
# > data filtering for diagnosis (IBD vs IBS)
allDataPhenos <- allDataPhenos[allDataPhenos$Diagnosis.3 %in% c("IBD","IBS"), ]
# > impute missing BMI to the mean (this is only few cases)
allDataPhenos$BMI[is.na(allDataPhenos$BMI)] <- median(allDataPhenos$BMI,na.rm = T)
# > fix variable types
phenoDataRdy <- allDataPhenos
phenoDataRdy$Diagnosis <- as.factor(as.character(phenoDataRdy$Diagnosis.3))
phenoDataRdy$Gender <- as.factor(as.character(phenoDataRdy$Gender))
# > prepare data for modelling
phenoDataRdyM <- phenoDataRdy[,colnames(phenoDataRdy) %in% c("Age","Gender","BMI","Diagnosis")]

# Build prediction models for Age+Sex+BMI (null model)
# =======================================================
# > ROC optimization
doMLModelling(outFolderName = "final_models_v5_MB_IBDvsIBS_null_roc",allDM=list(phenoDataRdyM),allDMn = c("Phenos"),
              responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric="ROC",
              tButLC = 10, tButXV = 10, tRep = 5,mtds = c("glm","svmRadial"),
              doDataPrep = F, doDataSep = T, doLC = T, doRFE = F,doRawMdls = T,doOptLC = F,doOptMdls = T,doDataPrepPreselection = F,
              smoothROCs = F,doRFEpreselection = F)

# > Kappa optimization
doMLModelling(outFolderName = "final_models_v5_MB_IBDvsIBS_null_kappa",allDM=list(phenoDataRdyM),allDMn = c("Phenos"),
              responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric="Kappa",
              tButLC = 10, tButXV = 10, tRep = 5,mtds = c("glm","svmRadial"),
              doDataPrep = F, doDataSep = T, doLC = T, doRFE = F,doRawMdls = T,doOptLC = F,doOptMdls = T,doDataPrepPreselection = F,
              smoothROCs = F,doRFEpreselection = F)

# LOAD NORMALIZED & CORRECTED DATA
# ==============================================================
allDataR <- read.table('data_step1_v4/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv',stringsAsFactors = F,sep=',',header = T)

# ==============================================================================
# ==============================================================================
# STEP 3: IBD modeling
# ==============================================================================
allDataRP_MB <- allDataR
allDataRP_MB <- allDataRP_MB[allDataRP_MB$Diagnosis.3 %in% c("IBD","IBS"), ]
allDataRP_MB$Diagnosis <- as.factor(as.character(allDataRP_MB$Diagnosis.3))

# make data models (BIOMARKERS)
mdlBmFc <- allDataRP_MB[,c("Diagnosis","Calprot")]
mdlBmHBD2 <- allDataRP_MB[,c("Diagnosis","HBD2")]
mdlBmChrA  <- allDataRP_MB[,c("Diagnosis","ChrA")]
mdlBmFcHBD2 <- allDataRP_MB[,c("Diagnosis","Calprot","HBD2")]
mdlBmFcChrA <- allDataRP_MB[,c("Diagnosis","Calprot","ChrA")]
mdlBmChrAHBD2 <- allDataRP_MB[,c("Diagnosis","ChrA","HBD2")]
mdlBmAll <- allDataRP_MB[,c("Diagnosis","ChrA","HBD2","Calprot")]

# make data models (MICROBIOME)
mMG_PWY <- allDataRP_MB[,c(grep('^Diagnosis$',colnames(allDataRP_MB)),grep('__',colnames(allDataRP_MB)),grep('PWY',colnames(allDataRP_MB)))]

# > make data models
# == taxa: phyla ==
mTaxP_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
#    + biomarkers
mTaxP_PWY_BM_N <- cbind.data.frame(mTaxP_PWY_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])
#    phyla only
mTaxP_N <- mTaxP_PWY_N[,-grep("PWY",colnames(mTaxP_PWY_N))]
#    phyla + BM
mTaxP_BM_N <- mTaxP_PWY_BM_N[,-grep("PWY",colnames(mTaxP_PWY_BM_N))]

# == taxa: genera ==
#    pathways & genera
mTaxG_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
#  + Cal
#mTaxG_PWY_Cal_N <- cbind.data.frame(mTaxG_PWY_N,Calprot=allDataRP_MB$Calprot)
#  + Biomarker
mTaxG_PWY_BM_N <- cbind.data.frame(mTaxG_PWY_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])
#    genera only
mTaxG_N <- mTaxG_PWY_N[,-grep("PWY",colnames(mTaxG_PWY_N))]
#  + Cal
#mTaxG_Cal_N <- cbind.data.frame(mTaxG_N,Calprot=allDataRP_MB$Calprot)
#  + Biomarker
mTaxG_BM_N <- cbind.data.frame(mTaxG_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])

# == taxa: species ==
#    pathways & species
mTaxS_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
#  + Cal
#mTaxS_PWY_Cal_N <- cbind.data.frame(mTaxS_PWY_N,Calprot=allDataRP_MB$Calprot)
#  + Biomarker
mTaxS_PWY_BM_N <- cbind.data.frame(mTaxS_PWY_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])
#    species only
mTaxS_N <- mTaxS_PWY_N[,-grep("PWY",colnames(mTaxS_PWY_N))]
#  + Cal
#mTaxS_Cal_N <- cbind.data.frame(mTaxS_N,Calprot=allDataRP_MB$Calprot)
#  + Biomarker
mTaxS_BM_N <- cbind.data.frame(mTaxS_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])

# === pathways ===
# pathways only
mPWY_N <- mTaxS_PWY_N[,c(grep('Diagnosis',colnames(mTaxS_PWY_N)),grep("PWY",colnames(mTaxS_PWY_N)))]
#  + Cal
#mPWY_Cal_N <- cbind.data.frame(mPWY_N,Calprot=allDataRP_MB$Calprot)
#  + Biomarker
mPWY_BM_N <- cbind.data.frame(mPWY_N,allDataRP_MB[,c("Calprot","HBD2","ChrA")])

## === test ===
# allDM <- list(mdlBmFc, mTaxG_N)
# allDMn <- list("Calprot", "tGenera")
# doMLModelling(outFolderName = "test_models_MB_IBDvsIBS_v5_clinical_roc",allDM=allDM,allDMn = allDMn,
#               responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric = "ROC",
#               tButLC = 10, tButXV = 5, tRep = 3,mtds = c("glm","svmRadial"),rfeR = 3,
#               doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T,doDataPrepPreselection = T,
#               smoothROCs = F,doRFEpreselection = T)
## END TEST 

# # > bundle data for modelling
# allDM <- list(mdlBmFc,mdlBmHBD2,mdlBmChrA,mdlBmFcHBD2,mdlBmFcChrA,mdlBmChrAHBD2,mdlBmAll,
#               mTaxP_N,     mTaxG_N,     mTaxS_N,     mPWY_N,     mMG_PWY_N,
#               mTaxP_Cal_N, mTaxG_Cal_N, mTaxS_Cal_N, mPWY_Cal_N, mMG_PWY_Cal_N,
#               mTaxP_BM_N,  mTaxG_BM_N,  mTaxS_BM_N,  mPWY_BM_N,  mMG_PWY_BM_N)
# allDMn <- list("Calprot", "HBD2",   "ChrA",   "Calprot_HBD2","Calprot_ChrA","ChrA_HBD2",  "Calprot_HBD2_ChrA",
#                "tPhyla","tGenera","tSpecies","pPWY", "PWY_Tax",
#                "tPhyla_Cal","tGenera_Cal","tSpecies_Cal","pPWY_Cal", "PWY_Tax_Cal",
#                "tPhyla_BM", "tGenera_BM", "tSpecies_BM", "pPWY_BM",  "PWY_Tax_BM")
# 

allDM <- list(mTaxP_N,mTaxP_BM_N)
allDMn <- list("Phyla","Phyla_BM")


#allDM <- list(mTaxG_N)
#allDMn <- list("tGenera")

# # > build models (GLM)
# doMLModelling(outFolderName = "final_models_MB_IBD_activity_v5_clinical_kappa",allDM=allDM,allDMn = allDMn,
#               responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric = "Kappa",
#               tButLC = 10, tButXV = 10, tRep = 5,mtds = c("glm","svmRadial"),rfeR = 10,
#               doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T,doDataPrepPreselection = T,
#               smoothROCs = F,doRFEpreselection = T)

doMLModelling(outFolderName = "final_models_MB_IBDvsIBS_kappa_extra",allDM=allDM,allDMn = allDMn,
              responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric = "Kappa",
              tButLC = 10, tButXV = 10, tRep = 5,mtds = c("glm","svmRadial"),rfeR = 10,
              doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T,doDataPrepPreselection = T,
              smoothROCs = F,doRFEpreselection = T)

# doMLModelling(outFolderName = "final_models_MB_IBDvsIBS_kappa_extra",allDM=allDM,allDMn = allDMn,
#               responseVar = "Diagnosis",posC = "IBD",keepTrSplit = T,optimiseMetric = "Kappa",
#               tButLC = 10, tButXV = 10, tRep = 5,mtds = c("glm","svmRadial"),rfeR = 10,
#               doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T,doDataPrepPreselection = T,
#               smoothROCs = F,doRFEpreselection = T)


# ==============================================================================
# ==============================================================================
# STEP 4: do comparison of models
# ==============================================================================
# ==============================================================================

# > load models
fitPhyla <- readRDS()
fitGenera <- readRDS()
fitSpecies <- readRDS()
fitPWY <- readRDS()
fitSpecPWY <- readRDS()

# > load testset
tSet <- read.table()

# > compare 

# COMPARE BIOMARKER MODELS
# ===========================================
# > load models (IBD vs IBS)
fitCalprot <- readRDS('')
fitChrA <- readRDS()
fitHBD2 <- readRDS()
fitCalprotChrA <- readRDS()
fitCalprotHBD2 <- readRDS()
fitCalprotChrAHBD2 <- readRDS()
# > load testset
tSet <- read.table()
# > compare 