# ==================================================
# By R.Gacesa, UMCG (2019)

#  Do analysis for flares
# ==================================================
library(caret)
library(vegan)

library(doSNOW)
cl <- makeCluster(12)
registerDoSNOW(cl)
#registerDoSEQ()

# load stuff
setwd("D:/UMCG/ML_IBD_v5")
source('../myLibs_v3/R_ML_scripts_v4.R')
source('../myLibs_v3/R_Microbiome_scripts.R')
source('../myLibs_v3/R_Misc.R')

# color blind palette
cbPalette <- c("#0072B2", "#D55E00", "#CC79A7")

# output folder
#outFolderName = "plots_final"

# CLINICAL PHENOS
# Age, Sex, BMI
# =====================================
allDataR <- read.table('data_step1_v5/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv',stringsAsFactors = F,sep=',',header = T)
allDataPhenos <- subsetMicrobiomeDF(inDF = allDataR,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,getDivs = F,verbose = T)
allDataPhenos <- allDataPhenos[allDataPhenos$Diagnosis %in% c("CD","UC","IBDU"),]
allDataPhenos$BMI[is.na(allDataPhenos$BMI)] <- median(allDataPhenos$BMI,na.rm = T)

# ===> Disease activity (Flares)
phenoDataRdy <- allDataPhenos
phenoDataRdy$Diagnosis <- "Inactive"
phenoDataRdy$Diagnosis[phenoDataRdy$IBD.flare.clinical == "in.flare"] <- "Active"
phenoDataRdy$IBD.active.clinical <- factor(phenoDataRdy$Diagnosis,levels=c("Inactive","Active"))

phenoDataRdy <- phenoDataRdy[phenoDataRdy$Diagnosis %in% c("Inactive","Active") & phenoDataRdy$Diagnosis.3 == "IBD",]
phenoDataRdy$Diagnosis <- as.factor(phenoDataRdy$Diagnosis)
phenoDataRdy$Gender <- as.factor(phenoDataRdy$Gender)
phenoDataRdyM <- phenoDataRdy[,colnames(phenoDataRdy) %in% c("Age","Gender","BMI","Diagnosis")]


# CALCULATE SUMMARY STATISTICS vs IBD activity
# =================================================
biomarkersSummary <- makeMultiSummaryByResponse(phenoDataRdy,
                                                ftrs = c("HBD2.nontransformed","ChrA.nontransformed","Calprot.nontransformed"),
                                                response = "Diagnosis",
                                                doTestVsControl = T,controlClass = "N",verbose = T,doSort = F,includeTotals = T)
#  > save it 
write.table(biomarkersSummary, 'final_plots_IBDactivity_clinical/flares_biomarkers_summary.csv', sep=',', row.names = F) 

# TEST BIOMARKERS (IBD activity)
# =====================================
t <- testOneFeature(dataIn = phenoDataRdy,saveFolder = 'final_plots_IBDactivity_clinical/flares_HBD2',feature = "HBD2.nontransformed",
                    doSave = T,display = "P", title = "Human beta defensin 2 (HBD2)",onlyShowSig = T,
                    yLab = "HBD2 (ng/g)",doViolin = F,doPlots = T,ylim = c(0,650),nrTests = 3,xLab = "IBD active",retPlot = T)
t[[2]] <- t[[2]] + scale_color_manual(values = cbPalette)
ggsave(plot = t[[2]],filename = 'final_plots_IBDactivity_clinical/flares_HBD2.png',width = 4,height = 6)
write.table(t,file = 'final_plots_IBDactivity_clinical/IBDactive_HBD2.csv',sep=',',row.names = F)

t <- testOneFeature(dataIn = phenoDataRdy,saveFolder = 'final_plots_IBDactivity_clinical/flares_CgA',feature = "ChrA.nontransformed",
                    doSave = T,display = "P", title = "Chromogranin A (CgA)",onlyShowSig = T,retPlot = T,
                    yLab = "CgA (nmol/g)",doViolin = F,doPlots = T,ylim = c(0,5.0),nrTests = 3,xLab = "IBD active")
t[[2]] <- t[[2]] + scale_color_manual(values = cbPalette)
ggsave(plot = t[[2]],filename = 'final_plots_IBDactivity_clinical/flares_CgA.png',width = 4,height = 6)
write.table(t,file = 'final_plots_IBDactivity_clinical/IBDactive_CgA.csv',sep=',',row.names = F)

t <- testOneFeature(dataIn = phenoDataRdy,saveFolder = 'final_plots_IBDactivity_clinical/flares_FCal',feature = "Calprot.nontransformed",
                    doSave = T,display = "P", title = "Fecal Calprotectin (FCal)",onlyShowSig = T,retPlot = T,
                    yLab = "FCal (ug/g)",doViolin = F,doPlots = T,ylim = c(0,1200),nrTests = 3,xLab = "IBD active")
t[[2]] <- t[[2]] + scale_color_manual(values = cbPalette)
ggsave(plot = t[[2]],filename = 'final_plots_IBDactivity_clinical/flares_FCal.png',width = 4,height = 6)
write.table(t,file = 'final_plots_IBDactivity_clinical/IBDactive_FCal.csv',sep=',',row.names = F)


# Build prediction models for Age+Sex+BMI (null model)
# =======================================================
doMLModelling(outFolderName = "models_v5_flares_clinical_null_ROC",allDM=list(phenoDataRdyM),allDMn = c("Phenos"),
              responseVar = "Diagnosis",posC = "Active",keepTrSplit = T,optimiseMetric="ROC",
              tButLC = 10, tButXV = 10, tRep = 3,mtds = c("svmRadial"),
              doDataPrep = F, doDataSep = T, doLC = T, doRFE = F,doRawMdls = T,doOptLC = F,doOptMdls = T,doDataPrepPreselection = F,
              smoothROCs = F,doRFEpreselection = F)

doMLModelling(outFolderName = "models_v5_flares_clinical_null_Kappa",allDM=list(phenoDataRdyM),allDMn = c("Phenos"),
              responseVar = "Diagnosis",posC = "Active",keepTrSplit = T,optimiseMetric="Kappa",
              tButLC = 10, tButXV = 10, tRep = 3,mtds = c("svmRadial"),
              doDataPrep = F, doDataSep = T, doLC = T, doRFE = F,doRawMdls = T,doOptLC = F,doOptMdls = T,doDataPrepPreselection = F,
              smoothROCs = F,doRFEpreselection = F)

# LOAD NORMALIZED & CORRECTED DATA
# ==============================================================
allDataR <- read.table('data_step1_v5/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv',stringsAsFactors = F,sep=',',header = T)
# allDataPhenos <- subsetMicrobiomeDF(inDF = allDataR,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,getDivs = F,verbose = T)
# allDataPhenos <- allDataPhenos[allDataPhenos$Diagnosis %in% c("CD","UC","IBDU"),]
# allDataPhenos$BMI[is.na(allDataPhenos$BMI)] <- median(allDataPhenos$BMI,na.rm = T)

# ==============================================================================
# ==============================================================================
# STEP 3: IBD activity modeling
# ==============================================================================
allDataRP_MB <- allDataR
allDataRP_MB <- allDataRP_MB[!is.na(allDataRP_MB$Diagnosis.3),]
allDataRP_MB <- allDataRP_MB[allDataRP_MB$Diagnosis.3 == "IBD",]
allDataRP_MB$Diagnosis <- "N"
allDataRP_MB$Diagnosis[allDataRP_MB$IBD.flare.clinical == "in.flare"] <- "Y"
allDataRP_MB <- allDataRP_MB[allDataRP_MB$Diagnosis %in% c("Y","N"),]

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
# doMLModelling(outFolderName = "tests_IBD_active",allDM=allDM,allDMn = allDMn,
#               responseVar = "Diagnosis",posC = "Y",keepTrSplit = T,
#               tBut = 2, tRep = 2,mtds = c("glm","svmRadial"),
#               doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T)
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

allDM <- list(mdlBmFc,mdlBmHBD2,mdlBmChrA,mdlBmFcHBD2,mdlBmFcChrA,mdlBmChrAHBD2,mdlBmAll,
              mTaxG_N,     mTaxS_N,     mPWY_N,
              mTaxG_BM_N,  mTaxS_BM_N,  mPWY_BM_N)
allDMn <- list("Calprot", "HBD2",   "ChrA",   "Calprot_HBD2","Calprot_ChrA","ChrA_HBD2",  "Calprot_HBD2_ChrA",
               "tGenera","tSpecies","pPWY",
               "tGenera_BM", "tSpecies_BM", "pPWY_BM")


#allDM <- list(mTaxG_N)
#allDMn <- list("tGenera")

# # > build models (ROC)
doMLModelling(outFolderName = "final_models_MB_IBD_activity_v5_clinical_ROC",allDM=allDM,allDMn = allDMn,
              responseVar = "Diagnosis",posC = "Y",keepTrSplit = T,optimiseMetric = "Kappa",
              tButLC = 3, tButXV = 10, tRep = 3,mtds = c("glm"),rfeR = 10,
              doDataPrep = T, doDataSep = T, doLC = T, doRFE = T,doRawMdls = T,doOptLC = T,doOptMdls = T,doDataPrepPreselection = T,
              smoothROCs = F,doRFEpreselection = T)

# # # > build models (ROC)
# doMLModelling(outFolderName = "final_models_MB_IBD_activity_v5_clinical_roc",allDM=allDM,allDMn = allDMn,
#               responseVar = "Diagnosis",posC = "Y",keepTrSplit = T,optimiseMetric = "ROC",
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