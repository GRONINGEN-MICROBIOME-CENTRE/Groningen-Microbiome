# =================================================
# By R.Gacesa, UMCG (2019)
# =================================================
# load libs
library(tidyr)
library(vegan)
#library(Rtsne)
#library(VennDiagram)
setwd("UMCG/ML_IBD_v4/")
source('../myLibs_v3/R_ML_scripts_v4.R')
source('../myLibs_v3/R_Microbiome_scripts.R')
#library(mice)

# ========================================================================================================
# ------------------ Exploratory data analysis and data preparation for step 2 ------------
# =========================================================================================================

# OPTIONS
# ===============
doAbundancies = T
outFolderName = "IBD_ML_exp_LLDMIBS_v5"
outFolderData = "data_step1_v5"

if (dir.exists(paste0(outFolderName,'/d3_microbiome'))) {unlink(paste0(outFolderName,'/d3_microbiome'),recursive = T,force = T)}
if (dir.exists(paste0(outFolderName,'/d3_pathways')))   {unlink(paste0(outFolderName,'/d3_pathways'),recursive = T,force = T)}
if (dir.exists(paste0(outFolderName,'/d4_microbiome'))) {unlink(paste0(outFolderName,'/d4_microbiome'),recursive = T,force = T)}
if (dir.exists(paste0(outFolderName,'/d4_pathways')))   {unlink(paste0(outFolderName,'/d4_pathways'),recursive = T,force = T)}

if (!dir.exists(paste0(outFolderName,'/d3_microbiome'))) {dir.create(paste0(outFolderName,'/d3_microbiome'))}
if (!dir.exists(paste0(outFolderName,'/d3_pathways')))   {dir.create(paste0(outFolderName,'/d3_pathways'))}
if (!dir.exists(paste0(outFolderName,'/d4_microbiome'))) {dir.create(paste0(outFolderName,'/d4_microbiome'))}
if (!dir.exists(paste0(outFolderName,'/d4_pathways')))   {dir.create(paste0(outFolderName,'/d4_pathways'))}

allDataPrep <- read.table('data_step1_v5/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv',sep=',',header=T)

# ================================
# prep models [IBS vs IBD vs HC]
# ================================
allDataIBSIBDHC <- allDataPrep
allDataIBSIBDHC$Diagnosis <- allDataIBSIBDHC$Diagnosis.3

mMG_PWY <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('k__',colnames(allDataIBSIBDHC)),grep('PWY',colnames(allDataIBSIBDHC)) )]
# == biomarkers ==
# calprotectin only
mBM_Cal  <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('Calprot',colnames(allDataIBSIBDHC)))]
# HBD2 only
mBM_HBD2 <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('HBD2',colnames(allDataIBSIBDHC)))]
# ChrA only
mBM_ChrA <-  allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('ChrA',colnames(allDataIBSIBDHC)))]
# 2 out of 3
mBM_CalHBD2 <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('Calprot',colnames(allDataIBSIBDHC)),grep('HBD2',colnames(allDataIBSIBDHC)))]
mBM_CalChrA <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('Calprot',colnames(allDataIBSIBDHC)),grep('ChrA',colnames(allDataIBSIBDHC)))]
mBM_HBD2ChrA <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('HBD2',colnames(allDataIBSIBDHC)),grep('ChrA',colnames(allDataIBSIBDHC)))]
# 3 biomarkers
mBM_CalHBD2ChrA <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('Calprot',colnames(allDataIBSIBDHC)),
                                 grep('HBD2',colnames(allDataIBSIBDHC)),grep('ChrA',colnames(allDataIBSIBDHC)))]

# == taxa: phyla ==
#    pathways & phyla
mTaxP_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#  phyla only
mTaxP_N <- mTaxP_PWY_N[,-grep("PWY",colnames(mTaxP_PWY_N))]

# == taxa: genera ==
#    pathways & genera
mTaxG_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    genera only
mTaxG_N <- mTaxG_PWY_N[,-grep("PWY",colnames(mTaxG_PWY_N))]

# == taxa: species ==
#    pathways & species
mTaxS_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    species only
mTaxS_N <- mTaxS_PWY_N[,-grep("PWY",colnames(mTaxS_PWY_N))]

# === pathways ===
# pathways only
mPWY_N <- mTaxS_PWY_N[,c(grep('Diagnosis',colnames(mTaxS_PWY_N)),grep("PWY",colnames(mTaxS_PWY_N)))]


# compare abundances of pathways:
#=====================================
mPWY_N <- mPWY_N[!is.na(mPWY_N$Diagnosis),]
tests <- testAbundances(dataIn = mPWY_N,saveFolder = paste(outFolderName,'/d3_pathways/pwy',sep=""),pathway = T,
                        onlyShowSig = T,doPlots = doAbundancies)

dataIn = mPWY_N;saveFolder = paste(outFolderName,'/d3_pathways/pwy',sep="");pathway = T;onlyShowSig = T;doPlots = doAbundancies

write.table(tests,paste0(outFolderName,'/d3_pathways/pwy_statistics.csv'),sep=',',row.names = F)

# compare abundances of PHYLA
# ==============================================================
mTaxP_N <- mTaxP_N[!is.na(mTaxP_N$Diagnosis),]
tests <- testAbundances(mTaxP_N,saveFolder = paste(outFolderName,'/d3_microbiome/phyla',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d3_microbiome/phyla_statistics.csv',sep=''),sep=',',row.names = F)

# comparae abundances of GENERA
# ====================================================================
mTaxG_N <- mTaxG_N[!is.na(mTaxG_N$Diagnosis),]
tests <- testAbundances(mTaxG_N,saveFolder = paste(outFolderName,'/d3_microbiome/genera',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d3_microbiome/genera_statistics.csv',sep=''),sep=',',row.names = F)

# compare abundances of SPECIES
# ====================================================================
mTaxS_N <- mTaxS_N[!is.na(mTaxS_N$Diagnosis),]
tests <- testAbundances(mTaxS_N,saveFolder = paste(outFolderName,'/d3_microbiome/specie',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d3_microbiome/species_statistics.csv',sep=''),sep=',',row.names = F)

# ===================================
# prep models [IBS vs UC vs CD vs HC]
# ====================================
# refactor diagnosis
allDataIBSUCCDHC <- allDataPrep
allDataIBSUCCDHC <- allDataIBSUCCDHC[allDataIBSUCCDHC$Diagnosis!="IBDU",]
allDataIBSUCCDHC <- allDataIBSUCCDHC[allDataIBSUCCDHC$Diagnosis!="IBS-SR",]
allDataIBSUCCDHC$Diagnosis <- as.factor(as.character(allDataIBSUCCDHC$Diagnosis))

mMG_PWY <- allDataIBSUCCDHC[,c(grep('Diagnosis',colnames(allDataIBSUCCDHC)),grep('k__',colnames(allDataIBSUCCDHC)),grep('PWY',colnames(allDataIBSUCCDHC)) )]
# == biomarkers ==
# calprotectin only
# == taxa: phyla ==
#    pathways & phyla
mTaxP_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#  phyla only
mTaxP_N <- mTaxP_PWY_N[,-grep("PWY",colnames(mTaxP_PWY_N))]

# == taxa: genera ==
#    pathways & genera
mTaxG_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    genera only
mTaxG_N <- mTaxG_PWY_N[,-grep("PWY",colnames(mTaxG_PWY_N))]

# == taxa: species ==
#    pathways & species
mTaxS_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    species only
mTaxS_N <- mTaxS_PWY_N[,-grep("PWY",colnames(mTaxS_PWY_N))]

# === pathways ===
# pathways only
mPWY_N <- mTaxS_PWY_N[,c(grep('Diagnosis',colnames(mTaxS_PWY_N)),grep("PWY",colnames(mTaxS_PWY_N)))]

# compare abundances of pathways:
#=====================================
tests <- testAbundances(mPWY_N,saveFolder = paste(outFolderName,'/d4_pathways/pwy',sep=""),pathway = T,onlyShowSig = T,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d4_pathways/pwy_statistics.csv',sep=""),sep=',',row.names = F)

# compare abundances of PHYLA
# ==============================================================
tests <- testAbundances(mTaxP_N,saveFolder = paste(outFolderName,'/d4_microbiome/phyla',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d4_microbiome/phyla_statistics.csv',sep=''),sep=',',row.names = F)

# comparae abundances of GENERA
# ====================================================================
tests <- testAbundances(mTaxG_N,saveFolder = paste(outFolderName,'/d4_microbiome/genera',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d4_microbiome/genera_statistics.csv',sep=''),sep=',',row.names = F)

# compare abundances of SPECIES
# ====================================================================
tests <- testAbundances(mTaxS_N,saveFolder = paste(outFolderName,'/d4_microbiome/specie',sep=''),pathway = F,doPlots = doAbundancies)
write.table(tests,paste(outFolderName,'/d4_microbiome/species_statistics.csv',sep=''),sep=',',row.names = F)

# =============================================================================
# =============================================================================
# =================== TESTS FOR IBD CASES ONLY ================================
# =============================================================================
# -> here we check: Disease location, Disease score & Flare state

# OPTIONS
# ===============
doAbundancies = T
outFolderName = "IBD_ML_exp_LLDMIBS_v5"
outFolderData = "data_step1_v5"
allData <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv'),sep=',',header=T)
# ==> output folder, create/clean as necessary

if (dir.exists(paste0(outFolderName,'/exploratory_flares')))   {unlink(paste0(outFolderName,'/exploratory_flares'),force = T,recursive = T)}
if (!dir.exists(paste(outFolderName,'/exploratory_flares'))) {dir.create(paste0(outFolderName,'/exploratory_flares'))}

# ======================================================================
# === exploratory data analysis: flare "time-line" ==========
# ======================================================================
allData$flareState <- allData$IBD.flare.clinical
allData <- allData[!is.na(allData$flareState),]
allData$IBD.closest.flare[allData$flareState=="pre.flare"] <- abs(allData$IBD.closest.flare[allData$flareState=="pre.flare"]) * -1
allData$flareState <- as.factor(allData$flareState)

# -- log cut; not very good --
#allData$flareBin <- log(abs(allData$IBD.closest.flare))
#allData$flareBin[allData$flareState == "In.Flare"] <- NA
#allData$flareBinCut <- cut(allData$flareBin,c(0,1,2,3,4,5,6,7,8,9,10))
#allData$flareBinCut <- as.character(allData$flareBinCut)
#allData$flareBinCut[allData$flareState == "In.Flare"] <- "(0)"
#allData$flareBinCut <- as.factor(allData$flareBinCut)

# -- 3 month cut; seems to work fine --
allData$flareBinCut <- cut(allData$IBD.closest.flare,c(-10000,-360,-90,-1,1,
                                                       90,360,10000))
allData[is.na(allData$flareBinCut),]$flareBinCut <- "(360,1e+04]"

#allData$flareCut2 <- cut(allData$IBD.closest.flare,c(-10000,-360,-90,-60,-30,-7,-1,1,
#                                                   7,30,60,90,180,360,10000))
allData <- allData[allData$Diagnosis!="IBDU" & allData$Diagnosis!="HC",]
allData$flareState <- factor(as.character(allData$flareState),levels = c("pre.flare","in.flare","post.flare"))

# ====================================================================
# ------------ calprot -----------------
# ====================================================================
bmData <- allData
bmData <- bmData[!is.na(bmData$Calprot),]
bmData <- bmData[bmData$Diagnosis.4!="HC",]
bmData$Diagnosis <- as.factor(as.character(bmData$Diagnosis.4))

ifUC <- bmData[bmData$flareState=="in.flare" & bmData$Diagnosis=="UC",]
nifUC <- bmData[bmData$flareState!="in.flare" & bmData$Diagnosis=="UC",]
ifCD <- bmData[bmData$flareState=="in.flare" & bmData$Diagnosis=="CD",]
nifCD <- bmData[bmData$flareState!="in.flare" & bmData$Diagnosis=="CD",]

# plot Calprot flarestate
# -- all
outFolderF = paste(outFolderName,'/exploratory_flares/',sep= '')
dodge = position_dodge(width=0.8)
jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
g <- ggplot(bmData,aes(col=flareState,x=Diagnosis,y=Calprot )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("Calprot binning, all IBD cases") + ylab("Calprot (ug/g),ln-transformed")
print(g)

wilcox.test(ifCD$Calprot,nifCD$Calprot)
wilcox.test(ifUC$Calprot,nifUC$Calprot)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_sIBD_fstate.png',sep=''))

g <- ggplot(bmData,aes(col=flareBinCut,x=Diagnosis,y=Calprot )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("Calprot binning, all IBD cases") + ylab("Calprot (ug/g),ln-transformed")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_sIBD_cut.png',sep=''))

g <- ggplot(bmData[bmData$Diagnosis=="UC",],aes(col=flareBinCut,x=flareBinCut,y=Calprot )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("Calprot binning, UC cases")+ ylab("Calprot (ug/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_sIBD_UC_cut.png',sep=''))

g <- ggplot(bmData[bmData$Diagnosis=="CD",],aes(col=flareBinCut,x=flareBinCut,y=Calprot )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("Calprot binning, CD cases")+ ylab("Calprot (ug/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_sIBD_CD_cut.png',sep=''))

#  - do tests -
# ALL IBD
g <- testOneFeature(bmData,feature = "Calprot",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprotectin vs IBD flare state",display = "Pn")[[2]]
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_IBDj_tests.png',sep=''))

# CD
g <- testOneFeature(bmData[bmData$Diagnosis=="CD",],feature = "Calprot",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprotectin vs CD flare state",display = "Pn")[[2]]
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_IBD_CD_tests.png',sep=''))

# UC
g <- testOneFeature(bmData[bmData$Diagnosis=="UC",],feature = "Calprot",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprotectin vs UC flare state",display = "Pn")[[2]]
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_IBD_UC_tests.png',sep=''))


# other stuff
# --> CD location
g <- testOneFeature(bmData[bmData$Diagnosis=="CD" & !is.na(bmData$IBD.Location),],feature = "Calprot",responseVar = "IBD.Location",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprot vs CD location",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_CDloc.png',sep=''))

# --> Mon S
g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.S.UC),],feature = "Calprot",responseVar = "IBD.Mon.S.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprot vs UC Montreal-S Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_MonS.png',sep=''))
# --> Mon E
g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.E.UC),],feature = "Calprot",responseVar = "IBD.Mon.E.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprot vs UC Montreal-E Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_MonE.png',sep=''))
# --> MOn B

g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.B.CD),],feature = "Calprot",responseVar = "IBD.Mon.B.CD",retPlot = T,doSave = F,doViolin = F,
                    yLab = "Calprot (ug/g) ,ln-transformed",title="Calprot vs CD Montreal-B Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_Calprot_MonB.png',sep=''))
#wilcox.test(bmData[!is.na(bmData$IBD.Mon.B.CD) & bmData$IBD.Mon.B.CD=="B3penetrating",]$Calprot,bmData[!is.na(bmData$IBD.Mon.B.CD) & bmData$IBD.Mon.B.CD!="B3penetrating",]$Calprot)

# ====================================================================
# ======================== ChrA ================
# ====================================================================
bmData <- allData
bmData <- bmData[!is.na(bmData$ChrA),]
bmData$ChrA <- log(bmData$ChrA)
bmData$ChrA.Raw <- NULL
bmData$Diagnosis <- as.factor(as.character(bmData$Diagnosis.4))

# plot ChrA flarestate
# -- all
outFolderF = paste(outFolderName,'/exploratory_flares/',sep= '')

dodge = position_dodge(width=0.8)
jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
g <- ggplot(bmData,aes(col=flareState,x=Diagnosis,y=ChrA )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("ChrA binning, all IBD cases") + ylab("ChrA (nmol/g),ln-transformed")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_sIBD_fstate.png',sep=''))

g <- ggplot(bmData,aes(col=flareBinCut,x=Diagnosis,y=ChrA )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("ChrA binning, all IBD cases") + ylab("ChrA (nmol/g),ln-transformed")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_sIBD_cut.png',sep=''))

g <- ggplot(bmData[bmData$Diagnosis=="UC",],aes(col=flareBinCut,x=flareBinCut,y=ChrA )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("ChrA binning, UC cases")+ ylab("ChrA (nmol/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_sIBD_UC_cut.png',sep=''))
wilcox.test(bmData[bmData$Diagnosis=="UC" & bmData$flareBinCut=="(1,90]",]$ChrA,
            bmData[bmData$Diagnosis=="UC" & bmData$flareBinCut!="(1,90]",]$ChrA)

g <- ggplot(bmData[bmData$Diagnosis=="CD",],aes(col=flareBinCut,x=flareBinCut,y=ChrA )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("ChrA binning, CD cases")+ ylab("ChrA (nmol/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_sIBD_CD_cut.png',sep=''))
wilcox.test(bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut=="(-90,-1]",]$ChrA,
            bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut!="(-90,-1]",]$ChrA)

#  - do tests -
# ALL IBD
g <- testOneFeature(bmData,feature = "ChrA",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs IBD flare state",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_IBDj_tests.png',sep=''))

# CD
g <- testOneFeature(bmData[bmData$Diagnosis=="CD",],feature = "ChrA",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs CD flare state",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_IBD_CD_tests.png',sep=''))

# UC
g <- testOneFeature(bmData[bmData$Diagnosis=="UC",],feature = "ChrA",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs UC flare state",display = "Pn")[[2]]
print(g)
wilcox.test(bmData[bmData$Diagnosis=="UC" & bmData$flareState=="In.Flare",]$ChrA,
            bmData[bmData$Diagnosis=="UC" & bmData$flareState!="In.Flare",]$ChrA)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_IBD_UC_tests.png',sep=''))

# other stuff
# --> CD location
g <- testOneFeature(bmData[bmData$Diagnosis=="CD" & !is.na(bmData$IBD.Location),],feature = "ChrA",responseVar = "IBD.Location",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs CD location",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_CDloc.png',sep=''))

# --> Mon S
g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.S.UC),],feature = "ChrA",responseVar = "IBD.Mon.S.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs UC Montreal-S Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_MonS.png',sep=''))

g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.E.UC),],feature = "ChrA",responseVar = "IBD.Mon.E.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs UC Montreal-E Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_MonE.png',sep=''))

g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.B.CD),],feature = "ChrA",responseVar = "IBD.Mon.B.CD",retPlot = T,doSave = F,doViolin = F,
                    yLab = "ChrA (nmol/g) ,ln-transformed",title="ChrA vs CD Montreal-B Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_ChrA_MonB.png',sep=''))
#wilcox.test(bmData[!is.na(bmData$IBD.Mon.B.CD) & bmData$IBD.Mon.B.CD=="B3penetrating",]$ChrA,bmData[!is.na(bmData$IBD.Mon.B.CD) & bmData$IBD.Mon.B.CD!="B3penetrating",]$ChrA)

# ====================================================================
# ======================== HBD2 ================
# ====================================================================
bmData <- allData
bmData <- bmData[!is.na(bmData$HBD2),]
bmData$Diagnosis <- as.factor(as.character(bmData$Diagnosis))

# plot HBD2 flarestate
# -- all
outFolderF = paste(outFolderName,'/exploratory_flares/',sep= '')

dodge = position_dodge(width=0.8)
jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
g <- ggplot(bmData,aes(col=flareState,x=Diagnosis,y=HBD2 )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("HBD2 binning, all IBD cases") + ylab("HBD2 (ng/g),ln-transformed")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_sIBD_fstate.png',sep=''))

g <- ggplot(bmData,aes(col=flareBinCut,x=Diagnosis,y=HBD2 )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("HBD2 binning, all IBD cases") + ylab("HBD2 (ng/g),ln-transformed")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_sIBD_cut.png',sep=''))

g <- ggplot(bmData[bmData$Diagnosis=="UC",],aes(col=flareBinCut,x=flareBinCut,y=HBD2 )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("HBD2 binning, UC cases")+ ylab("HBD2 (ng/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_sIBD_UC_cut.png',sep=''))
wilcox.test(bmData[bmData$Diagnosis=="UC" & bmData$flareBinCut=="(-1,1]",]$HBD2,
            bmData[bmData$Diagnosis=="UC" & bmData$flareBinCut!="(-1,1]",]$HBD2)

g <- ggplot(bmData[bmData$Diagnosis=="CD",],aes(col=flareBinCut,x=flareBinCut,y=HBD2 )) +
  geom_jitter(alpha=0.4,position=jdodge) +
  geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
  ggtitle("HBD2 binning, CD cases")+ ylab("HBD2 (ng/g),ln-transformed") + xlab("Time to Flare (days)")
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_sIBD_CD_cut.png',sep=''))
wilcox.test(bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut=="(-90,-1]",]$HBD2,
            bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut!="(-90,-1]",]$HBD2)
wilcox.test(bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut=="(-1,1]",]$HBD2,
            bmData[bmData$Diagnosis=="CD" & bmData$flareBinCut!="(-1,1]",]$HBD2)



#  - do tests -
# ALL IBD
g <- testOneFeature(bmData,feature = "HBD2",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs IBD flare state",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_IBDj_tests.png',sep=''))

# CD
g <- testOneFeature(bmData[bmData$Diagnosis=="CD",],feature = "HBD2",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs CD flare state",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_IBD_CD_tests.png',sep=''))

# UC
g <- testOneFeature(bmData[bmData$Diagnosis=="UC",],feature = "HBD2",responseVar = "flareState",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs UC flare state",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_IBD_UC_tests.png',sep=''))

# other stuff
# --> CD location
g <- testOneFeature(bmData[bmData$Diagnosis=="CD" & !is.na(bmData$IBD.Location),],feature = "HBD2",responseVar = "IBD.Location",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs CD location",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_CDloc.png',sep=''))

# --> Mon S
g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.S.UC),],feature = "HBD2",responseVar = "IBD.Mon.S.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs UC Montreal-S Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_MonS.png',sep=''))

g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.E.UC),],feature = "HBD2",responseVar = "IBD.Mon.E.UC",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs UC Montreal-E Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_MonE.png',sep=''))

g <- testOneFeature(bmData[!is.na(bmData$IBD.Mon.B.CD),],feature = "HBD2",responseVar = "IBD.Mon.B.CD",retPlot = T,doSave = F,doViolin = F,
                    yLab = "HBD2 (ng/g) ,ln-transformed",title="HBD2 vs CD Montreal-B Criterium",display = "Pn")[[2]]
print(g)
ggsave(g,filename = paste(outFolderF,'exp_flares_HBD2_MonB.png',sep=''))

# ======================================================================
#     ===================== MAKE SUMMARY TABLES ===============
# ======================================================================
source('../myLibs_v3//R_ML_scripts_v4.R')
source('../myLibs_v3/R_Microbiome_scripts.R')
source('../myLibs_v3/R_Misc.R')
# load data
outFolderData <- "data_step1_v5"
outFolderName = "IBD_ML_exp_LLDMIBS_v5"

if (!dir.exists(outFolderName)) {dir.create(outFolderName)}
if (!dir.exists(paste0(outFolderName,'/summary_tables'))) {dir.create(paste0(outFolderName,'/summary_tables'))}

inData <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv'),header=T,sep=',')

# make summary for taxa, DIAGNOSIS.3
# ==========================================
toTest <- colnames(subsetMicrobiomeDF(inDF = inData,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,getDivs = F))
resT <- makeMultiSummaryByResponse(inDF = inData,ftrs = toTest,doTestVsControl = T,controlClass = "HC",
                           response ="Diagnosis.3",doFDR = T,verbose = T)
resT$shortName <- NA
for (c in c(1:nrow(resT))) {
  resT$shortName[c] <- purgeMGNameOne(as.character(resT$Var[c]))
}
# > save it
write.table(resT, paste0(outFolderName,'/summary_tables/taxa_summary_diagnosis3.csv'),sep=',',row.names = F) 

# make summary for taxa, DIAGNOSIS.4
# ==========================================
resT4 <- makeMultiSummaryByResponse(inDF = inData,ftrs = toTest,doTestVsControl = T,controlClass = "HC",
                                   response ="Diagnosis.4",doFDR = T,verbose = T)
resT4$shortName <- NA
for (c in c(1:nrow(resT4))) {
  resT4$shortName[c] <- purgeMGNameOne(as.character(resT4$Var[c]))
}
#  > save it 
write.table(resT4, paste0(outFolderName,'/summary_tables/taxa_summary_diagnosis4.csv'),sep=',',row.names = F) 

# make summary for taxa, FLARES (Marjolein)
# ==========================================
inDataIBD <- inData
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
inDataIBD$IBD.Activity.Clinical <- "Inactive"
inDataIBD$IBD.Activity.Clinical[inDataIBD$IBD.flare.clinical=="in.flare"] <- "Active"
resT.Flare <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = toTest,doTestVsControl = T,controlClass = "Inactive",
                                    response ="IBD.Activity.Clinical",doFDR = T)
#  > save it 
write.table(resT.Flare, paste0(outFolderName,'/summary_tables/taxa_summary_IBDactivity.csv'),sep=',',row.names = F) 

# make summary for taxa, IBD-subtypes
# ==========================================
inDataIBD <- inData
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
inDataIBD$Diagnosis.4 <- as.factor(as.character(inDataIBD$Diagnosis.4))
resT.IBDtype <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = toTest,doTestVsControl = T,controlClass = "CD",
                                         response ="Diagnosis.4",doFDR = T)
#  > save it 
write.table(resT.IBDtype, paste0(outFolderName,'/summary_tables/taxa_summary_IBDtype.csv'),sep=',',row.names = F) 


# =========================================
# PWYS:
# =========================================

# make summary for PWYs, DIAGNOSIS.3
# ==========================================
totestPWY <- colnames(subsetMicrobiomeDF(inData,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,getDivs = F))
resPWYT3 <- makeMultiSummaryByResponse(inDF = inData,ftrs = totestPWY,doTestVsControl = T,controlClass = "HC",
                                       response ="Diagnosis.3",doFDR = T)
#  > save it 
write.table(resPWYT3, paste0(outFolderName,'/summary_tables/pwys_summary_diagnosis3.csv'),sep=',',row.names = F) 

# make summary for PWYs, DIAGNOSIS.4
# ==========================================
totestPWY <- colnames(subsetMicrobiomeDF(inData,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,getDivs = F))
resPWYT4 <- makeMultiSummaryByResponse(inDF = inData,ftrs = totestPWY,doTestVsControl = T,controlClass = "HC",
                                    response ="Diagnosis.4",doFDR = T)
#  > save it 
write.table(resT4, paste0(outFolderName,'/summary_tables/pwys_summary_diagnosis4.csv'),sep=',',row.names = F) 

# make summary for PWYs, FLARES (Marjolein)
# ==========================================
totestPWY <- colnames(subsetMicrobiomeDF(inData,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,getDivs = F))
inDataIBD <- inData
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
inDataIBD$IBD.Activity.Clinical <- "Inactive"
inDataIBD$IBD.Activity.Clinical[inDataIBD$IBD.flare.clinical=="in.flare"] <- "Active"
resT.Flare <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = totestPWY,doTestVsControl = T,controlClass = "Inactive",
                                         response ="IBD.Activity.Clinical",doFDR = T)
#  > save it 
write.table(resT.Flare, paste0(outFolderName,'/summary_tables/pwys_summary_IBDactivity.csv'),sep=',',row.names = F)

# make summary for PWYs, IBD type
# ==========================================
totestPWY <- colnames(subsetMicrobiomeDF(inData,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,getDivs = F))
inDataIBD <- inData
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD$Diagnosis.4 <- as.character(inDataIBD$Diagnosis.4)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
resT.IBDtype <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = totestPWY,doTestVsControl = T,controlClass = "CD",
                                         response ="Diagnosis.4",doFDR = T)
#  > save it 
write.table(resT.IBDtype, paste0(outFolderName,'/summary_tables/pwys_summary_IBDtype.csv'),sep=',',row.names = F)


# ==== make summary for phenotypes, FLARES (Marjolein)
ftrsToTest <- c("Age","BMI","Calprot.nontransformed","HBD2.nontransformed","ChrA.nontransformed")
ftrsToTestB <- c("Gender")
inDataIBD <- inData
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
inDataIBD$IBD.Activity.Clinical <- "Inactive"
inDataIBD$IBD.Activity.Clinical[inDataIBD$IBD.flare.clinical=="in.flare"] <- "Active"
resT.Flare <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = ftrsToTest,doTestVsControl = T,controlClass = "Inactive",
                                         response ="IBD.Activity.Clinical",doFDR = T,doSort = F)
write.table(resT.Flare, paste0(outFolderName,'/summary_tables/features_num_summary_IBDactivity.csv'),sep=',',row.names = F)
resT.Flare2 <- makeMultiSummaryByResponse(inDF = inDataIBD,ftrs = ftrsToTestB,doTestVsControl = T,controlClass = "Inactive",
                                         response ="IBD.Activity.Clinical",doFDR = T,doSort = F)
write.table(resT.Flare2, paste0(outFolderName,'/summary_tables/features_cat_summary_IBDactivity.csv'),sep=',',row.names = F)



# =================================================================
# =================================================================
# ===== make microbiome comparison, FLARES (Marjolein) =====
# 1) 3-category flares
# =================================================================
# =================================================================
inData <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv'),sep=',',header=T)
inDataIBD$Diagnosis.3 <- as.character(inDataIBD$Diagnosis.3)
inDataIBD <- inDataIBD[inDataIBD$Diagnosis.3=="IBD" & !is.na(inDataIBD$Diagnosis.3),]
inDataIBD$IBD.Activity.Clinical <- "Inactive"
inDataIBD$IBD.Activity.Clinical[inDataIBD$IBD.flare.clinical=="in.flare"] <- "Active"

if (dir.exists(paste0(outFolderName,'/flares_3cat_microbiome')))  {unlink(paste0(outFolderName,'/flares_3cat_microbiome'),recursive = T,force = T)}
if (dir.exists(paste0(outFolderName,'/flares_3cat_pathways')))  {unlink(paste0(outFolderName,'/flares_3cat_pathways'),recursive = T,force = T)}
if (!dir.exists(paste0(outFolderName,'/flares_3cat_microbiome'))) {dir.create(paste0(outFolderName,'/flares_3cat_microbiome'))}
if (!dir.exists(paste0(outFolderName,'/flares_3cat_pathways')))   {dir.create(paste0(outFolderName,'/flares_3cat_pathways'))}

# ================================
# prep models [FLARES]
# ================================
inDataIBD$Diagnosis <- as.factor(as.character(inDataIBD$IBD.flare.clinical))
mMG_PWY <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('k__',colnames(inDataIBD)),grep('PWY',colnames(inDataIBD)) )]
# == biomarkers ==
# calprotectin only
mBM_Cal  <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)))]
# HBD2 only
mBM_HBD2 <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)))]
# ChrA only
mBM_ChrA <-  inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
# 2 out of 3
mBM_CalHBD2 <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)))]
mBM_CalChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
mBM_HBD2ChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
# 3 biomarkers
mBM_CalHBD2ChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),
                                      grep('HBD2',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]

# == taxa: phyla ==
#    pathways & phyla
mTaxP_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#  phyla only
mTaxP_N <- mTaxP_PWY_N[,-grep("PWY",colnames(mTaxP_PWY_N))]

# == taxa: genera ==
#    pathways & genera
mTaxG_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    genera only
mTaxG_N <- mTaxG_PWY_N[,-grep("PWY",colnames(mTaxG_PWY_N))]

# == taxa: species ==
#    pathways & species
mTaxS_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    species only
mTaxS_N <- mTaxS_PWY_N[,-grep("PWY",colnames(mTaxS_PWY_N))]

# === pathways ===
# pathways only
mPWY_N <- mTaxS_PWY_N[,c(grep('Diagnosis',colnames(mTaxS_PWY_N)),grep("PWY",colnames(mTaxS_PWY_N)))]


# compare abundances of pathways:
#=====================================
mPWY_N <- mPWY_N[!is.na(mPWY_N$Diagnosis),]
tests <- testAbundances(dataIn = mPWY_N,saveFolder = paste0(outFolderName,'/flares_3cat_pathways/pwy'),pathway = T,
                        onlyShowSig = T,doPlots = doAbundancies)

dataIn = mPWY_N;saveFolder = paste0(outFolderName,'/flares_3cat_pathways/pwy');pathway = T;onlyShowSig = T;doPlots = doAbundancies

write.table(tests,paste0(outFolderName,'/flares_3cat_pathways/pwy_statistics.csv'),sep=',',row.names = F)

# compare abundances of PHYLA
# ==============================================================
mTaxP_N <- mTaxP_N[!is.na(mTaxP_N$Diagnosis),]
tests <- testAbundances(mTaxP_N,saveFolder = paste0(outFolderName,'/flares_3cat_microbiome/phyla'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_3cat_microbiome/phyla_statistics.csv'),sep=',',row.names = F)

# comparae abundances of GENERA
# ====================================================================
mTaxG_N <- mTaxG_N[!is.na(mTaxG_N$Diagnosis),]
tests <- testAbundances(mTaxG_N,saveFolder = paste0(outFolderName,'/flares_3cat_microbiome/genera'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_3cat_microbiome/genera_statistics.csv'),sep=',',row.names = F)

# compare abundances of SPECIES
# ====================================================================
mTaxS_N <- mTaxS_N[!is.na(mTaxS_N$Diagnosis),]
tests <- testAbundances(mTaxS_N,saveFolder = paste0(outFolderName,'/flares_3cat_microbiome/specie'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_3cat_microbiome/species_statistics.csv'),sep=',',row.names = F)

# =================================================================
# ===== make microbiome comparison, FLARES (Marjolein) =====
# 1) 2-category flares
# =================================================================
if (dir.exists(paste0(outFolderName,'/flares_yn_microbiome'))) {unlink(paste0(outFolderName,'/flares_yn_microbiome'),recursive = T,force = T)}
if (dir.exists(paste0(outFolderName,'/flares_yn_pathways')))   {unlink(paste0(outFolderName,'/flares_yn_pathways'),recursive = T,force = T)}
if (!dir.exists(paste0(outFolderName,'/flares_yn_microbiome'))) {dir.create(paste0(outFolderName,'/flares_yn_microbiome'))}
if (!dir.exists(paste0(outFolderName,'/flares_yn_pathways')))   {dir.create(paste0(outFolderName,'/flares_yn_pathways'))}

# ================================
# prep models [FLARES]
# ================================
inDataIBD$Diagnosis <- as.character(inDataIBD$IBD.flare.clinical)
inDataIBD$Diagnosis[inDataIBD$Diagnosis!="in.flare"] <- "inactive.IBD"
inDataIBD$Diagnosis[inDataIBD$Diagnosis=="in.flare"] <- "active.IBD"
inDataIBD$Diagnosis <- as.factor(inDataIBD$Diagnosis)

mMG_PWY <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('k__',colnames(inDataIBD)),grep('PWY',colnames(inDataIBD)) )]
# == biomarkers ==
# calprotectin only
mBM_Cal  <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)))]
# HBD2 only
mBM_HBD2 <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)))]
# ChrA only
mBM_ChrA <-  inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
# 2 out of 3
mBM_CalHBD2 <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)))]
mBM_CalChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
mBM_HBD2ChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('HBD2',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]
# 3 biomarkers
mBM_CalHBD2ChrA <- inDataIBD[,c(grep('Diagnosis',colnames(inDataIBD)),grep('Calprot',colnames(inDataIBD)),
                                grep('HBD2',colnames(inDataIBD)),grep('ChrA',colnames(inDataIBD)))]

# == taxa: phyla ==
#    pathways & phyla
mTaxP_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#  phyla only
mTaxP_N <- mTaxP_PWY_N[,-grep("PWY",colnames(mTaxP_PWY_N))]

# == taxa: genera ==
#    pathways & genera
mTaxG_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    genera only
mTaxG_N <- mTaxG_PWY_N[,-grep("PWY",colnames(mTaxG_PWY_N))]

# == taxa: species ==
#    pathways & species
mTaxS_PWY_N <- purgeMGNames(filterMetaGenomeDF(mMG_PWY,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = F))
#    species only
mTaxS_N <- mTaxS_PWY_N[,-grep("PWY",colnames(mTaxS_PWY_N))]

# === pathways ===
# pathways only
mPWY_N <- mTaxS_PWY_N[,c(grep('Diagnosis',colnames(mTaxS_PWY_N)),grep("PWY",colnames(mTaxS_PWY_N)))]


# compare abundances of pathways:
#=====================================
mPWY_N <- mPWY_N[!is.na(mPWY_N$Diagnosis),]
tests <- testAbundances(dataIn = mPWY_N,saveFolder = paste0(outFolderName,'/flares_yn_pathways/pwy'),pathway = T,
                        onlyShowSig = T,doPlots = doAbundancies)

dataIn = mPWY_N;saveFolder = paste0(outFolderName,'/flares_yn_pathways/pwy')

write.table(tests,paste0(outFolderName,'/flares_yn_pathways/pwy_statistics.csv'),row.names = F)

# compare abundances of PHYLA
# ==============================================================
mTaxP_N <- mTaxP_N[!is.na(mTaxP_N$Diagnosis),]
tests <- testAbundances(mTaxP_N,saveFolder = paste0(outFolderName,'/flares_yn_microbiome/phyla'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_yn_microbiome/phyla_statistics.csv'),sep=',',row.names = F)

# comparae abundances of GENERA
# ====================================================================
mTaxG_N <- mTaxG_N[!is.na(mTaxG_N$Diagnosis),]
tests <- testAbundances(mTaxG_N,saveFolder = paste0(outFolderName,'/flares_yn_microbiome/genera'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_yn_microbiome/genera_statistics.csv'),sep=',',row.names = F)

# compare abundances of SPECIES
# ====================================================================
mTaxS_N <- mTaxS_N[!is.na(mTaxS_N$Diagnosis),]
tests <- testAbundances(mTaxS_N,saveFolder = paste0(outFolderName,'/flares_yn_microbiome/specie'),pathway = F,doPlots = doAbundancies)
write.table(tests,paste0(outFolderName,'/flares_yn_microbiome/species_statistics.csv'),sep=',',row.names = F)