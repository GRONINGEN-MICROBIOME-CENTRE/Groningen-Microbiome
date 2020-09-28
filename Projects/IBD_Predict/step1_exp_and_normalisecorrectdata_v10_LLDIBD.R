# =================================================
# By R.Gacesa, UMCG (2019)
# =================================================
# load libs
library(tidyr)
library(vegan)
library(Rtsne)
#library(VennDiagram)
setwd("D:/UMCG/ML_IBD_v5/")
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

# ==> WD
#setwd("C:/Users/ranko/Documents/UMCG/ML_IBD_vs_IBS/")
allData <- read.table('data_clean/withBMDataMergedCleanRSPathways_v10.csv',stringsAsFactors = F,sep=',',header = T)

# ==> output folder, create/clean as necessary
doClean = T
if (doClean & dir.exists(outFolderName)) {unlink(outFolderName,recursive = T)}
#if (doClean & dir.exists(outFolderData)) {unlink(outFolderData,recursive = T)}
if (!dir.exists(outFolderName)) {dir.create(outFolderName)}
if (!dir.exists(outFolderData)) {dir.create(outFolderData)}
if (!dir.exists(paste(outFolderName,'/exploratory',sep=''))) {dir.create(paste(outFolderName,'/exploratory',sep=''))}
if (!dir.exists(paste(outFolderName,'/d3_microbiome',sep=''))) {dir.create(paste(outFolderName,'/d3_microbiome',sep=''))}
if (!dir.exists(paste(outFolderName,'/d3_pathways',sep='')))   {dir.create(paste(outFolderName,'/d3_pathways',sep=''))}
if (!dir.exists(paste(outFolderName,'/d4_microbiome',sep=''))) {dir.create(paste(outFolderName,'/d4_microbiome',sep=''))}
if (!dir.exists(paste(outFolderName,'/d4_pathways',sep='')))   {dir.create(paste(outFolderName,'/d4_pathways',sep=''))}


# select cohort
# ==========================
allData <- allData[allData$Cohort=="LLD_IBD" | allData$Cohort=="MIBS",]
# >>> CLEAN DATA:
# no diagnosis is bad:
allData <- allData[!is.na(allData$Diagnosis),]
# high percentage of unclassified taxa is also bad (basically means something went wrong during sequencing)
allData$unclassified[is.na(allData$unclassified)] <- 0.0
allDataT <- allData[as.numeric(as.character(allData$unclassified)) > 50.0,]
allData <- allData[as.numeric(as.character(allData$unclassified)) < 50.0,]
allData$unclassified <- NULL
rm(allDataT)
# - MC has only one case, we can't use this
allData <- allData[allData$Diagnosis != "MC",]
# fix factor for Diagnosis
allData$Diagnosis <- as.factor(as.character(allData$Diagnosis))
# fix factor for gender
allData$Gender <- as.factor(as.character(allData$Gender))
# make sure other variables are in right format
#allData$Calprot <- as.numeric(as.character(allData$CalprotD40))
allData$CalprotD40 <- NULL
allData$HBD2 <- as.numeric(as.character(allData$HBD2))
allData$Calprot <- as.numeric(as.character(allData$Calprot))
allData$ChrA <- as.numeric(as.character(allData$ChrA))
allData$Age <- as.numeric(allData$Age)
for (i in grep('k__',colnames(allData))) {
  allData[,i] <- as.numeric(as.character(allData[,i]))
}
# keep smokers (too many NAs [they are actually " " in data, but still NA])
allData$Smoking <- NULL

# ======================================================
# FILTERING !!!
# ======================================================

# # REMOVE people for which we don't have biomarkers
# allData <- allData[!is.na(allData$Calprot) & !is.na(allData$HBD2) & !is.na(allData$ChrA),]
# 
# # REMOVE stoma people - their metagenome is messed up
# allData$StomaPouchType[is.na(allData$StomaPouchType)] <- "None"
# allData$StomaPouch[is.na(allData$StomaPouch)] <- "N"
# allData$StomaPouch[allData$StomaPouchType != "None"] <- "Y"
# allData <- allData[allData$StomaPouch != "Y",]
# allData$StomaPouch <- NULL
# allData$StomaPouchType <- NULL

# assume NA means N drugs & smoking (but not biomarkers)
allData$Antidiarrhoea[is.na(allData$Antidiarrhoea)] <- "N"
allData$Antidiarrhoea <- as.factor(allData$Antidiarrhoea)
allData$PPI[is.na(allData$PPI)] <- "N"
allData$PPI <- as.factor(allData$PPI)
#allData$Smoking[is.na(allData$Smoking)] <- "N"
allData$NSAID[is.na(allData$NSAID)] <- "N"
allData$NSAID <- as.factor(allData$NSAID)

allData$SSRI[is.na(allData$SSRI)] <- "N"
allData$SSRI <- as.factor(allData$SSRI)

allData$Laxatives[is.na(allData$Laxatives)] <- "N"
allData$Laxatives <- as.factor(allData$Laxatives)

allData$Statins[is.na(allData$Statins)] <- "N"
allData$Statins <- as.factor(allData$Statins)

allData$Antibiotics[is.na(allData$Antibiotics)] <- "N"
allData$Antibiotics <- as.factor(allData$Antibiotics)

allData$Immunosup[is.na(allData$Immunosup)] <- "N"
allData$Immunosup <- as.factor(allData$Immunosup)

allData$Oral_Contraceptives[is.na(allData$Oral_Contraceptives)] <- "N"
allData$Oral_Contraceptives <- as.factor(allData$Oral_Contraceptives)

# fix below 40 values for Calprot
allData$Calprot[allData$Calprot < 40.0] <- 40.0
allData$Calprot[allData$Calprot > 4000.0] <- 4000

# refactor diagnosis (3-class: HC vs IBD vs IBS), IBS-SR is not included
allData$Diagnosis.3 <- as.character(allData$Diagnosis)
allData$Diagnosis.3[grep("IBS",allData$Diagnosis)] <- "IBS" 
allData$Diagnosis.3[allData$Diagnosis=="IBS-SR"] <- NA
allData$Diagnosis.3[allData$Diagnosis=="UC"] <- "IBD" 
allData$Diagnosis.3[allData$Diagnosis=="CD"] <- "IBD"
allData$Diagnosis.3[allData$Diagnosis=="IBDU"] <- "IBD" 

# refactor diagnosis (4-class: HC vs UC vs CD), IBS-SR is not included and IBDU is not included
allData$Diagnosis.4 <- as.character(allData$Diagnosis)
allData$Diagnosis.4[grep("IBS",allData$Diagnosis)] <- "IBS" # HC vs CD vs UC vs IBD
allData$Diagnosis.4[allData$Diagnosis=="IBS-SR"] <- NA
allData$Diagnosis.4[grep("IBDU",allData$Diagnosis)] <- NA
allData$ChrA.orig <- allData$ChrA


# ============= TEST BIOMARKERS [before fixing] ================
# ChrA
testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.4), ],feature = "ChrA",responseVar = "Diagnosis.4",
               saveFolder = paste(outFolderName,'/exploratory/e_ChrA_raw_04',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "ChrA Raw (ug/g)",doSave = T,doPlots = T,doViolin=F,title = "ChrA (raw measurements)")

testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.3), grep('ChrA.Raw',colnames(allData),invert = T)],feature = "ChrA",responseVar = "Diagnosis.3",
               saveFolder = paste(outFolderName,'/exploratory/e_ChrA_raw_03',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "ChrA Raw (ug/g)",doSave = T,doPlots = T,doViolin=F,title = "ChrA (raw measurements)")

# --------------------------------------
# FIX BAD LAB MEASUREMENTS:
# MIBS & LLD: ChrA values should be divided by 50
# --> LLD
allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="HC",]$ChrA <- allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="HC",]$ChrA / 50.0 
allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="IBS-SR",]$ChrA <- allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="IBS-SR",]$ChrA / 50.0 
# --> MIBS
allData[allData$Cohort=="MIBS",]$ChrA <- allData[allData$Cohort=="MIBS",]$ChrA / 50.0 
# IBD measurements: ChrA values should be divided by 1000
# --> IBD
allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis %in% c("IBDU","UC","CD","MC"),]$ChrA <- allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis %in% c("IBDU","UC","CD","MC"),]$ChrA / 1000.0

# ======= compare with new measruemtn ======
# ==========================================

newChr <- read.table('data_raw/chr_remeasurment.csv',sep=',',header = T,stringsAsFactors = F)
newChr <- newChr[,c("Original.ID","ChrA.nmol.g","Box","Box.pos")]
colnames(newChr) <- c("ID","ChrA.new","Box","Box.pos")
oldChr <- allData[,c("ID","ChrA","Diagnosis")]
colnames(oldChr) <- c("ID","ChrA.old","Diagnosis")
mergedChr <- merge(newChr,oldChr,by="ID")
chrLM <- lm(data=mergedChr,ChrA.new ~ ChrA.old + 0,method = "qr")
mergedChr$Diagnosis <- as.character(mergedChr$Diagnosis)
mergedChr$Diagnosis[mergedChr$Diagnosis=="IBDU"] <- "UC"
mergedChr$Diagnosis[mergedChr$Diagnosis=="IBS-SR"] <- "HC"
mergedChr$Diagnosis <- as.factor(mergedChr$Diagnosis)
mergedChr$Box <- as.factor(mergedChr$Box)

mergedChrL <- mergedChr[,c(1,3,4,2,5,6)]
#mergedChrL$ID <- NULL
mergedChrL <- gather(mergedChrL,'Dataset','ChrA.nmol.g',ChrA.new:ChrA.old)
mergedChrL$ALI <- "IBD.new"
mergedChrL$ALI[mergedChrL$Diagnosis == "HC" & mergedChrL$Box.pos <= 40] <- "LLD.old"
mergedChrL$ALI[mergedChrL$Diagnosis == "HC" & mergedChrL$Box.pos > 40] <- "LLD.new"
mergedChrL$ALI <- as.factor(mergedChrL$ALI)
mergedChrL$Dataset <- as.factor(mergedChrL$Dataset)

# summary
sX = 6; sY=4
g <- ggplot(mergedChrL,aes(x=ALI,y=ChrA.nmol.g,col=Dataset)) + geom_boxplot() + ggtitle('Summary boxplot, aliquots') + ylim(0,2.5)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_ali.png',sep=''),width = sX,height = sY)


mergedChrL$Diagnosis.2 <- as.character(mergedChrL$Diagnosis)
mergedChrL$Diagnosis.2[mergedChrL$Diagnosis.2!="HC"] <- "IBD"
mergedChrL$Diagnosis.2 <- as.factor(mergedChrL$Diagnosis.2)

g <- ggplot(mergedChrL,aes(x=Diagnosis.2,y=ChrA.nmol.g,col=Dataset)) + geom_boxplot() + ggtitle('Summary boxplot, diagnosis') + ylim(0,2.5)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_diag.png',sep=''),width = sX,height = sY)

t <- mergedChrL[mergedChrL$Dataset=="ChrA.new",]
g <- ggplot(t,aes(x=Diagnosis.2,y=ChrA.nmol.g,col=Diagnosis.2)) + geom_boxplot() + ggtitle('Summary boxplot, diagnosis') + ylim(0,2.5)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_diag_newonly.png',sep=''),width = sX,height = sY)
g <- ggplot(t,aes(x=ChrA.nmol.g,col=Diagnosis.2)) + geom_density(size=1.5) + ggtitle('Density (new ChrA), diagnosis') + xlim(0,2.5)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_diag_newonly_density.png',sep=''),width = sX,height = sY)
g <- ggplot(t,aes(x=ChrA.nmol.g,col=Diagnosis)) + geom_density(size=1.5) + ggtitle('Density (new ChrA), diagnosis') + xlim(0,2.5)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_diag3_newonly_density.png',sep=''),width = sX,height = sY)

#   > get distribution of new vs old measurements
tst <- mergedChrL[mergedChrL$Dataset=="ChrA.new",]
g <- ggplot(tst,aes(col=ALI,x=ChrA.nmol.g)) + geom_density(size=1.5) + ggtitle('ChrA new measurements')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_new_density.png',sep=''),width = sX,height = sY)

chrAnew_LLD <- mergedChrL[mergedChrL$Dataset == "ChrA.new" & mergedChrL$Diagnosis=="HC",]$ChrA.nmol.g
chrAold_LLD <- mergedChrL[mergedChrL$Dataset == "ChrA.old" & mergedChrL$Diagnosis=="HC",]$ChrA.nmol.g 
chrAoldAll_LLD <- allData[allData$Diagnosis=="HC" | allData$Diagnosis=="IBS-SR",]$ChrA; chrAoldAll_LLD <- chrAoldAll_LLD[!is.na(chrAoldAll_LLD)]
tst2 <- rbind.data.frame(data.frame(Data="New",ChrA=chrAnew_LLD),data.frame(Data="Old",ChrA=chrAold_LLD))
tst2 <- rbind.data.frame(tst2,data.frame(Data="Old_All",ChrA=chrAoldAll_LLD))
g <- ggplot(tst2,aes(col=Data,x=ChrA)) + geom_density(size=1.5) + xlim(0,3) + ggtitle('ChrA Old vs New')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_new_vs_old_density.png',sep=''),width = sX,height = sY)

# ======== correct HC =========
fixed <- tst2
chrAdelta <- abs(median(fixed$ChrA[fixed$Data=="Old"]) - median(fixed$ChrA[tst2$Data=="New"]))
#chrAdelta <- abs(median(fixed$ChrA[fixed$Data=="Old_All"]) - median(fixed$ChrA[tst2$Data=="New"]))
#chrAdelta <- mlv(fixed$ChrA[fixed$Data=="New"], method = "mfv")[[1]] -  mlv(fixed$ChrA[fixed$Data=="Old_All"], method = "mfv")[[1]]
fixed$ChrA[fixed$Data!='New'] <- fixed$ChrA[fixed$Data!='New'] + chrAdelta
g <- ggplot(fixed,aes(col=Data,x=ChrA)) + geom_density(size=1.5) + xlim(0,3) + ggtitle('ChrA corrected measurements')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_fixed_d.png',sep=''),width = sX,height = sY)
g <- ggplot(fixed,aes(col=Data,x=Data,y=ChrA)) + geom_boxplot() + ylim(0,3) + ggtitle('ChrA corrected measurements')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_fixed_b.png',sep=''),width = sX,height = sY)
#plot(fixed$ChrA[fixed$Data=="New"],fixed$ChrA[fixed$Data=="Old"]) 
lmON <- lm(fixed$ChrA[fixed$Data=="New"] ~ fixed$ChrA[fixed$Data=="Old"] + 0)
f2 <- data.frame(New=fixed$ChrA[fixed$Data=="New"],Old=fixed$ChrA[fixed$Data=="Old"])
g <- ggplot(f2,aes(x=Old,y=New)) + geom_point() + geom_smooth(method="lm",formula=y ~ x + 0) + ggtitle('Correlation for ChrA New vs Old')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_fixed_corr.png',sep=''),width = sX,height = sY)

f3 <- rbind.data.frame(fixed[fixed$Data=="Old_All",],data.frame(Data="IBD",ChrA=allData$ChrA[allData$Diagnosis.3=="IBD"]))

# ========= correct IBD ==============
tst2 <- mergedChrL[mergedChrL$Diagnosis != "HC",]
g <- ggplot(data = tst2,aes(x=Dataset,y=ChrA.nmol.g,col=Dataset)) + geom_boxplot() + ggtitle('ChrA IBD re-measurement')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_IBD_b.png',sep=''),width = sX,height = sY)
g <- ggplot(data = tst2,aes(x=ChrA.nmol.g,col=Dataset)) + geom_density(size=1.5) + xlim(0,3) + ggtitle('ChrA IBD re-measurement')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_IBD_d.png',sep=''),width = sX,height = sY)

chrAdeltaIBD <- median(tst2$ChrA.nmol.g[tst2$Data=="ChrA.new"]) - median(tst2$ChrA.nmol.g[tst2$Data=="ChrA.old"])
fixed <- tst2
fixed$ChrA.nmol.g[fixed$Data=="ChrA.old"] <- fixed$ChrA[fixed$Data=="ChrA.old"] + chrAdeltaIBD
g <- ggplot(data = fixed,aes(x=Dataset,y=ChrA.nmol.g,col=Dataset)) + geom_boxplot() + ggtitle('ChrA IBD 1st correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_fixed_IBD_b.png',sep=''),width = sX,height = sY)
g <- ggplot(data = fixed,aes(x=ChrA.nmol.g,col=Dataset)) + geom_density(size=1.5)+ xlim(0,3) + ggtitle('ChrA IBD 1st correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/remeasure_ChrA_fixed_IBD_b.png',sep=''),width = sX,height = sY)

# general
#ggplot(mergedChrL,aes(y=ChrA.nmol.g,x=Diagnosis,col=Dataset,position_dodge())) + geom_boxplot(outlier.alpha = 0.5)+ ggtitle('Summary boxplot, diagnosis')
#ggplot(mergedChrL,aes(y=ChrA.nmol.g,x=Box,col=Dataset,position_dodge())) + geom_boxplot(outlier.alpha = 0.5)+ ggtitle('Summary boxplot, boxes')

# IBD boxplots
#ggplot(mergedChrL[mergedChrL$Diagnosis != "HC",],aes(y=ChrA.nmol.g,x=Diagnosis,col=Dataset,position_dodge())) + geom_boxplot(outlier.alpha = 0.5)+ ggtitle('Summary boxplot, IBD ChrA')
#ggplot(mergedChrL[mergedChrL$Diagnosis != "HC",],aes(x=ChrA.nmol.g,col=Dataset)) + geom_density(size=1.5) + ggtitle('Summary boxplot, IBD ChrA')

# IBD only correlation
#ggplot(mergedChr[mergedChr$Diagnosis != "HC",],aes(ChrA.new,ChrA.old)) + geom_point() + geom_smooth(method = "lm",formula= y ~ x + 0) + xlim(0,2) + ylim(0,2)
#lmIBD <- lm(mergedChr[mergedChr$Diagnosis != "HC",], formula=ChrA.new ~ ChrA.old + 0)

# ================= DO EXTRA ADJUSTMENT FOR MEANS OF CHRA ========================
# move ChrA y chrAdelta
allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="HC",]$ChrA <- allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="HC",]$ChrA + chrAdelta
allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="IBS-SR",]$ChrA <- allData[allData$Cohort=="LLD_IBD" & allData$Diagnosis=="IBS-SR",]$ChrA + chrAdelta
# --> MIBS
allData[allData$Cohort=="MIBS",]$ChrA <- allData[allData$Cohort=="MIBS",]$ChrA + chrAdelta
# --> IBD
for (r in c(1:nrow(allData))) {
  if (!is.na(allData$ChrA[r]) & allData$Cohort[r] == "LLD_IBD" & !is.na(allData$Diagnosis.3[r] == "IBD") & allData$Diagnosis.3[r] == "IBD") {
    allData$ChrA[r] <- allData$ChrA[r] + chrAdeltaIBD
  }
}
# =================================================================================

# retest
g <- ggplot(allData[!is.na(allData$Diagnosis.3),],aes(y=ChrA,x=Diagnosis.3,col=Diagnosis.3)) + geom_boxplot() + ylim(0,4) + ggtitle('ChrA 2nd correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/e_ChrA_d3_fixed_b.png',sep=''),width = sX,height = sY)

g <- ggplot(allData[!is.na(allData$Diagnosis.3),],aes(x=ChrA,col=Diagnosis.3)) + geom_density(size=1.5) + xlim(0,4) + ggtitle('ChrA 2nd correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/e_ChrA_d3_fixed_d.png',sep=''),width = sX,height = sY)

g <- ggplot(allData[!is.na(allData$Diagnosis.4),],aes(y=ChrA,x=Diagnosis.4,col=Diagnosis.4)) + geom_boxplot() + ylim(0,4) + ggtitle('ChrA 2nd correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/e_ChrA_d4_fixed_b.png',sep=''),width = sX,height = sY)

g <- ggplot(allData[!is.na(allData$Diagnosis.4),],aes(x=ChrA,col=Diagnosis.4)) + geom_density(size=1.5) + xlim(0,4) = ggtitle('ChrA 2nd correction')
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/e_ChrA_d4_fixed_d.png',sep=''),width = sX,height = sY)

# do some stats
chrA.IBD <- allData$ChrA[allData$Diagnosis.3=="IBD"]; chrA.IBD <- chrA.IBD[!is.na(chrA.IBD)]
chrA.CD <- allData$ChrA[allData$Diagnosis.4=="CD"]; chrA.IBD <- chrA.IBD[!is.na(chrA.IBD)]
chrA.UC <- allData$ChrA[allData$Diagnosis.4=="UC"]; chrA.IBD <- chrA.IBD[!is.na(chrA.IBD)]
chrA.HC  <- allData$ChrA[allData$Diagnosis.3=="HC"]; chrA.HC <- chrA.IBD[!is.na(chrA.HC)]
chrA.IBS <- allData$ChrA[allData$Diagnosis.3=="IBS"]; chrA.IBS <- chrA.IBD[!is.na(chrA.IBS)]
# --> basically there is nothing there!
wilcox.test(chrA.CD,chrA.HC)
wilcox.test(chrA.UC,chrA.HC)
wilcox.test(chrA.IBD,chrA.HC)
wilcox.test(chrA.IBS,chrA.HC) #bit here?

# ============= TEST BIOMARKERS +================
# ChrA
# DEBUG
dataIn = allData[!is.na(allData$Diagnosis.4)];feature = "ChrA";responseVar = "Diagnosis.4";
               saveFolder = paste(outFolderName,'/exploratory/e_ChrA_fixed_04',sep="");display = "P";onlyShowSig = T;
               retPlot = T;yLab = "ChrA Fixed (ug/g)";doSave = T;doPlots = T;doViolin=F;title = "ChrA (fixed measurements)"

testOneFeature(dataIn = allData,feature = "ChrA",responseVar = "Diagnosis.4",ylim=c(0,4),
               saveFolder = paste(outFolderName,'/exploratory/e_ChrA_fixed_04',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "ChrA Fixed (ug/g)",doSave = T,doPlots = T,doViolin=F,title = "ChrA (fixed measurements)")


testOneFeature(dataIn = allData,feature = "ChrA",responseVar = "Diagnosis.3",ylim=c(0,4),
               saveFolder = paste(outFolderName,'/exploratory/e_ChrA_fixed_03',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "ChrA Fixed (ug/g)",doSave = T,doPlots = T,doViolin=F,title = "ChrA (fixed measurements)")

# HBD2
testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.4),],feature = "HBD2",responseVar = "Diagnosis.4",
               saveFolder = paste(outFolderName,'/exploratory/e_HBD2_04',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "HBD2 Raw (ng/g)",doSave = T,doPlots = T,doViolin=F,,ylim = c(0,500))

testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.3),],feature = "HBD2",responseVar = "Diagnosis.3",
               saveFolder = paste(outFolderName,'/exploratory/e_HBD2_03',sep=""),display = "P",onlyShowSig = T,
               retPlot = F,yLab = "HBD2 Raw (ng/g)",doSave = T,doPlots = T,doViolin=F,ylim = c(0,500))

# Calprot
testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.4),],feature = "Calprot",responseVar = "Diagnosis.4",
               saveFolder = paste(outFolderName,'/exploratory/e_Cal_04',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "Calprotectin Raw (nmol/g)",doSave = T,doPlots = T,doViolin=F,ylim = c(0,750))

testOneFeature(dataIn = allData[!is.na(allData$Diagnosis.3),],feature = "Calprot",responseVar = "Diagnosis.3",
               saveFolder = paste(outFolderName,'/exploratory/e_Cal_03',sep=""),display = "P",onlyShowSig = T,
               retPlot = T,yLab = "Calprotectin Raw (nmol/g)",doSave = T,doPlots = T,doViolin=F,ylim=c(0,750))

bmt <- allData[,c("ID","Diagnosis","Cohort","ChrA","Calprot","HBD2")]
write.table(bmt,paste0(outFolderData,'/LLDIBDMIBS_biomarkers.csv'),sep=',',row.names = F)


# clear IDs (useless)
#allData$ID <- NULL

#dataIn = allData[!is.na(allData$Diagnosis.4), grep('HBD2',colnames(allData),invert = T)];feature = "HBD2";responseVar = "Diagnosis.4";
#saveFolder = 'C:/Users/ranko/Documents/UMCG/ML_IBD_v2/test/test';display = "FDR";onlyShowSig = T;
#retPlot = T;yLab = "testor";doSave = F;doPlots = T;doViolin=F

# --------------------------------------


# - compare diagnosis vs biomarker
sX = 9
sY = 6
g <- ggplot(data=allData,aes(x=Diagnosis,y=HBD2,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4) + ylim(0,500)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_bar_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allData,aes(x=Diagnosis,y=ChrA,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4) + ylim(0,1500)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_bar_ChrA.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allData,aes(x=Diagnosis,y=Calprot,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4) + ylim(0,2000)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_bar_Calprot.png',sep=''),width = sX,height = sY)
# - pretty histograms
g <- ggplot(data=allData,aes(x=HBD2,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_hist_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allData,aes(x=ChrA,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_hist_ChrA.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allData,aes(x=Calprot,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/raw_diag5_hist_Calprot.png',sep=''),width = sX,height = sY)

# log-transform Biomarkers to 'smoothen' the distribution
allData$Calprot.nontransformed <- (allData$Calprot)
allData$HBD2.nontransformed <- (allData$HBD2)
allData$ChrA.nontransformed <- (allData$ChrA)

allData$Calprot <- log(allData$Calprot)
allData$HBD2 <- log(allData$HBD2)
allData$ChrA <- log(allData$ChrA)

# NOPE kill cases containing NA biomarkers
#allData <- allData[!is.na(allData$Calprot) & !is.na(allData$HBD2) & !is.na(allData$ChrA),]

# filter pathway data and taxonomy for abundance & stuff:
allData <- filterHumannDF(inDF = allData,rescale = F,verbose = T)
allData <- filterMetaGenomeDF(inDF=allData,keepDomains = c("Bacteria","Archaea"),rescaleTaxa = T,keepLevels = c("S","G","F","O","C","P"))
write.table(allData,file=paste0(outFolderData,"/allData_Taxa_PWY_cleaned_filtered_v10.csv"),sep=',',row.names = F)
# =========================================================================
# general parameters for data preparation 
# =========================================================================

# =================================================================================================
# =================================================================================================
# ------------------ ANALYSIS: clinical IBS & IBD (Self-reported IBS is excluded from datasets) -
# =================================================================================================
# =================================================================================================

# =========================================================================
#   PREP LLD and MIBS
# =========================================================================
# ===============
# SELECT COHORT
# ===============
allDataR <- allData #[allData$Cohort=="IBD" | allData$Cohort=="MIBS",]
# cleanup
#allDataR$Cohort <- NULL
# prep diagnosis (IBS & IBD, remove subtypes)
allDataR$Diagnosis <- as.character(allDataR$Diagnosis)
# only IBS vs IBD
#allDataR <- allDataR[!allDataR$Diagnosis=="HC",]
# correct for IBS
#allDataR <- allDataR[!allDataR$Diagnosis=="IBS-SR",]
allDataR$Diagnosis[grep("IBS-SR",allDataR$Diagnosis)] <- "TEMPTEMP"
allDataR$Diagnosis[grep("IBS",allDataR$Diagnosis)] <- "IBS"
allDataR$Diagnosis[grep("TEMPTEMP",allDataR$Diagnosis)] <- "IBS-SR"

allDataR$Diagnosis <- as.factor(as.character(allDataR$Diagnosis))

# do some exploratory data analysis
# ================================================

# - compare diagnosis vs biomarker
g <- ggplot(data=allDataR,aes(x=Diagnosis,y=HBD2,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bar_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataR,aes(x=Diagnosis,y=ChrA,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bar_ChrA.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataR,aes(x=Diagnosis,y=Calprot,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.4)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bar_Calprot.png',sep=''),width = sX,height = sY)
# - pretty histograms
g <- ggplot(data=allDataR,aes(x=HBD2,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_hist_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataR,aes(x=ChrA,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_hist_ChrA.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataR,aes(x=Calprot,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_hist_Calprot.png',sep=''),width = sX,height = sY)

# shannon
shannonRdy <- allDataR[,c(grep('Diagnosis',colnames(allDataR)),grep('s__',colnames(allDataR))) ]
shannonRdy$div.shannon <- diversity(shannonRdy[,-grep('Diagnosis',colnames(shannonRdy))],MARGIN = 1, index="shannon")
g <- ggplot(data=shannonRdy,aes(x=Diagnosis,y=div.shannon,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(alpha=0.2)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bar_shannon.png',sep=''),width = sX,height = sY)
# B-C
bcReady <- allDataR[,grep('s__',colnames(allDataR))]
bcurtis <- vegdist(bcReady,method = "bray")
r.pca <- prcomp(bcurtis, center = F,scale. = F)
pcas <- as.data.frame(r.pca$rotation[,1:5])
pcas$Diagnosis <- allDataR$Diagnosis
# PC 1-2
g <- ggplot(pcas ,aes(x = PC1,y=PC2,col=Diagnosis )) + geom_point(size=2.25,alpha=0.4) #+ scale_color_manual(values=c("#009933", "#CC0000", "#FF9933"))
centroids <- aggregate(cbind(PC1,PC2)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_12.png',sep=''),width = sX,height = sY)
# PC 2-3
g <- ggplot(pcas ,aes(x = PC2,y=PC3,col=Diagnosis )) + geom_point(size=2.25,alpha=0.4) 
centroids <- aggregate(cbind(PC2,PC3)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_23.png',sep=''),width = sX,height = sY)
# PC 3-4
g <- ggplot(pcas ,aes(x = PC3,y=PC4,col=Diagnosis )) + geom_point(size=2.25,alpha=0.4) 
centroids <- aggregate(cbind(PC3,PC4)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_34.png',sep=''),width = sX,height = sY)
# PC 4-5
g <- ggplot(pcas ,aes(x = PC4,y=PC5,col=Diagnosis )) + geom_point(size=2.25,alpha=0.4) 
centroids <- aggregate(cbind(PC4,PC5)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_45.png',sep=''),width = sX,height = sY)

# ===================================
#  refactor diagnosis
#=================================
allDataIBSIBDHC <- allDataR
allDataIBSIBDHC$Diagnosis <- as.character(allDataIBSIBDHC$Diagnosis)
allDataIBSIBDHC$Diagnosis[allDataIBSIBDHC$Diagnosis=="UC"] <- "IBD"
allDataIBSIBDHC$Diagnosis[allDataIBSIBDHC$Diagnosis=="CD"] <- "IBD"
allDataIBSIBDHC$Diagnosis[allDataIBSIBDHC$Diagnosis=="IBDU"] <- "IBD"
# get rid of self-reported IBS
allDataIBSIBDHC <- allDataIBSIBDHC[!allDataIBSIBDHC$Diagnosis=="IBS-SR",]
allDataIBSIBDHC$Diagnosis <- as.factor(as.character(allDataIBSIBDHC$Diagnosis))

# - correlate biomarkers
g <- ggplot(allDataIBSIBDHC,aes(x=Calprot,y=HBD2,col=Diagnosis)) + geom_point()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_cor_calprot_vs_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(allDataIBSIBDHC,aes(x=ChrA,y=HBD2,col=Diagnosis)) + geom_point()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_cor_ChrA_vs_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(allDataIBSIBDHC,aes(x=Calprot,y=ChrA,col=Diagnosis)) + geom_point()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_cor_ChrA_vs_calprot.png',sep=''),width = sX,height = sY)
# - replot with refactored diagnosis
g <- ggplot(data=allDataIBSIBDHC,aes(x=Diagnosis,y=HBD2,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bar_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataIBSIBDHC,aes(x=Diagnosis,y=ChrA,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bar_ChrA.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataIBSIBDHC,aes(x=Diagnosis,y=Calprot,col=Diagnosis)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter()
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bar_Calprot.png',sep=''),width = sX,height = sY)
#       - pretty histograms
g <- ggplot(data=allDataIBSIBDHC,aes(x=HBD2,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_hist_HBD2.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataIBSIBDHC,aes(x=Calprot,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_hist_Calprot.png',sep=''),width = sX,height = sY)
g <- ggplot(data=allDataIBSIBDHC,aes(x=ChrA,fill=Diagnosis)) + geom_histogram() + facet_wrap( ~ Diagnosis)
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_hist_ChrA.png',sep=''),width = sX,height = sY)
# - formal testing
g <- testOneFeature(dataIn=allDataIBSIBDHC,feature = 'HBD2', saveFolder = paste(outFolderName,'/exploratory/diag3_.png',sep=''),doSave = F,onlyShowSig = T,
               responseVar = "Diagnosis",nrTests = 3,yLab = "HBD2 (log transformed)",retPlot = T)[[2]]
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_fortest_HBD2.png',sep=''),width = sX,height = sY)
g <- testOneFeature(dataIn=allDataIBSIBDHC,feature = 'ChrA', saveFolder = paste(outFolderName,'/exploratory/diag3_.png',sep=''),doSave = F,onlyShowSig = T,
                    responseVar = "Diagnosis",nrTests = 3,yLab = "ChrA (log transformed)",retPlot = T)[[2]]
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_fortest_ChrA.png',sep=''),width = sX,height = sY)
g <- testOneFeature(dataIn=allDataIBSIBDHC,feature = 'Calprot', saveFolder = paste(outFolderName,'/exploratory/diag3_.png',sep=''),doSave = F,onlyShowSig = T,
                    responseVar = "Diagnosis",nrTests = 3,yLab = "Calprot (log transformed)",retPlot = T)[[2]]
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_fortest_Calprot.png',sep=''),width = sX,height = sY)

# prep Bray-curtis
bcReady <- allDataIBSIBDHC[,grep('s__',colnames(allDataIBSIBDHC))]
bcurtis <- vegdist(bcReady,method = "bray")
r.pca <- prcomp(bcurtis, center = F,scale. = F)
pcas <- as.data.frame(r.pca$rotation[,1:5])
pcas$Diagnosis <- allDataIBSIBDHC$Diagnosis
# PC 1-2
g <- ggplot(pcas ,aes(x = PC1,y=PC2,col=Diagnosis )) + geom_point(size=1.25,alpha=0.6) + scale_color_manual(values=c("#009933", "#CC0000", "#FF9933"))
centroids <- aggregate(cbind(PC1,PC2)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bc_12.png',sep=''),width = sX,height = sY)
# PC 2-3
g <- ggplot(pcas ,aes(x = PC2,y=PC3,col=Diagnosis )) + geom_point(size=1.75,alpha=0.6) + scale_color_manual(values=c("#009933", "#CC0000", "#FF9933"))
centroids <- aggregate(cbind(PC2,PC3)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bc_23.png',sep=''),width = sX,height = sY)
# PC 3-4
g <- ggplot(pcas ,aes(x = PC3,y=PC4,col=Diagnosis )) + geom_point(size=1.75,alpha=0.6) + scale_color_manual(values=c("#009933", "#CC0000", "#FF9933"))
centroids <- aggregate(cbind(PC3,PC4)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bc_34.png',sep=''),width = sX,height = sY)
# PC 4-5
g <- ggplot(pcas ,aes(x = PC4,y=PC5,col=Diagnosis )) + geom_point(size=1.75,alpha=0.6) + scale_color_manual(values=c("#009933", "#CC0000", "#FF9933"))
centroids <- aggregate(cbind(PC4,PC5)~Diagnosis,pcas,mean)
g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5))
print(g)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bc_45.png',sep=''),width = sX,height = sY)

# prep shannon
shannonRdy <- allDataIBSIBDHC[,c(grep('Diagnosis',colnames(allDataIBSIBDHC)),grep('s__',colnames(allDataIBSIBDHC))) ]
shannonRdy$div.shannon <- diversity(shannonRdy[,-grep('Diagnosis',colnames(shannonRdy))],MARGIN = 1, index="shannon")
g <- ggplot(shannonRdy,aes(x=Diagnosis,col=Diagnosis,y=div.shannon)) + geom_boxplot(outlier.alpha = 0.0) + geom_jitter(width = 0.25)
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_bar_Shannon.png',sep=''),width = sX,height = sY)
g <- testOneFeature(dataIn=shannonRdy,feature = 'div.shannon', saveFolder = paste(outFolderName,'/microbiome/shannon.png',sep=''),doSave = F,onlyShowSig = T,
               responseVar = "Diagnosis",nrTests = ncol(shannonRdy),yLab = "Shannon Diversity Index",retPlot = T)[[2]]
ggsave(plot=g,paste(outFolderName,'/exploratory/diag3_fortest_Shannon.png',sep=''),width = sX,height = sY)

# ========================
# PREP DATA FOR MODELLING
# ========================
allDataR <- allData #[allData$Cohort=="IBD" | allData$Cohort=="MIBS",]
# cleanup
#allDataR$Cohort <- NULL
# prep diagnosis (IBS & IBD, remove subtypes)
allDataR$Diagnosis <- as.character(allDataR$Diagnosis)
# only IBS vs IBD
#allDataR <- allDataR[!allDataR$Diagnosis=="HC",]
# correct for IBS
#allDataR <- allDataR[!allDataR$Diagnosis=="IBS-SR",]
allDataR$Diagnosis[grep("IBS-SR",allDataR$Diagnosis)] <- "TEMPTEMP"
allDataR$Diagnosis[grep("IBS",allDataR$Diagnosis)] <- "IBS"
allDataR$Diagnosis[grep("TEMPTEMP",allDataR$Diagnosis)] <- "IBS-SR"

allDataR$Diagnosis <- as.factor(as.character(allDataR$Diagnosis))
allDataPrep <- allDataR
# handle wrong PF_reads and PFReads
for (r in c(1:nrow(allDataPrep))) {
  if (!is.na(allDataPrep$PF_reads[r])) {
    allDataPrep$PFReads[r] <- allDataPrep$PF_reads[r]
  }
}
allDataPrep$PF_reads <- NULL

# impute BMI and PFReds to the median
allDataPrep$BMI[is.na(allDataPrep$BMI)] <- median(allDataPrep$BMI,na.rm = T)
allDataPrep$PFReads[is.na(allDataPrep$PFReads)] <- median(allDataPrep$BMI,na.rm = T)

# save 'ready data'
write.table(allDataPrep,paste0(outFolderData,'/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv'),sep=',',row.names = F)

allDataPrep <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv'),sep=',',header = T)

arcsinTransformAtStart = T # should be on
# Normalise & Scale data ('main data')
if (arcsinTransformAtStart) {allDataPrep <- asinSqrtNormalise(allDataPrep)}

corN <- c("Gender","Age","BMI","Laxatives","PPI","Antibiotics","PFReads")
# pre-normalization
write.table(allDataPrep,paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_uncorr_v10.csv'),sep=',',row.names = F)
# do linear correction for phenotypes
#iData = allDataPrep;corrNames = corN;correctZeros = F;corrMG = T;corrPWY = T
allDataPrepC <- linearCorrectMGPwy(iData = allDataPrep,corrNames = corN,correctZeros = F,corrMG = T,corrPWY = T)
# save normalized & corrected data
write.table(allDataPrepC,paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv'),sep=',',row.names = F)

# === test stuff ===
# ====================
# raw [OK]
tst1 <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_rdyfornormalisation.csv'),sep=',',header=T)
rowSums(tst1[,grep('k__',colnames(tst1))][1:5,])
# normalized [OK]
#tst2 <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_uncorr_v10.csv'),sep=',',header=T)
rowSums(tst2[,grep('s__',colnames(tst2))][1:5,])
# corrected [OK]
tst3 <- read.table(paste0(outFolderData,'/LLDIBD_Taxa_PWYs_data_nor_corr_v10.csv'),sep=',',header=T)
#rowSums(tst3[,grep('s__',colnames(tst3))][1:5,])

# ggplot(data=tst1,
#        aes(x=allDataPrep$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii,
#                      y=allDataPrepC$k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii)) + 
#          geom_point() + geom_smooth()
