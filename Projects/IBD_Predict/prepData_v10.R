# IBD prediction data preparation
# =========================================
#
# -> also includes pathway data!
library(plyr)

# version 9: now with final variant of Marjolein's data
# => also includes SCSI and HBI
# => also includes resections

# wd:
setwd('D:/UMCG/ML_IBD_v5/')
source('../myLibs_v3/R_Microbiome_scripts.R')

# prep MIBS
MIBSmeta <- read.csv('data_raw/MIBS_Phenotypes_v2.csv',stringsAsFactors = F) 
MIBSmeta <- MIBSmeta[!is.na(MIBSmeta$Diagnosis),]
MIBSmeta$Cohort <- "MIBS"
MIBSmetaphlan <- as.data.frame(t(read.csv('data_raw/MIBS_taxonomy_metaphlan2_092017.tsv',sep='\t',comment.char = '#',row.names = 1,stringsAsFactors = F)))
MIBSmetaphlan$ID <- rownames(MIBSmetaphlan)
MIBSmetaphlan$ID <- gsub('_metaphlan','',MIBSmetaphlan$ID)
MIBS <- merge(MIBSmeta,MIBSmetaphlan,by = "ID")
# some MIBS are in fact in LLD metaphlan file, add those as well
LLDmetaphlan <- as.data.frame(t(read.csv('data_raw/LLD_taxonomy_metaphlan2_092017.txt',sep='\t',comment.char = '#',row.names = 1,stringsAsFactors = F)))
LLDmetaphlan$ID <- rownames(LLDmetaphlan)
LLDmetaphlan$ID <- gsub('_metaphlan','',LLDmetaphlan$ID)
for (i in c(1:nrow(LLDmetaphlan))) {
  mID <- LLDmetaphlan$ID[i]
  if (substr(mID[1],1,1) == "X") {
    mID <- substring(mID,2)
  }
  LLDmetaphlan$ID[i] <- mID
}
MIBS2 <- merge(MIBSmeta,LLDmetaphlan,by="ID")
# merge MIBS & MIBS2
MIBS <- rbind.fill(MIBS,MIBS2)
MIBS2 <- NULL
# find missing IDs

# prep HMP2
HMP2meta <- read.csv('data_raw/hmp2_w1818/hmp2_metadata_w1818_filtered_v2.csv',stringsAsFactors = F) 
HMP2meta <- HMP2meta[!is.na(HMP2meta$Diagnosis),]
HMP2meta$Cohort <- "HMP2"
HMP2metaphlan <- as.data.frame(t(read.csv('data_raw/hmp2_w1818/hmp2_taxonomic_profiles_w1818.tsv',sep='\t',comment.char = "",row.names = 1,stringsAsFactors = F)))
HMP2metaphlan$ID <- rownames(HMP2metaphlan)
HMP2metaphlan$ID <- gsub('_taxonomic_profile','',HMP2metaphlan$ID)
HMP2 <- merge(HMP2meta,HMP2metaphlan,by = "ID")
#HMP2$SCCAI <- NULL
#HMP2$HBI <- NULL

# =======================================
# prep LLD/IBD
# =======================================
#  load IBD metaphlan (this is IBD only!)
# =============================================
IBDmetaphlan <- as.data.frame(t(read.csv('data_raw/IBD_taxonomy_metaphlan2_082017.txt',sep='\t',comment.char = '#',row.names = 1,stringsAsFactors = F)))
IBDmetaphlan$ID <- rownames(IBDmetaphlan)
IBDmetaphlan$ID <- gsub('_metaphlan','',IBDmetaphlan$ID)

# load IBD/LLD metadata $ 350 - 400 - 400
LLDIBDmeta <- read.csv('data_raw/IBD_LLD_metadata_filtered_wFEC_v3.csv',stringsAsFactors = F)

LLDIBDmeta$ID <- NA
LLDIBDmeta$Cohort <- "LLD_IBD"
for (i in c(1:nrow(LLDIBDmeta))) {
  mID <- LLDIBDmeta$IndividualWGSsampleID[i]
  if (!is.na(mID) & substr(mID,1,1) == "'") {
    mID <- substring(mID,2)
  }
  LLDIBDmeta$IndividualWGSsampleID[i] <- mID
}
# add self reported IBS
srIBS <- read.csv('data_raw/LLD_SR_IBS_IDs.csv',stringsAsFactors = F)
for (i in srIBS$ID) {
  gr <- grep(i,LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    LLDIBDmeta$Diagnosis[gr] <- "IBS-SR"
  }
  else {
    print (paste(i,'not found!'))
  }
}
# add extra calprotectin measurements
LLDcalProt <- read.csv("data_raw/calpro_lld.csv",stringsAsFactors = F)
for (i in c(1:nrow(LLDcalProt))) {
  gr <- grep(LLDcalProt$ID[i],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    LLDIBDmeta$Calprot[gr] <-  LLDcalProt$FecalCal40[i]
  }
  else {
    print (paste(i,'not found!'))
  }
}

# add extra ID links for new samples
LLDIBDmetaIDLinks <- read.csv('data_raw/LLD_IBD_idLinks.csv',stringsAsFactors = F)
for (col in c(1:9)) {
    LLDIBDmeta[,col] <- as.character(LLDIBDmeta[,col])
}
print ('linking new LLD/IBD IDs')
fndC <- 0
for (i in c(1:nrow(LLDIBDmetaIDLinks))) {
  #print(LLDIBDmetaIDLinks$Collaborator.Participant.ID[i])
  gr <- grep(LLDIBDmetaIDLinks$Collaborator.Participant.ID[i],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    fndC = fndC + 1
    #print (LLDIBDmetaIDLinks$Collaborator.Participant.ID[gr])
    LLDIBDmeta$WGSProjectID[gr] <- LLDIBDmetaIDLinks$G.[i]
  }
}
print (paste('Linked',fndC,'IDs; total = ',nrow(LLDIBDmetaIDLinks) ))
#write.csv(LLDIBDmeta,'data_raw/IBD_LLD_metadata_filtered.csv')

# handle 'problematic' IDs
print ('linking "problematic" IDs')
LLDprobIDs <- read.csv('data_raw/LLD_problematic_ids.txt',stringsAsFactors = F,sep='\t')
cF = 0
for (i in c(1:nrow(LLDprobIDs))) {
  gr <- grep(LLDprobIDs$ID.1[i],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    cF = cF + 1
    LLDIBDmeta$WGSProjectID[gr] <- LLDprobIDs$ID[i]
  }
  else {
    print(LLDprobIDs[i,])
  }
}
print (paste('linked ',cF,'IDs'))
#IndividualWGSsampleID

print (' >> connecting stoma/pouch phenotypes (this take a while - bad code)')
# add stoma/pouch to LLD/IBD data
stomaPouch <- read.csv('data_raw/Metadata_Pouch_Stoma.csv',stringsAsFactors = F)
# kubj ri ptohinla ofd
idColsSP <- c(1,2,3,4,6)
idColsLLDIBDmeta <- c(1:9)
LLDIBDmeta$StomaPouchType <- NA
LLDIBDmeta <- LLDIBDmeta[!is.na(LLDIBDmeta$Diagnosis),]
fnd = c()
nfnd <- c()
for (rSP in c(1:nrow(stomaPouch))) {
  print (rSP)
  g <- grepl(stomaPouch$UMCGIBDResearchIDorLLDeepID[rSP],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (sum(g) > 0) {
    LLDIBDmeta$StomaPouchType[g] <- stomaPouch$CurrentStomaOrPouchType[rSP]
    fnd <- c(fnd,rSP)
  } else {
    nfnd <- c(nfnd,rSP)
  }
}
LLDIBDmeta$StomaPouchType[LLDIBDmeta$StomaPouchType=="none"] <- NA
print (paste('stoma/pounch type: found',length(unique(fnd))))

# -- add SSCI and HBI --
sschbi <- read.table('data_raw/IBD_HBI_SSCAI.csv',sep=',',header = T)
idColsLLDIBDmeta <- c(1:9)
LLDIBDmeta$HBI <- NA
LLDIBDmeta$HBI <- as.character(LLDIBDmeta$HBI)
LLDIBDmeta$SCCAI <- NA
LLDIBDmeta$SCCAI <- as.character(LLDIBDmeta$SCCAI)
#LLDIBDmeta$IBD.active <- NA
#LLDIBDmeta$IBD.active <- as.character(LLDIBDmeta$IBD.active)
fnd = c()
nfnd = c()
sschbi$ResearchID <- as.character(sschbi$ResearchID)
for (rSP in c(1:nrow(sschbi))) {
  gr <- grepl(sschbi$ResearchID[rSP],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (sum(gr) > 0) {
    fnd <- c(fnd,sschbi$ResearchID[rSP])
    LLDIBDmeta$HBI[gr] <- sschbi$HBI[rSP]
    LLDIBDmeta$SCCAI[gr] <- sschbi$SSCAI[rSP]
    #LLDIBDmeta$IBD.active[gr] <- sschbi$IBD.active[rSP]
  } else {
    nfnd <- c(nfnd,sschbi$ResearchID[rSP])
  }
}

print ('connecting IBD metadata to metagenome IDs')
#  - basically first 9 rows are IDs
LLDIBDmeta$ID <- NA
IBDmetaphlan$ID <- as.character(IBDmetaphlan$ID)
fndC = 0
nFnd = c()
for (i in c(1:nrow(IBDmetaphlan))) {
  #print (IBDmetaphlan$ID[i])
  fnd = F
  for (col in c(1:9)) {
    if (!fnd) {
    #LLDIBDmeta[,col] <- as.character(LLDIBDmeta[,col])
      g <- grep(IBDmetaphlan$ID[i],LLDIBDmeta[,col])
      if (length(g) > 0 & !fnd)
      {
        fID <- LLDIBDmeta[g,col][1]
        #print (paste('found',IBDmetaphlan$ID[i],fID))
        fnd = T
        fndC = fndC + 1
        LLDIBDmeta$ID[g] <- fID
      }
    }
  }
  if (!fnd) {
    nFnd <- c(nFnd,IBDmetaphlan$ID[i])
  }
}
print(paste('IDs total',nrow(IBDmetaphlan),'; found:',fndC,'; missed:',nrow(IBDmetaphlan)-fndC))

# now link to faecal calprotectin & other biomarkers (new measurements, 02/05/2018 !!
# ===================================================================================
newHBDCHR <- read.table('data_raw/IBD_HBD_CHRA.csv',header=T,stringsAsFactors = F,sep=',')
# >>>>>>>>>>> CAL <<<<<<<<<<<<<<<
newCal <- newHBDCHR[newHBDCHR$PARAMETER_NAME=="calprotectine (faeces)",]
newCal$PATIENT_ID <- NULL; newCal$PARAMETER_NAME <- NULL; newCal$UNIT <- NULL
colnames(newCal) <- c("RID","Calprot.ug.g")
newCal$Calprot.ug.g[newCal$Calprot.ug.g=="<19.5"] <- 19.5
newCal$Calprot.ug.g[newCal$Calprot.ug.g==">800"] <- 800.0
newCal$Calprot.ug.g <- as.numeric(as.character(newCal$Calprot.ug.g))

#newCal$RID %in% LLDIBDmeta$UMCGIBDResearchIDorLLDeepID
# put calprotectin where it belongs
fnd <- 0
for (c in c(1:nrow(newCal))) {
  gr <- grep(newCal$RID[c],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (gr > 0) {
    fnd = fnd + length(gr)
    LLDIBDmeta$Calprot[gr] <- newCal$Calprot.ug.g[c]
  }
}
print (paste(' -> linked',fnd,' Calprotectins! out of ',nrow(newCal)))

# >>>>>>>>>>>>>>>> CHRA <<<<<<<<<<<<<<<<<<<
newChrA <- newHBDCHR[newHBDCHR$PARAMETER_NAME=="Chromogranine-A (faeces)",]
newChrA$PATIENT_ID <- NULL; newChrA$PARAMETER_NAME <- NULL; newChrA$UNIT <- NULL
colnames(newChrA) <- c("RID","ChrA")
newChrA$ChrA <- as.numeric(as.character(newChrA$ChrA))
#newChrA$RID %in% LLDIBDmeta$UMCGIBDResearchIDorLLDeepID
LLDIBDmeta$ChrA <- NA
fnd <- 0
for (c in c(1:nrow(newChrA))) {
  gr <- grep(newChrA$RID[c],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    fnd = fnd + length(gr)
    LLDIBDmeta$ChrA[gr] <- newChrA$ChrA[c]
  }
}
print (paste(' -> linked',fnd,' ChrA! out of ',nrow(newChrA)))

# >>>>>>>>>>>>>>>> HBD <<<<<<<<<<<<<<<<<<<
newHBD<- newHBDCHR[newHBDCHR$PARAMETER_NAME=="beta-defensin2",]
newHBD$PATIENT_ID <- NULL; newHBD$PARAMETER_NAME <- NULL; newHBD$UNIT <- NULL
colnames(newHBD) <- c("RID","HBD2")
# fix non-numeric
newHBD$HBD2[newHBD$HBD2==">750"] <- "750.0"
newHBD$HBD2[newHBD$HBD2==">800"] <- "800.0"
#
newHBD$HBD2 <- as.numeric(as.character(newHBD$HBD2))
#newHBD$RID %in% LLDIBDmeta$UMCGIBDResearchIDorLLDeepID
LLDIBDmeta$HBD2 <- NA
fnd <- 0
for (c in c(1:nrow(newHBD))) {
  gr <- grep(newHBD$RID[c],LLDIBDmeta$UMCGIBDResearchIDorLLDeepID)
  if (length(gr) > 0) {
    fnd = fnd + length(gr)
    LLDIBDmeta$HBD2[gr] <- newHBD$HBD2[c]
  }
}
print (paste(' -> linked',fnd,' HBDs! out of ',nrow(newHBD)))

# now link to faecal calprotectin & other biomarkers (new measurements for LLD, 08/05/2018 !!
# ===================================================================================
newBM <- read.table('data_raw/BioMarkers_valerie.csv',header=T,stringsAsFactors = F,sep=',')
newBM$FecalCalprotectin[newBM$FecalCalprotectin=="under40"] <- 40.0
newBM$FecalCalprotectin[newBM$FecalCalprotectin=="<40"] <- 40.0
newBM$FecalCalprotectin <- as.numeric(as.character(newBM$FecalCalprotectin))
newBM$FecalCalprotectin[newBM$FecalCalprotectin < 40.0] <- 40.0

print (' --> before Valerie data:')
print(paste(' CALPROT: has data:',sum(!is.na(LLDIBDmeta$Calprot)),'; NA:',sum(is.na(LLDIBDmeta$Calprot))  ))
print(paste(' HBD2: has data:',sum(!is.na(LLDIBDmeta$HBD2)),'; NA:',sum(is.na(LLDIBDmeta$HBD2))  ))
print(paste(' CHRA: has data:',sum(!is.na(LLDIBDmeta$ChrA)),'; NA:',sum(is.na(LLDIBDmeta$ChrA))  ))

fnd = rep(F,nrow(newBM))
for (c in c(1:nrow(newBM))) {
  # iterate over all ID rows
  newBMid <- newBM$Ids[c]
  for (idr in c(1:9,26)) {
    gr <- grep(newBMid,LLDIBDmeta[[idr]])
    if (length(gr) > 0) {
      if (!is.na(newBM$Chromogranin_A[c])) {LLDIBDmeta$ChrA[gr] <- newBM$Chromogranin_A[c]}
      if (!is.na(newBM$FecalCalprotectin[c])) {LLDIBDmeta$Calprot[gr] <- newBM$FecalCalprotectin[c]}
        if (!is.na(newBM$HBD_2[c])) {LLDIBDmeta$HBD2[gr] <- newBM$HBD_2[c]}
      fnd[c] <- T
    }
  }
}
print ('  --> after valerie data')
print(paste(' CALPROT: has data:',sum(!is.na(LLDIBDmeta$Calprot)),'; NA:',sum(is.na(LLDIBDmeta$Calprot))  ))
print(paste(' HBD2: has data:',sum(!is.na(LLDIBDmeta$HBD2)),'; NA:',sum(is.na(LLDIBDmeta$HBD2))  ))
print(paste(' CHRA: has data:',sum(!is.na(LLDIBDmeta$ChrA)),'; NA:',sum(is.na(LLDIBDmeta$ChrA))  ))
# test for debug
#newBM[!fnd & !grepl('IBS',newBM$Ids) & !grepl('MICO',newBM$Ids),]
#write.table(LLDIBDmeta,'mici.csv',sep=',')

# ===================================
# ADD MARJOLEIN DATA (flares, final)
# ===================================
flareData <- read.table('data_raw/IBD_diseaseactivity_final_marjolein.txt',sep='\t',header = T)
flareData$UMCGNoFromZIC <- as.character(flareData$UMCGNoFromZIC)
flareData$UMCGIBDDNAID <- as.character(flareData$UMCGIBDDNAID)
#table(flareData$Disease_activity_Categorical)
flareData$Disease_activity_Categorical <- as.character(flareData$Disease_activity_Categorical)
flareData$Disease_activity_Categorical[flareData$Disease_activity_Categorical=="after a flare"] <- "post.flare"
flareData$Disease_activity_Categorical[flareData$Disease_activity_Categorical=="before a flare"] <- "pre.flare"
flareData$Disease_activity_Categorical[flareData$Disease_activity_Categorical=="during a flare"] <- "in.flare"

flareData$Days_Since_Last_Flare[flareData$Days_Since_Last_Flare=="None"] <- NA
flareData$IBD.time.from.flare <- as.numeric(as.character(flareData$Days_Since_Last_Flare))
flareData$Days_Since_Last_Flare <- NULL
flareData$IBD.time.to.flare <- as.numeric(flareData$Days_Until_Next_Flare)
flareData$Days_Until_Next_Flare <- NULL

flareData$IBD.closest.flare <- pmin(abs(flareData$IBD.time.from.flare),abs(flareData$IBD.time.to.flare),na.rm = T)

#flareData <- flareData[!is.na(flareData$ResearchID),]

LLDIBDmeta$IBD.time.to.flare <- NA
LLDIBDmeta$IBD.time.from.flare <- NA
LLDIBDmeta$IBD.closest.flare <- NA
LLDIBDmeta$IBD.flare.clinical <- NA

fnd = rep(F,nrow(flareData))
for (c in c(1:nrow(flareData))) {
  # iterate over all ID rows
  flareId <- flareData$UMCGNoFromZIC[c]
  flareId2 <- flareData$UMCGIBDDNAID[c]
  for (idr in c(1:9,26)) {
    gr <- grep(flareId,LLDIBDmeta[[idr]])
    if (length(gr) > 0) {
      #print (gr)
      LLDIBDmeta$IBD.time.from.flare[gr] <- flareData$IBD.time.from.flare[c]
      LLDIBDmeta$IBD.time.to.flare[gr]   <- flareData$IBD.time.to.flare[c]
      LLDIBDmeta$IBD.closest.flare[gr]   <- flareData$IBD.closest.flare[c]
      LLDIBDmeta$IBD.flare.clinical[gr] <- flareData$Disease_activity_Categorical[c]
      fnd[c] <- T
    } else {
      gr <- grep(flareId2,LLDIBDmeta[[idr]])
      LLDIBDmeta$IBD.time.from.flare[gr] <- flareData$IBD.time.from.flare[c]
      LLDIBDmeta$IBD.time.to.flare[gr]   <- flareData$IBD.time.to.flare[c]
      LLDIBDmeta$IBD.closest.flare[gr]   <- flareData$IBD.closest.flare[c]
      LLDIBDmeta$IBD.flare.clinical[gr] <- flareData$Disease_activity_Categorical[c]
      fnd[c] <- T
    }
  }
}
print (paste('linked',sum(fnd),'flare IDs / ',nrow(flareData),' to IBD ids'))

#TODO ================ ADD MAP =================
#TODO
#mapData <- read.table(file = 'data_raw/MAP_MGS_ids.csv',sep=',',header=T,quote = '"',na.strings = 'NA',dec = '.')
#TODO ==========================================
#mapData$ID <- NA
#for (c in c(1:nrow(LLDIBDmeta))) {
#  print(c)
#  for (cc in c(1:nrow(mapData))) {
#    if (!is.na(LLDIBDmeta$UMCGNoFromZIC[[c]]) & !is.na(mapData$SampleID[cc]) & (LLDIBDmeta$UMCGNoFromZIC[[c]] == mapData$SampleID[cc] )) {
#      mapData$ID[cc] <- LLDIBDmeta$ID[c]
#    }
#  }
#}
#mapData$SampleID <- NULL
#mapData$SAMPLE <- NULL
#mapData <- mapData[!is.na(mapData$ID),]
# merge with MAP
#LLDIBDmeta <- merge(LLDIBDmeta,mapData,by="ID",all = T)

#TODO ========= END OF ADD MAP ===================

# merge and cleanup
IBD <- merge(LLDIBDmeta,IBDmetaphlan,by="ID")

IBD[,c(2:10)] <- NULL

# LLD metaphlan
# =============================
print ('linking LLD metaphlan to LLD/IBD metadata')
LLDmetaphlan <- as.data.frame(t(read.csv('data_raw/LLD_taxonomy_metaphlan2_092017.txt',sep='\t',comment.char = '#',row.names = 1,stringsAsFactors = F)))
LLDmetaphlan$ID <- rownames(LLDmetaphlan)
LLDmetaphlan$ID <- gsub('_metaphlan','',LLDmetaphlan$ID)
#  - basically first 9 rows are IDs
LLDmetaphlan$ID <- as.character(LLDmetaphlan$ID)
fndC = 0
nFnd = c()
for (i in c(1:nrow(LLDmetaphlan))) {
  #print (LLDmetaphlan$ID[i])
  fnd = F
  for (col in c(1:9)) {
    if (!fnd) {
      #LLDIBDmeta[,col] <- as.character(LLDIBDmeta[,col])
      mID <- LLDmetaphlan$ID[i]
      if (substr(mID[1],1,1) == "X") {
        mID <- substring(mID,2)
      }
      LLDmetaphlan$ID[i] <- mID
      g <- grep(LLDmetaphlan$ID[i],LLDIBDmeta[,col])
      if (length(g) > 0 & !fnd)
      {
        fID <- LLDIBDmeta[g,col][1]
        #print (paste('found',LLDmetaphlan$ID[i],fID))
        fnd = T
        fndC = fndC + 1
        LLDIBDmeta$ID[g] <- fID
      }
    }
  }
  if (!fnd) {
    nFnd <- c(nFnd,LLDmetaphlan$ID[i])
  }
}
print(paste('IDs total',nrow(LLDmetaphlan),'; found:',fndC,'; missed:',nrow(LLDmetaphlan)-fndC))
# merge and cleanup
LLD <- merge(LLDIBDmeta,LLDmetaphlan,by="ID")
#write.csv(LLDIBDmeta,'testor.csv')
LLD[,c(2:10)] <- NULL

# final join
allData <- rbind.fill(LLD,MIBS,IBD,HMP2)
# fix calprot, introduce calprot40
allData$Calprot[allData$Calprot=="under40"] <- "40.0"
allData$Calprot[allData$Calprot=="<40"] <- "40.0"
allData$Calprot[allData$Calprot=="#VALUE!"] <- NA
allData$Calprot <- as.numeric(allData$Calprot)

allData$CalprotD40 <- allData$Calprot
allData$CalprotD40[allData$CalprotD40=="under40"] <- "40"
allData$CalprotD40[allData$CalprotD40=="<40"] <- "40"
allData$CalprotD40 <- as.numeric(allData$CalprotD40)
allData$CalprotD40 <- ceiling(allData$CalprotD40/40.0)
allDataNonMG <- allData[,grep('k__',colnames(allData),invert=T)]

#write.csv(allDataNonMG,'testor.csv')
write.csv(allData,'data_clean/allDataMergedClean_v10.csv',row.names = F)

# ========================================
# PATHWAYS 
# ========================================
# clean data & prep for merging
prepPathways <- function(inFile) {
  humann <- read.csv(inFile,stringsAsFactors = F,sep="\t")
  humannPur <- humann[grep('\\|',humann$X..Pathway,invert = T),]
  colnames(humannPur) <- gsub('_kneaddata_merged_Abundance','',colnames(humannPur))
  colnames(humannPur) <- gsub('_Abundance','',colnames(humannPur))
  rownames(humannPur) <- humannPur$X..Pathway
  humannPur$X..Pathway <- NULL
  humannPur <- as.data.frame(t(humannPur))
  humannPur$UNINTEGRATED <- NULL
  humannPur$UNMAPPED <- NULL
  for (c in grep('PWY',colnames(humannPur),invert = T)) {
    colnames(humannPur)[c] <- paste("PWY:",colnames(humannPur)[c])
  }
  rn <- c()
  for (r in rownames(humannPur)) {
    if (substr(r,1,1) == "X") {
      r <- substring(r,2)
    }
    rn <- c(rn,r)
  }
  rownames(humannPur) <- rn
  humannPur
}

MIBSpath <- prepPathways('data_raw/MIBS_humann2_uniref90_092017.tsv')
LLDpath <- prepPathways('data_raw/LLD_humann2_Uniref90_092017.txt') # X entered names at start!
IBDpath <- prepPathways('data_raw/IBD_humann2_pathways_uniref90_082017.txt') # X entered names at start!
HMP2path <- prepPathways('data_raw/hmp2_w1818/hmp2_pathabundance_w1818.tsv') # X entered names at start!
# merge with metadata
MIBSpath$ID <- rownames(MIBSpath)
MIBSall <- merge(allData,MIBSpath,by="ID")
LLDpath$ID <- rownames(LLDpath)
LLDall <- merge(allData,LLDpath,by="ID")
IBDpath$ID <- rownames(IBDpath)
IBDall <- merge(allData,IBDpath,by="ID")
HMP2path$ID <- rownames(HMP2path)
HMP2all <- merge(allData,HMP2path,by="ID")

hmp2indID <- read.table('data_raw/hmp2_w1818/hmp2_metadata_w1818_ID_participants_key.csv',sep=',',header = T)
hmp2indID$ProjectID <- NULL
HMP2all_pID <- merge(HMP2all,hmp2indID,by="ID",all = T)

allDataAll <- rbind.fill(LLDall,MIBSall,IBDall,HMP2all_pID)
allDataAllRescaled <- filterHumannDF(inDF = allDataAll,rescale = T,verbose = T,presPerc = 0.0,minMRelAb = 0.0,minMedRelAb = 0.0)
# fixes 
allDataAllRescaled$Participant.ID <- NULL
# calculate activity of IBD
allDataAllRescaled$SCCAI <- as.numeric(as.character(allDataAllRescaled$SCCAI))
allDataAllRescaled$HBI <- as.numeric(as.character(allDataAllRescaled$HBI))
allDataAllRescaled$IBD.active.HBISCCAI <- "N"
allDataAllRescaled$IBD.active.HBISCCAI[allDataAllRescaled$HBI >= 5 | allDataAllRescaled$SCCAI >= 3] <- "Y"

write.csv(allDataAllRescaled,'data_clean/withBMDataMergedCleanRSPathways_v10.csv',row.names = F)

allDataAllPheno <- allDataAllRescaled[,-grep('__',colnames(allDataAllRescaled))]
allDataAllPheno <- allDataAllPheno[,-grep('PWY',colnames(allDataAllPheno))]