# DMP heritability study-wide FDR calculation
# ============================================
# > 1) collects data from all jobs for heritability analysis
# > 2) loads permutation results
# > 3) loads results from heritability analysis and
#   constructs p-values for each taxon / pwy 
#   by using all perm. results for taxa / pathways
# ============================================
library(plyr)

# 1/a) collect results for individual taxa / pwy runs
# ===================================================
# NOTE: set the path to appropriate location
# example: inF <- "D:/Vbox/shared/dag/git_14_05/DMP/heritability_analysis_v2/"
inF <- "."
setwd(inF)

inFolder <- paste0(inF,'/mock_data_heritability_results/')
print(paste0('> collecting taxonomy jobs from ',inFolder))
# collect taxa jobs
inFiles = list.files(inFolder,pattern="*.csv")
res <- NULL
for (i in inFiles) {
  if (grepl('taxa',i)) {
    print(paste0(' > collecting tax job ',i))
    inDF <- read.table(paste0(inFolder,'/',i),sep=',',header=T,stringsAsFactors=F)
    res <- rbind.fill(res,inDF)
  }
}
if ("Pv_ID" %in% colnames(res)) {
  res <- res[order(res$Pv_ID),]
}
write.table(res,paste0('results_mockdata_taxa.csv'),sep=',',row.names=F)

# 1/b) collect pwy jobs
# ====================================
# NOTE: only collect bacterial pathways (drop eukaryia, plants...)
inBacPath <- read.table('bacpaths.txt',header=F,stringsAsFactors = F)

print(paste0('> collecting pathway jobs from ',inFolder))
inFiles = list.files(inFolder,pattern="*.csv")
res <- NULL
for (i in inFiles) {
  if (grepl('pwy',i)) {
    print(paste0(' > collecting pwy job ',i))
    inDF <- read.table(paste0(inFolder,'/',i),sep=',',header=T,stringsAsFactors=F)
    #print(inDF$Trait)
    if (sum(grepl(inDF$Trait,inBacPath$V1) > 0)) {
      res <- rbind.fill(res,inDF)
    } else {
      print (paste0(' >> ',inDF$Trait,' is not bacterial pathway, skipping!'))
    }
  }
}
if ("Pv_ID" %in% colnames(res)) {
  res <- res[order(res$Pv_ID),]
}
write.table(res,paste0('results_pwys_mockdata.csv'),sep=',',row.names=F)
print (' ====================== DONE ====================')

# 2) collect permutation runs
# ==============================================
inFolder <- paste0(inF,'/mock_data_permutation_runs')
print(paste0('> collecting permutation results from ',inFolder))
# collect taxa jobs
inFiles = list.files(inFolder,pattern="*.RDS")
resPwy <- NULL
resTax <- NULL
for (i in inFiles) {
  if (grepl('tax',i)) {
    print(paste0(' > collecting run ',i))
    inPerm <- readRDS(paste0(inFolder,'/',i))
    inPerm2 <- inPerm[["permutations"]]
    resTax <- rbind.data.frame(resTax,inPerm2)
  } else if (grepl('pwy',i)) {
    print(paste0(' > collecting run ',i))
    inPerm <- readRDS(paste0(inFolder,'/',i))
    inPerm2 <- inPerm[["permutations"]]
    resPwy <- rbind.data.frame(resPwy,inPerm2)
  }
}
resAll <- rbind.data.frame(resPwy,resTax)

# 3/a) load results & calculate empirical study-wise p-values for Taxa
# ==================================================
# > main results
inRes <- read.table('results_mockdata_taxa.csv',sep=',',header = T)
inRes$SW_PV_ID <- NA
inRes$SW_PV_COHOUSING.ID_DMP <- NA
inRes$SW_PV_famID <- NA
inRes$SW_PV_Residual <- NA
# > model selection results
# calculate study-wise permutation p-values (using taxa runs only)
resSel <- resTax
for (rn in c(1:nrow(inRes))) {
  # H2
  inRes$SW_PV_ID[rn] <- sum(inRes$VE_ID[rn] - resSel$ID <= 0)/nrow(resSel)
  # COHOUSING
  inRes$SW_PV_COHOUSING.ID_DMP[rn] <- sum(inRes$VE_COHOUSING.ID_DMP [rn] - resSel$COHOUSING.ID_DMP <= 0)/nrow(resSel)
  # Family ID
  inRes$SW_PV_famID[rn] <- sum(inRes$VE_famID[rn] - resSel$famID <= 0)/nrow(resSel)
}
inRes2 <- inRes[order(inRes$Pv_ID),]
inRes2$SW_PV_Residual <- NULL
# FDRs (permutations)
inRes2$SW_FDR_ID <- p.adjust(inRes2$SW_PV_ID)
inRes2$SW_FDR_famID <- p.adjust(inRes2$SW_PV_famID)
inRes2$SW_FDR_COHOUSING.ID_DMP <- p.adjust(inRes2$SW_PV_COHOUSING.ID_DMP)
# FDRs (LOG-test)
inRes2$FDR_ID <- p.adjust(inRes2$Pv_ID)
inRes2$FDR_famID <- p.adjust(inRes2$Pv_famID)
inRes2$FDR_COHOUSING.ID_DMP <- p.adjust(inRes2$Pv_COHOUSING.ID_DMP)
# sort
inRes2 <- inRes2[order(inRes2$Pv_ID),]
# test output
inRes2[,c("Trait.short","VE_ID","CI_ID","empPV_ID","Pv_ID","FDR_ID","SW_PV_ID","SW_FDR_ID","Pv_COHOUSING.ID_DMP",
          "SW_PV_COHOUSING.ID_DMP",
          "Pv_famID","SW_PV_famID")][1:10,]

# fix Inf in CIs (Inf = 1)
for (cc in grep('CI',colnames(inRes2))) {
  inRes2[[cc]] <- as.character(inRes2[[cc]])
  inRes2[[cc]] <- gsub('Inf','1.00',inRes2[[cc]])
}

inRes3 <- inRes2[,c("Trait","Trait.short","N","Nz","Prev",
                    "Model",
                    "VE_ID","CI_ID","LogL_ID","Pv_ID","FDR_ID","SW_PV_ID","SW_FDR_ID",
                    "VE_COHOUSING.ID_DMP","CI_COHOUSING.ID_DMP","LogL_COHOUSING.ID_DMP","Pv_COHOUSING.ID_DMP","FDR_COHOUSING.ID_DMP","SW_PV_COHOUSING.ID_DMP","SW_FDR_COHOUSING.ID_DMP",
                    "VE_famID","CI_famID","LogL_famID","Pv_famID","FDR_famID","SW_PV_famID","SW_FDR_famID",
                    "VE_Residual","CI_Residual","LogL_Residual",
                    "Est_Age","SE_Age","P.Chisq_Age",
                    "Est_Sex","SE_Sex","P.Chisq_Sex",
                    "Est_META.DNA.postclean.reads","SE_META.DNA.postclean.reads","P.Chisq_META.DNA.postclean.reads",
                    "Est_META.POOP.Freq","SE_META.POOP.Freq","P.Chisq_META.POOP.Freq"
                    )] 

write.table(x = inRes3,file = 'results_mockdata_withFDRs_and_CIs_taxa.csv',sep=',',row.names = F)

# =======================================================
# =======================================================
# load results (PWYs)
# =======================================================
# =======================================================
inRes <- read.table('results_pwys_mockdata.csv',sep=',',header=T)
inRes$SW_PV_ID <- NA
inRes$SW_PV_COHOUSING.ID_DMP <- NA
inRes$SW_PV_famID <- NA
inRes$SW_PV_Residual <- NA
# calculate study-wise permutation p-values (using pathway runs only)
resSel <- resPwy
for (rn in c(1:nrow(inRes))) {
  # H2
  inRes$SW_PV_ID[rn] <- sum(inRes$VE_ID[rn] - resSel$ID <= 0)/nrow(resSel)
  # COHOUSING
  inRes$SW_PV_COHOUSING.ID_DMP[rn] <- sum(inRes$VE_COHOUSING.ID_DMP [rn] - resSel$COHOUSING.ID_DMP <= 0)/nrow(resSel)
  # Family ID
  inRes$SW_PV_famID[rn] <- sum(inRes$VE_famID[rn] - resSel$famID <= 0)/nrow(resSel)
}
inRes2 <- inRes[order(inRes$Pv_ID),]
inRes2$SW_PV_Residual <- NULL
# FDRs (permutations)
inRes2$SW_FDR_ID <- p.adjust(inRes2$SW_PV_ID)
inRes2$SW_FDR_famID <- p.adjust(inRes2$SW_PV_famID)
inRes2$SW_FDR_COHOUSING.ID_DMP <- p.adjust(inRes2$SW_PV_COHOUSING.ID_DMP)
# FDRs (LOG-test)
inRes2$FDR_ID <- p.adjust(inRes2$Pv_ID)
inRes2$FDR_famID <- p.adjust(inRes2$Pv_famID)
inRes2$FDR_COHOUSING.ID_DMP <- p.adjust(inRes2$Pv_COHOUSING.ID_DMP)
# sort
inRes2 <- inRes2[order(inRes2$Pv_ID),]
# example output
inRes2[,c("Trait","VE_ID","CI_ID","CI_famID","CI_COHOUSING.ID_DMP","empPV_ID","Pv_ID","FDR_ID","SW_PV_ID","SW_FDR_ID","Pv_COHOUSING.ID_DMP",
          "SW_PV_COHOUSING.ID_DMP",
          "Pv_famID","SW_PV_famID")][1:10,]

# fix Inf in CIs [Inf = 1]
for (cc in grep('CI',colnames(inRes2))) {
  inRes2[[cc]] <- as.character(inRes2[[cc]])
  inRes2[[cc]] <- gsub('Inf','1.00',inRes2[[cc]])
}

inRes3 <- inRes2[,c("Trait","Trait.short","N","Nz","Prev",
                    "Model",
                    "VE_ID","CI_ID","LogL_ID","Pv_ID","FDR_ID","SW_PV_ID","SW_FDR_ID",
                    "VE_COHOUSING.ID_DMP","CI_COHOUSING.ID_DMP","LogL_COHOUSING.ID_DMP","Pv_COHOUSING.ID_DMP","FDR_COHOUSING.ID_DMP","SW_PV_COHOUSING.ID_DMP","SW_FDR_COHOUSING.ID_DMP",
                    "VE_famID","CI_famID","LogL_famID","Pv_famID","FDR_famID","SW_PV_famID","SW_FDR_famID",
                    "VE_Residual","CI_Residual","LogL_Residual",
                    "Est_Age","SE_Age","P.Chisq_Age",
                    "Est_Sex","SE_Sex","P.Chisq_Sex",
                    "Est_META.DNA.postclean.reads","SE_META.DNA.postclean.reads","P.Chisq_META.DNA.postclean.reads",
                    "Est_META.POOP.Freq","SE_META.POOP.Freq","P.Chisq_META.POOP.Freq"
)] 
write.table(x = inRes3,file = 'results_mockdata_withFDRs_and_CIs_pwys.csv',sep=',',row.names = F)
