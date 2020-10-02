# =======================================================================
# By R.Gacesa, UMCG (2019)

#  Make pretty plots for model comparisons (after models are ready)
# =======================================================================
library(caret)
library(vegan)

# set working directory and load extra libraries
#setwd("D:/UMCG/ML_IBD_v5/")
source('../myLibs_v3/R_ML_scripts_v4.R')
source('../myLibs_v3/R_Microbiome_scripts.R')

# =========================================================
# COMPARISOR FUNCTION 
# =========================================================
runComparisons <- function(mName,mNameF,inFolder,botsFac = 0.1) {

  bootsXV = 10*botsFac
  bootsTest = 100*botsFac
  outFolder = paste0('final_plots/model_comparisons_',mNameF,'_',mtd)
  # cleanup
  print(' >> cleaning up')
  if (dir.exists(outFolder)) {
    unlink(outFolder,recursive = T)
  }
  dir.create(outFolder)
  
  print(' >> loading null model')
  # > load "null model" and training set
  # ===========================================
  fitPhenosASB <- readRDS(paste0('final_models_null/',inFolder,'/raw_fittedmdls/mdl_raw_glm_Phenos.rds'))
  fitPhenosASBtstSet <- read.table(paste0('final_models_null/',inFolder,'/inputData/Phenos_tstSet.csv'),sep=',',header = T)
  
  print(' >> loading optimised models')
  # > load optimised models
  # ===========================================
  fitFCal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_FCal_v2.rds'))
  fitCgA <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_CgA_v2.rds'))
  fitHBD2 <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_HBD2_v2.rds'))
  
  fitFCalCgA <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_FCal_CgA_v2.rds'))
  fitFCalHBD2 <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_FCal_HBD2_v2.rds'))
  fitCgAHBD2 <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_CgA_HBD2_v2.rds'))
  fitFCalCgAHBD2 <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_FCal_HBD2_CgA_v2.rds'))
  
  # MB
  fitMB_Phy <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tPhyla_v2.rds'))
  fitMB_Gen <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tGenera_v2.rds'))
  fitMB_Spec <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tSpecies_v2.rds'))
  fitMB_PWYs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_pPWY_vMax.rds'))
  fitMB_TaxaPWYs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_PWY_Tax_vMax.rds'))
  
  # MB + Calprotectin
  fitMB_Phy_Cal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tPhyla_Cal_vMax.rds'))
  fitMB_Gen_Cal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tGenera_Cal_vMax.rds'))
  fitMB_Spec_Cal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tSpecies_Cal_vMax.rds'))
  fitMB_PWYs_Cal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_pPWY_Cal_vMax.rds'))
  fitMB_TaxaPWYs_Cal <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_PWY_Tax_Cal_vMax.rds'))
  
  # MB + All BMs
  fitMB_Phy_BMs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tPhyla_BM_vMax.rds'))
  fitMB_Gen_BMs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tGenera_BM_vMax.rds'))
  fitMB_Spec_BMs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_tSpecies_BM_vMax.rds'))
  fitMB_PWYs_BMs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_pPWY_BM_vMax.rds'))
  fitMB_TaxaPWYs_BMs <- readRDS(paste0(inFolder,'/opt_fittedmdls/mdl_opt_',mtd,'_PWY_Tax_BM_vMax.rds'))
  
  modelsOpt <- list(fitFCal,fitCgA,fitHBD2,fitFCalCgA,fitFCalHBD2,fitCgAHBD2,fitFCalCgAHBD2,
                    fitMB_Phy,fitMB_Gen,fitMB_Spec,fitMB_PWYs,fitMB_TaxaPWYs,
                    fitMB_Phy_Cal,fitMB_Gen_Cal,fitMB_Spec_Cal,fitMB_PWYs_Cal,fitMB_TaxaPWYs_Cal,
                    fitMB_Phy_BMs,fitMB_Gen_BMs,fitMB_Spec_BMs,fitMB_PWYs_BMs,fitMB_TaxaPWYs_BMs)
  modelsOptN <- list("o.FCal","o.CgA","o.HBD2","o.FCal_CgA","o.FCal_HBD2","o.CgA_HBD2","o.BMs",
                     "o.Phyla","o.Genera","o.Species","o.PWYs","o.TaxaPWYs",
                     "o.Phyla.FCal","o.Genera.FCal","o.Species.FCal","o.PWYs.FCal","o.TaxaPWYs.FCal",
                     "o.Phyla.BM","o.Genera.BM","o.Species.BM","o.PWYs.BM","o.TaxaPWYs.BM")
  
  print(' >> loading raw models')
  # > load "raw" models
  # ===========================================
  fitRawFCal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_FCal.rds'))
  fitRawCgA <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_CgA.rds'))
  fitRawHBD2 <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_HBD2.rds'))
  
  fitRawFCalCgA <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_FCal_CgA.rds'))
  fitRawFCalHBD2 <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_FCal_HBD2.rds'))
  fitRawCgAHBD2 <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_CgA_HBD2.rds'))
  fitRawFCalCgAHBD2 <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_FCal_HBD2_CgA.rds'))
  
  # MB
  fitRawMB_Phy <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tPhyla.rds'))
  fitRawMB_Gen <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tGenera.rds'))
  fitRawMB_Spec <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tSpecies.rds'))
  fitRawMB_PWYs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_pPWY.rds'))
  fitRawMB_TaxaPWYs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_PWY_Tax.rds'))
  
  # MB + Calprotectin
  fitRawMB_Phy_Cal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tPhyla_Cal.rds'))
  fitRawMB_Gen_Cal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tGenera_Cal.rds'))
  fitRawMB_Spec_Cal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tSpecies_Cal.rds'))
  fitRawMB_PWYs_Cal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_pPWY_Cal.rds'))
  fitRawMB_TaxaPWYs_Cal <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_PWY_Tax_Cal.rds'))
  
  # MB + All BMs
  fitRawMB_Phy_BMs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tPhyla_BM.rds'))
  fitRawMB_Gen_BMs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tGenera_BM.rds'))
  fitRawMB_Spec_BMs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_tSpecies_BM.rds'))
  fitRawMB_PWYs_BMs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_pPWY_BM.rds'))
  fitRawMB_TaxaPWYs_BMs <- readRDS(paste0(inFolder,'/raw_fittedmdls/mdl_raw_',mtd,'_PWY_Tax_BM.rds'))
  
  modelsRaw <- list(fitRawFCal,fitRawCgA,fitRawHBD2,fitRawFCalCgA,fitRawFCalHBD2,fitRawCgAHBD2,fitRawFCalCgAHBD2,
                    fitRawMB_Phy,fitRawMB_Gen,fitRawMB_Spec,fitRawMB_PWYs,fitRawMB_TaxaPWYs,
                    fitRawMB_Phy_Cal,fitRawMB_Gen_Cal,fitRawMB_Spec_Cal,fitRawMB_PWYs_Cal,fitRawMB_TaxaPWYs_Cal,
                    fitRawMB_Phy_BMs,fitRawMB_Gen_BMs,fitRawMB_Spec_BMs,fitRawMB_PWYs_BMs,fitRawMB_TaxaPWYs_BMs)
  modelsRawN <- list("r.FCal","r.CgA","r.HBD2","r.FCal_CgA","r.FCal_HBD2","r.CgA_HBD2","r.BMs",
                     "r.Phyla","r.Genera","r.Species","r.PWYs","r.TaxaPWYs",
                     "r.Phyla.FCal","r.Genera.FCal","r.Species.FCal","r.PWYs.FCal","r.TaxaPWYs.FCal",
                     "r.Phyla.BM","r.Genera.BM","r.Species.BM","r.PWYs.BM","r.TaxaPWYs.BM")
  
  print(' >> loading test sets')
  # > load testset
  tSetRaw <- read.table(paste0(inFolder,'/inputData/PWY_Tax_BM_tstSet.csv'),sep=',',header=T)
  # > load transformer
  prepper <- readRDS(paste0(inFolder,'/inputData/FCal_HBD2_CgA_trSet_prepMdl.RData'))
  #  >> transform
  tSetRdy <- predict(prepper,tSetRaw)
  
  print(' >> doing optimised VS raw comparison')
  # ============= optimised vs RAW =====================
  # ===============================================================
  # > compare (X-validation)
  # for (c in c(1:length(modelsRawN))) {
  #   g <- compareModelsTrainingCV(fittedMdls = list(modelsRaw[[c]],modelsOpt[[c]]),
  #                                modelNames = c(modelsRawN[[c]],modelsOptN[[c]]),
  #                                mtd = mtd,
  #                                posClass = pClass,
  #                                roc.smooth = T,
  #                                annotateAUConly = T,
  #                                roc.conf.boot = bootsXV,
  #                                tit = ''
  #   )
  #   nrRaw <- length(predictors(modelsRaw[[c]]))
  #   nrOpt <- length(predictors(modelsOpt[[c]]))
  #   g <- g + theme(legend.position = "top") + ggtitle(paste0(gsub('r\\.','',modelsRawN[[c]]),", opt <",nrOpt,"> vs raw <",nrRaw,"> (x-validation)"))
  #   print(g)
  #   ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/_raw_vs_opt_',gsub('r\\.','',modelsRawN[[c]]),'_XV_',mtd,'.png'))
  # }
  
  tSetRdyPN <- purgeMGNames(tSetRdy)
  # > compare (test-set)
  # for (c in c(1:length(modelsRawN))) {
  #   gg <- compareMdlsDatasets(mdls = list(modelsRaw[[c]],modelsOpt[[c]]),
  #                             dataSets = list(tSetRdyPN),
  #                             mdNames = c(modelsRawN[[c]],modelsOptN[[c]]),
  #                             response = "Diagnosis",
  #                             roc.smooth = T,
  #                             specSensAnnot = T,
  #                             target = F,
  #                             removeLegend = F,
  #                             posClass = pClass,
  #                             roc.conf.boot = bootsTest,
  #                             tit = ""
  #   )
  #   nrRaw <- length(predictors(modelsRaw[[c]]))
  #   nrOpt <- length(predictors(modelsOpt[[c]]))
  #   gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(gsub('r\\.','',modelsRawN[[c]]),", opt <",nrOpt,"> vs raw <",nrRaw,"> (Test)"))
  #   print(gg[[1]])
  #   ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/_raw_vs_opt_',gsub('r\\.','',modelsRawN[[c]]),'_test_',mtd,'.png'))
  # }
  
  print(' >> comparing biomarkers')
  bootsXV = 100*botsFac
  bootsTest = 1000*botsFac
  # ============= biomarkers individually ==============
  # ===============================================================
  # > compare (X-validation)
  g <- compareModelsTrainingCV(fittedMdls = list(fitPhenosASB,fitFCal,fitCgA,fitHBD2,fitFCalCgAHBD2),
                               modelNames = c("Phenos","FCal","CgA","HBD2","BMs"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               roc.conf.boot = bootsXV,
                               tit = paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Biomarkers_individual_XV_',mtd,'.png'))
  
  #> compare (test-set)
  gg <- compareMdlsDatasets(mdls = list(fitPhenosASB,fitFCal,fitCgA,fitHBD2,fitFCalCgAHBD2),
                            dataSets = list(fitPhenosASBtstSet,tSetRdy,tSetRdy,tSetRdy,tSetRdy),
                            mdNames = c("Phenos","FCal","CgA","HBD2","BMs"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = paste0(mName," [",mtd,"] (Biomarkers, test set)")
  )
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Biomarkers, test set )"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_IBD_vs_IBS_Biomarkers_individual_Test_',mtd,'.png'))
  write.table(gg[[3]],file= paste0(outFolder,'/',mNameF,'_Biomarkers_individual_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  print(' >> comparing combined biomarkers')
  # ========== biomarkers combined ==========
  # ===============================================================
  # > compare (X-validation)
  g <- compareModelsTrainingCV(fittedMdls = list(fitFCalCgA,fitFCalHBD2,fitCgAHBD2,fitFCalCgAHBD2),
                               modelNames = c("FCal+CgA","FCal+HBD2","CgA+HBD2","All"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Biomarkers_combos_XV_',mtd,'.png'))
  
  # > compare (test-set)
  gg <- compareMdlsDatasets(mdls = list(fitFCalCgA,fitFCalHBD2,fitCgAHBD2,fitFCalCgAHBD2),
                            dataSets = list(tSetRdy),
                            mdNames = c("FCal+CgA","FCal+HBD2","CgA+HBD2","All"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = paste0(mName," [",mtd,"] (Biomarkers, test set)")
  )
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Biomarkers, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Biomarkers_combos_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Biomarkers_combos_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  print(' >> comparing microbiome')
  # ============= microbiome (individual levels) ==================
  # ===============================================================
  # XV
  g <- compareModelsTrainingCV(fittedMdls = list(fitMB_Phy,fitMB_Gen,fitMB_Spec,fitMB_PWYs,fitMB_TaxaPWYs),
                               modelNames = c("Phyla","Genera","Species","PWYs","PWYs+Taxa"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_levels_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  gg <- compareMdlsDatasets(mdls = list(fitMB_Phy,fitMB_Gen,fitMB_Spec,fitMB_PWYs,fitMB_TaxaPWYs),
                            dataSets = list(tSetRdyPN),
                            mdNames = c("Phyla","Genera","Species","PWYs","PWYs+Taxa"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = paste0(mName," [",mtd,"] (Biomarkers, test set)")
  )
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_levels_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_levels_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  print(' >> comparing microbiome + calprot')
  # ================== MICROBIOME + CALPROT =====================
  # ===============================================================
  g <- compareModelsTrainingCV(fittedMdls = list(fitMB_Phy_Cal,fitMB_Gen_Cal,fitMB_Spec_Cal,fitMB_PWYs_Cal,fitMB_TaxaPWYs_Cal),
                               modelNames = c("FC.Ph","FC.Gen","FC.Sp","FC.PWY","FC.PTax"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = paste0(mName," [",mtd,"] (Biomarkers, x-validation)"))
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + FCal, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_FCal_levels_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  gg <- compareMdlsDatasets(mdls = list(fitMB_Phy_Cal,fitMB_Gen_Cal,fitMB_Spec_Cal,fitMB_PWYs_Cal,fitMB_TaxaPWYs_Cal),
                            dataSets = list(tSetRdyPN,tSetRdyPN,tSetRdyPN,tSetRdyPN,tSetRdy),
                            mdNames = c("FC.Ph","FC.Gen","FC.Sp","FC.PWY","FC.PTax"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = '')
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + FCal, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_FCal_levels_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_FCal_levels_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  print(' >> comparing microbiome + biomarkers')
  # ================== MICROBIOME + BM =====================
  # ===============================================================
  g <- compareModelsTrainingCV(fittedMdls = list(fitMB_Phy_BMs,fitMB_Gen_BMs,fitMB_Spec_BMs,fitMB_PWYs_BMs,fitMB_TaxaPWYs_BMs),
                               modelNames = c("BM.Ph","BM.Gen","BM.Sp","BM.PWY","BM.PTax"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = '')
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_BM_levels_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  gg <- compareMdlsDatasets(mdls = list(fitMB_Phy_BMs,fitMB_Gen_BMs,fitMB_Spec_BMs,fitMB_PWYs_BMs,fitMB_TaxaPWYs_BMs),
                            dataSets = list(tSetRdyPN,tSetRdyPN,tSetRdyPN,tSetRdyPN,tSetRdy),
                            mdNames = c("BM.Ph","BM.Gen","BM.Sp","BM.PWY","BM.PTax"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = '')
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_FCal_levels_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_BM_levels_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  print(' >> doing extra plots')
  
  #   ================== FINAL HAPPY PLOT 1 [species] =====================
  # ===============================================================
  # FCal, all 3, Microbiome+FCal, Microbiome+3 BM
  
  g <- compareModelsTrainingCV(fittedMdls = list(fitPhenosASB,fitMB_Spec,fitMB_Spec_BMs),
                               modelNames = c("Phenos","MB","MB+BMark"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = '')
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBsBM_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  gg <- compareMdlsDatasets(mdls = list(fitPhenosASB,fitMB_Spec,fitMB_Spec_BMs),
                            dataSets = list(fitPhenosASBtstSet,tSetRdyPN,tSetRdyPN),
                            mdNames = c("Phenos","MB","MB+BMark"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = '')
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBsBM_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_vs_MBsBM_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  #   ================== FINAL HAPPY PLOT 1 [genera] =====================
  # ===============================================================
  # FCal, all 3, Microbiome+FCal, Microbiome+3 BM
  
  g <- compareModelsTrainingCV(fittedMdls = list(fitPhenosASB,fitMB_Gen,fitMB_Gen_BMs),
                               modelNames = c("Phenos","MB","MB+BMark"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = '')
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBsBM_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  gg <- compareMdlsDatasets(mdls = list(fitPhenosASB,fitMB_Gen,fitMB_Gen_BMs),
                            dataSets = list(fitPhenosASBtstSet,tSetRdyPN,tSetRdyPN),
                            mdNames = c("Phenos","MB","MB+BMark"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = '')
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBsBM_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_vs_MBgBM_Test_',mtd,'.csv'),row.names = F,sep=',')
  
  
  #   ================== FINAL HAPPY PLOT 2 =====================
  # ===============================================================
  # FCal, all 3, Microbiome+FCal, Microbiome+3 BM
  
  g <- compareModelsTrainingCV(fittedMdls = list(fitPhenosASB,fitFCal,fitFCalCgAHBD2,fitMB_Gen,fitMB_Gen_BMs,fitMB_Spec,fitMB_Spec_BMs),
                               modelNames = c("Ph","FCal","BM","Gen","Gen+BM","Spec","Spec+BM"),
                               mtd = "glm",
                               posClass = pClass,
                               roc.smooth = T,
                               annotateAUConly = T,
                               removeLegend = F,
                               roc.conf.boot = bootsXV,
                               tit = '')
  g <- g + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, x-validation)"))
  print(g)
  ggsave(plot=g,scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBzBM2_XV_',mtd,'.png'))
  
  # test set
  tSetRdyPN <- purgeMGNames(tSetRdy)
  # compareMdlsDatasets(mdls = list(fitPhenosASB,fitMB_Spec_BMs,fitMB_Spec,fitMB_Gen,fitMB_Spec),
  #                     dataSets = list(fitPhenosASBtstSet,tSetRdyPN,tSetRdyPN,tSetRdyPN,tSetRdyPN),
  #                     mdNames = c("Ph","Spec.BM","Gen.BM","Spec","fitMB_Gen","BM"),
  #                     response = "Diagnosis",
  #                     roc.smooth = T,
  #                     specSensAnnot = T,
  #                     target = F,
  #                     removeLegend = F,
  #                     posClass = pClass,
  #                     roc.conf.boot = bootsTest,
  #                     tit = '')
  
  gg <- compareMdlsDatasets(mdls =     list(fitPhenosASB,      fitFCal,  fitFCalCgAHBD2,fitMB_Gen,fitMB_Gen_BMs,fitMB_Spec,fitMB_Spec_BMs),
                            dataSets = list(fitPhenosASBtstSet,tSetRdyPN,tSetRdyPN,     tSetRdyPN,tSetRdyPN,    tSetRdyPN, tSetRdyPN),
                            mdNames = c("Ph","FCal","BM","Gen","Gen+BM","Spec","Spec+BM"),
                            response = "Diagnosis",
                            roc.smooth = T,
                            specSensAnnot = T,
                            target = F,
                            removeLegend = F,
                            posClass = pClass,
                            roc.conf.boot = bootsTest,
                            tit = '')
  gg[[1]] <- gg[[1]] + theme(legend.position = "top") + ggtitle(paste0(mName," [",mtd,"] (Microbiome + BM, test set)"))
  print(gg[[1]])
  ggsave(plot=gg[[1]],scale = scl,filename = paste0(outFolder,'/',mNameF,'_Microbome_vs_MBzBM2_Test_',mtd,'.png'))
  write.table(gg[3],file= paste0(outFolder,'/',mNameF,'_Microbome_vs_MBzBM2_Test_',mtd,'.csv'),row.names = F,sep=',')
}

# =======================================================================================
# comparison function 2 (only makes table for test and training set, no fancy plots)
# ==========================================================================================

setwd('~/UMCG/ML_IBD_v5/')
outFolder <- "model_results_IBD_vs_IBS_happy"
inFolderModels <- "models_happy_MB_IBDvsIBS_kappa/opt_fittedmdls/"
testSetPath <- "models_happy_MB_IBDvsIBS_kappa/inputData/PWY_Tax_BM_tstSet.csv"
testSetPrepperPath <- "models_happy_MB_IBDvsIBS_kappa/inputData/PWY_Tax_BM_trSet_prepMdl.RData"
posClass = "IBD"
responseVar = "Diagnosis"

runComparisons <- function(outFolder,inFolderModels,testSetPath,testSetPrepperPath,responseVar,posClass) {
  
  outFolder = paste0(outFolder)
  # cleanup
  print(' >> cleaning up')
  if (dir.exists(outFolder)) {
    unlink(outFolder,recursive = T)
  }
  dir.create(outFolder)
  
  # load test set and transformer
  print ('loading test set and preprocessor')
  testSetRaw <- read.table(testSetPath,sep=',',header=T)
  testSetPrepper <- readRDS(testSetPrepperPath)
  print (' > pre-processing test set')
  testSetPrepped <- predict(testSetPrepper,testSetRaw)
  
  # assess X-validation of training set
  # =========================================
  # iterate over all models in inFolder
  print('>> assessing X-validation performance of models')
  resultsCV = NULL
  mdls <- dir(inFolderModels)
  mdls <- mdls[grep('_v2\\.',mdls)]
  for (md in mdls) {
    fittedMdl <- readRDS(paste0(inFolderModels,'/',md))
    print(paste0(' >> analyzing x-validation performance for ',md))
    confM <- confusionMatrix(fittedMdl$pred$pred,fittedMdl$pred$obs)
    # acc, kappa
    acc <- confM[[3]][1]
    kappa <- confM[[3]][2]
    # sens spec ppv npv prec recall F1 prev bacc
    sens <- confM[[4]][1]
    spec <- confM[[4]][2]
    ppv <- confM[[4]][3]
    npv <- confM[[4]][4]
    f1 <- confM[[4]][7]
    bacc <- confM[[4]][11]
    # get auc
    r <- roc(predictor = fittedMdl$pred[[posClass]], response = fittedMdl$pred$obs, auc=T,percent=T,smooth=F)
    auc <- as.numeric(r$auc)/100
    nrVars <- length(fittedMdl$coefnames)
    oneRow <- data.frame(Model=md,Method=fittedMdl$method,NrVars=nrVars,ACC=acc,Kappa=kappa,Sensitivity=sens,Specificity=spec,PPV=ppv,NPV=npv,B.ACC=bacc,F1=f1,AUC=auc)
    resultsCV <- rbind.data.frame(resultsCV,oneRow)
    row.names(resultsCV) <- NULL
  }
  # assess test set preformance
  # =======================================
  # iterate over all models in inFolder
  resultsTest = NULL
  mdls <- dir(inFolderModels)
  mdls <- mdls[grep('_v2\\.',mdls)]
  for (md in mdls) {
    fittedMdl <- readRDS(paste0(inFolderModels,'/',md))
    print(paste0(' >> analyzing test set peformance for ',md))
    predictedTest <- predict(fittedMdl,testSetPrepped)
    confM <- confusionMatrix(testSetPrepped[[responseVar]],predictedTest)
    # acc, kappa
    acc <- confM[[3]][1]
    kappa <- confM[[3]][2]
    # sens spec ppv npv prec recall F1 prev bacc
    sens <- confM[[4]][1]
    spec <- confM[[4]][2]
    ppv <- confM[[4]][3]
    npv <- confM[[4]][4]
    f1 <- confM[[4]][7]
    bacc <- confM[[4]][11]
    # get auc
    r <- roc(predictor = fittedMdl$pred[[posClass]], response = fittedMdl$pred$obs, auc=T,percent=T,smooth=F)
    auc <- as.numeric(r$auc)/100
    nrVars <- length(fittedMdl$coefnames)
    oneRow <- data.frame(Model=md,Method=fittedMdl$method,NrVars=nrVars,ACC=acc,Kappa=kappa,Sensitivity=sens,Specificity=spec,PPV=ppv,NPV=npv,B.ACC=bacc,F1=f1,AUC=auc)
    resultsTest <- rbind.data.frame(resultsTest,oneRow)
    row.names(resultsTest) <- NULL
  }
  ret = list(CV=resultsCV,test=resultsTest)
  ret
}

write.table(ret[["CV"]],'model_results_IBD_vs_IBS_happy/results_XV_v2.csv',sep=',',row.names = F)
write.table(ret[["test"]],'model_results_IBD_vs_IBS_happy/results_test_v2.csv',sep=',',row.names = F)

# ==============================================================================
# ==============================================================================
# STEP 4: do comparison of models
# ==============================================================================
# ==============================================================================

# # ======================== IBD VS HC ==============================
# scl = 1.25
# mName = "IBD vs HC"
# mNameF = "IBD_vs_HC"
# mtd = 'glm'
# inFolder = 'final_models_MB_IBD_vs_HC'
# pClass = "IBD"
# runComparisons(mName,mNameF,inFolder,botsFac = 0.1)
# # > SVM
# scl = 1.25
# mName = "IBD vs HC"
# mNameF = "IBD_vs_HC"
# mtd = 'svmRadial'
# inFolder = 'final_models_MB_IBD_vs_HC'
# pClass = "IBD"
# runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# ======================== IBD VS IBS ==============================
scl = 1.25
mName = "IBD vs IBS"
mNameF = "IBD_vs_IBS"
mtd = 'glm'
pClass = "IBD"
#outFolder = paste0('final_plots/model_comparisons_',mNameF,'_',mtd)
inFolder = 'final_models_MB_IBD_vs_IBS'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# ======================== IBD_UC_vs_CD ==============================
scl = 1.25
mName = "UC vs CD"
mNameF = "IBD_UC_vs_CD"
#mtd = 'glm'
pClass = "CD"
inFolder = 'final_models_MB_IBD_UC_vs_CD'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# ======================== IBD location ==============================
scl = 1.25
mName = "Colonic IBD"
mNameF = "IBD_loc_colon_vs_ileum"
mtd = 'glm'
pClass = "colon"
inFolder = 'final_models_MB_IBD_loc_colon_vs_ileum'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# ======================== Flares (Marjolein) ========================
scl = 1.25
mName = "IBD activity"
mNameF = "MB_IBD_flares_Marjolein"
#mtd = 'glm'
pClass = "Y"
inFolder = 'final_models_MB_IBD_flares_Marjolein'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)


# ======================== Flares (HBI-SCCAI) ========================
scl = 1.25
mName = "IBD activity"
mNameF = "MB_IBD_flares_HBISCCAI"
#mtd = 'glm'
pClass = "Y"
inFolder = 'final_models_MB_IBD_flares_HBISCCAI'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# ======================== Flares (Marjolein) ========================
scl = 1.25
mName = "IBD activity"
mNameF = "MB_IBD_flares_Marjolein"
mtd = 'glm'
pClass = "Y"
inFolder = 'final_models_MB_IBD_flares_Marjolein'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

# 
# # ======================== Flares (HBI-SCCAI) ========================
# scl = 1.25
# mName = "IBD activity"
# mNameF = "MB_IBD_flares_HBISCCAI"
# #mtd = 'glm'
# pClass = "Y"
# inFolder = 'final_models_MB_IBD_flares_HBISCCAI'

runComparisons(mName,mNameF,inFolder,botsFac = 0.1)

