# ===========================================
# DAG3 ANALYSIS, diseases ...
# >> version 4 (control group is whole population)
# >> include only some diseases
# ===========================================

# load libaries
library(arsenal)
require(knitr)
require(survival)
library(MatchIt)

setwd('~/UMCG/DAG3_stats/')
source('../myLibs/R_Microbiome_scripts.R')
source('../myLibs/R_ML_scripts_v3b.R')
source('../myLibs/R_Misc.R')

# =========================================
# MAIN 
# =========================================
options("optmatch_max_problem_size" = Inf)
# load data (filtered - asin(sqrt) normalised - rescaled - corrected for technicals)
#inDF <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor_corrtech.csv',sep=',',header=T)
inDF <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor_corrtechASB.csv',sep=',',header=T)

inDFprev <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled_PA.csv',sep=',',header=T)
inDFprev <- purgeMGNames(inDFprev)

# load diversity data
# NOTE: diversity is calculated from non-normalised data
inDFdiv <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_diversitiesonly.csv',sep=',',header=T)
# clean names of Taxa
inDFf <- purgeMGNames(inDF)
# select taxa for testing
toTest <- colnames(subsetMicrobiomeDF(inDFf,getTaxa = T,getPhenos = F,verbose = F,getVFs = F,getCARDs = F,getPWYs = F))
# PWYs
inDFpwy <- subsetMicrobiomeDF(inDFf,getTaxa = F,getPhenos = T,verbose = F,getVFs = F,getCARDs = F,getPWYs = T)
inDFpwy <- inDFpwy[,c(1,23:285)]
# PHENOS
inDFphenos <- subsetMicrobiomeDF(inDFf,getTaxa = F,getPhenos = T,verbose = F,getVFs = F,getCARDs = F,getPWYs = F)

# SETTINGS
# ==========================================
outFolder <- "output_healthDiseases_v4"
plotNominalPV <- F # if true, uses nominal P-value cut-off for plotting, otherwise uses multiple-testing corrected
taxRetestThresh <- 0.005
pwyRetestThresh <- 0.0005
doPlots <- F  # if F, make no plots!
plotTopX <- 10 # if > 0, plots top X regardless of stats
minAge <- 20 # throw out kids!
diagnosticPlots <- T
# ==========================================

# load phenotypes
# ==========================================
phenos <- read.table('../DAG3_meta/v2/health/curated/DAG3_T1T2_Health.csv',sep=',',header=T,quote = '"')
colnames(phenos) <- tolower(colnames(phenos))
#  -- select phenotypes to test --
selCols <- colnames(phenos)[c(2:46,57,58,65:67,69:77,80,81)]
phenNames  <- selCols

# define controls
# ==============
specificControls = F
controlIDs <- ""
sexAgeBMIMatch <- F
maxCtrlRatio <- 10
# =============

# CLEANUP 
# ==========================================
setwd('~/UMCG/DAG3_stats/')
unlink(outFolder,recursive = T,force = T)
dir.create(outFolder)
dir.create(paste0(outFolder,'/plots_PWYs'))
dir.create(paste0(outFolder,'/plots_prev'))
dir.create(paste0(outFolder,'/plots_ab'))
dir.create(paste0(outFolder,'/plots_div'))
dir.create(paste0(outFolder,'/tables'))
if (diagnosticPlots) {
  dir.create(paste0(outFolder,'/summaries'))
}

# run analyses
# > age limit (if required)
inDFphenos <- inDFphenos[inDFphenos$age > minAge,]
# > counter
phenNR <- 0
for (phen in selCols) {
  print(phen)
  phDF <- phenos[!is.na(phenos[[phen]]),c("dag3_sampleid",phen)]
  # merge with phenos for testing
  colnames(phDF)[1] <- "ID"
  cnt <- 0
  # > results dataframe (diversity)
  resultsDiv <- data.frame()
  # > advance phen counter
  phenNR <- phenNR + 1
  # > grind ...
  print(paste0(' >> GRINDING ',phen))
  phDFm <- merge(phDF,inDFdiv,by="ID")
  phDFm2 <- phDFm
  # merge with other phenos
  phDFm2 <- merge(phDFm2,inDFphenos,by.x="ID",by.y="DAG3_sampleID")
  # merge with taxa
  inDFft <- inDFf[,c("DAG3_sampleID",toTest)]
  phDFm2 <- merge(phDFm2,inDFft,by.x="ID",by.y="DAG3_sampleID")
  # merge with taxa (prevalence)
  toTestP <- paste0("p.",toTest)
  inDFftp <- inDFprev[,c("DAG3_sampleID",toTest)]
  colnames(inDFftp) <- paste0('p.',colnames(inDFftp))
  phDFm2 <- merge(phDFm2,inDFftp,by.x="ID",by.y="p.DAG3_sampleID")
  # merge with pwys
  phDFm2 <- merge(phDFm2,inDFpwy,by.x="ID",by.y="DAG3_sampleID")

  # DEFINE CONTROLS
  # > use positives and control list
  if (specificControls) {
    toKeep <- factor(c(as.character(phDFm2$ID[phDFm2[[phen]] == "Y"]),as.character(controlIDs)))
    phDFm2[[phen]] <- as.character(phDFm2[[phen]])
    phDFm2[[phen]][phDFm2$ID %in% controlIDs] <- "N"
    phDFm2 <- phDFm2[!is.na(phDFm2[[phen]]),]
    phDFm2 <- phDFm2[!phDFm2[[phen]]=="",]
    phDFm2 <- phDFm2[phDFm2$ID %in% toKeep,]
  } else {
    phDFm2[[phen]] <- as.character(phDFm2[[phen]])
    phDFm2[[phen]][is.na(phDFm2[[phen]])] <- "N"
    phDFm2[[phen]][phDFm2[[phen]]==""] <- "N"
  }
  # fix factors
  if (length(unique(phDFm2[[phen]])) == 2) {
    if ("y" %in% unique(phDFm2[[phen]]) | "n" %in% unique(phDFm2[[phen]])) {
      phDFm2[[phen]] <- factor(phDFm2[[phen]],levels = c("n","y"))
    } else if ("Y" %in% unique(phDFm2[[phen]]) | "N" %in% unique(phDFm2[[phen]])) {
      phDFm2[[phen]] <- factor(phDFm2[[phen]],levels = c("N","Y"))
    } else if ("0" %in% unique(phDFm2[[phen]]) | "1" %in% unique(phDFm2[[phen]])) {
      phDFm2[[phen]] <- factor(phDFm2[[phen]],levels = c("0","1"))
    }
  }
  # do case-control matching
  if (sexAgeBMIMatch) {
    if (length(levels(phDFm2[[phen]])) == 2) {
      frm <- reformulate(response = phen,termlabels = c("age","gender"))
      # mcs <- matchControls(frm,data=phDFm2,contlabel = "N",caselabel = "Y")
      # phDFm2 <- phDFm2[!is.na(mcs$factor),]
      drk <- phDFm2[,c("ID",phen,"gender","age")]
      drk[[phen]] <- as.character(drk[[phen]])
      drk[[phen]][drk[[phen]] == "Y"] <- 1
      drk[[phen]][drk[[phen]] == "N"] <- 0
      drk[[phen]] <- as.integer(drk[[phen]])
      options("optmatch_max_problem_size" = Inf)
      m.out <- matchit(frm, data = drk, method = "optimal",ratio=round(min(maxCtrlRatio,max(1,2000/sum(drk[[phen]] == "1"))),0))
      m.data <- match.data(m.out)
      phDFm2 <- phDFm2[phDFm2$ID %in% m.data$ID,]
    }
  }
  # diagnostics plots
  if (diagnosticPlots) {
    g <- ggplot(phDFm2,aes_string(x="age",col=phen)) + geom_density(size=2)
    ggsave(plot=g,filename = paste0(outFolder,'/summaries/phen_',phenNR,'_age.png'))
    g <- ggplot(phDFm2,aes_string(x="BMI",col=phen)) + geom_density(size=2)
    ggsave(plot=g,filename = paste0(outFolder,'/summaries/phen_',phenNR,'_BMI.png'))
    g <- ggplot(phDFm2,aes_string(x="gender",col=phen,fill=phen)) + geom_bar(position = "fill")
    ggsave(plot=g,filename = paste0(outFolder,'/summaries/phen_',phenNR,'_gender.png'))
    frm <- reformulate(response = phen,termlabels = c("age","gender","BMI"))
    write.table(x=summary(tableby(frm,data=phDFm2),text=T),
                file = paste0(outFolder,'/summaries/phen_',phenNR,'_summary.csv'),sep=',',row.names = F)
  }
  print("   >>> TESTING DIVERSITY <<< ")
  # SHANNON (diversity)
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.S.shannon",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  cnt <- cnt + 1
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Shannon diversity [Species]") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_sdivS_',cnt,'.png'))
  }
  # SHANNON (genera)
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.G.shannon",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  cnt <- cnt + 1
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Shannon diversity [Genera]") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_sdivG_',cnt,'.png'))
  }
  # SHANNON (pathways)
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.PWY.shannon",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  cnt <- cnt + 1
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Shannon diversity [PWYs]") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_sdivPWY_',cnt,'.png'))
  }
  # NUMBER OF SPECIES
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.nrS",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  cnt <- cnt + 1
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Number of Species") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_nrSpec_',cnt,'.png'))
  }
  # NUMBER OF GENERA
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.nrG",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  cnt <- cnt + 1
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Number of Genera") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_nrGen_',cnt,'.png'))
  }
  # NUMBER OF Virulence factors
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.nrVF",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  cnt <- cnt + 1
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Number of Virulence factors") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_nrVF_',cnt,'.png'))
  }
  # NUMBER OF CARDs
  t <- testOneFeature(dataIn=phDFm2,feature = "DIV.nrCARD",doSave = F,responseVar = phen,retPlot = T,display = "Pn")
  cnt <- cnt + 1
  resultsDiv <- rbind.data.frame(resultsDiv,t[[1]])
  if (doPlots) {
    g <- t[[2]] + theme(legend.position = "none") + ylab("Number of Antibiotic resistance genes") + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
    #print(g)
    ggsave(plot=g,filename = paste0(outFolder,'/plots_div/phen_',phenNR,'_plot_nrCARD_',cnt,'.png'))
  }
  resultsDiv$PHENOTYPE <- phen
  write.table(resultsDiv,paste0(outFolder,'/tables/phen_',phenNR,'_diversities.csv'),sep=',',row.names = F)
  
  # grind through species to see if there is anything in there
  # ======= ABUNDANCE =========
  print("   >>> TESTING ABUNDANCE <<< ")
  cnt = 0
  allTests <- data.frame()
  # > run pair-wise mann-u-whitney tests for relative abundance
  for (t in toTest) {
    #print(t)
    oneTest <- testOneFeature(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = F,discardZeros = F)
    oneTest$pValue[is.na(oneTest$pValue)] <- 1.0
    oneTest$FDR <- p.adjust(oneTest$pValue)
    allTests <- rbind.data.frame(allTests,oneTest)
  }
  allTests <- allTests[order(allTests$FDR),]
  allTests$PHENOTYPE <- phen
  allTests$Multitest.FDR <- p.adjust(allTests$pValue)
  write.table(allTests,paste0(outFolder,'/tables/phen_',phenNR,'_taxAbundances.csv'),sep=',',row.names = F)
  
  if (plotNominalPV) {
    toRetest <- unique(allTests$feature[allTests$FDR < taxRetestThresh])
  } else {
    toRetest <- unique(allTests$feature[allTests$Multitest.FDR < taxRetestThresh])
  }
  if (length(toRetest) < plotTopX) {toRetest <- unique(allTests$feature)[c(1:plotTopX)]}
  cnt <- 0
  if (length(toRetest) > 0 & doPlots) {
    for (t in toRetest) {
      #print(t)
      oneTest <- testOneFeature(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = T,discardZeros = F,retPlot = T,doSave = F)
      g <- oneTest[[2]] + theme(legend.position = "none") + ylab(t) + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
      #print(g)
      cnt <- cnt + 1
      ggsave(plot=g,filename = paste0(outFolder,'/plots_ab/phen_',phenNR,'_plot_abund_',cnt,'.png'))
    }
  }
  
  
  # ============ PREVALENCE OF TAXA ===============
  print("   >>> TESTING PREVALENCE <<< ")
  # > run pair-wise mann-u-whitney tests for relative abundance in CD vs HC
  allTests <- data.frame()
  phDFm2[[phen]] <- factor(phDFm2[[phen]])
  for (t in toTestP) {
    #print(t)
    #phDFm2[[t]][phDFm2[[t]] > 0 ] <- 1
    phDFm2[[t]] <- factor( phDFm2[[t]],levels=c("0","1"))
    if (length(unique(phDFm2[[t]])) > 1) {
      oneTest <- testOneFeaturePrevalence(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = F,retPlot = F)
      oneTest <- oneTest
      oneTest$pValue[is.na(oneTest$pValue)] <- 1.0
      oneTest$FDR <- p.adjust(oneTest$pValue)
      allTests <- rbind.data.frame(allTests,oneTest)
    }
  }
  allTests <- allTests[order(allTests$FDR),]
  allTests$PHENOTYPE <- phen
  allTests$Multitest.FDR <- p.adjust(allTests$pValue)
  write.table(allTests,paste0(outFolder,'/tables/phen_',phenNR,'_taxPrevalences.csv'),sep=',',row.names = F)
  
  if (plotNominalPV) {
    toRetest <- unique(allTests$feature[allTests$FDR < taxRetestThresh])
  } else {
    toRetest <- unique(allTests$feature[allTests$Multitest.FDR < taxRetestThresh])
  }
  if (length(toRetest) < plotTopX) {toRetest <- unique(allTests$feature)[c(1:plotTopX)]}
  cnt <- 0
  if (length(toRetest) > 0 & doPlots) {
    for (t in toRetest) {
      #print(t)
      
      # -- DEBUG --
      #dataIn = phDFm2; saveFolder = F; feature=t; doSave=F; display="P"; onlyShowSig=T; yLab=NA; title=NA; na.rem=T
      #responseVar=phen; stTest="chisquared"; doPlots = T; nrTests=-1; xLab=NA; retPlot=T; ylim=NA
      # ---
      
      oneTest <- testOneFeaturePrevalence(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = T,retPlot = T,doSave = F,addNRs=T)
      g <- oneTest[[2]] + theme(legend.position = "none") + ylab(t) + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
      #print(g)
      cnt <- cnt + 1
      ggsave(plot=g,filename = paste0(outFolder,'/plots_prev/phen_',phenNR,'_plot_prev_',cnt,'.png'))
    }
  }
  
  # ============ ABUNDANCE OF PATHWAYS ===========
  print(' >>> testing pathways <<<')
  allTests <- data.frame()
  # > run pair-wise mann-u-whitney tests for relative abundance in CD vs HC
  toRetest <- c()
  for (t in colnames(phDFm2)[grep('PWY\\.',colnames(phDFm2))] ) {
    #print(t)
    oneTest <- testOneFeature(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = F,discardZeros = F)
    oneTest$pValue[is.na(oneTest$pValue)] <- 1.0
    oneTest$FDR <- p.adjust(oneTest$pValue)
    if (min(oneTest$FDR) < 0.05) {toRetest <- c(toRetest,t)}
    allTests <- rbind.data.frame(allTests,oneTest)
  }
  allTests <- allTests[order(allTests$FDR),]
  allTests$PHENOTYPE <- phen
  allTests$Multitest.FDR <- p.adjust(allTests$pValue)
  write.table(allTests,paste0(outFolder,'/tables/phen_',phenNR,'_pwyAbundance.csv'),sep=',',row.names = F)
  if (length(toRetest) < plotTopX) {toRetest <- unique(allTests$feature)[c(1:plotTopX)]}
  cnt <- 0
  if (length(toRetest) > 0 & doPlots) {
    for (t in toRetest) {
      print(t)
      oneTest <- testOneFeature(phDFm2,saveFolder = F,feature = t,responseVar = phen,doPlots = T,discardZeros = F,retPlot = T,doSave = F)
      g <- oneTest[[2]] + theme(legend.position = "none") + ylab(t) + ggtitle('') + xlab(phenNames[[phenNR]]) + theme(text = element_text(size=15))
      #print(g)
      cnt <- cnt + 1
      ggsave(plot=g,filename = paste0(outFolder,'/plots_PWYs/phen_',phenNR,'_plot_PWY_abund_',cnt,'.png'))
    }
  }
}

# do pretty (and not very useful) heatmaps
hmFolder <- paste0(outFolder,'/plots_summary')
unlink(paste0(hmFolder),recursive = T,force = T)
dir.create(paste0(hmFolder))

tblFolder <- paste0(outFolder,"/tables/")
# other options 
mT = T
FDRc = 0.005
ssp = ""

# DIVERSITIES
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("diversities"),nrToPlot = 10,outFile = paste0(hmFolder,"/heatmap_diversities_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,clusterPhenos = T,subsetPhenos = ssp,metric="statSig")
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("diversities"),nrToPlot = 10,outFile = paste0(hmFolder,"/heatmap_diversities_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,clusterPhenos = T,subsetPhenos = ssp,metric="meanEffectRatio")

# GENERA
#  > prevalence
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxPrevalences"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_g_prevalence_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="g__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="statSig")
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxPrevalences"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_g_prevalence_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="g__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="meanEffectRatio")
#  > abundance
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxAbundances"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_g_abundances_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="g__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="statSig")
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxAbundances"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_g_abundances_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="g__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="meanEffectRatio")

# SPECIES
#  > prevalences
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxPrevalences"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_s_prevalence_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="s__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="statSig")
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxPrevalences"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_s_prevalence_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="s__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="meanEffectRatio")
#  > abundances
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxAbundances"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_s_abundances_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="s__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="statSig")
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("taxAbundances"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_taxa_s_abundances_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,subsetFeatures="s__",clusterPhenos = T,subsetPhenos = ssp,
                          metric="meanEffectRatio")

# PATHWAYS
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("pwyAbundance"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_pwys_stats.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,clusterPhenos = T,subsetPhenos = ssp,addXextra=600)
p <- aggregatePlotHeatmap(tblFolder,doFeatures=c("pwyAbundance"),nrToPlot = 30,outFile = paste0(hmFolder,"/heatmap_pwys_effects.png"),
                          multiTest=mT,pVcutoff = FDRc,enrichmentCutoff = 1.005,clusterPhenos = T,subsetPhenos = ssp,addXextra=600,
                          metric="meanEffectRatio")


# LOGISTIC REGRESSION ANALYSIS
# ====================================
# ====================================
# load data (filtered - asin(sqrt) normalised - rescaled - corrected for technicals)
inDF <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor.csv',sep=',',header=T)
# clean names of Taxa
inDFf <- purgeMGNames(inDF)

#inDFprev <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled_PA.csv',sep=',',header=T)
#inDFprev <- purgeMGNames(inDFprev)

# load diversity data
inDFdiv <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_diversitiesonly.csv',sep=',',header=T)
colnames(inDFdiv)[[1]] <- "DAG3_sampleID"
# load phenotypes
phenos <- read.table('../DAG3_meta/v2/health/curated/DAG3_T1T2_Health.csv',sep=',',header=T,quote = '"')
colnames(phenos) <- tolower(colnames(phenos));colnames(phenos)[1] <- "DAG3_sampleID"

#  -- select phenotypes to test --
selCols <- colnames(phenos)[c(2:46,57,58,65:67,69:77,80,81)]
selCols <- c(selCols,"age","gender","BMI")
phenNames  <- selCols



# >> merge 
phDFm2 <- merge(inDFf,inDFdiv,by="DAG3_sampleID")
phDFm2 <- merge(phDFm2,phenos,by="DAG3_sampleID")
# > age limit
phDFm2 <- phDFm2[phDFm2$age > minAge,]

#  -- do correlation matrix --
# phenSel <- as.data.frame(phDFm2[,selCols])
# phenSel <- apply(phenSel,MARGIN = 2,as.character)
# phenSel[phenSel=="Y"] <- 1
# phenSel[phenSel=="N"] <- 0
# phenSel[phenSel=="M"] <- 1
# phenSel[phenSel=="F"] <- 0
# phenSel <- apply(phenSel,MARGIN = 2,as.numeric)
# phenCor <- cor(phenSel,use = "pairwise.complete.obs",method = "spearman")
# heatmap(phenCor)

# test test test
#toTest <- subsetMicrobiomeDF(phDFm2,)
#for (t in toTest) {
#}

phDFm2b <- phDFm2
phDFm2b$age <- (phDFm2b$age - mean(phDFm2b$age))/sd(phDFm2b$age)
phDFm2b$BMI <- (phDFm2b$BMI - mean(phDFm2b$BMI))/sd(phDFm2b$BMI)
phDFm2b[,selCols][is.na(phDFm2b[,selCols])] <- "N"

frm <- reformulate(response = "g__Faecalibacterium", termlabels = selCols)
#frm <- reformulate(response = "g__Faecalibacterium", termlabels = c("age","gender","BMI","ibd.any"))
test <- lm(frm, data=phDFm2b )
drk <- as.data.frame(summary(test)$coefficients)
drk <- drk[order(drk$Estimate),]
ggplot(phDFm2b,aes(x=BMI,y=g__Faecalibacterium,col=ibd.any)) + geom_point() + geom_smooth(method="lm")

