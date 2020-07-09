# =============================
#
# 09/07/2020, R.Gacesa (UMCG)
# DAG3 association plots
#
# =============================
trimPhenos <- function(inPh) {
  inPh <- gsub('META\\.','',inPh)
  inPh <- gsub('ANTHRO\\.','',inPh)
  inPh <- gsub('SOCIOEC\\.','',inPh)
  inPh <- gsub('MED\\.','',inPh)
  inPh <- gsub('MEDS\\.','',inPh)
  inPh <- gsub('SOCIOECONOMIC\\.','',inPh)
  inPh <- gsub('INCOME\\.','',inPh)
  inPh <- gsub('WORK\\.','',inPh)
  inPh <- gsub('DISEASES\\.','',inPh)
  inPh <- gsub('Cardiovascular\\.','',inPh)
  inPh <- gsub('BLOOD\\.','',inPh)
  inPh <- gsub('EXP\\.','',inPh)
  inPh <- gsub('EARLYLIFE\\.','',inPh)
  inPh <- gsub('DIET\\.','',inPh)
  inPh <- gsub('Scores\\.','',inPh)
  inPh <- gsub('BIRTH\\.','',inPh)
  inPh <- gsub('DEMOGRAPHICS\\.','',inPh)
  inPh <- gsub('Scores\\.','',inPh)
  inPh <- gsub('HEALTH\\.','',inPh)
  inPh <- gsub('RAND\\.','',inPh)
  inPh <- gsub('Autoimmune\\.','',inPh)
  inPh
}

# > set wd
setwd('D:/UMCG/DAG3_statistics_V2/')
#  > plotting scripts live in R_Misc.R
source('../myLibs_v3/R_Misc.R')
#  > other microbiome-related scripts live in R_Microbiome_scripts.R
#    NOTE: these are required
source('../myLibs_v3/R_Microbiome_scripts.R')

# load Alex's last version of results for DAG3
# ==================================================
inData <- read.table('results_DAG3_June17/taxa_Jun17.txt',sep='\t',header=T,stringsAsFactors = F)

# plots per adonis group
# ===============================================
inAdonis <- read.table('D:/UMCG/DAG3_results_v26/adonis/adonis_taxa_v26_rdy_grouped.csv',sep=',',header=T,stringsAsFactors = T)
inAdonis$Significant <- as.character(inAdonis$Significant)
inAdonis$Significant[inAdonis$FDR.BH > 0.05 & inAdonis$pval < 0.05] <- "Nominal"
#inAdonis$Variable <- factor(as.character(inAdonis$Variable),levels = inAdonis$Variable[order(inAdonis$R2,decreasing = F)])
inAdonis$Variable2 <- factor(as.character(inAdonis$Variable2),levels = inAdonis$Variable2[order(inAdonis$R2,decreasing = F)])

inDF <- inData
inDF$phenotype <- trimPhenos(inDF$phenotype)

# plot config
#cntplot <- 0
# do plots per group
for (grp in unique(inAdonis$Group.smaller)) {
#  cntplot <- cntplot + 1
  # GENERA
  phenosToPlot <- inAdonis$Var[inAdonis$Group.smaller==grp]
  phenosToPlot <- trimPhenos(phenosToPlot)
  
  print(paste0(' >> PLOTTING ',grp,' vs Genera'))
  nrFeatures = 30
  nrPhenos = length(phenosToPlot)
  phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                       phenosToPlot = phenosToPlot, # what is plotted
                                       featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                       nrFeaturesToPlot = nrFeatures,
                                       clusterPhenos = T,clusterFeatures = T) # how many top shared features to plot
  save_pheatmap_png(phm,filename = paste0('D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_',grp,'_Genera.png'),
                    width = 1200+70*nrPhenos,height = 70*nrFeatures+1100,res = 300)
  
  # SPECIES
  print(paste0(' >> PLOTTING ',grp,' vs Species'))
  nrFeatures = 50
  nrPhenos = length(phenosToPlot)
  phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                       phenosToPlot = phenosToPlot, # what is plotted
                                       featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                       nrFeaturesToPlot = nrFeatures,
                                       clusterPhenos = T,clusterFeatures = T) # how many top shared features to plot
  save_pheatmap_png(phm,filename = paste0('D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_',grp,'_Species.png'),
                    width = 1200+70*nrPhenos,height = 70*nrFeatures+1100,res = 300)
}

# =======================================
# "healthy microbiome" patterns
# =======================================
phenosToPlot <- c("MED.DISEASES.Cancer.AnyNonBasal",
                  "MED.DISEASES.Cancer.Colon",
                  "MED.DISEASES.Gastrointestinal.FGID.Bloating",
                  "MED.DISEASES.Gastrointestinal.FGID.Diarrhea",
                  "MED.DISEASES.Gastrointestinal.FGID.Constipation",
                  "MED.DISEASES.Gastrointestinal.Rome3_IBS.Any",
                  "MED.DISEASES.Hepatologic.Gallstones",
                  "MED.DISEASES.Mental.Depression",
                  "MED.INDICES.FattyLiverIndex.T1.Class",
                  "MED.INDICES.NAFLD.T1.CLASS",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.Celiac",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.UC",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD",
                  "MED.DISEASES.Endocrine.DiabetesT2",
                  "MED.DISEASES.Cardiovascular.Colesterol.high",
                  "MED.DISEASES.Cardiovascular.Hypertension",
                  "MED.HEALTH.RAND.Health.General",
                  "MED.HEALTH.RAND.Health.Change.1y",
                  "MED.DISEASES.None.No.Diseases"
                  )
phenos2 <- c("MED.HEALTH.RAND.Health.General",
             "MED.HEALTH.RAND.Health.Change.1y",
             "MED.DISEASES.None.No.Diseases")

phenosToPlot <- trimPhenos(phenosToPlot)
inDF$phenotype <- trimPhenos(inDF$phenotype)
phenosToPlot2 <- trimPhenos(phenos2)

nrFeatures = 30
print(paste0(' >> PLOTTING Healthy microbiome vs Genera'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot2, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51, # increase number of colors for gradient
                                     doRefactor = T,   # transform multi-class variables into multiple columns
                                     clusterPhenos = F,
                                     refactorIncludeAllClasses = T) # do not cluster phenotypes

save_pheatmap_png(phm,filename = 'D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_health_metrics_effectsize_Genera.png',width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)

nrFeatures = 50
print(paste0(' >> PLOTTING Healthy microbiome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot2, # what is plotted
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51, # increase number of colors for gradient
                                     doRefactor = T,   # transform multi-class variables into multiple columns
                                     clusterPhenos = F,
                                     refactorIncludeAllClasses = T) # do not cluster phenotypes

save_pheatmap_png(phm,filename = 'D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_health_metrics_effectsize_Species.png',width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)

# === variant 2 (diseases & metrics) ===
nrFeatures = 30
print(paste0(' >> PLOTTING Healthy microbiome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51, # increase number of colors for gradient
                                     doRefactor = T,   # transform multi-class variables into multiple columns
                                     clusterPhenos = F) # do not cluster phenotypes

save_pheatmap_png(phm,filename = 'D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_health_and_disease_effectsize_Genera.png',width = 1400+70*(nrPhenos*2+7),
                  height = 70*nrFeatures+1800,res = 300)


print(paste0(' >> PLOTTING Healthy microbiome vs Species'))
nrFeatures = 50
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51, # increase number of colors for gradient
                                     doRefactor = T,   # transform multi-class variables into multiple columns
                                     clusterPhenos = F) # do not cluster phenotypes

save_pheatmap_png(phm,filename = 'D:/UMCG/DAG3_results_v26/plots/associations/taxonomy_health_and_disease_effectsize_Species.png',width = 1400+70*(nrPhenos*2+7),
                  height = 70*nrFeatures+1800,res = 300)
