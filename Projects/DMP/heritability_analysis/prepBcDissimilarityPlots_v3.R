# ==================================
# By: Weersma Group, UMCG (2020)
# ==================================

# ============================================================
# ============================================================
#  >> Pair-wise Bray-Curtis dissimilarity plots
# ============================================================
# ============================================================

setwd('# placeholder, replace with the folder with this file')
source('../myLibs_v3/R_Misc.R')
source('../myLibs_v3/R_Microbiome_scripts.R')
source('../myLibs_v3/R_ML_scripts_v4.R')

library(ggplot2)
library(dplyr)
library(plyr)
library(vegan)

# input file: pair-wise Bray-Curtis dissimilarity between participants (generated using vegan:vegdist function with method="bray")
# NOTE: this file contains participant personal information and is therefore not provided on the github repo,
#       it is avaliable on request, subject to approval by UMCG ethics comittee and Lifelines Biobank

allResAnnot <- read.table('./data/df_allpairs_annotated.csv',sep=',',header = T,stringsAsFactors = F)

toPlotVar <- c("BC_Spec","BC_PWY","BC_VF","BC_CARD") 
toPlotYs <- c("Bray-Curtis distance of Species",
              "Bray-Curtis distance of Pathways",
              "Bray-Curtis distance of Virulence factors",
              "Bray-Curtis distance of Antibiotic resistance")

cbPalette <- c("#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
txtSize = 16

# >>> children vs siblings vs rnd vs partners, cohabitating only
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("PARTNERS","RND.PAIR","PARENT_CHILD","SIBLINGS"),]
dfToPlot <- dfToPlot[dfToPlot$COHAB | dfToPlot$RELATIONSHIP.0=="RND.PAIR",]
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","PARTNERS","PARENT_CHILD","SIBLINGS"))
dfToPlot$RELATIONSHIP.0 <- mapvalues(dfToPlot$RELATIONSHIP.0 , from = c("PARENT_CHILD", "SIBLINGS"), to = c("PAR_CH", "SIBL"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005,ylim = c(0.0,1.2))
  g <- t[[2]] + xlab("Relationship") + ggtitle('') + ylab(toPlotYs[c]) + theme_classic() +
     scale_color_manual(values = cbPalette) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + ylim(0.0,1.2)
  print(g)
  ggsave(plot=g,filename = paste0('./Plots/microbiome_similarity_relationship_0_cohab_',tpv,'.png'),width =5.5,height=6,dpi = 320)
}

# >>> children vs siblings vs rnd vs partners, non-cohabitating only
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","PARENT_CHILD","SIBLINGS"),]
dfToPlot <- dfToPlot[!dfToPlot$COHAB | dfToPlot$RELATIONSHIP.0=="RND.PAIR" | dfToPlot$RELATIONSHIP.0=="PARTNERS",]
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","PARENT_CHILD","SIBLINGS"))
dfToPlot$RELATIONSHIP.0 <- mapvalues(dfToPlot$RELATIONSHIP.0 , from = c("PARENT_CHILD", "SIBLINGS"), to = c("PAR_CH", "SIBL"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',ylim = c(0.0,1),
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005)
  
  g <- t[[2]] + xlab("Relationship (cohabitating pairs)") + ggtitle('') + ylab(toPlotYs[c]) + theme_classic() +
     scale_color_manual(values = cbPalette) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + ylim(c(0.0,1.0))
  print(g)
  ggsave(plot=g,filename = paste0('./Plots/heritability/microbiome_similarity_relationship_0_noncohab_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}

# >>> children (cohab) vs children (sep) vs random
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","PARENT_CHILD"),]
dfToPlot$RELATIONSHIP.0[dfToPlot$RELATIONSHIP.0=="PARENT_CHILD" & dfToPlot$COHAB] <- "PAR_CH_COH"
dfToPlot$RELATIONSHIP.0[dfToPlot$RELATIONSHIP.0=="PARENT_CHILD" & !dfToPlot$COHAB] <- "PAR_CH_SEP"
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005)
  g <- t[[2]] + xlab("Relationship (cohabitating pairs)") + ggtitle('') + ylab(toPlotYs[c]) + theme_classic() +
     scale_color_manual(values = cbPalette) + ylim(0.2,1) + theme(legend.position="none") + theme(text = element_text(size=txtSize))
  print(g)
  ggsave(plot=g,filename = paste0('./Plots/heritability/microbiome_similarity_relationship_parchild_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}

# >>> siblings (cohab) vs siblings (sep) vs random
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","SIBLINGS"),]
dfToPlot$RELATIONSHIP.0[dfToPlot$RELATIONSHIP.0=="SIBLINGS" & dfToPlot$COHAB] <- "SIB.COH"
dfToPlot$RELATIONSHIP.0[dfToPlot$RELATIONSHIP.0=="SIBLINGS" & !dfToPlot$COHAB] <- "SIB.SEP"
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("SIB.COH","SIB.SEP","RND.PAIR"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005)
  g <- t[[2]] + xlab("Relationship (cohabitating pairs)") + 
    ggtitle('') + ylab(toPlotYs[c]) + theme_classic() + scale_color_manual(values = cbPalette) + ylim(0.0,1) + theme(legend.position="none") + theme(text = element_text(size=15))
  print(g)
  ggsave(plot=g,filename = paste0('D:/UMCG/DAG3_results_v26/plots/heritability/microbiome_similarity_relationship_siblings_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}

# >>> 1st-degree (cohab) vs 1st-degree (sep) vs random
cbPalette2 <- c("#E69F00", "#555555","#D55E00","#56B4E9", "#009E73", "#CC79A7","#F0E442","#0072B2")

dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","SIBLINGS","PARENT_CHILD"),]
dfToPlot$RELATIONSHIP.0[(dfToPlot$RELATIONSHIP.0=="SIBLINGS" | dfToPlot$RELATIONSHIP.0=="PARENT_CHILD") & dfToPlot$COHAB] <- "1stDEG.COH"
dfToPlot$RELATIONSHIP.0[(dfToPlot$RELATIONSHIP.0=="SIBLINGS" | dfToPlot$RELATIONSHIP.0=="PARENT_CHILD") & !dfToPlot$COHAB] <- "1stDEG.SEP"
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","1stDEG.SEP","1stDEG.COH"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display='',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.05,ylim = c(0.0,1))
  g <- t[[2]] + xlab("Relationship (cohabitating vs non-cohabitating pairs)") + 
    ggtitle('') + ylab(toPlotYs[c]) + theme_classic() + scale_color_manual(values = cbPalette2) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + ylim(0.0,1)
  print(g)
  ggsave(plot=g,filename = paste0('./Plots/microbiome_similarity_relationship_1stDeg_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}