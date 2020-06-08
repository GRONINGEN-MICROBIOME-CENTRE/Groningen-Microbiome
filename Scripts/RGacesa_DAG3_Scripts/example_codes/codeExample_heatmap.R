# ===================================
library(ggplot2)
library(dplyr)
library(pheatmap)

# set WD and load libraries
setwd('D:/UMCG/Codes')
#  > plotting scripts live in R_Misc.R
source('../myLibs_v3/R_Misc.R')
#  > other microbiome-related scripts live in R_Microbiome_scripts.R
#    NOTE: these are required
source('../myLibs_v3/R_Microbiome_scripts.R')

# # DEBUG
# == this is for internal debugging, leave it commented unless debugging the code ==
# nrFeaturesToPlot = 20
# featuresToPlot = "G"
# clusterPhenos = T
# clusterFeatures = T
# modelToUse = "FullModel"
# sigStat = "Padj"
# sigStatNominal = "Pvalue"
# sigValue = 0.05
# sigNominal = 0.05
# statToPlot = "Padj"
# addText = "Direction"
# plotValuesSig = "all"
# plotTextSig = "all"
# directionValue = "effect.mu"
# transformValues=T
# logLimit = 10.0
# sigAll = 0.9
# textSize=20


# EXAMPLE Plot 1: Top X bugs associated with mental disorders
# =================================================
dataPath = 'D:/UMCG/DAG3_phenotype_associations/pheno.Taxa.results.03May.txt'

inData <- read.table(dataPath,sep='\t',header=T,stringsAsFactors = F)
# recode names to nice names
inData$phenotype <- gsub("ANTHRO\\.BMI","BMI",inData$phenotype)
inData$phenotype[inData$phenotype =="MED.DISEASES.None.No.Diseases"] <- "No_Disease"
inData$phenotype <- gsub("MED.DISEASES.","",inData$phenotype)
inData$phenotype <- gsub("IBS.Any","IBS",inData$phenotype)
inData$phenotype <- gsub("Mental.","",inData$phenotype)
inData$phenotype <- gsub("Neurological.","",inData$phenotype)
inData$phenotype <- gsub("Rome3_","",inData$phenotype)
inData$phenotype <- gsub("Gastrointestinal.","",inData$phenotype)
inData$phenotype <- gsub("^Any","All Psychiatric Disorders",inData$phenotype)
inData$phenotype <- gsub("^Social\\.phobia","Social Phobia",inData$phenotype)
inData$phenotype <- gsub("^Fatigue\\.Chronic","Chronic Fatigue Syndrome",inData$phenotype)
inData$phenotype <- gsub("^Other\\.anxiety","Anxiety Disorders",inData$phenotype)
inData$phenotype <- gsub('Panic\\.disorder',"Panic Disorder",inData$phenotype)

phenosToCheck <- c("BMI",
                   "All Psychiatric Disorders","Depression","Social Phobia",
                   "IBS","Fibromyalgia","Chronic Fatigue Syndrome",
                   "Anxiety Disorders","Panic Disorder","No_Disease"
                   )

# variant 1: default values:
#  > colors correspond to -log(p-value) of full model, direction is plotted based on mu
#  > always plots colors (even for non-significant results)
#  > puts direction notes (+/-) only on significant results (@FDR < 0.05)
#  > phenotypes & taxa are both clustered
# ============================================
nrFeatures = 10
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures)
save_pheatmap_png(phm,filename = 'happy_plot_1.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 1b: few phenotypes & taxa to demonstrate plot size scaling
#    Note: save_pheatmap_png is used to scale plot size due to limitations of pheatmap function
#          this works well for numbers of features and phenotypes, but requires some manual adjustment
#          if phenotype names or taxa names become very long (as with species - see Variant 1d)
# ================================================================================================
nrFeatures = 5
phenosToCheck2 <- phenosToCheck[1:5]
nrPhenos = length(phenosToCheck2)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck2,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures)
save_pheatmap_png(phm,filename = 'happy_plot_1b.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 1c: more taxa to demonstrate plot size scaling
# =============================================================
nrFeatures = 30
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures)
save_pheatmap_png(phm,filename = 'happy_plot_1c.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)


# variant 1d: species instead of genera
# NOTE: we need to increase height of plot as species have longer names, 
#       this is not automagical (and might never be)
# ============================================
nrFeatures = 30
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "S",
                                   nrFeaturesToPlot = nrFeatures)
save_pheatmap_png(phm,filename = 'happy_plot_1d.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1200,res = 350)

# variant 1e: phyla instead of genera
# NOTE: function is resistant to putting more features in then there are in total, in this case, it will just plot everything there is
#       this does make plot wider then necessary though, so don't do it
# ===========================================================================
nrFeatures = 20 # <- more then total number of phyla
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "P",
                                   nrFeaturesToPlot = nrFeatures)
save_pheatmap_png(phm,filename = 'happy_plot_1e.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)


# variant 2: default values, but plots significance instead of direction
#  > colors correspond to -log(p-value) of full model, direction is plotted based on mu
#  > always plots colors (even for non-significant results)
#  > puts * note only on significant results (@FDR < 0.05)
#  > phenotypes & taxa are both clustered
# ============================================
nrFeatures = 10
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures,addText="Significance")
save_pheatmap_png(phm,filename = 'happy_plot_2.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 4: plot Chisq of full model
#
# NOTES:
#  > color gradient here can get messy due to assymetric values, increasing color number (nrColors) helps
# ============================================================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures,statToPlot = "Chisq",nrColors = 21)
save_pheatmap_png(phm,filename = 'happy_plot_4.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 5: plot R2 of full model
#
# NOTES:
#  > color gradient here can get messy due to assymetric values, increasing color number (nrColors) helps
# ============================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures,statToPlot = "R2",nrColors = 21)
save_pheatmap_png(phm,filename = 'happy_plot_5.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 6: plot effect size of full model
#
# NOTES:
#  > this one throws an exception and dies with controlled error:
# "ERROR: stat "effect" + model "FullModel" does not exist in data, try different modelToUse or statToPlot!"
#   >> this is working as intended as we don't have effect size of full model
# ============================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",
                                   nrFeaturesToPlot = nrFeatures,statToPlot = "effect",nrColors = 21)
#save_pheatmap_png(phm,filename = 'happy_plot_5.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# PLOTS THAT REQUIRE FURTHER TESTING
# > NOTES:
#   - plots for modelToUse != "FullModel" are currently experimental (there might be bugs in there, further testing is required)
# ===========================================================

# variant 7: plot mu instead of full model
# ===================================================================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",modelToUse = "mu",
                                   nrFeaturesToPlot = nrFeatures,nrColors = 11)
save_pheatmap_png(phm,filename = 'happy_plot_7.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# variant 7b: plot nu instead of full model
# ===================================================================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",modelToUse = "nu",
                                   nrFeaturesToPlot = nrFeatures,nrColors = 11)
save_pheatmap_png(phm,filename = 'happy_plot_7b.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)

# WORK IN PROGRESS / EXPERIMENTAL PLOTS THAT MIGHT NOT BE WORKING CORRECTLY
# > NOTES:
#   - plots for effectsize are currently experimental (there might be bugs in there, further testing is required)
# ===============================================================================================================

# variant 8: plot effect size of mu
# NOTES:
#  > known bug/feature: very small positive values can be significant, so we get white colors 
# ===================================================================================
nrFeatures = 20
nrPhenos = length(phenosToCheck)
phm <- plotTaxaAssociationsHeatmap(phenosToPlot = phenosToCheck,inData = inData,featuresToPlot = "G",modelToUse = "mu",
                                   nrFeaturesToPlot = nrFeatures,statToPlot = "effect",nrColors = 21)
save_pheatmap_png(phm,filename = 'happy_plot_8.png',width = 1000+70*nrFeatures,height = 70*nrPhenos+1000,res = 350)