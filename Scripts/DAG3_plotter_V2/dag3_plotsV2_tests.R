# =============================
# testing new DAG3 plotter
# =============================
# > set wd
setwd('D:/UMCG/DAG3_statistics_V2/')
#  > plotting scripts live in R_Misc.R
source('../myLibs_v3/R_Misc.R')
#  > other microbiome-related scripts live in R_Microbiome_scripts.R
#    NOTE: these are required
source('../myLibs_v3/R_Microbiome_scripts.R')

#  > LOAD DATA
inData <- read.table('results_DAG3_June17/taxa_Jun17.txt',sep='\t',header=T,stringsAsFactors = F)


# # debug
# # ===========================
# phenosToPlot = c("ANTHRO.BMI","ANTHRO.Sex","ANTHRO.AGE","MED.HEALTH.RAND.Health.General","MED.HEALTH.RAND.Health.Change.1y","MED.DISEASES.None.No.Diseases")
# # inData = inData
# # nrFeaturesToPlot = 20
# # featuresToPlot = "G"
# # clusterPhenos = F
# # clusterFeatures = T
# # sigStat = "PadjBH"
# # sigStatNominal = "Pvalue"
# # sigValue = 0.05 
# # sigNominal = 0.05 
# # statToPlot = "effect.size" # Pvalue, PadjBH, PadjBonferroni, R2, F.stat, effect.size
# # addText = "Direction"
# # plotValuesSig = "padjSig" # all, nomSig, padjSig
# # plotTextSig = "padjSig" 
# # transformValues=T
# # logLimit = 10.0
# # sigAll = 0.9 
# # textSize=15 
# # colLow = "Orange" 
# # colHigh = "Darkblue"
# # colTextSize = 12
# # rowTextSize = 12
# # legendTextSize = 10
# # cellWidth=15 
# # cellHeight=15
# # nrColors = 13
# # flipCoords = T
# # doRefactor = F
# 
# p1 <- plotAssociationsDag3HeatmapV2(inData = inData,phenosToPlot = phenosToPlot,nrFeaturesToPlot=20,plotValuesSig = "nomSig",
#                               statToPlot="effect.size",clusterPhenos = F,clusterFeatures = T,doRefactor = T,nrColors = 59,
#                               refactorIncludeAllClasses = F)



# load Alex's last version of results for DAG3
# ==================================================
inData <- read.table('results_DAG3_June17/taxa_Jun17.txt',sep='\t',header=T,stringsAsFactors = F)

# define what to plot:
phenosToPlot = c("ANTHRO.BMI","ANTHRO.Sex","ANTHRO.AGE","MED.HEALTH.RAND.Health.General","MED.HEALTH.RAND.Health.Change.1y","MED.DISEASES.None.No.Diseases")
phenosToPlot2 = c("ANTHRO.BMI","ANTHRO.Sex","ANTHRO.AGE","MED.HEALTH.RAND.Health.General","MED.HEALTH.RAND.Health.Change.1y","MED.DISEASES.None.No.Diseases",
                  unique(inData$phenotype)[grep('MED\\.DISEASES\\.Gastrointestinal\\.',unique(inData$phenotype))])

# Variant 1: default values:
# ============================================================================
#  > colors correspond to -log(p-value) of model, direction is plotted based on effect size
#  > plots colors for nominal significance
#  > puts direction notes (+/-) only on significant results (@FDR < 0.05)
#  > phenotypes & taxa are both clustered
#  > uses "model average" effect size for categorical variables
# ========================================================================
nrFeatures = 10
nrPhenos = length(phenosToPlot)
phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures) # how many top shared features to plot
                                     
save_pheatmap_png(phm,filename = 'happy_plot_1.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# variant 1b: less taxa and phenotypes to demonstrate plot size scaling
#    Note: save_pheatmap_png is used to scale plot size due to limitations of pheatmap function
#          this works well for numbers of features and phenotypes, but requires some manual adjustment
#          if phenotype names or taxa names become very long (as with species - see Variant 1d)
# ================================================================================================
nrFeatures = 5
nrPhenos = length(phenosToPlot)
phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot[1:3], # what is plotted
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures) # how many top shared features to plot

save_pheatmap_png(phm,filename = 'happy_plot_1b.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# variant 1c: more taxa & more phenotypes to demonstrate plot size scaling
# =============================================================
nrFeatures = 50
nrPhenos = length(phenosToPlot2)
phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot2, # what is plotted
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures) # how many top shared features to plot

save_pheatmap_png(phm,filename = 'happy_plot_1c.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# variant 1e: phyla instead of genera
# NOTE: function is resistant to putting more features in then there are in total, in this case, it will just plot everything there is
#       this does make plot bigger then necessary though
# ===========================================================================
nrFeatures = 20 # <- more then total number of phyla
nrPhenos = length(phenosToPlot)
phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "P", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures) # how many top shared features to plot

save_pheatmap_png(phm,filename = 'happy_plot_1e.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# =========================================================================================================
# Variant 2: default values, but plots significance instead of direction
# =========================================================================================================
#  > colors correspond to -log(p-value) of full model, direction is plotted based on mu
#  > always plots colors (even for non-significant results)
#  > puts * note only on significant results (@FDR < 0.05)
#  > phenotypes & taxa are both clustered
# ============================================
nrFeatures = 10
nrPhenos = length(phenosToPlot)

phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     addText="Significance") # put * instead of +/- on significant results (FDR significant)

save_pheatmap_png(phm,filename = 'happy_plot_2.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# =========================================================================================================
# Variant 4: plot effect size instead of p-values
# =========================================================================================================
# NOTES:
#  > color gradient here can get messy due to small but significant effect sizes, increasing color number (nrColors) helps
#  > in case of categorical variables, it uses effect size of non-control for binary (it is always 0 for control)
#     if variable has > 2 classes, it plots effect.size.asInteger value
# =========================================================================================================
nrFeatures = 20
nrPhenos = length(phenosToPlot)

phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51) # increase number of colors for gradient

save_pheatmap_png(phm,filename = 'happy_plot_4.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# =========================================================================================================
# Variant 4/b: plot R2 instead of p-values
# =========================================================================================================
# NOTES:
#  > color gradient here can get messy due to small but significant effect sizes, increasing color number (nrColors) helps
#  > in case of categorical variables, it uses effect size of non-control for binary (it is always 0 for control)
#     if variable has > 2 classes, it plots effect.size.asInteger value
# =========================================================================================================
nrFeatures = 20
nrPhenos = length(phenosToPlot)

phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "R2", # plot effect size
                                     nrColors = 51) # increase number of colors for gradient

save_pheatmap_png(phm,filename = 'happy_plot_4b.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)


# =========================================================================================================
# Variant 4/c: plot F statistic instead of p-values
# =========================================================================================================
# NOTES:
#  > color gradient here can get messy due to small but significant effect sizes, increasing color number (nrColors) helps
#  > in case of categorical variables, it uses effect size of non-control for binary (it is always 0 for control)
#     if variable has > 2 classes, it plots effect.size.asInteger value
# =========================================================================================================
nrFeatures = 20
nrPhenos = length(phenosToPlot)

phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "F.stat", # plot effect size
                                     nrColors = 51) # increase number of colors for gradient

save_pheatmap_png(phm,filename = 'happy_plot_4c.png',width = 1000+70*nrPhenos,height = 70*nrFeatures+1000,res = 300)

# ==============================================================================================================================
# Variant 5: plot effect size instead of p-values, but expand multi-class variables and plot each case-control effect size
# =============================================================================================================================
# NOTES:
#  > color gradient here can get messy due to small but significant effect sizes, increasing color number (nrColors) helps
#  > if doing this, it is good idea to turn off clustering of phenotypes (otherwise columns corresponding to same phenotype can get shuffled)
#  > in case of categorical variables, scripts discards control and transforms classes into multiple columns
#    (ex: MED.HEALTH.RAND.Health.General (levels: poor:mediocre:good:very good:excellent) becomes:
#           > column 1: MED.HEALTH.RAND.Health.General.C1.mediocre
#           > column 2: MED.HEALTH.RAND.Health.General.C2.good
#           > column 3: MED.HEALTH.RAND.Health.General.C3.very good
#           > column 4: MED.HEALTH.RAND.Health.General.C4.excellent
#           > poor value (first class) is not plotted as it is control and all effect sizes are 0 for it
#  > this takes longer as code needs to grind through big table
# =========================================================================================================
nrFeatures = 20
nrPhenos = length(phenosToPlot)

phm <- plotAssociationsDag3HeatmapV2(inData = inData, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     featuresToPlot = "G", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     statToPlot = "effect.size", # plot effect size
                                     nrColors = 51, # increase number of colors for gradient
                                     doRefactor = T,   # transform multi-class variables into multiple columns
                                     clusterPhenos = F # do not cluster phenotypes
                                     ) 
save_pheatmap_png(phm,filename = 'happy_plot_5.png',width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)


