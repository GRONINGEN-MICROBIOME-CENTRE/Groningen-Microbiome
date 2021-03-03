# ===============================================================================
# ===============================================================================
# R. Gacesa, UMCG (2021)
#
# Microbiome data analysis student workshop
# ===============================================================================
# ===============================================================================

# ======================================================
# ======================================================
# ======================================================
#
# MICROBIOME WORKSHOP, parts 3 & 4:
#  PART III:
#  > examination of multiple microbiome samples
#  > calculation of Beta diversity (between-microbiome distance)
#  
#  PART IV:
#  > case-control analysis of microbiome taxa between healthy and disease
#  
#
# ======================================================
# ======================================================
# ======================================================

setwd('D:/UMCG/2021_microbiome_data_lecture/')

# load libraries
# ===========================
library(vegan)
library(ggplot2)
library(tidyr)

# =====================================================================================================================
# =====================================================================================================================
# PART 3: exploratory data analysis, multiple samples
# =====================================================================================================================
# =====================================================================================================================
# load microbiome data: 
inMB <- read.table('p2_microbiomes_data.csv',sep=',',header=T,stringsAsFactors = F)
# load metadata: 
inMeta <- read.table('p2_microbiomes_meta.csv',sep=',',header=T,stringsAsFactors = F)
# merge microbiome with metadata
inDFm <- merge(inMB,inMeta,by="ID")

# examine data
inMB[1:5,1:5]
head(inMeta)

# i) plot samples, phylum-level
inDFm2 <- inDFm
rownames(inDFm2) <- inDFm2$ID
inDFm2 <- inDFm2[,-grep('c__',colnames(inDFm2))]
inDFm2 <- inDFm2[,grep('p__',colnames(inDFm2))]
inDFm2 <- inDFm2[,colSums(inDFm2) != 0,]
inDFm2$ID <- row.names(inDFm2)
#   > shorten names of phyla for prettier plot
#     > iterate over phyla columns, shorten each
for (rn in grep('p__',colnames(inDFm2))) { colnames(inDFm2)[rn] <- strsplit(colnames(inDFm2)[rn],split = '\\.')[[1]][2] }

#  > data has to be in long format (each sample + phylum combination in separate row) for ggplot
inDFm2long <- gather(inDFm2, Phylum, Rel.Abundance,p__Actinobacteria:p__Euryarchaeota, factor_key=TRUE)
#   > examine long data
head(inDFm2long)
#   > plot phylum-level composition per sample
ggplot(inDFm2long,aes(x=ID,y=Rel.Abundance,col=Phylum,fill=Phylum)) + geom_col() + theme(legend.position = "bottom")

# plot variant 2: heatmap of phyla VS samples
inDFm3 <- inDFm2
inDFm3$ID <- NULL
# plain heatmap
heatmap(log(as.matrix(inDFm3)+0.00001),Rowv = NA, Colv = NA)
# with hierarchial clustering of taxa and samples
heatmap(log(as.matrix(inDFm3)+0.00001))

# ii) Diversity:
# ============================

# a) plot shannon diversity of species for each sample
colsSpecies <- grep('s__',colnames(inDFm))
inDFm$DIV.species <- vegan::diversity(inDFm[,colsSpecies],index = "shannon")
ggplot(inDFm, aes(y=DIV.species,x="All.Samples")) + geom_boxplot(outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

# b) plot it split by healthy / disease
table(inDFm$Diagnosis)
ggplot(inDFm, aes(y=DIV.species,x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

# c) is this difference statistically significant?
wilcox.test(inDFm$DIV.species[inDFm$Diagnosis=="Disease"],inDFm$DIV.species[inDFm$Diagnosis=="HC"])

# d) how do distributions look, are there any outliers / suspicious samples?
ggplot(inDFm, aes(x=DIV.species,col=Diagnosis)) + geom_density(size=2)

# e) what happens if we examine different diversity metric (inverse simpson)
# - why?
inDFm$DIV.species.is <- vegan::diversity(inDFm[,colsSpecies],index = "invsimpson")
ggplot(inDFm, aes(y=DIV.species.is,x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)
wilcox.test(inDFm$DIV.species.is[inDFm$Diagnosis=="Disease"],inDFm$DIV.species.is[inDFm$Diagnosis=="HC"])

# what about other metadata?
# > biomarker
ggplot(inDFm, aes(y=DIV.species,x=Biomarker)) + geom_point() + geom_smooth(method = "lm")
# > is it linked to disease?
ggplot(inDFm, aes(y=Biomarker,x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)
# > BMI
ggplot(inDFm, aes(y=DIV.species,x=BMI)) + geom_point() + geom_smooth(method = "lm")
ggplot(inDFm, aes(y=DIV.species,x=BMI,col=Diagnosis)) + geom_point() + geom_smooth(method = "lm")
# > AGE
ggplot(inDFm, aes(y=DIV.species,x=Age)) + geom_point() + geom_smooth(method = "lm")
ggplot(inDFm, aes(y=DIV.species,x=Age,col=Diagnosis)) + geom_point() + geom_smooth(method = "lm")

# iii) Beta diversity
# ============================

# 0) calculate B/C dissimilarity between first two samples
rownames(inDFm) <- inDFm$ID
colsSpecies <- grep('s__',colnames(inDFm))
bc12 <- vegan::vegdist(inDFm[1:2,colsSpecies],method = "bray")
print(bc12)


# a) calculate distance matrix of Bray-Curtis dissimilarity between samples on species level
# ======================================================
rownames(inDFm) <- inDFm$ID
colsSpecies <- grep('s__',colnames(inDFm))
bcSpec <- vegan::vegdist(inDFm[,colsSpecies],method = "bray")
# b) perform PCoA on it
pcoaMat <- cmdscale(bcSpec, eig = T, k=3)
variance <- head(eigenvals(pcoaMat)/sum(eigenvals(pcoaMat)))
var_pc1 <- as.integer(variance[1]*100)
var_pc2 <- as.integer(variance[2]*100)
var_pc3 <- as.integer(variance[3]*100)
# c) convert to dataframe for happiness and plotting
pcoaDF <- as.data.frame(pcoaMat$points)
# d) annotate dataframe
colnames(pcoaDF) <- c("PC1","PC2","PC3")
#    merge with metadata
pcoaDF <- merge(pcoaDF,inDFm,by='row.names')    

# e) plot first two principal coordinates
# ======================================================
centroids <- aggregate(cbind(PC1,PC2) ~ Diagnosis,pcoaDF,mean)

ggplot(pcoaDF,aes(x=PC1,y=PC2,col=Diagnosis)) + 
        geom_point(size=2) +
        theme_classic() + 
        geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC1,y=centroids$PC2),alpha=1) + 
        xlab(paste0('PC 1,',' (',var_pc1,'% variance)')) + 
        ylab(paste0('PC 2,',' (',var_pc2,'% variance)'))

# f) plot 2nd and 3rd principal coordinates
# ======================================================
centroids <- aggregate(cbind(PC2,PC3) ~ Diagnosis,pcoaDF,mean)
ggplot(pcoaDF,aes(x=PC2,y=PC3,col=Diagnosis)) + 
  geom_point(size=2) +
  theme_classic() + 
  geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC2,y=centroids$PC3),alpha=1) + 
  xlab(paste0('PC 2,',' (',var_pc2,'% variance)')) + 
  ylab(paste0('PC 3,',' (',var_pc3,'% variance)'))


# iv) quantifying the variance explained by disease variable
# (Permutational Multivariate Analysis of Variance Using Distance Matrices)
# ======================================================
adonis(bcSpec ~ inDFm$Diagnosis,permutations = 1000)


# v) advanced/optional: biplot: quantification of factors shaping PCoA plots
# ===================================================================
# a) select species and calculate bray-curtis distance matrix
rownames(inDFm) <- inDFm$ID
colsSpecies <- grep('s__',colnames(inDFm))
bcSpec <- vegan::vegdist(inDFm[,colsSpecies],method = "bray")
# b) perform PCoA on it
pcoaMat <- cmdscale(bcSpec, eig = T, k=3)
var_pc1 <- as.integer(variance[1]*100)
var_pc2 <- as.integer(variance[2]*100)
var_pc3 <- as.integer(variance[3]*100)
# c) convert to dataframe for happiness and plotting
pcoaDF <- as.data.frame(pcoaMat$points)
# d) annotate dataframe
colnames(pcoaDF) <- c("PC1","PC2","PC3")
#    merge with metadata
pcoaDF <- merge(pcoaDF,inDFm,by='row.names')
#    calculate centers of groups
centroids <- aggregate(cbind(PC1,PC2) ~ Diagnosis,pcoaDF,mean)

# e) calculate biplot weights
sp <- inDFm[,grep('s__',colnames(inDFm))]
colnames(sp)=sapply(strsplit(colnames(sp), ".", fixed=TRUE), tail, 1)
colnames(sp)=gsub(".*s__","",colnames(sp))
pcoaMat2 <- cmdscale(bcSpec)
en=envfit(pcoaMat2,inDFm[,c("Age","Diagnosis","Biomarker","BMI")], na.rm=T, permutations=10000)

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_factors = as.data.frame(scores(en, display = "factors"))
en_coords_merged <- rbind.data.frame(en_coord_cont,en_coord_factors[rownames(en_coord_factors) %in% c("GenderM","DiagnosisDisease"),])
coord = en_coord_cont/2

# f) plot results
ggplot(pcoaDF,aes(x=PC1,y=PC2,col=Diagnosis)) +
  geom_point(size=2) +
  theme_classic() +
  geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC1,y=centroids$PC2),alpha=1) +
  xlab(paste0('PC 1,',' (',var_pc1,'% variance)')) +
  ylab(paste0('PC 2,',' (',var_pc2,'% variance)')) +
  geom_segment(data=en_coords_merged, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), size =1, colour = "firebrick",arrow = arrow()) +
  geom_text(data = en_coords_merged, aes(x = Dim1, y = Dim2-0.05), colour = "firebrick", fontface = "bold", label = row.names(en_coords_merged)) + 
  theme(legend.position = "bottom")

# ==============================================================
# ==============================================================
# Part IV) compare microbiome taxa between healthy and disease
# ==============================================================
# ==============================================================
inDFm2 <- inDFm
# find all taxa (for comparison of healthy vs disease)
taxaC <- grep('^k__',colnames(inDFm))
# define dataframe for results
res <- NULL
# shorten names for convenience
for (taxC in taxaC) {
  taxName <- colnames(inDFm2)[taxC]
  taxNameS <- strsplit(colnames(inDFm2)[taxC],split = '\\.')[[1]]
  taxShortName = taxNameS[length(taxNameS)]
  colnames(inDFm2)[taxC] <- taxShortName
}
# test features one by one
for (taxC in taxaC) {
  # select column for test
  taxName <- colnames(inDFm2)[taxC]
  # do non-parametric test of healthy samples vs disease samples
  oneTest <- wilcox.test(inDFm2[[taxC]][inDFm2$Diagnosis=="HC"],
                         inDFm2[[taxC]][inDFm2$Diagnosis=="Disease"])
  # parse test results and add them to dataframe
  res <- rbind.data.frame(res,data.frame(taxon=taxName, # name of taxon
                                         abHC=mean(inDFm2[[taxC]][inDFm2$Diagnosis=="HC"]), # mean of healthy samples
                                         abD=mean(inDFm2[[taxC]][inDFm2$Diagnosis=="Disease"]), # mean of disease samples
                                         pV=oneTest$p.value)) # test p-value
}
# order results by p-value
res <- res[order(res$pV),]
# remove rows that include NAs - these appear when for 0 vs 0 comparisons
res <- res[complete.cases(res),]
# do multiple testing correction for p-values
res$FDR <- p.adjust(res$pV)
# show results wiht false discovery rate < 20%
res[res$FDR < 0.2,]

# plot example:
#  > raw
ggplot(inDFm2,aes(y=g__Escherichia,
                 x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

# > log-transformed relative abundances: note we add arbitrary small number* to abundances to not lose 0s due to log transform
# *: 1*10^-6 is approximately our detection threshold for relative abundance so it is ~equal to 0 in our data
ggplot(inDFm2,aes(y=log(g__Escherichia+0.000001),
                  x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

