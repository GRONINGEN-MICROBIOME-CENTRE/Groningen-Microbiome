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
# MICROBIOME WORKSHOP
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
# PART 1: exploratory data analysis, one sample
# =====================================================================================================================
# =====================================================================================================================

# > 0) load microbiome data (sample02)
# ======================================
inDF <- read.table(file = 'sample02_metaphlan.txt',sep='\t',header=F)
colnames(inDF) <- c("Taxon","TaxID","Rel.Abundance","Coverage","Aligned.Reads")
head(inDF)

## > i) find relative abundances of unknown & kingdoms (archaea, bacteria, viruses)
# ======================================================================================================================
##  > we grep for taxon names starting with k__, including any number of chars and |; and then we remove it
inDFk <- inDF[-grep('^k__.*\\|',inDF$Taxon),]

# > ii) find relative abundances of phyla and unknown
# ======================================================================================================================
## two steps: 
## 1)  we grep for taxon names starting with k__, including any number of chars and |c__; and then we remove it
##     this leaves all kingdoms, unknown and all phyla
inDFp <- inDF[-grep('^k__.*\\|c__',inDF$Taxon),]
## 2) then remove kingdoms by grepping for k__<anything>|p__ OR "UNKNOWN"
inDFp <- inDFp[grep('UNKNOWN|^k__.*\\|p__',inDFp$Taxon),]

# > iii) extract all species, remove unclassified reads and re-calculate relative abundances
# ======================================================================================================================
#  1) extract species (by grepping for |s__); Note: this also removes UNKNOWN 
inDFf <- inDF[grep('\\|s__',inDF$Taxon),]
#  2) re-normalize to 1 by dividing all values by sum of species
inDFf$Rel.Abundance <- inDFf$Rel.Abundance/sum(inDFf$Rel.Abundance)

# > iv) plot phylum-level composition of microbiome sample
# ======================================================================================================================
inDFp <- inDF[-grep('^k__.*\\|c__',inDF$Taxon),]
inDFp <- inDFp[grep('UNKNOWN|^k__.*\\|p__',inDFp$Taxon),]
inDFp$Taxon <- factor(inDFp$Taxon,levels = inDFp$Taxon[order(inDFp$Rel.Abundance)])
ggplot(inDFp,aes(x=Taxon,y=Rel.Abundance,col=Taxon,fill=Taxon)) + geom_col() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) 
ggplot(inDFp,aes(x="",y=Rel.Abundance,col=Taxon,fill=Taxon)) + geom_col() + coord_polar("y", start=0)

# > v) remove "unknown taxa" & re-normalize, then plot top-ten most abundant species and rest of data grouped as "Other"
# ======================================================================================================================
#  1) extract species (by grepping for |s__); Note: this also removes UNKNOWN 
inDFs <- inDF[grep('\\|s__',inDF$Taxon),]
inDFs$Taxon <- as.character(inDFs$Taxon)
#  2) re-normalize to 1 by dividing all values by sum of species
inDFs$Rel.Abundance <- inDFs$Rel.Abundance/sum(inDFs$Rel.Abundance)
#  3) sort
inDFs <- inDFs[order(inDFs$Rel.Abundance,decreasing = T),]
inDFs2 <- inDFs[1:10,]
inDFs3 <- rbind.data.frame(inDFs2,
                 data.frame(Taxon="OTHER",TaxID=-1,Rel.Abundance=sum(inDFs$Rel.Abundance[11:nrow(inDFs)]),Coverage=-1,Aligned.Reads=-1))
#  4) clean-up species names
for (rr in c(1:nrow(inDFs3))) {
     tsplit <- strsplit(inDFs3$Taxon[rr],"\\|")[[1]]
     inDFs3$Taxon[rr] <- tsplit[length(tsplit)]
  }

inDFs3$Taxon <- factor(inDFs3$Taxon,levels = inDFs3$Taxon[order(inDFs3$Rel.Abundance)])
ggplot(inDFs3,aes(x=Taxon,y=Rel.Abundance,col=Taxon,fill=Taxon)) + geom_col() + theme(axis.text.x = element_text(angle = 90,hjust = 1)) 
ggplot(inDFs3,aes(x="",y=Rel.Abundance,col=Taxon,fill=Taxon)) + geom_col() + coord_polar("y", start=0)


# > vi) calculate microbiome diversity
# a) species richness (total number of species)
# b) genera richness (total number of genera)
# c) shannon diversity of species
# d) inverse simpson diversity of species
# question: does re-normalization after removal of "UNKNOWN" reads change the diversity?

# ======================================================================================================================
# > a) extract species (by grepping for |s__); Note: this also removes UNKNOWN 
inDFs <- inDF[grep('\\|s__',inDF$Taxon),]
specRichness <- nrow(inDFs)
print(specRichness)
# > b) extract genera: subset for genera then remove species
inDFg <- inDF[grep('\\|g__',inDF$Taxon),]
inDFg <- inDFg[-grep('\\|s__',inDFg$Taxon),]
genRichness <- nrow(inDFg)
print(genRichness)
# > c) shannon diversity: we use vegan:diversity to calculate it
inDFs <- inDF[grep('\\|s__',inDF$Taxon),]
#inDFs$Rel.Abundance <- inDFs$Rel.Abundance/sum(inDFs$Rel.Abundance)
diversity(inDFs$Rel.Abundance,index = "shannon")

# > d) inverse simpson diversity (vegan::diversity can calculate it as well)
inDFs <- inDF[grep('\\|s__',inDF$Taxon),]
diversity(inDFs$Rel.Abundance,index = "invsimpson")

# > vii) examine our three samples:
#   a) at phylum level
#   b) at level of species diversity
#   c) is there anything unusual in these samples?
# =====================================================
examineSample <- function(inFile) {
  inS1 <- read.table(file = inFile,sep='\t',header=F)
  colnames(inS1) <- c("Taxon","TaxID","Rel.Abundance","Coverage","Aligned.Reads")
  head(inS1)
  # phyla plot
  inS1p <- inS1[-grep('^k__.*\\|c__',inS1$Taxon),]
  inS1p <- inS1p[grep('UNKNOWN|^k__.*\\|p__',inS1p$Taxon),]
  inS1p$Taxon <- factor(inS1p$Taxon,levels = inS1p$Taxon[order(inS1p$Rel.Abundance)])
  print(ggplot(inS1p,aes(x=Taxon,y=Rel.Abundance,col=Taxon,fill=Taxon)) + geom_col() + theme(axis.text.x = element_text(angle = 90,hjust = 1)))
  inS1s <- inS1[grep('\\|s__',inS1$Taxon),]
  print(diversity(inS1s$Rel.Abundance,index = "shannon"))
}

examineSample('sample01_metaphlan.txt')
examineSample('sample02_metaphlan.txt')
examineSample('sample03_metaphlan.txt')

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

# 0) plot samples, phylum-level
inDFm2 <- inDFm
rownames(inDFm2) <- inDFm2$ID
inDFm2 <- inDFm2[,-grep('c__',colnames(inDFm2))]
inDFm2 <- inDFm2[,grep('p__',colnames(inDFm2))]
inDFm2 <- inDFm2[,colSums(inDFm2) != 0,]
inDFm2$ID <- row.names(inDFm2)


inDFm2long <- gather(inDFm2, Phylum, Rel.Abundance,k__Bacteria.p__Actinobacteria:k__Archaea.p__Euryarchaeota, factor_key=TRUE)

ggplot(inDFm2long,aes(x=ID,y=Rel.Abundance,col=Phylum,fill=Phylum)) + geom_col()

inDFm2long2 <- inDFm2long
ggplot(inDFm2long,aes(x=ID,y=Rel.Abundance,col=Phylum,fill=Phylum)) + geom_col()

# or ...
inDFm3 <- inDFm2
inDFm3$ID <- NULL
heatmap(log(as.matrix(inDFm3)+0.00001),Rowv = NA, Colv = NA)


# i) Diversity:
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

# ii) Beta diversity
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


# iii) quantifying the variance explained
# (Permutational Multivariate Analysis of Variance Using Distance Matrices)
# ======================================================
adonis(bcSpec ~ inDFm$Diagnosis,permutations = 1000)


# # iv) advanced/optional: biplot
# # ======================================================
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

# calculate biplot weights
sp <- inDFm[,grep('s__',colnames(inDFm))]
colnames(sp)=sapply(strsplit(colnames(sp), ".", fixed=TRUE), tail, 1)
colnames(sp)=gsub(".*s__","",colnames(sp))
pcoaMat2 <- cmdscale(bcSpec)
en=envfit(pcoaMat2,inDFm[,c("Age","Diagnosis","Biomarker","BMI")], na.rm=T, permutations=10000)

en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
en_coord_factors = as.data.frame(scores(en, display = "factors"))
en_coords_merged <- rbind.data.frame(en_coord_cont,en_coord_factors[rownames(en_coord_factors) %in% c("GenderM","DiagnosisDisease"),])
coord = en_coord_cont/2

ggplot(pcoaDF,aes(x=PC1,y=PC2,col=Diagnosis)) +
  geom_point(size=2) +
  theme_classic() +
  geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC1,y=centroids$PC2),alpha=1) +
  xlab(paste0('PC 1,',' (',var_pc1,'% variance)')) +
  ylab(paste0('PC 2,',' (',var_pc2,'% variance)')) +
  geom_segment(data=en_coords_merged, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), size =1, colour = "firebrick",arrow = arrow()) +
  geom_text(data = en_coords_merged, aes(x = Dim1, y = Dim2-0.05), colour = "firebrick", fontface = "bold", label = row.names(en_coords_merged))

# ==============================================================
# ==============================================================
# Part IV) compare microbiome taxa between healthy and disease
# ==============================================================
# ==============================================================
inDFm2 <- inDFm
taxaC <- grep('^k__',colnames(inDFm))
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
  taxName <- colnames(inDFm2)[taxC]
  oneTest <- wilcox.test(inDFm2[[taxC]][inDFm2$Diagnosis=="HC"],inDFm2[[taxC]][inDFm2$Diagnosis=="Disease"])
  res <- rbind.data.frame(res,data.frame(taxon=taxName,
                                         abHC=mean(inDFm2[[taxC]][inDFm2$Diagnosis=="HC"]),
                                         abD=mean(inDFm2[[taxC]][inDFm2$Diagnosis=="Disease"]),
                                         pV=oneTest$p.value))
}
res <- res[order(res$pV),]
res <- res[complete.cases(res),]
res$FDR <- p.adjust(res$pV)
res[res$FDR < 0.2,]

ggplot(inDFm2,aes(y=g__Escherichia,
                 x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

ggplot(inDFm2,aes(y=log(g__Escherichia+0.000001),
                  x=Diagnosis,col=Diagnosis)) + geom_violin() + geom_boxplot(width=0.2,outlier.alpha = F) + geom_jitter(width = 0.1,height = 0)

