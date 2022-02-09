# ======================================================
# ======================================================
# ======================================================
# MICROBIOME WORKSHOP, exercise 2
# ======================================================
# ======================================================
# ======================================================

setwd('~/2022_Teaching/codes/')

# load libraries
# ===========================
library(vegan)
library(ggplot2)
library(tidyr)

# > 0) load microbiome data (sample02)
# ======================================
inDF <- read.table(file = 'sample02_metaphlan.txt',sep='\t',header=F)
colnames(inDF) <- c("Taxon","TaxID","Rel.Abundance","Coverage","Aligned.Reads")
head(inDF)

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
length(inDFs$Taxon)
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
