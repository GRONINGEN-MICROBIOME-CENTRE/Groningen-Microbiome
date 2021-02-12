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
