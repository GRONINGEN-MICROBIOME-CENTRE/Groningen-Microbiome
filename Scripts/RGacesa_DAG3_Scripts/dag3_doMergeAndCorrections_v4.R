# library(rgdal)

library(geojsonio) # for geocoding
library(ggmap) # basic stuff for maps
library(leaflet) # for fancy maps
library(htmltools) # required for leaflet pretty maps
library(sf) # for reverse geocoding
library(dplyr) # for 'group by'
library(vegan)
library(Rtsne)
library(nortest)
library(compositions)



# =================================================================
# =================================================================
# MAIN 
# =================================================================
# =================================================================
setwd('C:/Users/ranko/Documents/UMCG/DAG3_stats/')
source('leafletPatch.R')
source('../myLibs/R_Microbiome_scripts.R')
source('../myLibs/R_ML_scripts_v3b.R')
source('../myLibs/R_Misc.R')

# ======================================================
# ======================================================
# load data and prep it
# ======================================================
# ======================================================
# > load data and prep it
metaTownsJ <- read.table('data_mappen/dutch_stats_2015_cleaned.tsv',sep='\t',header=T)
dutchTowns2 <- readRDS(file="data_mappen/dutchTownsMapRdy.Rdata")
taxRdy <- read.table('data_processed/DAG3_metaphlan_unfiltered_cleaned_rescaled_nontransformed.tsv',header=T,sep='\t',stringsAsFactors = F)
smplsAnnot <- read.table('data_processed/DAG3_reverse_geocoded.tsv',sep='\t',header = T)
smplsAnnot$PSEUDOIDEXT <- NULL
metaTownsJf <- metaTownsJ[,c(1,5,8,9,10,11)]
colnames(metaTownsJf) <- c("municipaly","POP.KM2","HOUSES.KM2","URBAN.INDEX","ADDRESS.DENSITY","PROVINCE")
# > prep microbiome
taxRdyF <- filterMetaGenomeDF(taxRdy,keepLevels = c("P","C","O","F","G","S"), rescaleTaxa = T,presPerc = 0.01,minMRelAb = 0.00001)
taxRdyFann <- merge(taxRdyF,smplsAnnot,by="DAG3_sampleID")
taxRdyFann <- merge(taxRdyFann,metaTownsJf,by="municipaly")
taxRdyFann$URBAN.INDEX <- as.factor(taxRdyFann$URBAN.INDEX)

# MERGE WITH PHENOTYPES
taxRdyFann$LOG.POP.KM2 <- log(taxRdyFann$POP.KM2)
taxRdyFann$LOG.HOUSES.KM2 <- log(taxRdyFann$HOUSES.KM2)
taxRdyFann$LOG.ADDRESS.DENSITY <- log(taxRdyFann$ADDRESS.DENSITY)

qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
qc1$status <- NULL
tt <- merge(taxRdyFann,qc1,by.x="DAG3_sampleID",by.y="sample")
qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned_v2.csv',sep=',',header=T)
qc2$id.external <- NULL
qc2$sample.nr <- NULL
qc2$nano.260.280 <- NULL
tt <- merge(tt,qc2,by.x="DAG3_sampleID",by.y="sample")
qs <- read.table('data/DAG3_qBasic.csv',sep=',',header=T)
tt <- merge(tt,qs,by="DAG3_sampleID")
tt <- reorderMicrobiomeDF(tt)

# address NAs in phenotypes
phenos <- colnames(subsetMicrobiomeDF(tt,getPhenos=T,getTaxa=F))
for (p in phenos) {
  if (sum(is.na(tt[[p]])) > 0) {print(paste0(' > feature ',p,' has ',sum(is.na(tt[[p]])), ' NAs '))}
}
tt$BMI[is.na(tt$BMI)] <- median(tt$BMI,na.rm = T)

write.table(tt,'data_processed/profiling/DAG3_metaphlan_filtered_cleaned_rescaled_nontransformed.tsv',sep=',',row.names=F)

# MERGE WITH pathways
pws <- read.table('data/DAG3_humann_pathabundances_merged_clean.tsv',sep='\t',header=T,quote = "",stringsAsFactors = F)
rownames(pws) <- pws$ID
pws$ID <- NULL
pws <- as.data.frame(t(pws))
pws$DAG3_sampleID <- as.character(rownames(pws))
rownames(pws) <- NULL
pws$DAG3_sampleID <- gsub('_kneaddata_merged_pathabundance','',pws$DAG3_sampleID)

pwsF <- filterHumannDF(pws, rescale = T,presPerc = 0.00001,minMRelAb = 0.0000001)
write.table(pwsF,'data_processed/profiling/DAG3_pathways_unfiltered_rescaled.csv',sep=',',row.names = F)
pwsF <- filterHumannDF(pws, rescale = T,presPerc = 0.05,minMRelAb = 0.0001)
write.table(pwsF,'data_processed/profiling//DAG3_pathways_filtered_rescaled.csv',sep=',',row.names = F)

ttt <- merge(tt,pwsF,by="DAG3_sampleID")
ttt <- reorderMicrobiomeDF(ttt)
write.table(ttt,'data_processed/profiling/DAG3_taxa_pwys_filtered_rescaled.csv',sep=',',row.names = F)

# MERGE WITH CARD
card <- read.table('data/DAG3_CARD_2018_merged.txt',sep='\t',header=T,quote='"',stringsAsFactors = F)
row.names(card) <- card$ID
card$ID <- NULL
card <- as.data.frame(t(card))
card$DAG3_sampleID <- as.character(rownames(card))
rownames(card) <- NULL

write.table(card,'data_processed/profiling/DAG3_CARD_unfiltered.csv',sep=',',row.names = F)
cardF <- filterCardDF(card, rescale = T,presPerc = 0.00001,minMRelAb = 0.0000001)
write.table(card,'data_processed/profiling/DAG3_CARD_unfiltered_rescaled.csv',sep=',',row.names = F)
cardF <- filterCardDF(card, rescale = T,presPerc = 0.01,minMRelAb = 0.0001)
write.table(card,'data_processed/profiling/DAG3_CARD_filtered_rescaled.csv',sep=',',row.names = F)
colnames(cardF)[2:ncol(cardF)] <- paste0('CARD_',colnames(cardF)[2:ncol(cardF)])

ttt <- merge(ttt,cardF,by="DAG3_sampleID")
ttt <- reorderMicrobiomeDF(ttt)
write.table(ttt,'data_processed/profiling/DAG3_taxa_pwys_CARD_filtered_rescaled.csv',sep=',',row.names = F)

# MERGE WITH VFDB
vfdb <- read.table('data/DAG3_VFDB_SB_vfs_merged.tsv',sep='\t',stringsAsFactors = F,header=T)
row.names(vfdb) <- vfdb$ID
vfdb$ID <- NULL
vfdb <- as.data.frame(t(vfdb))
vfdb$DAG3_sampleID <- as.character(rownames(vfdb))
rownames(vfdb) <- NULL
write.table(vfdb,'data_processed/profiling/DAG3_VFDB_unfiltered.csv',sep=',',row.names = F)
vfdb <- filterVfdbDF(vfdb, rescale = T,presPerc = 0.00001,minMRelAb = 0.0000001)
write.table(card,'data_processed/profiling/DAG3_VFDB_unfiltered_rescaled.csv',sep=',',row.names = F)
vfdb <- filterVfdbDF(vfdb, rescale = T,presPerc = 0.01,minMRelAb = 0.0001)
write.table(card,'data_processed/profiling/DAG3_VFDB_filtered_rescaled.csv',sep=',',row.names = F)

ttt <- merge(ttt,vfdb,by="DAG3_sampleID")
ttt <- reorderMicrobiomeDF(ttt)

write.table(ttt,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled.csv',sep=',',row.names = F)

# PRESENCE - ABSENCE TABLE
tttColnames <- colnames(subsetMicrobiomeDF(ttt,getPhenos = F,getTaxa = T,getPWYs = T,getCARDs = T,getVFs = T))
tttPA <- ttt
tttPA[,tttColnames][tttPA[,tttColnames] > 0] <- 1
tttPA[,tttColnames][tttPA[,tttColnames] <= 0] <- 0
write.table(tttPA,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled_PA.csv',sep=',',row.names = F)

# ==================================
# NORMALISE
# ==================================
# do normality tests on RAW (on non-zeros)
res = data.frame()
inDF <- subsetMicrobiomeDF(ttt,getPhenos = F,getTaxa = T,getPWYs = T,getCARDs = T,getVFs = T)
# 
# for (c in colnames(inDF)) {
#   #print (paste0(' > testing ',c))
#   nz <- inDF[[c]][inDF[[c]] > 0]
#   ptest <- pearson.test(nz)
#   res <- rbind.data.frame(res,data.frame(Feature=c,Pearson_PtoDF=ptest[[1]]/ptest[[6]],Pearson_Pv=ptest[[2]])) 
# }
# res$PCorr <- p.adjust(res$Pearson_Pv,method = "bonferroni")
# res <- res[order(res$Pearson_PtoDF),]
# print (paste0(' >> ',sum(res$PCorr < 0.05),' / ',nrow(res),' Features are *non-normal* !' ))
# summary(res$Pearson_PtoDF)

# normalise (using asin-sqrt)
toNormalise <- colnames(inDF)
dfNorm <- ttt
for (c in toNormalise) {
  dfNorm[[c]] <- asin(sqrt(dfNorm[[c]]))
}
dfNormCLR <- ttt
for (c in toNormalise) {
  dfNormCLR[[c]] <- clr(dfNorm[[c]])
}

# do normality tests on NORMALISED (on non-zeros)
res2 = data.frame()
inDF <- subsetMicrobiomeDF(dfNorm,getPhenos = F,getTaxa = T,getPWYs = T,getCARDs = T,getVFs = T)
for (c in colnames(inDF)) {
  nz <- inDF[[c]][inDF[[c]] > 0]
  ptest <- pearson.test(nz)
  res2 <- rbind.data.frame(res2,data.frame(Feature=c,Pearson_PtoDF=ptest[[1]]/ptest[[6]],Pearson_Pv=ptest[[2]]))
}
res2$PCorr <- p.adjust(res2$Pearson_Pv,method = "bonferroni")
res2 <- res2[order(res2$Pearson_PtoDF),]
print (paste0(' >> ',sum(res2$PCorr < 0.05),' / ',nrow(res2),' Features are *non-normal* !' ))
# summary(res2$Pearson_PtoDF)
# example(s)
# hist(dfNorm[[res2$Feature[1] ]][dfNorm[[res2$Feature[1]]] > 0])

# do normality tests on NORMALISED (CLR)
res3 = data.frame()
inDF <- subsetMicrobiomeDF(dfNormCLR,getPhenos = F,getTaxa = T,getPWYs = T,getCARDs = T,getVFs = T)
for (c in colnames(inDF)) {
  nz <- inDF[[c]][inDF[[c]] > 0]
  ptest <- pearson.test(nz)
  res3 <- rbind.data.frame(res3,data.frame(Feature=c,Pearson_PtoDF=ptest[[1]]/ptest[[6]],Pearson_Pv=ptest[[2]]))
}
res3$PCorr <- p.adjust(res3$Pearson_Pv,method = "bonferroni")
res3 <- res3[order(res3$Pearson_PtoDF),]
print (paste0(' >> ',sum(res3$PCorr < 0.05),' / ',nrow(res3),' Features are *non-normal* !' ))

# save it
write.table(dfNorm,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor.csv',sep=',',row.names = F)
write.table(dfNormCLR,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_CLRnor.csv',sep=',',row.names = F)

# ==================================
# CALCULATE DIVERSITIES and other metrics
# ==================================
divMatrix <- calcDIVMetrics(ttt,ID="DAG3_sampleID")
write.table(divMatrix,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_diversitiesonly.csv',sep=',',row.names = F)

# ==================================
# CORRECT (NORMALISED DATA)
# ==================================
# > correct for technical parameters ("vol.ul","postclean.reads","conc.ng.ul")
dfNormLcTech <- linearCorrectMGPwy(dfNorm,corrNames=c("vol.ul","postclean.reads","conc.ng.ul"),negIsZero=F,debug = F,removeCorCol = F)
write.table(dfNormLcTech,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor_corrtech.csv',sep=',',row.names = F)
# > correct for age, sex, BMI
dfNormLcASB <- linearCorrectMGPwy(dfNormLcTech,corrNames=c("age","gender","BMI"),negIsZero=F,debug = T,removeCorCol = F)
write.table(dfNormLcASB,'data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor_corrtechASB.csv',sep=',',row.names = F)

# ======================================================================
# ======================================================================
# ======================================================================
# ======================================================================


