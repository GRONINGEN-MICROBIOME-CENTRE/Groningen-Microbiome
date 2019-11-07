library(geojsonio) # for geocoding
library(ggmap) # basic stuff for maps
library(leaflet) # for fancy maps
library(htmltools) # required for leaflet pretty maps
library(sf) # for reverse geocoding
library(dplyr) # for 'group by'
library(vegan)
library(Rtsne)

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

# =================================================================
# =================================================================
# DO STATS ON FILTERED, CORRECTED, DAG3
# =================================================================
# =================================================================

taxRdyFannCorrTechnical <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled.csv',sep=',',header=T)
#taxRdyFannCorrTechnASB <- read.table('data_processed/DAG3_taxa_annotated_filtered_correctedforTechAgeGenderBMI.csv',sep=',',header=T)

# =================================================================
# do species-based PCA (filtered species)
# =================================================================
if (file.exists('data_processed/DAG3_Species_PCAs_annotated.tsv')) {
  pcasCorrTech <- read.table('data_processed/DAG3_Species_PCAs_annotated.tsv',sep='\t',header=T)
} else {
  taxMeta <- subsetMicrobiomeDF(taxRdyFannCorrTechnical,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,getDivs = F)
  taxRdyForPCA <- subsetMicrobiomeDF(taxRdyFannCorrTechnical,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,getDivs = F)
  taxRdyForPCA <- taxRdyForPCA[grep('s__',colnames(taxRdyForPCA))]
  #taxRdyForPCA$DAG3_sampleID <- taxMeta$DAG3_sampleID
  bcurtis <- vegdist(taxRdyForPCA,method = "bray")
  r.pca <- prcomp(bcurtis, center = T,scale. = F)
  pcas <- as.data.frame(r.pca$rotation[,1:10])
  pcas$DAG3_sampleID <- taxRdyFannCorrTechnical$DAG3_sampleID
  shannon <- diversity(taxRdyForPCA)
  pcas$shannon <- shannon
  pcasM <- merge(pcas,taxMeta,by="DAG3_sampleID")
  write.table(pcasM,file = 'data_processed/DAG3_Species_PCAs_annotated.tsv',sep='\t',row.names = F)
  pcasCorrTech <- read.table(file = 'data_processed/DAG3_Species_PCAs_annotated.tsv',sep='\t',header=T)
}
# =================================================================
# do genera-based PCA (filtered genera)
# =================================================================
if (file.exists('data_processed/DAG3_Genera_PCAs_annotated.tsv')) {
  pcasCorrTech <- read.table('data_processed/DAG3_Genera_PCAs_annotated.tsv',sep='\t',header=T)
} else {
  taxMeta <- subsetMicrobiomeDF(taxRdyFannCorrTechnical,getPWYs = F,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = T,getDivs = F)
  taxRdyForPCA <- subsetMicrobiomeDF(taxRdyFannCorrTechnical,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,getDivs = F)
  taxRdyForPCA <- taxRdyForPCA[grep('s__',colnames(taxRdyForPCA))]
  #taxRdyForPCA$DAG3_sampleID <- taxMeta$DAG3_sampleID
  bcurtis <- vegdist(taxRdyForPCA,method = "bray")
  r.pca <- prcomp(bcurtis, center = T,scale. = F)
  pcas <- as.data.frame(r.pca$rotation[,1:10])
  pcas$DAG3_sampleID <- taxRdyFannCorrTechnical$DAG3_sampleID
  shannon <- diversity(taxRdyForPCA)
  pcas$shannon <- shannon
  pcasM <- merge(pcas,taxMeta,by="DAG3_sampleID")
  write.table(pcasM,file = 'data_processed/DAG3_Genera_PCAs_annotated.tsv',sep='\t',row.names = F)
  pcasCorrTech <- read.table(file = 'data_processed/DAG3_Genera_PCAs_annotated.tsv',sep='\t',header=T)
}

# =================================================================
# do species-based PCA (filtered species, corrected for tehnicals + BMI + Age + Sex)
# =================================================================
if (file.exists('data_processed/DAG3_FilCorrTechASB_PCAs_annotated.tsv')) {
  pcasCorrTechASB <- read.table('data_processed/DAG3_FilCorrTechASB_PCAs_annotated.tsv',sep='\t',header=T)
} else {
  taxRdyForPCA <- taxRdyFannCorrTechnASB[,grep('s__',colnames(taxRdyFannCorrTechnical))]
  taxMeta <- taxRdyFannCorrTechnASB[,-grep('__',colnames(taxRdyFannCorrTechnical))]
  bcurtis <- vegdist(taxRdyForPCA,method = "bray")
  r.pca <- prcomp(bcurtis, center = T,scale. = F)
  pcas <- as.data.frame(r.pca$rotation[,1:10])
  pcas$DAG3_sampleID <- taxRdyFannCorrTechnASB$DAG3_sampleID
  shannon <- diversity(taxRdyForPCA)
  pcas$shannon <- shannon
  pcasM <- merge(pcas,taxMeta,by="DAG3_sampleID")
  write.table(pcasM,file = 'data_processed/DAG3_FilCorrTechASB_PCAs_annotated.tsv',sep='\t',row.names = F)
  pcasCorrTechASB <- read.table(file = 'data_processed/DAG3_FilCorrTechASB_PCAs_annotated.tsv',sep='\t',header=T)
}

# =================================================================
# do adonis (takes forever)
# =================================================================
# ===================================================================
adonisRdy <- taxRdyFann
colnames(adonisRdy)[2] <- "sample"
qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned.csv',sep=',',header=T)
qc2$id.external <- NULL
qc2$sample.nr <- NULL
adonisRdy <- merge(adonisRdy,qc2,by="sample")
qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
adonisRdy <- merge(adonisRdy,qc1,by="sample")
adonisRdy$sample <- NULL
adonisRdy$status <- NULL
adonisRdy$isolation.person <- NULL
adonisRdy$isolation.date <- NULL
adonisRdy$POSTCODE <- NULL

adonisSpec <- adonisRdy[,grep('s__',colnames(adonisRdy))]
adonisVars <- adonisRdy[,-grep('__',colnames(adonisRdy))]
bcurtis <- vegdist(adonisSpec,method = "bray") 

#adonisVarsTouse <- adonisVars[,c(4,5,8,10,11,12,13,15,16,17,20,21)]
adonisVarsTouse <- adonisVars[,c(4,5,8,10,11,12,13,14,15,16,19,21)]#,10,11,12,13,15,16,17,20,21)]

#adonis <- adonis(bcurtis ~ .,data=adonisVarsTouse,permutations=100,parallel=F)

# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  ad<-adonis(bcurtis ~ adonisVars[,i],permutations=2500)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
adonis_meta <- adonis_meta[order(adonis_meta$`FDR(BH)`),]
rownames(adonis_meta) = colnames(adonisVarsTouse)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
ggsave(g,"uncorr_adonis_results.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.table(adonis_meta,"adonis_meta.csv",sep=",",row.names=F)


# repeat adonis for good measure
# =========================================
taxRdyFannCorr$nano.260.230 <- NULL
adonisSpec <- taxRdyFannCorr[,grep('s__',colnames(taxRdyFannCorr))]
adonisSpec[adonisSpec < 0] <- 0.0
taxRdyFannCorr[,grep('s__',colnames(adonisRdy))] <- adonisSpec

taxRdyFannCorr3P <- taxRdyFannCorr[taxRdyFannCorr$PROVINCE %in% c("Drenthe","Friesland","Groningen"),]
adonisSpec <- taxRdyFannCorr3P[,grep('s__',colnames(taxRdyFannCorr3P))]
adonisVars <- taxRdyFannCorr3P[,-grep('__',colnames(taxRdyFannCorr3P))]
bcurtis <- vegdist(adonisSpec,method = "bray") 

adonisVarsTouse <- adonisVars[,c(1,4,5,8,10,11,12,13,14,15)]
adonisVarsTouse$municipaly <- as.factor(as.character(adonisVarsTouse$municipaly))
adonisVarsTouse$URBAN.INDEX.N <- as.numeric(as.character(adonisVarsTouse$URBAN.INDEX))

# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  ad<-adonis(bcurtis ~ adonisVars[,i],permutations=2500)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
adonis_meta <- adonis_meta[order(adonis_meta$`FDR(BH)`),]
rownames(adonis_meta) = colnames(adonisVarsTouse)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
print(g)
ggsave(plot=g,filename = "corr_adonis_meta.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.table(adonis_meta,"corr_adonis_meta.csv",sep=",",row.names=F)
write.table(taxRdyFannCorr,"corr_DAG3.tsv",sep="\t",row.names=F)

# ======================================================================================
#      ================= DO STATISTICS ON TECH CORRECTED METAPHLAN =====================
# ======================================================================================
# > PCs & gender vs gender, age, BMI

# SEX
responseVar <- "gender"

plotPCAs(pcasCorrTech,responseVar = responseVar)

testOneFeature(dataIn = pcasCorrTech,feature = "PC1",responseVar = responseVar,doSave = F,retPlot = T,
                yLab = "Principal component 1", xLab = "Sex",title = "Sex VS Microbiome beta-diversity")

testOneFeature(dataIn = pcasCorrTech,feature = "PC2",responseVar = responseVar,doSave = F,retPlot = T,
               yLab = "Principal component 1", xLab = "Sex",title = "Sex VS Microbiome beta-diversity")

testOneFeature(dataIn = pcasCorrTech,feature = "shannon",responseVar = responseVar,doSave = F,retPlot = T,
               yLab = "Shannon diversity", xLab = "Sex",title = "Sex VS Microbiome alpha-diversity")

# AGE
responseVar <- "age"
frm <- reformulate(responseVar,"PC1")
mdl <- lm(frm,pcasCorrTech)
ggplotRegression(mdl)

frm <- reformulate(responseVar,"shannon")
mdl <- lm(frm,pcasCorrTech)
ggplotRegression(mdl)

mdl <- lm(frm,pcasCorrTech[pcasCorrTech$age > 20,])
ggplotRegression(mdl)
mdl <- lm(frm,pcasCorrTech[pcasCorrTech$age < 20,])
ggplotRegression(mdl)

mdl <- lm(frm,pcasCorrTech)
ggplotRegression(mdl)

# BMI
responseVar <- "BMI"
frm <- reformulate(responseVar,"PC1")
mdl <- lm(frm,pcasCorrTech[pcasCorrTech$BMI < 40,])
ggplotRegression(mdl)

frm <- reformulate(responseVar,"shannon")
mdl <- lm(frm,pcasCorrTech[pcasCorrTech$BMI < 40,])
ggplotRegression(mdl)

# tests for individual bugs ...
responseVar <- "URBAN.INDEX"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  krus <- kruskal.test(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,Krus.CSq=krus$statistic[[1]],Krus.df=krus$parameter[[1]],Krus.p=krus$p.value))
}

# 
# # > do tests
# #    > taxa vs Urban index
# responseVar <- "URBAN.INDEX"
# res <- data.frame()
# for (c in grep('__',colnames(taxRdyFann))) {
#   tax <- colnames(taxRdyFann)[[c]]
#   print (paste(' > testing',tax,'VS',responseVar))
#   frm <- reformulate(responseVar,tax)
#   krus <- kruskal.test(frm,taxRdyFann)
#   res <- rbind.data.frame(res,data.frame(TAXON=tax,Krus.CSq=krus$statistic[[1]],Krus.df=krus$parameter[[1]],Krus.p=krus$p.value))
# }
# resUrbanIndex <- res[order(res$Krus.p),]
# resUrbanIndex$FDR <- p.adjust(resUrbanIndex$Krus.p,method="fdr")
# resUrbanIndex[c(1:10),]
# #tax = "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii"
# for (tax in resUrbanIndex$TAXON[resUrbanIndex$FDR <= 0.005]) {
#   print(ggplot(taxRdyFann,aes_string(x="URBAN.INDEX",y=tax,col="URBAN.INDEX")) + geom_violin() + geom_boxplot(width=0.2))
# }
# =================================================

# ===========================================================
# STUFF after linear correction for Age + Sex + BMI
# ===========================================================
# URBAN INDEX
responseVar <- "URBAN.INDEX"

plotPCAs(pcasCorrTechASB,responseVar = responseVar)
pcasCorrTechASB$URBAN.INDEX <- as.factor(as.character(pcasCorrTechASB$URBAN.INDEX))

testOneFeature(dataIn = pcasCorrTechASB,feature = "PC1",responseVar = responseVar,doSave = F,retPlot = T,
               yLab = "Principal component 1", xLab = "URBAN.INDEX",title = "URBAN.INDEX VS Microbiome beta-diversity")

testOneFeature(dataIn = pcasCorrTechASB,feature = "PC2",responseVar = responseVar,doSave = F,retPlot = T,
               yLab = "Principal component 2", xLab = "URBAN.INDEX",title = "URBAN.INDEX VS Microbiome beta-diversity")

testOneFeature(dataIn = pcasCorrTechASB,feature = "shannon",responseVar = responseVar,doSave = F,retPlot = T,
               yLab = "Shannon diversity", xLab = "URBAN.INDEX",title = "URBAN.INDEX VS Microbiome alpha-diversity",ylim = c(2.5,4))

pcasCorrTechASB$URBAN.INDEX <- as.numeric(as.character(pcasCorrTechASB$URBAN.INDEX))
frm <- reformulate(responseVar,"shannon")
mdl <- lm(frm,pcasCorrTechASB$URBAN.INDEX)
ggplotRegression(mdl)

# lat (3 provinces only)
ggplotRegression(lm(data=pcasCorrTechASB[pcasCorrTechASB$PROVINCE %in% c("Drenthe","Groningen","Friesland"),], shannon ~ lat))
ggplotRegression(lm(data=pcasCorrTechASB[pcasCorrTechASB$PROVINCE %in% c("Drenthe","Groningen","Friesland"),], PC1 ~ lat))
ggplotRegression(lm(data=pcasCorrTechASB[pcasCorrTechASB$PROVINCE %in% c("Drenthe","Groningen","Friesland"),], PC2 ~ lat))



testOneFeature(dataIn=pcasCorrTechASB[pcasCorrTechASB$PROVINCE %in% c("Drenthe","Groningen","Friesland"),], 
               feature = "shannon",responseVar = "PROVINCE",doSave = F,retPlot = T,
               yLab = "Shannon-diversity", xLab = "PROVINCE",title = "PROVINCE VS Microbiome alpha-diversity"
               )

# lon

# do example bug test
# (numeric)
responseVar <- "URBAN.INDEX"
taxRdyFannCorrTechnASB[[responseVar]] <- as.numeric(as.character(taxRdyFannCorrTechnASB[[responseVar]]))
res <- data.frame()
for (c in grep('__',colnames(taxRdyFannCorrTechnASB))) {
  tax <- colnames(taxRdyFannCorrTechnASB)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  lmod <- lm(frm,taxRdyFannCorrTechnASB)
  smry <- summary(lmod)
  # labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
  #                    #"Intercept =",signif(fit$coef[[1]],2 ),
  #                    " Coef =",signif(fit$coef[[2]], 2),
  #                    " P =",signif(summary(fit)$coef[2,4], 2)))
  
  res <- rbind.data.frame(res,data.frame(TAXON=tax,Krus.CSq=krus$statistic[[1]],Krus.df=krus$parameter[[1]],Krus.p=krus$p.value))
}

# do example bug test
# (categorical)
responseVar <- "URBAN.INDEX"
res <- data.frame()
taxRdyFannCorrTechnASB[[responseVar]] <- as.factor(as.character(taxRdyFannCorrTechnASB[[responseVar]]))
for (c in grep('__',colnames(taxRdyFannCorrTechnASB))) {
  tax <- colnames(taxRdyFannCorrTechnASB)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  krus <- kruskal.test(frm,taxRdyFannCorrTechnASB)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,Krus.CSq=krus$statistic[[1]],Krus.df=krus$parameter[[1]],Krus.p=krus$p.value))
}
res$FDR <- p.adjust(res$Krus.p,method = "fdr")
res <- res[order(res$Krus.p,decreasing = F),]

taxRdyFannCorrTechnASB$URBAN.INDEX.C <- as.factor(as.character(taxRdyFannCorrTechnASB$URBAN.INDEX))

testOneFeature(dataIn = taxRdyFannCorrTechnASB,feature = "s__Sutterella_wadsworthensis",responseVar = "URBAN.INDEX.C",
                doSave = F,retPlot = T)

taxRdyFannCorrTechnASB$URBAN.INDEX.N <- as.numeric(as.character(taxRdyFannCorrTechnASB$URBAN.INDEX))
lmTest <- lm(taxRdyFannCorrTechnASB$s__Sutterella_wadsworthensis ~ taxRdyFannCorrTechnASB$URBAN.INDEX)
ggplotRegression(lmTest)
