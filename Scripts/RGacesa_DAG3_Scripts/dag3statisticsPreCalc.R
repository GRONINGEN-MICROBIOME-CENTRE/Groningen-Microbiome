library(vegan)
library(plyr)
library(ggplot2)
library(venn)
library(scales)
registerDoSEQ()
library("fastqcr")

# SET WD AND LOAD HELPER FUNCTIONS
setwd('~/UMCG/DAG3_stats/')
source('dag3helper.R')

# ====================================================================
# ====================================================================
# ===========                  MAIN                       ============
# ====================================================================
# ====================================================================
# grind through QC (note: we run this only once and save results)

# -- rerun
#qc.dir <- system.file("DAG3_rerun/qc_post/", package = "fastqcr")
# qcPostclean <- qc_aggregate("data/qc_postclean/")
# # --- load QC (post-clean) ---
# qcPostclean$seq.length <- NULL
# qcPostclean <- qcPostclean[qcPostclean$module=="Basic Statistics",]
# qcPostclean$module <- NULL
# qcPostclean <- as.data.frame(qcPostclean)
# qcPostclean <- qcPostclean[order(qcPostclean$tot.seq),]
# qcPostclean$tot.seq <- as.numeric(qcPostclean$tot.seq)
# qcPostclean$pct.gc <- as.numeric(qcPostclean$pct.gc)
# qcPostclean$pct.dup<- as.numeric(qcPostclean$pct.dup)/100.0
# qcPostclean$reads.dedup <- qcPostclean$tot.seq*(1.0-qcPostclean$pct.dup)
# qcPostclean$sample <- gsub("_kneaddata_merged","",qcPostclean$sample)
# qcPostclean$status <- NULL
# colnames(qcPostclean) <- c("sample","postclean.reads","postclean.gc","postclean.dup","postclean.reads.dedup")
# 
# # --- load raw metaphlan
inFile <- "data/DAG3_metaphlan_merged.txt"
taxRaw <- as.data.frame(t(read.table(inFile,sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"')))
#taxRaw <- taxRaw[,-grep('c__',colnames(taxRaw))]
taxRaw$sample <- rownames(taxRaw)
taxRaw$sample <- gsub("_metaphlan","",taxRaw$sample)
rownames(taxRaw) <- NULL
# 
# # --- merge ---
 taxRawQC <- merge.data.frame(taxRaw,qcPostclean,by="sample",all = T)
 taxRawQC.k.rr <- taxRawQC[,-grep('p__',colnames(taxRawQC))]
 taxRawQC.k.rr <- taxRawQC.k[order(taxRawQC.k$k__Bacteria,decreasing = F),]
# #taxRawQC.k <- taxRawQC.k[taxRawQC.k$k__Bacteria < 95,]
# write.table(taxRawQC.k.rr,file = "rerun_stats.csv",sep=',',row.names = F)

# ========= ALL DATA ==========
# ===========================

#qcPreclean <- qc_aggregate("data/qc_preclean/")
#saveRDS(qcPreclean,file='data/qcPreClean.RData')

#qcPostclean <- qc_aggregate("data/qc_postclean/")
#saveRDS(qcPostclean,file='data/qcPostClean.RData')

# --- load QC (post-clean) ---
qcPostclean <- readRDS('data/qcPostClean.RData')
qcPostclean$seq.length <- NULL
qcPostclean <- qcPostclean[qcPostclean$module=="Basic Statistics",]
qcPostclean$module <- NULL
qcPostclean <- as.data.frame(qcPostclean)
qcPostclean <- qcPostclean[order(qcPostclean$tot.seq),]
qcPostclean$tot.seq <- as.numeric(qcPostclean$tot.seq)
qcPostclean$pct.gc <- as.numeric(qcPostclean$pct.gc)
qcPostclean$pct.dup<- as.numeric(qcPostclean$pct.dup)/100.0
qcPostclean$reads.dedup <- qcPostclean$tot.seq*(1.0-qcPostclean$pct.dup)
qcPostclean$sample <- gsub("_kneaddata_merged","",qcPostclean$sample)
qcPostclean$status <- NULL
colnames(qcPostclean) <- c("sample","postclean.reads","postclean.gc","postclean.dup","postclean.reads.dedup")
write.table(qcPostclean,file = "DAG3_QC_postknead.csv",row.names = F,sep = ',')

# --- load raw metaphlan
#taxRaw <- prepData('data/DAG3_metaphlan_merged.txt',getLevel="P",minBac = -1.0,doDropZero = F,doRescale = F)
inFile <- 'data/DAG3_metaphlan_merged.txt'
taxRaw <- as.data.frame(t(read.table(inFile,sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"')))
taxRaw <- taxRaw[,-grep('c__',colnames(taxRaw))]
taxRaw$sample <- rownames(taxRaw)
taxRaw$sample <- gsub("_metaphlan","",taxRaw$sample)
rownames(taxRaw) <- NULL

# --- merge ---
taxRawQC <- merge.data.frame(taxRaw,qcPostclean,by="sample",all = T)
#taxRawQC <- merge.data.frame(taxRawQC,qcPreclean,by='sample',all=T)
taxRawQC.k <- taxRawQC[,-grep('p__',colnames(taxRawQC))]
taxRawQC.k <- taxRawQC.k[order(taxRawQC.k$k__Bacteria,decreasing = F),]

#taxRawQC.k <- taxRawQC.k[taxRawQC.k$k__Bacteria < 95,]

#taxRdy <- prepData('data/DAG3_metaphlan_merged.txt',getLevel="P",minBac = -1.0,doDropZero = T,doRescale = F)

write.table(taxRawQC.k,file = "DAG3_stats.csv",sep=',',row.names = F)

# ============ load data ==============
# which tax level are we lookig at
lvls = c("P","F","G","S")
lvls = c("P")
# ================================
for (lvl in lvls)
{
  # full level
  taxonName <- NULL
  if        (lvl == 'S') {taxonName <- "Species"
  } else if (lvl == 'G') {taxonName <- "Genera"
  } else if (lvl == 'F') {taxonName <- "Families"
  } else if (lvl == 'C') {taxonName <- "Classes"
  } else if (lvl == 'O') {taxonName <- "Orders"
  } else if (lvl == 'P') {taxonName <- "Phyla"
  } 
  
  # prep DAG3
  dag3TaxF <- prepData('data/DAG3_metaphlan_merged.txt',getLevel=lvl)
  dag3TaxF$Cohort <- "DAG3"
  dag3TaxF$Diagnosis <- "HC"
  dag3TaxF$ID <- gsub('_metaphlan','',rownames(dag3TaxF))

  # merged
  allCohTaxF <- dag3TaxF
  allCohTaxF[is.na(allCohTaxF)] <- 0.0
  allCohTaxF <- shortenNames(allCohTaxF)
  allCohTaxF$Diagnosis <- as.factor(allCohTaxF$Diagnosis)
  
  allCohTaxF <- dag3TaxF
  allCohTaxF$ID <- NULL
  
  # filter by prevalence & abundance
  allCohTaxF_001 <- filterAbundPrev(allCohTaxF,minAbundance = 0.001,minPrevalence=0.001,filterPerCohort=F)
  allCohTaxF_01 <- filterAbundPrev(allCohTaxF,minAbundance = 0.001,minPrevalence=0.005,filterPerCohort=F)
  allCohTaxF_1 <- filterAbundPrev(allCohTaxF,minAbundance = 0.001,minPrevalence=0.01,filterPerCohort=F)
  allCohTaxF_10 <- filterAbundPrev(allCohTaxF,minAbundance = 0.001,minPrevalence=0.1,filterPerCohort=F)
  
  # ============ do diversities ===============
  # plots
  pw = 7.5
  ph = 6
  
  # diversity (shannon)
  shannon <- diversity(subset(allCohTaxF,select=-c(Cohort,Diagnosis)),index = "shannon",MARGIN = 1)
  diversityDF <- cbind.data.frame(div.shannon=shannon,cohort=allCohTaxF$Cohort,diagnosis=allCohTaxF$Diagnosis)
  simpson <- diversity(subset(allCohTaxF,select=-c(Cohort,Diagnosis)),index = "simpson",MARGIN = 1)
  diversityDF <- cbind.data.frame(diversityDF,simpson=simpson,diagnosis=allCohTaxF$Diagnosis)
  specNr <- specnumber(subset(allCohTaxF,select=-c(Cohort,Diagnosis)),MARGIN = 1)
  diversityDF <- cbind.data.frame(diversityDF,species.richness=specNr)
  
  #hist(shannon)
  g <- ggplot(diversityDF,aes(x=cohort,y=shannon,col=cohort)) + geom_violin() + geom_jitter(width=0.2,height=0.001,alpha=0.2) + 
    geom_boxplot(alpha = 0.66,width=0.25,outlier.alpha = 0) + ggtitle(paste("Shannon diversity (",lvl,")",sep='')) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_shannon_comparison.png',sep=''),width = pw,height = ph)
  g <- ggplot(diversityDF,aes(x=shannon,col=cohort)) + geom_density(size=1.25) + 
    ggtitle(paste("Shannon div. density (",lvl,")",sep='')) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_shannon_comparison_density.png',sep=''),width = pw*3/4,height = ph)
  
  # diversity (simpson)
  g <- ggplot(diversityDF,aes(x=cohort,y=simpson,col=cohort)) + geom_violin() + geom_jitter(width=0.2,height=0.001,alpha=0.2) + 
    geom_boxplot(alpha = 0.66,width=0.25) + ggtitle(paste("Simpson diversity (",lvl,")",sep='')) + ylim(c(0.5,1)) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_simpson_comparison.png',sep=''),width = pw,height = ph)
  g <- ggplot(diversityDF,aes(x=simpson,col=cohort)) + geom_density(size=1.25) + 
    ggtitle(paste("Simpson div. density (",lvl,")",sep='')) + xlim(c(0.5,1)) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_simpson_comparison_density.png',sep=''),width = pw*3/4,height = ph)
  
  # NR of species
  g <- ggplot(diversityDF,aes(x=cohort,y=species.richness,col=cohort)) + geom_violin() + geom_jitter(width=0.2,height=0.1,alpha=0.2) + 
    geom_boxplot(alpha = 0.66,width=0.25) + ggtitle(paste(taxonName,"richness")) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_richness_comparison.png',sep=''),width = pw,height = ph)
  
  g <- ggplot(diversityDF,aes(x=species.richness,col=cohort)) + geom_density(size=1.25) + ggtitle("Species richness density (IBD)")
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_richness_comparison_density.png',sep=''),width = pw*3/4,height = ph)
  
  g <- ggplot(diversityDF[diversityDF$cohort=="IBD",],aes(x=species.richness,col=diagnosis)) + geom_density(size=1.25) + 
    ggtitle(ggtitle(paste(taxonName,"richness density (IBD)"))) + theme(text = element_text(size = 18))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_richness_comparison_density_ibd.png',sep=''),width = pw*3/4,height = ph)
  
  # ==== rarefication curves =====
  # ===============================
  nrBoots = 100
  steps <- c(10,20,30,40,50,75,100,150,200,300,400,503,750,1000,1083,1270,1500,2000,2500,3000,3500,4000,4220,5000,6000,7085,8000)
  
  inDFt <- allCohTaxF;
  inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
  rcurve <- doRarefaction(inDF = inDFt,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  g <- ggplot(rcurve,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.25) +
    theme(text = element_text(size = 18)) + xlab('Number of samples') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Rarefication curve (',lvl,', all)',sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_rare_0000.png',sep=''),width = pw*1.1,height = ph)
  
  inDFt <- allCohTaxF_001
  inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
  rcurve <- doRarefaction(inDF = inDFt,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  g <- ggplot(rcurve,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.25) +
    theme(text = element_text(size = 18)) + xlab('Number of samples') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Rarefication curve (',lvl,', prev 0.01%)',sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_rare_0001.png',sep=''),width = pw*1.1,height = ph)
  
  inDFt <- allCohTaxF_01
  inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
  rcurve <- doRarefaction(inDF = inDFt,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  g <- ggplot(rcurve,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.25) +
    theme(text = element_text(size = 18)) + xlab('Number of samples') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Rarefication curve (',lvl,', prev 0.1%)',sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_rare_001.png',sep=''),width = pw*1.1,height = ph)
  
  # prevalence 1
  inDFt <- allCohTaxF_1
  inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
  rcurve <- doRarefaction(inDF = inDFt,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  g <- ggplot(rcurve,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.25) +
    theme(text = element_text(size = 18)) + xlab('Number of samples') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Rarefication curve (',lvl,', prev 1%)',sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_rare_01.png',sep=''),width = pw*1.1,height = ph)
  
  # prevalence 10
  inDFt <- allCohTaxF_10
  inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
  rcurve <- doRarefaction(inDF = inDFt,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  g <- ggplot(rcurve,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.25) +
    theme(text = element_text(size = 18)) + xlab('Number of samples') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Rarefication curve (',lvl,', prev 10%)',sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_rare_1.png',sep=''),width = pw*1.1,height = ph)
  
  # do average taxa plots
  # ======= pretty pretty pie charts =========
  # - prevalence 0.01
  pie <- makePieChart(inDFt = allCohTaxF_01,cutLabel=0.05)
  write.table(allCohTaxF_01,file = paste('plots/',lvl,'_pie_table_01.csv',sep=''))
  ggsave(plot = pie,filename = paste('plots/',lvl,'_pie_01.png',sep=''),width = pw*1.1,height = ph)
  pie <- makePieChart(inDFt = allCohTaxF_1,cutLabel=0.05)
  ggsave(plot = pie,filename = paste('plots/',lvl,'_pie_1.png',sep=''),width = pw*1.1,height = ph)
  write.table(allCohTaxF_01,file = paste('plots/',lvl,'_pie_table_1.csv',sep=''))
  pie <- makePieChart(inDFt = allCohTaxF_10,cutLabel=0.05)
  ggsave(plot = pie,filename = paste('plots/',lvl,'_pie_10.png',sep=''),width = pw*1.1,height = ph)
  write.table(allCohTaxF_01,file = paste('plots/',lvl,'_pie_table_10.csv',sep=''))
  
    # find most prevalent taxa by cohort
  # plot prevalence - NR
  # =====================================
  psDf <- data.frame()
  for (prev.min in c( seq(0,0.1,0.001), seq(0.1,1,0.01)) ) {
    prevL <- getTaxPrev(allCohTaxF,responseVar = "Cohort",prevalenceMin = prev.min)
    for (coh in names(prevL)) {
      psDf <- rbind.data.frame(psDf,data.frame(Cohort=coh,min.prev=prev.min,NrTaxa=length(prevL[[coh]])))
    }
  }
  g <- ggplot(psDf,aes(x=min.prev,y=NrTaxa,col=Cohort)) + geom_line(size=1.25) + 
    geom_vline(xintercept=0.01,linetype="dashed") + geom_vline(xintercept = 0.05,linetype="dashed") + 
    theme(text = element_text(size = 18)) + xlab('Taxon prevalence') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Prevalence vs ',taxonName,sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_prevalence.png',sep=''),width = pw*1.1,height = ph)
  
  g <- ggplot(psDf,aes(x=log(min.prev),y=NrTaxa,col=Cohort)) + geom_line(size=1.25) + 
    geom_vline(xintercept=log(0.01),linetype="dashed") + geom_vline(xintercept = log(0.05),linetype="dashed") + 
    theme(text = element_text(size = 18)) + xlab('ln (prevalence)') + ylab(paste('Number of',taxonName)) + 
    ggtitle(paste('Prevalence vs ',taxonName,sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_prevalence_log.png',sep=''),width = pw*1.1,height = ph)
  
  # make pretty data frame and save it
  prevDF <- data.frame()
  cutOff <- 0.05
  prevL <- getTaxPrev(allCohTaxF,responseVar = "Cohort",prevalenceMin = cutOff)
  for (i in names(prevL)) {
    abundL <- getTaxAvgAbundance(allCohTaxF[allCohTaxF$Cohort==i,],abundanceMin = -1,responseVar = "Cohort")
    for (tx in names(prevL[[i]])) {
      txPrev <- prevL[[i]][tx]
      txMAb  <- abundL[[i]][tx]
      prevDF <- rbind.data.frame(prevDF,data.frame(cohort=i,taxon=tx,prevalence=txPrev,mean.abundance=txMAb) )
    }
  }
  write.table(prevDF,paste('data_out/',lvl,'_prev_abun.txt',sep=''),sep = '\t',row.names = F)
  
}
  # do venn diagrams
  c = 0
  plot_names <- c("00000","00001","0001","001","01")
  for (min.prev in c(0.0000001,0.0001,0.001,0.01,0.1)) {
    c = c + 1
    min.abundance = 0.0001
    inDFt <- filterAbundPrev(allCohTaxF,minAbundance = min.abundance,minPrevalence=min.prev,filterPerCohort=F)
    inDFt$Cohort[inDFt$Cohort=="HMP2.HC" | inDFt$Cohort=="HMP2.IBD"] <- "HMP2"
    prevL <- getTaxPrev(inDFt,responseVar = "Cohort",prevalenceMin = min.prev)
    
    nrNames <- 10000
    nmList <- list()
    for (n in names(prevL)) {
      nmList[[n]] <- names(prevL[[n]][1:nrNames])
    }
    plotName = paste('plots/',lvl,'_overlap_4_',plot_names[c],'.png',sep='')
    png(filename=plotName, type="cairo", units="in", width=10, height=7.5, pointsize=25, res=150)
    venn(nmList,zcol = c("green", "blue","red","orange"),title='burac',line=-1)
    title(main=paste(taxonName,' overlap (prev=',min.prev,')',sep=''),line=-1.5)
    dev.off()
    
    nmList3 <- nmList; nmList3[["IBD"]] <- NULL
    plotName = paste('plots/',lvl,'_overlap_3_',plot_names[c],'.png',sep='')
    png(filename=plotName, type="cairo", units="in", width=10, height=7.5, pointsize=25, res=150)
    venn(nmList3,zcol = c("green", "blue","red","orange"),title='burac',line=-1)
    title(main=paste(taxonName,' overlap (prev=',min.prev,')',sep=''),line=-1.5)
    dev.off()
    
    nmList2 <- nmList3; nmList2[["HMP2"]] <- NULL
    plotName = paste('plots/',lvl,'_overlap_2_',plot_names[c],'.png',sep='')
    png(filename=plotName, type="cairo", units="in", width=10, height=7.5, pointsize=25, res=150)
    venn(nmList2,zcol = c("green", "blue","red","orange"),title='burac',line=-1)
    title(main=paste(taxonName,' overlap (prev=',min.prev,')',sep=''),line=-1.5)
    dev.off()
    
  }
  
  # plot mean realtive abundance - NR
  abDf <- data.frame()
  for (ab.min in c(0,0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,seq(0.001,0.1,0.001), seq(0.1,0.15,0.01)) ) {
    prevL <- getTaxAvgAbundance(allCohTaxF,responseVar = "Cohort",abundanceMin = ab.min)
    for (coh in names(prevL)) {
      abDf <- rbind.data.frame(abDf,data.frame(Cohort=coh,min.abundance=ab.min,NrTaxa=length(prevL[[coh]])))
    }
  }
  nT <- max(abDf$NrTaxa)
  g <- ggplot(abDf,aes(x=min.abundance,y=NrTaxa,col=Cohort)) + geom_line(size=1.25) + ylim(0,min(100,nT)) +
    geom_vline(xintercept=0.001,linetype="dashed") + geom_vline(xintercept = 0.01,linetype="dashed") + theme(text = element_text(size = 18)) + 
    ggtitle(paste(taxonName,' vs rel. abundance',sep='')) + xlab('Relative abundance') + ylab(paste('NR of ',taxonName,sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_relabundance.png',sep=''),width = pw*1.1,height = ph)
  
  g <- ggplot(abDf,aes(x=log(min.abundance),y=NrTaxa,col=Cohort)) + geom_line(size=1.25) + theme(text = element_text(size = 18)) + 
    ggtitle(paste(taxonName,' vs rel. abundance',sep='')) + xlab('ln (relative abundance)') + ylab(paste('NR of ',taxonName,sep=''))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_relabundance_log.png',sep=''),width = pw*1.1,height = ph)
  
  # ============ do Bray-curtis ===========================
  # diversity (nr species)
  # do bray curtis
  fn <- paste('data_out/',lvl,'_pca_centered_scaled.txt',sep='')
  if (!file.exists(fn)) {
    print (' > Doing PCA analysis of B-C distance, this will take a while...')
    bcurtis <- vegdist(subset(allCohTaxF[,grep('__',colnames(allCohTaxF))],method = "bray"))
    # PCa scaled, centered
    r.pca.cs <- prcomp(bcurtis, center = T,scale. = F)
    r.pca.cs.10 <- as.data.frame(r.pca.cs$rotation[,1:10])
    write.table(r.pca.cs.10,fn,sep = '\t')
    print ('  >> PCA done')
  } else {
    r.pca.cs.10 <- read.table(fn)
  }
  
  # r.pca.c <- prcomp(bcurtis, center = T,scale. = F)
  # r.pca.c.10 <- as.data.frame(r.pca.c$rotation[,1:10])
  # write.table(r.pca.c.10,paste('data_out/',lvl,'_pca_centered.txt',sep=''),sep = '\t')
  # 
  # r.pca.s <- prcomp(bcurtis, center = F,scale. = T)
  # r.pca.s.10 <- as.data.frame(r.pca.s$rotation[,1:10])
  # write.table(r.pca.s.10,paste('data_out/',lvl,'_pca_scaled.txt',sep=''),sep = '\t')
  # 
  # r.pca.raw <- prcomp(bcurtis, center = F,scale. = F)
  # r.pca.raw.10 <- as.data.frame(r.pca.raw$rotation[,1:10])
  # write.table(r.pca.raw.10,paste('data_out/',lvl,'_pca_raw.txt',sep=''),sep = '\t')
  
  # test rpca
  # -- seems to work really really slow :()
  #r.pca.rpca <- rrpca(bcurtis)
  
  # plot them
  pcas <- r.pca.cs.10
  #pcas <- r.pca.raw.10
  pcas$Cohort <- allCohTaxF$Cohort
  
  centroids <- aggregate(cbind(PC1,PC2) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC1,y=PC2,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC1,y=centroids$PC2),alpha=1) + 
    ggtitle(paste(taxonName," Bray-Curtis PCA",sep='')) + theme(text = element_text(size = 17))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_BC_PCA_12.png',sep=''),width = pw,height = ph*1.1)
  
  centroids <- aggregate(cbind(PC2,PC3) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC2,y=PC3,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC2,y=centroids$PC3),alpha=1) + 
    ggtitle(paste(taxonName," Bray-Curtis PCA",sep='')) + theme(text = element_text(size = 17))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_BC_PCA_23.png',sep=''),width = pw,height = ph*1.1)
  
  centroids <- aggregate(cbind(PC3,PC4) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC3,y=PC4,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC3,y=centroids$PC4),alpha=1) + 
    ggtitle(paste(taxonName," Bray-Curtis PCA",sep='')) + theme(text = element_text(size = 17))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_BC_PCA_34.png',sep=''),width = pw,height = ph*1.1)
  
  centroids <- aggregate(cbind(PC4,PC5) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC4,y=PC5,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3.5,size=6,aes(x=centroids$PC4,y=centroids$PC5),alpha=1) + 
    ggtitle(paste(taxonName," Bray-Curtis PCA",sep='')) + theme(text = element_text(size = 17))
  print(g)
  ggsave(plot = g,filename = paste('plots/',lvl,'_BC_PCA_45.png',sep=''),width = pw,height = ph*1.1)
}



# make data comply!
# ===================================

# tmp
inDF <- allCohTaxF

transformCohortData <- function(inDF, toTrans="DAG3",lvl,referenceCoh="LLD",doTransformAS=T,extraTransforms=1,
                                outTables='transformations',outPlots='transformations/plots',sx=8,sy=5,doPlots=F,
                                doRawBC=F,bcCenter=F,bcScale=F) 
{
  origDF <- inDF
  # 0) do sqrt/arcsin/sqrt transformation (works better then just arcsin(sqrt))
  if (doTransformAS) {
    inDF[,grep('__',colnames(inDF))] <- asin(sqrt(inDF[,grep('__',colnames(inDF))]))
  }
  if (extraTransforms > 0) {
    for (r in c(1:extraTransforms)) {
      inDF[,grep('__',colnames(inDF))] <- sqrt(inDF[,grep('__',colnames(inDF))])
    }
  }
  inDF <- inDF[inDF$Cohort %in% c(toTrans,referenceCoh),]
  # a) take only taxa in common to alldatasets
  # filter for prevalence (> 0.05)
  inDF <- filterAbundPrev(inDF,minPrevalence = 0.05,filterPerCohort = F,minAbundance = -1,dropCols = T)
  
  inDFref <- inDF[inDF$Cohort == referenceCoh,]
  inDFref <- inDFref[,grep('__',colnames(inDFref))]
  inDFtoTrans <- inDF[inDF$Cohort == toTrans,]
  inDFtoTrans <- inDFtoTrans[,grep('__',colnames(inDFtoTrans))]
  # remove all-zero columns
  inDFref <- inDFref[,colSums(inDFref) > 0]
  inDFtoTrans <- inDFtoTrans[,colSums(inDFtoTrans) > 0]
  # keep only shared columns
  taxRef <- colnames(inDFref)
  taxToTrans <- colnames(inDFtoTrans)
  toKeep <- intersect(taxRef,taxToTrans)
  inDFtoTrans <- inDFtoTrans[,toKeep]
  inDFref <- inDFref[,toKeep]
  toKeep <- c(colnames(inDF)[grep('__',colnames(inDF),invert = T)],toKeep)
  inDF <- inDF[,toKeep]
  # go through, check concordance
  res <- data.frame()
  for (c in toKeep[grep('__',toKeep)]) {
    t.wilcox <- wilcox.test(inDFref[[c]][inDFref[[c]]!=0],inDFtoTrans[[c]][inDFtoTrans[[c]]!=0])
    t.nr.ref = sum(inDFref[[c]]!=0)
    t.prev.ref = t.nr.ref/nrow(inDFref)
    t.nr.toTrans = sum(inDFtoTrans[[c]]!=0)
    t.prev.toTrans = t.nr.toTrans/nrow(inDFtoTrans)
    res <- rbind.data.frame(res,data.frame(taxon=c,nrRef=t.nr.ref,prevRef=t.prev.ref,nrTar=t.nr.toTrans,prevTar=t.prev.toTrans,
                                           wilcox.p.abund=t.wilcox$p.value))
  }
  res <- res[order(res$wilcox.p,decreasing = F),]
  write.table(res,file = paste(outTables,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_raw_stats.csv',sep=''),sep=',',row.names = F)
  
  # do "raw" Bray-curtis
  if (doRawBC) {
    # PCa scaled, centered
    inDFt <- inDF[,grep('__',colnames(inDF))]
    print (' >> doing PCA of raw data')
    bcurtis.RAW <- vegdist(subset(inDFt,method = "bray"))
    # PCa scaled, centered
    r.pca.cs.RAW <- prcomp(bcurtis.RAW, center = bcCenter,scale. = bcScale)
    r.pca.cs.RAW.10 <- as.data.frame(r.pca.cs.RAW$rotation[,1:10])
    print ('  >>> PCA DONE <<<')
    pcas <- r.pca.cs.RAW.10
    pcas$Cohort <- inDF$Cohort
    centroids <- aggregate(cbind(PC1,PC2) ~ Cohort,pcas,mean)
    g <- ggplot(pcas,aes(x=PC1,y=PC2,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)+
      ggtitle(paste('B-C PCA (',toTrans,' vs ',referenceCoh,', ',length(res$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
    print(g)
    ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_raw_12.png',sep=''),width = sx*3/4,height = sy)
    
    centroids <- aggregate(cbind(PC2,PC3) ~ Cohort,pcas,mean)
    g <- ggplot(pcas,aes(x=PC2,y=PC3,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3),alpha=1)+
      ggtitle(paste('B-C PCA (',toTrans,' vs ',referenceCoh,', ',length(res$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
    print(g)
    ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_raw_23.png',sep=''),width = sx*3/4,height = sy)
    
    centroids <- aggregate(cbind(PC3,PC4) ~ Cohort,pcas,mean)
    g <- ggplot(pcas,aes(x=PC3,y=PC4,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4),alpha=1)+
      ggtitle(paste('B-C PCA (',toTrans,' vs ',referenceCoh,', ',length(res$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
    print(g)
    ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_raw_34.png',sep=''),width = sx*3/4,height = sy)
    
    centroids <- aggregate(cbind(PC4,PC5) ~ Cohort,pcas,mean)
    g <- ggplot(pcas,aes(x=PC4,y=PC5,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5),alpha=1)+
      ggtitle(paste('B-C PCA (',toTrans,' vs ',referenceCoh,', ',length(res$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
    print(g)
    ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_raw_45.png',sep=''),width = sx*3/4,height = sy)
  }
  
  # do transformations
  inDF2 <- inDF
  res2 <- data.frame()
  cnt = 0
  transformRuns <- 100
  for (c in res$taxon[c(1:length(res$taxon))]) {
  #for (c in res$taxon[c(1:5)]) {  # debug!
    print(paste(' > transforming abundances of ',c))
    cnt <- cnt + 1
    # plot (pre-normalisation)
    if (doPlots) {
      g <- ggplot(origDF[origDF$Cohort==toTrans | origDF$Cohort==referenceCoh,],aes_string(x=c,col="Cohort")) + geom_density(size=1.25) + 
        ggtitle(paste(c,"(Raw)")) + theme(text = element_text(size = 16))
      print(g)
      ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_',cnt,'_a.png',sep=''),width = sx,height = sy)
    }
    # plot (post normalisation)
    if (doPlots) {
      g <- ggplot(inDF,aes_string(x=c,col="Cohort")) + geom_density(size=1.25) +
        ggtitle(paste(c,"(Normalised)")) + theme(text = element_text(size = 16))
      print(g)
      ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_',cnt,'_b.png',sep=''),width = sx,height = sy)
    }
    
    t.wilcox <- wilcox.test(inDFref[[c]][inDFref[[c]]!=0],inDFtoTrans[[c]][inDFtoTrans[[c]]!=0])
    t.nr.ref = sum(inDFref[[c]]!=0)
    t.prev.ref = t.nr.ref/nrow(inDFref)
    t.nr.toTrans = sum(inDFtoTrans[[c]]!=0)
    t.prev.toTrans = t.nr.toTrans/nrow(inDFtoTrans)  
    for (r in c(1:transformRuns)) {
      # rescale
      transfTax <- NULL
      sdShift <- sd(inDF2[[c]][inDF2[[c]]!=0 & inDF2$Cohort==referenceCoh]) / sd(inDF2[[c]][inDF2[[c]]!=0 & inDF2$Cohort==toTrans])
      transfTax <- inDF2[[c]][inDF2$Cohort==toTrans]
      transfTax[transfTax!=0] <- (transfTax[transfTax!=0]-median(transfTax[transfTax!=0]))*sdShift+median(transfTax[transfTax!=0])
      # move
      tarShift <- median(inDF2[[c]][inDF2[[c]]!=0 & inDF2$Cohort==referenceCoh]) - median(transfTax[transfTax!=0])
      transfTax[transfTax!=0] <- transfTax[transfTax!=0] + tarShift
      inDF2[[c]][inDF2$Cohort==toTrans] <- transfTax
    }
    inDF2[[c]][inDF2[[c]] < 0] <- 0
    # plot (transformed)
    if (doPlots) {
      g <- ggplot(inDF2,aes_string(x=c,col="Cohort")) + geom_density(size=1.25) + ggtitle(paste(c,"(Transformed)")) + theme(text = element_text(size = 16))
      print(g)
      ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_',cnt,'_c.png',sep=''),width = sx,height = sy)
    }
    
    t.nr.post <- sum(inDF2[[c]]!=0 & inDF2$Cohort==toTrans)
    t.prev.post <- t.nr.post / sum(inDF2$Cohort==toTrans)
    wc.post <- wilcox.test(inDF2[[c]][inDF2[[c]]!=0 & inDF2$Cohort==referenceCoh],inDF2[[c]][inDF2[[c]]!=0 & inDF2$Cohort==toTrans])
    # do proportion test for prevalence
    v1 <- c(sum(inDF2[[c]][inDF2$Cohort=="LLD"] > 0),sum(inDF2[[c]][inDF2$Cohort=="DAG3"] > 0))
    v2 <- c(sum(inDF2$Cohort=="LLD"),sum(inDF2$Cohort=="DAG3"))
    propt <- prop.test(v1,v2)
    res2 <- rbind.data.frame(res2,data.frame(taxon=c,prev.REF=t.prev.ref,prev.POST=t.prev.post,
                                             OR.POST=max(t.prev.ref/t.prev.post,t.prev.post/t.prev.ref),
                                             wilcox.p.abund.POST=wc.post$p.value,
                                             prev.POST.prop.pv=propt$p.value))
  }
  rownames(res2) <- NULL
  
  res2$FDR.abund.POST <- p.adjust(res2$wilcox.p.abund.POST)
  res2$FDR.prev.POST <- p.adjust(res2$prev.POST.prop.pv)
  res2 <- res2[order(res2$wilcox.p.abund.POST),]
  res3 <- res2[order(res2$prev.POST.prop.pv*res2$wilcox.p.abund.POST,decreasing = T),]
  rownames(res3) <- NULL
  
  # cutoff
  fdrCutPrev <- 1.0e-50 # 0.00005 seems to work fine!
  fdrCutAbun <- 0.01 # 0.05 seems to work fine
  res4 <- res3[res3$FDR.prev.POST > fdrCutPrev & res3$FDR.abund.POST > fdrCutAbun,]
  write.table(res3,file = paste(outTables,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_transformed_stats.csv',sep=''),sep=',',row.names = F)
  write.table(res4,file = paste(outTables,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_selected_stats.csv',sep=''),sep=',',row.names = F)
  
  # now try B-C with corrected data
  inDF3 <- cbind.data.frame(inDF2[,grep('__',colnames(inDF2),invert = T)],inDF2[,as.character(res4$taxon)])
  # output
  inDF3t <- inDF3
  write.table(inDF3t,file = paste(outTables,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_selected_taxonomy.csv',sep=''),sep=',',row.names = F)
  # do BC (Transformed)
  inDF3 <- inDF3[,grep('__',colnames(inDF3))]
  print (' >> doing PCA of transformed data')
  bcurtis.POST <- vegdist(subset(inDF3,method = "bray"))
  # PCa scaled, centered
  r.pca.cs.POST <- prcomp(bcurtis.POST, center = bcCenter,scale. = bcScale)
  r.pca.cs.POST.10 <- as.data.frame(r.pca.cs.POST$rotation[,1:10])
  print ('  >>> PCA DONE <<<')
  
  pcas <- r.pca.cs.POST.10
  pcas$Cohort <- inDF3t$Cohort
  
  centroids <- aggregate(cbind(PC1,PC2) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC1,y=PC2,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)+
    ggtitle(paste('B-C PCA (',toTrans,' -> ',referenceCoh,', ',length(res4$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
  print(g)
  ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_trans_',length(res4$taxon),'_12.png',sep=''),width = sx*3/4,height = sy)
  
  centroids <- aggregate(cbind(PC2,PC3) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC2,y=PC3,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3),alpha=1)+
    ggtitle(paste('B-C PCA (',toTrans,' -> ',referenceCoh,', ',length(res4$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
  print(g)
  ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_trans_',length(res4$taxon),'_23.png',sep=''),width = sx*3/4,height = sy)
  
  centroids <- aggregate(cbind(PC3,PC4) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC3,y=PC4,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4),alpha=1)+
    ggtitle(paste('B-C PCA (',toTrans,' -> ',referenceCoh,', ',length(res4$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
  print(g)
  ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_trans_',length(res4$taxon),'_34.png',sep=''),width = sx*3/4,height = sy)
  
  centroids <- aggregate(cbind(PC4,PC5) ~ Cohort,pcas,mean)
  g <- ggplot(pcas,aes(x=PC4,y=PC5,col=Cohort)) + geom_point(alpha=0.33)+ geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5),alpha=1)+
    ggtitle(paste('B-C PCA (',toTrans,' -> ',referenceCoh,', ',length(res4$taxon),' taxa)',sep='')) + theme(text = element_text(size = 16))
  print(g)
  ggsave(plot=g,filename = paste(outPlots,'/',lvl,'_',toTrans,'_to_',referenceCoh,'_BC_trans_',length(res4$taxon),'_45.png',sep=''),width = sx*3/4,height = sy)
  
  inDF3t
}
# retest B/C plots

