library(vegan)
library(plyr)
library(ggplot2)
#library(venn)
library(scales)
library(plotly)
source('dag3helper.R')
source('../myLibs/R_Microbiome_scripts.R')
setwd('~/UMCG/DAG3_stats/')

inDF <- read.table('data_processed/DAG3_metaphlan_unfiltered_cleaned_rescaled_nontransformed.tsv',header=T,sep='\t')
inDFa <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled.csv',header=T,sep=',')
phenos <- subsetMicrobiomeDF(inDFa,getPhenos = T,getPWYs = F,getTaxa = F,getVFs = F,getCARDs = F)
inDF <- merge(inDF,phenos,by="DAG3_sampleID")
inDFtaxa <- subsetMicrobiomeDF(inDF,getPhenos = T,getPWYs = F,getTaxa = T,getVFs = F,getCARDs = F)
inDFtaxa$Cohort <- as.character(inDFtaxa$PROVINCE)
inDFtaxa <- inDFtaxa[inDFtaxa$PROVINCE %in% c("Drenthe","Groningen","Friesland"),]
inDFtaxa$Cohort <- as.character(inDFtaxa$Cohort)

# iterate and plot
lvls <- c("Species","Genera","Phyla")
lvlsg <- c("S","G","P")
for (i in c(1:length(lvls))) {
  inDFsub <- filterMetaGenomeDF(inDFtaxa,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,keepLevels = lvlsg[i])
  
  rarAll <- doRarefaction(inDF = inDFsub,doAll = T,bootstraps = 100)
  g <- ggplot(rarAll,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.1) + xlab("Number of sequenced samples") + geom_point(size=2) + 
    ylab(paste0(lvls[i]," identified")) + geom_errorbar(aes(ymin=spec.nr.mn-spec.nr.sd, ymax=spec.nr.mn+spec.nr.sd),size=1,width=1.5) + 
    theme(text = element_text(size=15)) + labs(col = "Province")
  ggsave(plot=g,paste0('plots_general/rarefaction_',lvls[i],'_all.png'),width = 8,height = 6)
  
  rarAllOnly <- rarAll[rarAll$cohort=="All",]
  g <- ggplot(rarAllOnly,aes(x=nr,y=spec.nr.mn)) + geom_line(size=1.1,col="blue") + xlab("Number of sequenced samples") + geom_point(size=2) + 
    ylab(paste0(lvls[i]," identified")) + geom_errorbar(aes(ymin=spec.nr.mn-spec.nr.sd, ymax=spec.nr.mn+spec.nr.sd),size=1,width=1.5) + 
    theme(text = element_text(size=15))
  ggsave(plot=g,paste0('plots_general/rarefaction_',lvls[i],'_allonly.png'),width = 8,height = 6)

  rarProvs <- doRarefaction(inDF = inDFsub,doAll = F,bootstraps = 100)
  g <- ggplot(rarProvs,aes(x=nr,y=spec.nr.mn,col=cohort)) + geom_line(size=1.1) + xlab("Number of sequenced samples") + geom_point(size=2) + 
    ylab(paste0(lvls[i]," identified")) + geom_errorbar(aes(ymin=spec.nr.mn-spec.nr.sd, ymax=spec.nr.mn+spec.nr.sd),size=1,width=1.5) + 
    theme(text = element_text(size=15)) + labs(col = "Province")
  ggsave(plot=g,paste0('plots_general/rarefaction_',lvls[i],'_provonly.png'),width = 8,height = 6)
  
}

# get core microbiome
lvls <- c("Species","Genera","Phyla")
lvlsg <- c("S","G","P")
for (i in c(1:length(lvls))) {
  inDFsub <- filterMetaGenomeDF(inDFtaxa,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,keepLevels = lvlsg[i])
  inDFsub$PROVINCE = "all"
  drk <- getAbundancesForProvince(inDF = inDFsub,getTopX = -1,minAbundance = 0.0000001,province = "all")
  # all abundance
  drk$Taxon <- purgeMGNamesVector(as.character(drk$Taxon))
  drk$Taxon <- as.factor(drk$Taxon)
  drk$Taxon <- factor(drk$Taxon, levels = unique(drk$Taxon[order(drk$Mean.Abundance,decreasing = T)]))
  drk <- drk[order(drk$Mean.Abundance,decreasing = T),]
  ggplot(drk,aes(x=Taxon,y=Mean.Abundance,col=Province,fill=Province)) + geom_bar(stat="identity",position="dodge") 
  
  ggplot(drk,aes(x=Province,y=asin(sqrt(Mean.Abundance)),col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge",col="gray90") + 
    theme(legend.position = "none")
  
  # core abundance (prevalence > 80%)
  drkCore <- drk[drk$Prevalence >= 0.8,]
  g <- ggplot(drkCore,aes(x=Taxon,y=Mean.Abundance,col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge",col="gray90") + 
    theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + theme(legend.position = "none")
  ggsave(plot=g,paste0('plots_general/abundance_core_',lvls[i],'.png'),width = 8,height = 6)
  
  #geom_errorbar(position="dodge",aes(ymin=Mean.Abundance,ymax=Mean.Abundance+SD.Abundance),col="gray40")
  # all prevalence
  drk <- drk[order(drk$Prevalence,decreasing = T),]
  drk$Taxon <- factor(drk$Taxon, levels = unique(drk$Taxon[order(drk$Prevalence,decreasing = T)]))  
  g <- ggplot(drk,aes(x=Province,y=Prevalence,col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge") + 
    theme(legend.position = "none") + xlab("Species") + theme(text = element_text(size=15))
  ggsave(plot=g,paste0('plots_general/prevalence_alltaxa_',lvls[i],'.png'),width = 8,height = 6)
  
  #   core prevalence (prevalence > 80%)
  drkCore <- drk[drk$Prevalence >= 0.8,]
  g <- ggplot(drkCore,aes(x=Taxon,y=Prevalence,col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge",col="gray90") + 
    theme(legend.position = "none") + xlab("Species") + theme(text = element_text(size=12)) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(plot=g,paste0('plots_general/prevalence_core_',lvls[i],'.png'),width = 8,height = 6)

  makePieChart(inDFt = inDFsub[,drkCore$Taxon])
    
  
}

# get core per province
lvls <- c("Species","Genera","Phyla")
lvlsg <- c("S","G","P")
for (i in c(1:length(lvls))) {
  inDFsub <- filterMetaGenomeDF(inDFtaxa,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,keepLevels = lvlsg[i])
  drk <- getAbundancesForProvince(inDF = inDFsub,getTopX = -1,minAbundance = 0.0000001)
  # top abundance
  drk$Taxon <- purgeMGNamesVector(as.character(drk$Taxon))
  drk$Taxon <- as.factor(drk$Taxon)
  drk$Taxon <- factor(drk$Taxon, levels = unique(drk$Taxon[order(drk$Mean.Abundance,decreasing = T)]))
  ggplot(drk,aes(x=Taxon,y=Mean.Abundance,col=Province,fill=Province)) + geom_bar(stat="identity",position="dodge")
  
  ggplot(drk,aes(x=Province,y=Mean.Abundance,col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge",col="gray90")# + 
    #geom_errorbar(position="dodge",aes(ymin=Mean.Abundance,ymax=Mean.Abundance+SD.Abundance),col="gray40")
  # top prevalence
  drk$Taxon <- factor(drk$Taxon, levels = unique(drk$Taxon[order(drk$Prevalence,decreasing = T)]))  
  ggplot(drk,aes(x=Province,y=Prevalence,col=Taxon,fill=Taxon)) + geom_bar(stat="identity",position="dodge",col="gray90")
  
}

