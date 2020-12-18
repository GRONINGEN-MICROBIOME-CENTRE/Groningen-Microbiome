# ========================================================================
# By: R.Gacesa, Weersma Group, UMCG (2020)
#
# Script plots results of heritability analysis
# (Panels A & B in Figure 2 in main DMP manuscript )

# =========================================================================
library(ggplot2)
library(tidyr)

# NOTES:
# ======================
# > requires prepFilesForGraphphlan.R to be run on parsed file and results file from poly run
# > requires set of R functions developed by Weersma Group, UMCG (R_Misc.R and R_Microbiome_scripts.R)
#    - these are avaliable from 

# load log & results file
# set working directory: NOTE this has to be set to appropriate path!
setwd('# placeholder path DAG3_heritability - set this one to path with this file! #')

source('../r_scripts_library/R_Misc.R')
source('../r_scripts_library/R_Microbiome_scripts.R')

if (!dir.exists('Plots')) {
  dir.create('Plots')
}

# input: heritability results
inDFm <- read.table('data/Heritablity_results_table_taxaonly.csv',sep=',',header=T,quote = '"',fill = T,stringsAsFactors = F)
# input: ped and dat files for POLY input
inDataHdr <- read.table('data/_dag3_clrnorm_ziz_rew_trimmed.dat',sep='\t',header=F)
inDataHdr$V1 <- gsub('T  ','',inDataHdr$V1)
inDataHdr$V1 <- gsub('C  ','',inDataHdr$V1)
inDataHdr$V1 <- gsub(' ','',inDataHdr$V1)
inData <- read.table('data/_dag3_clrnorm_ziz_rew_trimmed.ped',sep='\t',header=F)
colnames(inData) <- c('Family_ID','Subject_ID','Father_ID','Mother_ID','Sex',inDataHdr$V1)
inData$`EEND-OF-DATA` <- NULL

summs <- NULL
for (cn in colnames(inData)[grep('k__',colnames(inData))]) {
  cnp <- purgeMGNameOne(cn)
  print(cnp)
  sumOne <- makeOneSummary(inData[[cn]])
  sumOne$Var <- cn
  summs <- rbind.data.frame(summs,sumOne)
}
summs$Var <- gsub('__','_',summs$Var)

inDFm <- merge(inDFm,summs,by.x="Taxon",by.y="Var")

inDFm$FTYPE <- NA
inDFm$FTYPE[grep('^s_',inDFm$Taxon_Shortname)] <-"Taxon.S"
inDFm$FTYPE[grep('^g_',inDFm$Taxon_Shortname)] <-"Taxon.G"
inDFm$FTYPE[grep('^f_',inDFm$Taxon_Shortname)] <-"Taxon.F"
inDFm$FTYPE[grep('^c_',inDFm$Taxon_Shortname)] <-"Taxon.C"
inDFm$FTYPE[grep('^p_',inDFm$Taxon_Shortname)] <-"Taxon.P"
inDFm$FTYPE[grep('^k_',inDFm$Taxon_Shortname)] <-"Taxon.K"
inDFm$FTYPE[grep('^o_',inDFm$Taxon_Shortname)] <-"Taxon.O"

# debug output
#print(inDFm$FEATURE[is.na(inDFm$FTYPE)])

fTypes <- c("Taxon.S","Taxon.G","Taxon.F","Taxon.C","Taxon.O","Taxon.P")

cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for (oneType in fTypes) {
  inDFs <- inDFm[inDFm$FTYPE==oneType,]
  inDFs <- inDFs[inDFs$Persons >= 300,]
  inDFs$FDR <- p.adjust(inDFs$pvalue)
  inDFs$Taxon <- factor(as.character(inDFs$Taxon),levels=inDFs$Taxon[order(inDFs$H2)])
  inDFs <- inDFs[order(inDFs$H2),]
  inDFs$LBL = ""
  inDFs$LBL[inDFs$pvalue <= 0.05] <- "*"
  inDFs$LBL[inDFs$FDR <= 0.1] <- "**"
  if (grepl('Taxon',oneType)) {
    inDFs$Taxon_Shortname <- as.character(inDFs$Taxon_Shortname)
    inDFs$Taxon_Shortname <- paste0(inDFs$Taxon_Shortname,' [N=',inDFs$NrNonZero,']',sep='')
    inDFs$Taxon_Shortname <- factor(as.character(inDFs$Taxon_Shortname),levels=inDFs$Taxon_Shortname[order(inDFs$H2)])
  }
  inDFtoPlot <- inDFs[,c("Taxon_Shortname","LBL","H2")]
  inDFtoPlot$VarEnv <- 1-inDFtoPlot$H2
  inDFtoPlot$SEH2 <- inDFs$SEH2
  inDFtoPlot$LBL <- inDFs$LBL
  inDFtoPlotL <- gather(inDFtoPlot,"Var.Exp","Var.Exp.NR", H2:VarEnv,factor_key = T)
  
  inDFtoPlotL$Var.Exp <- as.character(inDFtoPlotL$Var.Exp)
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "H2"] <- "Genetics"
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VarEnv"] <- "Environment"
  
  inDFtoPlotL$Var.Exp <- factor(as.character(inDFtoPlotL$Var.Exp),level = c("Environment","Genetics"))
  
  inDFtoPlotL$LBL[inDFtoPlotL$Var.Exp=="Environment"] <- ""
  
  inDFtoPlotL$se.x.min <- 0
  inDFtoPlotL$se.x.min[inDFtoPlotL$Var.Exp=="Genetics"] <- pmax(0,inDFtoPlotL$Var.Exp.NR[inDFtoPlotL$Var.Exp=="Genetics"]-inDFtoPlotL$SEH2[inDFtoPlotL$Var.Exp=="Genetics"])
  
  inDFtoPlotL$se.x.max <- 0
  inDFtoPlotL$se.x.max[inDFtoPlotL$Var.Exp=="Genetics"] <- pmin(1,inDFtoPlotL$Var.Exp.NR[inDFtoPlotL$Var.Exp=="Genetics"]+inDFtoPlotL$SEH2[inDFtoPlotL$Var.Exp=="Genetics"])
  
  g <- ggplot(inDFtoPlotL,aes(x=Taxon_Shortname,y=Var.Exp.NR,fill=Var.Exp)) + geom_col(col="black", width=1,size=0.75) + 
    theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,1) +
    theme(axis.text.y = element_text(face="bold")) +
    geom_errorbar(ymin=inDFtoPlotL$se.x.min,ymax=inDFtoPlotL$se.x.max,width=0.25, linetype='solid') + 
    geom_text(data = inDFtoPlotL,
              aes(x = Taxon_Shortname, y=se.x.max, 
                  label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
              color="black", vjust=+0.75, angle = 0, hjust=-1,size=6) + ylim(0,1.0) + 
    ylab('Heritability (H2)') + xlab('') + 
    theme(legend.position="bottom") + 
    theme(text = element_text(size = 14)) + 
    coord_flip()
  print(g)
  ggsave(paste0('Plots/heritability_v2_',oneType,'.png'),height = 1.25+8/50*nrow(inDFs),width = 9,limitsize = F)
}