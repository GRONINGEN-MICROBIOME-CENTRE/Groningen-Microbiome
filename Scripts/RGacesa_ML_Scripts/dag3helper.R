# DAG3 statistics helper functions for various stuff
# =============================================
library(vegan)

# pie chart maker
# ============================================
makePieChart <- function(inDFt,cutLabel=0.01,doLabels=F) {
  avg_DF <- data.frame()
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  for (c in colnames(inDFt)[grep('__',colnames(inDFt))]) {
    avg_DF <- rbind.data.frame(avg_DF, data.frame(Taxon=c,mean=mean(inDFt[[c]]),sd=sd(inDFt[[c]]) ))
  }
  avg_DF$Taxon <- shortenNames2(avg_DF$Taxon)
  avg_DF$lbl <- avg_DF$mean
  avg_DF$NR <- 0
  avg_DF$lbl[avg_DF$lbl < cutLabel] <- 0
  avg_DF$lbl <- percent(avg_DF$lbl)
  avg_DF$lbl[avg_DF$lbl == "0.0%"] <- ""
  avg_DF <- avg_DF[order(avg_DF$mean,decreasing = T),]
  avg_DF$Taxon <- factor(avg_DF$Taxon, levels = avg_DF$Taxon[order(avg_DF$mean,decreasing = T)])
  c = 1
  for (r in c(1:nrow(avg_DF))) {
    if (avg_DF$lbl[r] != "") {
      avg_DF$NR[r] = c
      c = c + 1
    }
  }
  pie <- ggplot(avg_DF,aes(x="",y=mean,col=Taxon,fill=Taxon)) + geom_col() + coord_polar("y", start=0)
  pie <- pie + blank_theme + theme(axis.text.x=element_blank()) 
  if (doLabels) pie <- pie + geom_text(aes(y = mean,x=0.9 + NR %% 2 * 0.3,label = (lbl)), size=5,col="black")
  pie
}


#  =============== getTaxPrevalence ===================================
# -> small function for calculating prevalence of strains in dataset
getTaxPrev <- function(inDF,responseVar,prevalenceMin = 0.0) 
{
  resL = list()
  for (u in unique(inDF[[responseVar]])) {
    #print(u)
    tDF <- inDF[inDF[[responseVar]]==u,grep('__',colnames(inDF))]
    tDF[tDF > 0] <- 1.0
    resL[[u]] <- sort(apply(tDF,MARGIN = 2,FUN=sum),decreasing = T)
    resL[[u]] <- resL[[u]] / sum(inDF[[responseVar]]==u)
    resL[[u]] <- resL[[u]][resL[[u]] > prevalenceMin]
  }
  resL
}

getTaxAvgAbundance <- function(inDF,responseVar,abundanceMin = 0.0) 
{
  resL = list()
  for (u in unique(inDF[[responseVar]])) {
    tDF <- inDF[inDF[[responseVar]]==u,grep('__',colnames(inDF))]
    resL[[u]] <- sort(apply(tDF,MARGIN = 2,FUN=mean),decreasing = T)
    #resL[[u]] <- resL[[u]] / length(resL[[u]])
    resL[[u]] <- resL[[u]][resL[[u]] > abundanceMin]
  }
  resL
}

filterAbundPrev <- function(inDF,minAbundance = 0.001,minPrevalence=0.01,filterPerCohort=F,dropCols = T,abundDropZ=T) {
  # taxonomy
  tDF <- inDF[,grep('__',colnames(inDF))]
  # metadata
  mDF <- inDF[,grep('__',colnames(inDF),invert = T)]
  # prevalence filter
  prev <- colSums(tDF > 0,na.rm = T) / nrow(tDF)
  toKeep <- names(prev[prev >= minPrevalence])
  toDrop <- names(prev[prev < minPrevalence])
  #    drop below prevalence
  if (dropCols) {
    oDF <- tDF[,toKeep]
  } else {
    tDF[,toDrop] <- 0
    oDF <- tDF
  }
  
  # abundance filter
  #tDF <- oDF[,grep('__',colnames(oDF))]
  if (abundDropZ) {
    tDF[tDF==0] <- NA
  }
  #tDF <- oDF
  #abun <- colMeans(tDF,na.rm = T)
  abun <- apply(X = tDF,MARGIN = 2,FUN = function(x) mean(x,na.rm = T) )
  toKeep <- names(abun[abun >= minAbundance])
  toDrop <- names(abun[abun < minAbundance])
  if (abundDropZ) {
    tDF[is.na(tDF)] <- 0.0
  }
  if (dropCols) {
    #oDF <- subset(tDF, select = -c(toDrop))
    #oDF <- tDF[,-c(toDrop)]
    tDF[,c(toDrop)] <- list(NULL)
    oDF <- tDF
  } else {
    tDF[,toDrop] <- 0
    oDF <- tDF
  }
  oDF <- cbind.data.frame(mDF,oDF)
  oDF
}

#rarify (min abundance = 0.001)
doRarefaction <- function(inDF,replacements=F,doAll=F,
                          steps = c(10,20,30,40,50,75,100,150,200,300,400,500,750,1000,1500,2000,2500,3000,3500,4000,5000,6000,7000,8000,9000,10000),
                          bootstraps=10,extrapolate=F)
{
  allCohTaxFt <- inDF[,grep('__',colnames(inDF))]
  allCohTaxFt$Cohort <- inDF$Cohort
  boots <- bootstraps
  rareDf2 <- data.frame()
  cohorts = unique(inDF$Cohort)
  if (doAll) {cohorts <- c("All",cohorts)}
  #for (n in seq(10,nrow(allCohTaxF),100)) {
  for (coh in cohorts) {
    print (paste(" >> Rarefying",coh))
    for (n in steps) {
      if (coh != "All" & n > sum(allCohTaxFt$Cohort==coh) & !extrapolate) {
        break()
      } else if (!extrapolate & n > length(allCohTaxFt$Cohort)) {
        break()
      }
      print(n)
      mns <- c()
      for (b in c(1:boots)) {
        if (coh != "All") {
          t <- allCohTaxFt[allCohTaxFt$Cohort==coh,]
          if (extrapolate & n > length(t$Cohort)) {
            smp <- t[sample(nrow(t), n,replace = T), ]
          } else {
            smp <- t[sample(nrow(t), n,replace = replacements), ]
          }
        } else {
          if (extrapolate & n > length(t$Cohort)) {
            smp <- t[sample(nrow(t), n,replace = T), ]
          }
          else {
            smp <- t[sample(nrow(t), n,replace = replacements), ]
          }
        }
        smp$Cohort <- NULL
        csums <- colSums(smp)
        nrSpec <- sum(csums>0)
        mns <- c(mns,nrSpec)
      }
      mn <- mean(mns)
      sdev <- sd(mns)
      rareDf2 <- rbind.data.frame(rareDf2,data.frame(nr=n,spec.nr.mn=mn,spec.nr.sd=sdev,cohort=coh))
    }
  }
  rareDf2
}

# example:
# k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii
prepData <- function(inFile,getLevel="S",isHMP2=F,minBac=95,doDropZero=T,doRescale=T) {
  if (isHMP2) {taxRaw <- as.data.frame(read.table(inFile,sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"'))
  } else {taxRaw <- as.data.frame(t(read.table(inFile,sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"')))}
  if (doDropZero) {
    taxF <- taxRaw[rowSums(taxRaw) > 0,]
  } else {
    taxF <- taxRaw
  }
  taxF <- taxF[taxRaw$k__Bacteria >= minBac,]
  # take only bacteria
  taxF <- taxF[,grep('k__Bacteria',colnames(taxF))]
  # take only species
  taxF <- taxF[,grep('t__',colnames(taxF),invert = T)]
  taxF.spec <- taxF[,grep('s__',colnames(taxF))]
  # take only genera
  taxF <- taxF[,grep('s__',colnames(taxF),invert = T)]
  taxF.gen <- taxF[,grep('g__',colnames(taxF))]
  # take only families
  taxF <- taxF[,grep('g__',colnames(taxF),invert = T)]
  taxF.family <- taxF[,grep('f__',colnames(taxF))]
  # take only orders
  taxF <- taxF[,grep('f__',colnames(taxF),invert = T)]
  taxF.order <- taxF[,grep('o__',colnames(taxF))]
  # take only classes
  taxF <- taxF[,grep('o__',colnames(taxF),invert = T)]
  taxF.class <- taxF[,grep('c__',colnames(taxF))]
  # take only phyla
  taxF <- taxF[,grep('c__',colnames(taxF),invert = T)]
  taxF.phyla <- taxF[,grep('p__',colnames(taxF))]
  # rescale
  if (getLevel=="S") {taxF.final <- taxF.spec
  } else if (getLevel=="G") {taxF.final <- taxF.gen
  } else if (getLevel=="F") {taxF.final <- taxF.family
  } else if (getLevel=="O") {taxF.final <- taxF.order
  } else if (getLevel=="C") {taxF.final <- taxF.class
  } else if (getLevel=="P") {taxF.final <- taxF.phyla}
  if (doRescale) {
    taxF.final.rescaled <- taxF.final/rowSums(taxF.final)
  } else {
    taxF.final.rescaled <- taxF.final
  }
  # return
  colnames(taxF.final.rescaled) <- gsub('\\|','.',colnames(taxF.final.rescaled))
  taxF.final.rescaled
}
# prepDataHMP2 <- function(inFile) {
#   taxRaw <- as.data.frame(read.table(inFile,sep='\t',header = T,row.names = 1,stringsAsFactors = F))
#   taxF <- taxRaw[rowSums(taxRaw) > 0 & taxRaw$k__Bacteria >= 95,]
#   taxF <- taxF[,grep('t__',colnames(taxF),invert = T)]
#   taxF <- taxF[,colSums(taxF) > 0]
# 
#   # take only bacteria
#   taxF <- taxF[,grep('k__Bacteria',colnames(taxF))]
#   # take only species
#   taxF.spec <- taxF[,grep('s__',colnames(taxF))]
#   # rescale
#   taxF.spec.rescaled <- taxF.spec/rowSums(taxF.spec)
#   # return
#   taxF.spec.rescaled
# }
shortenNames <- function(inDF) {
  for (c in grep('__',colnames(inDF))) {
    #print (c)
    cn <- colnames(inDF)[c]
    #cns <- unlist(strsplit(cn,'\\|'))
    cns <- unlist(strsplit(cn,'\\.'))
    cnsf <- cns[length(cns)]
    colnames(inDF)[c] <- cnsf
    #print (cnsf)
  }
  inDF
}
shortenNames2 <- function(names) {
  outnames <- c()
  names <- as.character(names)
  cns <- (strsplit(names,'\\.?__'))
  for (c in cns) {
    outnames <- c(outnames,c[length(c)])
  }
  outnames
}


# takes annotated taxonomy and returns
# <getTopX> taxa for <province> (or all Provinces if province = "all" or "")
# ================================================================================
getAbundancesForProvince <- function(inDF, province="all",getTopX=-1) {
  tst <- inDF
  tst$NR <- 1.0
  tst$LOCATION <- NULL
  tstAgg <- aggregate(x = tst[,grep('__',colnames(tst))],by = list(tst$PROVINCE),FUN = function(x) {mean(x,na.rm = T)} )
  tstCnt <- aggregate(x = tst[,"NR"],by = list(tst$PROVINCE),FUN = function(x) {sum(x,na.rm = T)} )
  tstSD <- aggregate(x = tst[,grep('__',colnames(tst))],by = list(tst$PROVINCE),FUN = function(x) {sd(x,na.rm = T)} )
  
  tstAgg$PROVINCE <- tstAgg$Group.1
  tstAgg$Group.1 <- NULL
  tstAgg$NR <- tstCnt$x
  tstSD$PROVINCE <- tstSD$Group.1
  
  tstAggL <- data.frame()
  for (p in tstAgg$PROVINCE) {
    for (c in colnames(tstAgg)[grep('__',colnames(tstAgg))]) {
      tstAggL <- rbind.data.frame(tstAggL,data.frame(Province=p,Taxon=c,
                                                     Mean.Abundance=tstAgg[tstAgg$PROVINCE==p,c],
                                                     SD.Abundance=tstSD[tstSD$PROVINCE==p,c]) )
    }
  }
  
  # filter
  tstAggL$Mean.Abundance <- round(tstAggL$Mean.Abundance,digits = 4)
  tstAggL <- tstAggL[tstAggL$Mean.Abundance > 0,]
  if (!(province == "" | province=="all")) {
    tstAggL <- tstAggL[tstAggL$Province==province,]
  } 
  
  # get only top X
  if (getTopX > 0) {
    res <- data.frame()
    for (g in unique(tstAggL$Province)) {
      tt <- tstAggL[tstAggL$Province==g,]
      tt <- tt[order(tt$Mean.Abundance,decreasing = T)[1:min(nrow(tt),getTopX)],]
      res <- rbind.data.frame(res,tt)
    }
    tstAggL <- res
  }
  
  tstAggL
}


# takes annotated taxonomy and returns
# <getTopX> taxa for <gementee> (or all gementeen if gementee = "all" or "")
# ================================================================================
getAbundancesForGementee <- function(inDF, gementee="all",getTopX=-1,minSamples=3) {
  tst <- inDF
  tst$NR <- 1.0
  tstAgg <- aggregate(x = tst[,grep('__',colnames(tst))],by = list(tst$municipaly),FUN = function(x) {mean(x,na.rm = T)} )
  tstCnt <- aggregate(x = tst[,"NR"],by = list(tst$municipaly),FUN = function(x) {sum(x,na.rm = T)} )
  
  tstAgg$GEMENTEE <- tstAgg$Group.1
  tstAgg$Group.1 <- NULL
  tstAgg$NR <- tstCnt$x
  
  tstSD <- aggregate(x = tst[,grep('__',colnames(tst))],by = list(tst$municipaly),FUN = function(x) {sd(x,na.rm = T)} )
  tstSD$GEMENTEE <- tstSD$Group.1
  
  tstAggL <- data.frame()
  
  for (p in tstAgg$GEMENTEE) {
    print(paste(' >> processing',p))
    for (c in colnames(tstAgg)[grep('__',colnames(tstAgg))]) {
      if (tstAgg$NR[tstAgg$GEMENTEE==p] >= minSamples) {
        tstAggL <- rbind.data.frame(tstAggL,data.frame(Gementee=p,NR=tstAgg$NR[tstAgg$GEMENTEE==p],Taxon=c,
                                                       Mean.Abundance=tstAgg[tstAgg$GEMENTEE==p,c],
                                                       SD.Abundance=tstSD[tstSD$GEMENTEE==p,c]) )
      }
    }
  }
  
  # filter
  tstAggL$Mean.Abundance <- round(tstAggL$Mean.Abundance,digits = 4)
  tstAggL <- tstAggL[tstAggL$Mean.Abundance > 0,]
  if (!(gementee == "" | gementee=="all")) {
    tstAggL <- tstAggL[tstAggL$Gementee==gementee,]
  } 
  
  # get only top X
  if (getTopX > 0) {
    res <- data.frame()
    for (g in unique(tstAggL$Gementee)) {
      tt <- tstAggL[tstAggL$Gementee==g,]
      tt <- tt[order(tt$Mean.Abundance,decreasing = T)[1:min(nrow(tt),getTopX)],]
      res <- rbind.data.frame(res,tt)
    }
    tstAggL <- res
  }
  
  tstAggL
}

# takes annotated taxonomy and returns
# <getTopX> taxa for <gementee> (or all gementeen if gementee = "all" or "")
# ================================================================================
getDiversityForGementee <- function(inDF, gementee="all",divIndex="shannon",minSamples=3) {
  tst <- inDF
  tst$NR <- 1.0
  tst$diversity <- diversity(tst[,grep('__',colnames(tst))],index = divIndex)
  
  tstAgg <- aggregate(x = tst[,"diversity"],by = list(tst$municipaly),FUN = function(x) {mean(x,na.rm = T)} )
  tstCnt <- aggregate(x = tst[,"NR"],by = list(tst$municipaly),FUN = function(x) {sum(x,na.rm = T)} )
  
  tstAgg$GEMENTEE <- tstAgg$Group.1
  tstAgg$Group.1 <- NULL
  tstAgg$NR <- tstCnt$x
  
  tstSD <- aggregate(x = tst[,"diversity"],by = list(tst$municipaly),FUN = function(x) {sd(x,na.rm = T)} )
  tstSD$GEMENTEE <- tstSD$Group.1
  
  res <- tstAgg
  res$diversity.sd <- tstSD$x
  colnames(res) <- c("Diversity","Gementee","Nr.Samples","Diversity.SD")
  res <- res[res$Nr.Samples > minSamples,]
  res
}
# get most abundant taxon

# formatting stuff
format2D <- function(x) {formatC(x, digits = 2, format = "f")}
formatNtoM <- function(x) {
  formatC(x/1000000, digits = 1, format = "f") 
}
formatNtoPerc <- function(x) {paste(formatC(x*100, digits = 1, format = "f"),"%",sep='')} 

# random stuff
addNoise <- function(mtx) {
  if (!is.matrix(mtx)) mtx <- matrix(mtx, byrow = TRUE, nrow = 1)
  random.stuff <- matrix(runif(prod(dim(mtx)), min = -0.01, max = 0.01), nrow = dim(mtx)[1])
  random.stuff + mtx
}
