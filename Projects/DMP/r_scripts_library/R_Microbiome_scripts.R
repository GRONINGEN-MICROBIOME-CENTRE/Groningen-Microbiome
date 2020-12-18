library(coin)
library(vegan)

# Microbiome scripts
# by: Ranko Gacesa, UMCG
# ==============================================
doMicrobiomeAdonis <- function(inDF,
                               idCol = "ID",
                               permNR = 10000,
                               adonisVarsTouse,
                               adonisVarsToSkip = c("ID"),
                               targetData = "Taxa", 
                               pwyType = "MetaCyc",
                               verbose = T,
                               nrThreads = 4,
                               measureTime = F) {
  if (!(targetData %in% c("Taxa","PWYs"))) {
    exit("ERROR: targetData must be 'Taxa' or 'PWYs' ")
  }
  adonisResults <- NULL
  adonisVarsTouse <- adonisVarsTouse[!(adonisVarsTouse %in% adonisVarsToSkip)]
  inPhenos <- inDF
  if (verbose) {
    if (targetData=="Taxa") {
      spDF <- subsetMicrobiomeDF(inDF = inDF,verbose = verbose,getPWYs = F,getVFs = F,getTaxa = T,getCARDs = F,getPhenos = F,
                                 getDivs = F,pwyType = pwyType,getID = T,idName = idCol,idToRowName = F)
    } else if (targetData=="PWYs") {
      spDF <- subsetMicrobiomeDF(inDF = inDF,verbose = verbose,getPWYs = T,getVFs = F,getTaxa = F,getCARDs = F,getPhenos = F,
                                 getDivs = F,pwyType = pwyType,getID = T,idName = idCol,idToRowName = F)
    }
  }
  for (i in adonisVarsTouse) {
    if (verbose) {
      print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
      print ('  >> collecting complete cases')
    }
    if (measureTime) {t1 <- Sys.time()}
    inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,idCol)]
    allDF <- merge(x=inPhenosOneVarID,by.x=idCol,y=spDF,by.y=idCol)
    rownames(allDF) <- allDF$DAG3_sampleID
    allDF[[idCol]] <- NULL
    allDF <- allDF[complete.cases(allDF),]
    av <- allDF[[i]]
    allDF[[i]] <- NULL
    if (verbose) {print ('  >> calculating B/C distance')}
    inBC <- vegdist(allDF,method = "bray")
    if (verbose) {print ('  >> doing adonis')}
    nrRows <- length(av)
    if (length(av) < 3 | length(unique(av)) < 2) {
      if (verbose) {print(paste0(' >> WARNING: ',i,' has no useful data!!'))}
    } else {
      #print(paste0(' NR NAs: ',sum(is.na(av))))
      ad <- adonis(inBC ~ av,permutations=permNR,parallel = nrThreads)
      aov_table <- ad$aov.tab
      # accumulate results
      oneRow <- data.frame(Var=i,
                           NR_nonNA=nrRows,
                           DF=aov_table[1,1],
                           SumsOfSqs=aov_table[1,2],
                           MeanSqs=aov_table[1,3],
                           FModel=aov_table[1,4],
                           R2=aov_table[1,5],
                           pval=aov_table[1,6],
                           FDR.BH=NA,
                           Significant=NA)
      print(oneRow)
      adonisResults <- rbind.data.frame(adonisResults,oneRow)
      if (verbose) {print (paste0('--- ',i,' DONE! ---'))}
      if (measureTime) {t2 <- Sys.time()}
      if (verbose & measureTime) {print(t2-t1)}
    }
  }
  rownames(adonisResults) = adonisResults$Var
  adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
  adonisResults$Significant="No"
  adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"
  adonisResults <- adonisResults[order(adonisResults$pval),]
  adonisResults #write.table(adonisResults,'phenotypes_adonis_univar.csv',sep=',',row.names = F)
  
  # g <- ggplot(adonisResults, aes(reorder(row.names(adonisResults), R2), R2, fill=Significant)) + 
  #   geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  #   ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
  # ggsave(plot = g,filename = 'phenotypes_adonis_univar.png')
}


# PREP DATA (CEDNN)
# ===============================================================
# > load CeDNN microbiome (metaphlan "raw" output / taxa)
# ===============================================================
prepCleanMetaphlan <- function(inPath,
                               unclassifiedPercDrop=95,
                               presenceFilter = -1,
                               minRelativeAbundance = -1,
                               rescaleTaxa = T,
                               keepDomains = c("Archaea","Bacteria"),
                               keepLevels=c("P","C","O","F","G","S"),
                               novogeneIdsClean = T
) {
  inTaxa <- read.table(inPath,sep='\t',header=T)
  # reshuffle id column, transpose, turn to data trame
  rownames(inTaxa) <- inTaxa$ID
  inTaxa$ID <- NULL
  inTaxa <- as.data.frame(t.data.frame(inTaxa))
  inTaxa$ID <- rownames(inTaxa)
  print(paste('Loaded Metaphlan (',nrow(inTaxa),'samples)'))
  # transform sample ID names to real-world ones
  # (Novogene adds lots of stuff to names, we need to get rid of it)
  if (novogeneIdsClean) {
    print ('NOTE: doing cleaning of Novogene IDs (keeping format <AB>_<CD>_<EFG...>)')
    for (i in c(1:nrow(inTaxa))) {
      ss <- strsplit(inTaxa$ID[i],'_')[[1]]
      sss <- paste0(ss[1],'_',ss[2],'_',ss[3])
      inTaxa$ID[i] <- sss
    }
  }
  # get rid of duplicates (Novogene randomly sequenced some samples mores then once)
  if (sum(duplicated(inTaxa$ID) > 0)) {
    print(paste('WARNING: found ',sum(duplicated(inTaxa$ID) > 0),'duplicates, dropping them!'))
  }
  inTaxa <- inTaxa[!duplicated(inTaxa$ID),]
  rownames(inTaxa) <- NULL
  # drop failed samples (ideally these would not exist)
  if (sum(inTaxa$unclassified > unclassifiedPercDrop)) {
    print(paste('WARNING: found ',sum(inTaxa$unclassified > unclassifiedPercDrop),'bad samples (with unclassified reads >',unclassifiedPercDrop,'%), dropping them!'))
  }
  inTaxa <- inTaxa[inTaxa$unclassified < unclassifiedPercDrop,]
  # get rid of "nothing" taxon
  inTaxa$unclassified <- NULL
  # more fixing of ID column
  rownames(inTaxa) <- inTaxa$ID
  inTaxa$ID <- NULL
  # filter metaphlan: keep Bacteria and Archaea, drop strains, rescale relative abundances to make sure each level sums to 1
  inTaxaF <- filterMetaGenomeDF(inTaxa,presPerc = presenceFilter,minMRelAb=minRelativeAbundance,minMedRelAb = -1,rescaleTaxa = rescaleTaxa,verbose = T,keepDomains = keepDomains,
                                keepLevels = keepLevels)
  print(paste('Done, returning ',nrow(inTaxaF),'samples'))
  inTaxaF$ID <- rownames(inTaxaF)
  inTaxaF
}

# > load CeDNN microbiome (humann pathways, MetaCyc)
#TODO: this should also be done for other humann data layers (GO, InfoGO, KEGG, Pfam, ...)
# ===============================================================
prepCleanHumann <- function(inPath,
                            dropUnintegrated = T,
                            dropUnmapped = T,
                            dropTaxonSpecific = T,
                            presenceFilter = -1,
                            minRelativeAbundance = -1,
                            rescaleTaxa = T,
                            novogeneIdsClean = T
) {
  inDF <- read.table(inPath,sep='\t',header=T,quote ='',comment.char = '')
  #fix pathway ID (these tend to be weird coming out of humann)
  rownames(inDF) <- inDF$X..Pathway
  inDF$X..Pathway <- NULL
  #drop "junk" from humann (unintegrated/unmapped data & taxon-specific pathways)
  if (dropTaxonSpecific) {
    inDF <- inDF[grep('\\|',rownames(inDF),invert = T),]
  }
  if (dropUnintegrated) {
    inDF <- inDF[grep('UNINTEGRATED',rownames(inDF),invert = T),]
  }
  if (dropUnmapped) {
    inDF <- inDF[grep('UNMAPPED',rownames(inDF),invert = T),]
  }
  rownames(inDF)[grep('PWY',rownames(inDF),invert = T)] <- paste0('PWY_',rownames(inDF)[grep('PWY',rownames(inDF),invert = T)])
  inDF <- as.data.frame(t.data.frame(inDF))
  # fix sample IDs, remove duplicates
  inDF$ID <- rownames(inDF)
  if (novogeneIdsClean) {
    print ('NOTE: doing cleaning of Novogene IDs (keeping format <AB>_<CD>_<EFG...>)')
    for (i in c(1:nrow(inDF))) {
      ss <- strsplit(inDF$ID[i],'_')[[1]]
      sss <- paste0(ss[1],'_',ss[2],'_',ss[3])
      inDF$ID[i] <- sss
    }
  }
  if (sum(duplicated(inDF$ID) > 0)) {
    print(paste('WARNING: found ',sum(duplicated(inDF$ID) > 0),'duplicates, dropping them!'))
  }
  inDF <- inDF[!duplicated(inDF$ID),]
  rownames(inDF) <- inDF$ID
  inDF$ID <- NULL
  # make sure columns are actually numbers (otherwise filter dies; NOTE: this should not be necessary, but... )
  for (c in colnames(inDF)) {inDF[[c]] <- as.numeric(inDF[[c]])}
  # clean, rescale and save
  inDFt2 <- filterHumannDF(inDF,presPerc = presenceFilter,minMRelAb = minRelativeAbundance,minMedRelAb = -1,rescale = T,minSum = 1,verbose = T)
  inDFt2 <- inDFt2[,colSums(inDFt2)!=0]
  inDFt2$ID <- row.names(inDFt2)
  print(paste('Done, returning ',nrow(inDFt2),'samples'))
  inDFt2
}


purgeMGNames <- function(dF) {
  cn <- colnames(dF)
  cnt = 1
  for (i in cn) {
    if (grepl('k__',i)) {
      cn[cnt] <- strsplit(x = i,split = '\\.')[[1]][length(strsplit(x = i,split = '\\.')[[1]])]
    } else {
      cn[cnt] <- i
    }
    cnt = cnt + 1
  }
  colnames(dF) <- cn
  dF
}

purgeMGNameOne <- function(oneName) {
  if (grepl('\\.',oneName)) {
    cn <- strsplit(x = oneName,split = '\\.')[[1]][length(strsplit(x = oneName,split = '\\.')[[1]])]
  } else if (grepl('\\|',oneName)) {
    cn <- strsplit(x = oneName,split = '\\|')[[1]][length(strsplit(x = oneName,split = '\\|')[[1]])]
  } else {cn <- oneName}
  cn
}

purgeMGNamesRaw <- function(dF) {
  cn <- colnames(dF)
  cnt = 1
  for (i in cn) {
    if (grepl('k__',i)) {
      cn[cnt] <- strsplit(x = i,split = '\\|')[[1]][length(strsplit(x = i,split = '\\|')[[1]])]
    } else {
      cn[cnt] <- i
    }
    cnt = cnt + 1
  }
  colnames(dF) <- cn
  dF
}


# =========== parser for cor function, returns dataframe =============
# ====================================================================
findCorrelationsR <- function(corMat,cutoff=0.5) {
  corsDF <- data.frame()
  for (i in 1:nrow(corMat)){
    cors <- which(corMat[i,] > cutoff & corMat[i,] != 1.0)
    if (length(cors) > 0) {
      for (j in cors) {
        corsDF <- rbind.data.frame(corsDF,data.frame(A=i,B=j,AName=colnames(corMat)[i],
                                                     BName=rownames(corMat)[j],cor=corMat[i,j]))
      }
    }
  }
  corsDF <- corsDF[order(corsDF$cor,decreasing = T),]
  corsDF
}
# ================================================================================
# function for filtering pathway results (HUMANN)
# - takes dataframe (any)
# - converts NAs to 0
# - filters for relative abundances (median & mean)
#
#
# ================================================================================
filterHumannDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,minMedRelAb=0.0,minSum=90.0, rescale=T,verbose=T,type='MetaCyc') {
  
  nonPWYpwys <- c("ARG+POLYAMINE-SYN: superpathway of arginine and polyamine biosynthesis",
     "CHLOROPHYLL-SYN: chlorophyllide a biosynthesis I (aerobic, light-dependent)",
     "GLYCOLYSIS-E-D: superpathway of glycolysis and Entner-Doudoroff",
     "GLYCOLYSIS-TCA-GLYOX-BYPASS: superpathway of glycolysis, pyruvate dehydrogenase, TCA, and glyoxylate bypass",
     "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
     "GLYOXYLATE-BYPASS: glyoxylate cycle",
     "HEME-BIOSYNTHESIS-II: heme biosynthesis I (aerobic)",
     "MANNOSYL-CHITO-DOLICHOL-BIOSYNTHESIS: protein N-glycosylation (eukaryotic, high mannose)",
     "NAD-BIOSYNTHESIS-II: NAD salvage pathway II",                  
     "REDCITCYC: TCA cycle VIII (helicobacter)",
     "TCA-GLYOX-BYPASS: superpathway of glyoxylate bypass and TCA",
     "TCA: TCA cycle I (prokaryotic)")
  
  colnames(inDF)[colnames(inDF) %in% nonPWYpwys] <- paste0('PWY_',colnames(inDF)[colnames(inDF) %in% nonPWYpwys])
  
  if (type=='MetaCyc') {
    nonPWYdf <- as.data.frame(inDF[,-grep('PWY',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^EC_',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    nonPWYdf <- as.data.frame(inDF[,-grep('RXN',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^PF[01]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^GO',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    nonPWYdf <- as.data.frame(inDF[,-grep('^K[012]',colnames(inDF))])
    cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^K[012]',colnames(inDF))] ])
  }
  colnames(nonPWYdf) <- cnsNonPWYdf
  if (type=='MetaCyc') {
    yesPWYdf <- as.data.frame(inDF[,grep('PWY',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('PWY',colnames(inDF))] ])
  } else if (type=='EC') {
    yesPWYdf <- as.data.frame(inDF[,grep('^EC_',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^EC_',colnames(inDF))] ])
  } else if (type=='RXN') {
    yesPWYdf <- as.data.frame(inDF[,grep('RXN',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('RXN',colnames(inDF))] ])
  } else if (type=='PFAM') {
    yesPWYdf <- as.data.frame(inDF[,grep('^PF[01]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^PF[01]',colnames(inDF))] ])
  } else if (type=='GO') {
    yesPWYdf <- as.data.frame(inDF[,grep('^GO',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^GO',colnames(inDF))] ])
  } else if (type=='KEGG') {
    yesPWYdf <- as.data.frame(inDF[,grep('^K[012]',colnames(inDF))])
    cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^K[012]',colnames(inDF))] ])
  }
  
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }

  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }

  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'pathways!, ',length(colnames(yesPWYdf)),'pathways left!')); }
  
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

# ================================================================================
# function for filtering CARD results
# - takes dataframe (any)
# - converts NAs to 0
# - filters for relative abundances (median & mean) and prevalence (0 vs non-0)
#================================================================================
filterCardDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,grepPat="'gb\\|'",minMedRelAb=0.0,rescale=T,verbose=T) {
  nonPWYdf <- as.data.frame(inDF[,-grep(grepPat,colnames(inDF))])
  cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep(grepPat,colnames(inDF))] ])
  colnames(nonPWYdf) <- cnsNonPWYdf
  yesPWYdf <- as.data.frame(inDF[,grep(grepPat,colnames(inDF))])
  cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep(grepPat,colnames(inDF))] ])
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'CARDs!, ',length(colnames(yesPWYdf)),'CARDs left!')); }
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'CARDs!, ',length(colnames(yesPWYdf)),'CARDs left!')); }
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'CARDs!, ',length(colnames(yesPWYdf)),'CARDs left!')); }
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

# ================================================================================
# function for filtering VFDB results
# - takes dataframe (any)
# - converts NAs to 0
# - filters for relative abundances (median & mean) and prevalence (0 vs non-0)
#================================================================================
filterVfdbDF <- function(inDF,presPerc = 0.05,minMRelAb = 0.001,minMedRelAb=0.0,rescale=T,verbose=T) {
  nonPWYdf <- as.data.frame(inDF[,-grep('^VF',colnames(inDF))])
  cnsNonPWYdf <- colnames(inDF[colnames(inDF)[-grep('^VF',colnames(inDF))] ])
  colnames(nonPWYdf) <- cnsNonPWYdf
  yesPWYdf <- as.data.frame(inDF[,grep('^VF',colnames(inDF))])
  cnsYesPWYdf <- colnames(inDF[colnames(inDF)[grep('^VF',colnames(inDF))] ])
  # replaces NAs with 0s
  for (c in colnames(yesPWYdf)) {
    yesPWYdf[,c][is.na(yesPWYdf[,c])] <- 0.0
  }
  # rescale to rel ab (if rescale = T)
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    nrnZ = as.numeric(sum(yesPWYdf[,c]!=0.0))
    if (nrnZ/as.numeric(nrow(yesPWYdf)) < presPerc) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'VFs!, ',length(colnames(yesPWYdf)),'VFs left!')); }
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = mean(yesPWYdf[,c])
    if ( mn < minMRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'VFs!, ',length(colnames(yesPWYdf)),'VFs left!')); }
  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in colnames(yesPWYdf)) {
    mn = median(yesPWYdf[,c])
    if ( mn < minMedRelAb) {
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    yesPWYdf <- yesPWYdf[,!(colnames(yesPWYdf) %in% toRemove)]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'CARDs!, ',length(colnames(yesPWYdf)),'VFs left!')); }
  # do final rescale
  if (rescale==T) {
    if (verbose) {print ('  >> rescaling')}
    rsums <- rowSums(yesPWYdf)
    rsums[rsums==0] <- 1.0
    yesPWYdf <- yesPWYdf/rsums
  }
  inDF <- cbind.data.frame(nonPWYdf,yesPWYdf)
  if (verbose) {print ('> DONE')}
  inDF
}

# ====================================================
# > reorders phenotypes VS data (includes: taxa, pwys, CARD, VFs)
# ====================================================
reorderMicrobiomeDF <- function(inDF,verbose=T) {
  if (verbose) {print ("RE-ORDERING PHENOTYPES VS DATA ...")}
  nonPhenoCols <- c()
  nonPhenoColsTaxa <- c(grep("^[dgtspcfko]__",colnames(inDF)))
  if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsTaxa),' TAXA'))}
  nonPhenoColsPWYs <- c(grep("PWY",colnames(inDF)))
  if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsPWYs),' PWYs'))}
  nonPhenoColsVFs <- c(grep("^VF",colnames(inDF)))
  if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsVFs),' VFs'))}
  nonPhenoColsCARDs <- c(grep("CARD",colnames(inDF)))
  if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsCARDs),' CARDs'))}
  nonPhenoCols <- c(nonPhenoCols,nonPhenoColsTaxa,nonPhenoColsPWYs,nonPhenoColsVFs,nonPhenoColsCARDs)
  phenoCols <- c(1:ncol(inDF)); phenoCols <- phenoCols[!(phenoCols %in% nonPhenoCols)]
  if (verbose) {print(paste0(' > FOUND ',length(phenoCols),' *Phenotypes*'))}
  inDF <- inDF[, c(phenoCols,nonPhenoCols)]
}

# ====================================================
# > Subsets stuff from merged microbiome dataframe
# ====================================================
subsetMicrobiomeDF <- function(inDF,verbose=T,getPWYs=F,getVFs=F,getTaxa=T,getCARDs=F,getPhenos=T,getDivs=F,pwyType="MetaCyc",getID=T,idName="ID",idToRowName=F) {
  if (!(pwyType %in% c('MetaCyc','EC','RXN','PFAM','GO','KEGG'))) {
    exit("ERROR: pwyType must be one of 'MetaCyc','EC','RXN','PFAM','GO','KEGG' ")
  }
  if (verbose) {print ("SUBSETTING Microbiome Dataframe  ...")
    if (getPhenos) {print("  > getting phenotypes")}
    if (getTaxa) {print("  > getting taxa")}
    if (getPWYs) {print(paste0("  > getting pathways [",pwyType,"]"))}
    if (getCARDs) {print("  > getting CARDs")}
    if (getVFs) {print("  > getting VFs")}
    if (getDivs) {print("  > getting Diversities")}
    if (getID) {print(paste0("  > getting ID column = ",idName))}
  }
  nonPhenoCols <- c()
  nonPhenoColsTaxa <- c(grep("^[dgtspcfko]__",colnames(inDF)))
  toGet <- c()
  if (getTaxa) {
    toGet <- nonPhenoColsTaxa
    if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsTaxa),' TAXA'))}
  }
  nonPhenoColsPWYs <- c()
  if (pwyType=='MetaCyc' | pwyType=='All' | pwyType=='Metacyc') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("PWY",colnames(inDF)))
    nonPhenoColsPWYs <- nonPhenoColsPWYs[!(nonPhenoColsPWYs %in% grep('^DIV\\.',colnames(inDF)))]
  }
  if (pwyType=='EC' | pwyType=='All') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("^EC_",colnames(inDF)))
  }
  if (pwyType=='RXN' | pwyType=='All') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("RXN",colnames(inDF)))
  }
  if (pwyType=='PFAM' | pwyType=='All' | pwyType=='Pfam') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("^PF[01]",colnames(inDF)))
  }
  if (pwyType=='GO' | pwyType=='All') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("^GO",colnames(inDF)))
  }
  if (pwyType=='KEGG' | pwyType=='All' | pwyType=='Kegg') {
    nonPhenoColsPWYs <- c(nonPhenoColsPWYs,grep("^K[012]",colnames(inDF)))
  }
  nonPhenoColsPWYs <- unique(nonPhenoColsPWYs)
  if (getPWYs) {
    toGet <- c(toGet,nonPhenoColsPWYs)
    if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsPWYs),' PWYs'))}
  }
  nonPhenoColsVFs <- c(grep("^VF",colnames(inDF)))
  if (getVFs) {
    toGet <- c(toGet,nonPhenoColsVFs)
    if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsVFs),' VFs'))}
  }
  nonPhenoColsCARDs <- c(grep("^CARD",colnames(inDF)))
  if (getCARDs) {
    toGet <- c(toGet,nonPhenoColsCARDs)
    if (verbose) {print(paste0(' > FOUND ',length(nonPhenoColsCARDs),' CARDs'))}
  }
  divCols <- grep("^DIV\\.",colnames(inDF))
  if (getDivs) {
    if (verbose) {print(paste0(' > FOUND ',length(divCols),' *Diversity metrics*'))}
    toGet <- c(toGet,divCols)
  }
  idCols <- grep(paste0("^",idName,"$"),colnames(inDF))
  if (getID) {
    toGet <- unique(c(toGet,idCols))
    if (verbose) {print(paste0(' > FOUND ',length(idCols),' *ID COLUMNS*'))}
  }
  nonPhenoCols <- c(nonPhenoCols,nonPhenoColsTaxa,nonPhenoColsPWYs,nonPhenoColsVFs,nonPhenoColsCARDs,divCols,idCols)
  phenoCols <- c(1:ncol(inDF)); phenoCols <- phenoCols[!(phenoCols %in% nonPhenoCols)]
  if (getPhenos) {
    toGet <- c(phenoCols,toGet)
    if (verbose) {print(paste0(' > FOUND ',length(phenoCols),' *Phenotypes*'))}
  }
  inDF <- inDF[, unique(toGet)]
  if (idToRowName) {
    rownames(inDF) <- inDF[[idName]]
    inDF[[idName]] <- NULL
  }
  inDF
}

# ================================================================================
# calculator for diversity metrics
# ================================================================================
calcDIVMetrics <- function(inDF,IDcol="RowNames",pwyType="MetaCyc",
                           metrics=c("shannon","simpson","invsimpson","richness"),
                           DIVlvls=c("taxS","taxG","taxF","taxO","taxC","taxP","PWY","CARD","VF") ) {
  # select IDs column
  if (IDcol == "RowNames") {
    DIVMatrix <- data.frame(RN=rownames(inDF))
  } else {
    DIVMatrix <- data.frame(IDcol=inDF[[IDcol]])
    colnames(DIVMatrix)[1] <- IDcol
  }
  # iterate over stuff to calculate
  #  > diversity levels (taxa per taxon lvl, pathways, VFs, CARDs)
  for (l in DIVlvls) {
    # TAXA
    # ================================
    if (grepl('^tax.$',l)) {
      toUse <- gsub('^tax','',l)
      inDFf <- subsetMicrobiomeDF(inDF,getTaxa = T,getPhenos = F,getVFs = F,getPWYs = F,getCARDs = F,pwyType = pwyType,getID = T)
      if (IDcol=="RowNames") {
      } else {
        rownames(inDFf) <- inDFf[[IDcol]]
        inDFf[[IDcol]] <- NULL
      }
      inDFf <- filterMetaGenomeDF(inDFf,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                  keepDomains = "All",keepLevels = c(toUse))
      # iterate over metrics, calculate each
      # NOTE: richness is not implemented in vegan, requires special treatment
      for (m in metrics) {
        print(paste0('  > calculating ',m,'[',l,']'))
        if (m=="richness") {
          inDFpa <- inDFf
          inDFpa[inDFpa > 0] <- 1
          dv <- rowSums(inDFpa)
        } else {
          dv <- diversity(inDFf,index = m)
        }
        DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
        colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.',toUse,'.',m)
      }
      # PWYs
      # =======================
    } else if (l=="PWY") {
      inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = F,getPWYs = T,getCARDs = F,getID = T,pwyType = pwyType)
      if (IDcol=="RowNames") {
      } else {
        rownames(inDFf) <- inDFf[[IDcol]]
        inDFf[[IDcol]] <- NULL
      }
      for (m in metrics) {
        print(paste0('  > calculating ',m, '[PWYs,',pwyType,']'))
        if (m=="richness") {
          inDFpa <- inDFf
          inDFpa[inDFpa > 0] <- 1
          dv <- rowSums(inDFpa)
        } else {
          dv <- diversity(inDFf,index = m)
        }
        DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
        colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.PWY.',m)
      }
      # CARD
      # ==========================
    } else if (l=="CARD") {
      inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = F,getPWYs = F,getCARDs = T,getID = T)
      if (IDcol=="RowNames") {
      } else {
        rownames(inDFf) <- inDFf[[IDcol]]
        inDFf[[IDcol]] <- NULL
      }
      for (m in metrics) {
        print(paste0('  > calculating ',m, '[CARDs]'))
        if (m=="richness") {
          inDFpa <- inDFf
          inDFpa[inDFpa > 0] <- 1
          dv <- rowSums(inDFpa)
        } else {
          dv <- diversity(inDFf,index = m)
        }
        DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
        colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.CARD.',m)
      }
      # VFs
      # =======================
    } else if (l=="VF") {
      inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = T,getPWYs = F,getCARDs = F,getID = T)
      if (IDcol=="RowNames") {
      } else {
        rownames(inDFf) <- inDFf[[IDcol]]
        inDFf[[IDcol]] <- NULL
      }
      for (m in metrics) {
        print(paste0('  > calculating ',m, '[VFs]'))
        if (m=="richness") {
          inDFpa <- inDFf
          inDFpa[inDFpa > 0] <- 1
          dv <- rowSums(inDFpa)
        } else {
          dv <- diversity(inDFf,index = m)
        }
        DIVMatrix <- cbind.data.frame(DIVMatrix,dv)
        colnames(DIVMatrix)[ncol(DIVMatrix)] <- paste0('DIV.VF.',m)
      }
    } 
  }
  return(DIVMatrix)
}

# ================================================================================
# ================================================================================
# calculator for bray-curtis metrix
# ================================================================================
calcBCMatrix <- function(inDF,IDcol="RowNames",pwyType="MetaCyc",
                           div="taxS" ) {
  # find and calculate bray-curtis
  #  > diversity levels (taxa per taxon lvl, pathways, VFs, CARDs)
  # TAXA
  # ================================
  l <- div
  if (grepl('^tax.$',l)) {
    toUse <- gsub('^tax','',l)
    inDFf <- subsetMicrobiomeDF(inDF,getTaxa = T,getPhenos = F,getVFs = F,getPWYs = F,getCARDs = F,pwyType = pwyType,getID = T)
    if (IDcol=="RowNames") {
    } else {
      rownames(inDFf) <- inDFf[[IDcol]]
      inDFf[[IDcol]] <- NULL
    }
    inDFf <- filterMetaGenomeDF(inDFf,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,rescaleTaxa = T,verbose = T,
                                keepDomains = "All",keepLevels = c(toUse))
    # PWYs
    # =======================
  } else if (l=="PWY") {
    inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = F,getPWYs = T,getCARDs = F,getID = T,pwyType = pwyType)
    if (IDcol=="RowNames") {
    } else {
      rownames(inDFf) <- inDFf[[IDcol]]
      inDFf[[IDcol]] <- NULL
    }
  # CARD
  # ==========================
  } else if (l=="CARD") {
    inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = F,getPWYs = F,getCARDs = T,getID = T)
    if (IDcol=="RowNames") {
    } else {
      rownames(inDFf) <- inDFf[[IDcol]]
      inDFf[[IDcol]] <- NULL
    }
    # VFs
    # =======================
  } else if (l=="VF") {
    inDFf <- subsetMicrobiomeDF(inDF,getTaxa = F,getPhenos = F,getVFs = T,getPWYs = F,getCARDs = F,getID = T)
    if (IDcol=="RowNames") {
    } else {
      rownames(inDFf) <- inDFf[[IDcol]]
      inDFf[[IDcol]] <- NULL
    }
  }
  # actually calc B/C
  dv <- vegdist(inDFf,method = "bray")
  return(dv)
}
# ================================================================================


# ================================================================================
# ================================================================================
# function for filtering microbiome (metaphlan) results
# - takes dataframe (any)
# - replaces NAs with 0.0
# - filters OUT taxa present in less then presPerc samples
# - filters OUT taxa with mean releative abundance < minMRelAb
# - filters OUT taxa with median relative abundance < minMedRelAb
# - filters OUT taxonomic levels not in keepLevels ()
#   -> keepLevels should be vector of following: T = strain, S = species, G = Genera, F = families
#                                                O = Orders, C = classes, P = phyla, K = Kingdoms
#   -> example: keepLevels=c('S','G') keeps only species and genera
# - filters OUT domains not in keepDomains: 
#   -> input is vector of following: Eukaryota, Bacteria, Viruses, Archaea
#   -> example: keepDomains=c('B') keeps only bacteria, removes viruses, archaea and eukarya
# - if rescaleTaxa = True, rescales relative abundancies after filtering
# returns modified dataframe
# NOTES:
# - assumes metaphlan encoding (k__<kingdom>.o__<order> ... ); it will mess stuff
# up if non-metagenome rows are encoded like this!
# DEFAULTS: 
# - removes non-bacteria
# - keeps all except kingdoms
# ================================================================================

#TODO: implement keepDomains
filterMetaGenomeDF <- function(inDF,presPerc = 0.1,minMRelAb = 0.01,minMedRelAb=0.0, rescaleTaxa=F,verbose=T,
                             keepDomains=c('Bacteria','Archaea'),
                             keepLevels=c('T','S','G','F','O','C','P')) {
  
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  tColsNMG = grep('k__',colnames(inDF),invert = T) # colums with microbiome
  # replaces NAs in microbiome with 0s
  for (c in tCols) {
    inDF[,c][is.na(inDF[,c])] <- 0.0
  }
  # filter for presence
  # -----------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    nrnZ = as.numeric(sum(inDF[,c]!=0.0))
    if ( (nrnZ/as.numeric(nrow(inDF))) < presPerc) {
      # if (verbose) {
      #   print (paste('col',c,': ',colnames(inDF)[c],'; nr non Zero:',nrnZ,'=',nrnZ/as.numeric(nrow(inDF)),'>> Column removed!'))
      # }
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('^[dgtspcfko]__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > presence filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }
  
  # filter for abundance (mean)
  # ---------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = mean(inDF[,c])
    if ( mn < minMRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; mean rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  tCols = grep('k__',colnames(inDF)) # colums with microbiome
  if (verbose) {print (paste(' > mean abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }

  # filter for abundance (median)
  # -----------------------------
  nrRemoved = 0
  toRemove = c()
  for (c in tCols) {
    mn = median(inDF[,c])
    if ( mn < minMedRelAb) {
      #if (verbose) {
      #  print (paste('col',c,': ',colnames(inDF)[c],'; median rel abundance:',mn,' >> Column removed!'))
      #}
      nrRemoved = nrRemoved + 1
      toRemove <- c(toRemove,c)
    }
  }
  if (length(toRemove) > 0) {
    inDF <- inDF[,-toRemove]
  }
  if (verbose) {print (paste(' > median abundance filter: Removed',nrRemoved,'taxa!, ',length(tCols),'taxa left!')); }  
  
  # keep domains
  # -----------------------------
  if (keepDomains!='All' & keepDomains !="") {
  toKeep = grep('k__',colnames(inDF),invert = T)
  for (d in keepDomains) {
    #d print(d)
    toKeep = c(toKeep,grep(paste('k__',d,sep=''),colnames(inDF)))
  }
  inDF <- inDF[,toKeep]
  }

  # remove taxonomic levels
  # -----------------------------
  inDFnonTaxa <- as.data.frame(inDF[,grep('k__',colnames(inDF),invert=T)])
  colnames(inDFnonTaxa) <- colnames(inDF)[grep('k__',colnames(inDF),invert=T)]
  inDF2 <- inDF[,grep('k__',colnames(inDF),invert=F)]
  # pick strains (T)
  taxaTCols <- grep('t__',colnames(inDF2))
  taxaT <- inDF2[,taxaTCols]
  if (length(taxaT) > 0) {inDF2 <- inDF2[,-taxaTCols]}
  # pick species (S)
  taxaSCols <- grep('s__',colnames(inDF2))
  taxaS <- inDF2[,taxaSCols]
  if (length(taxaS) > 0) {inDF2 <- inDF2[,-taxaSCols]}
  # pick genera (G)
  taxaGCols <- grep('g__',colnames(inDF2))
  taxaG <- inDF2[,taxaGCols]
  if (length(taxaG) > 0) {inDF2 <- inDF2[,-taxaGCols]}
  # pick families (F)
  taxaFCols <- grep('f__',colnames(inDF2))
  taxaF <- inDF2[,taxaFCols]
  if (length(taxaF) > 0) {inDF2 <- inDF2[,-taxaFCols]}
  # pick orders (O)
  taxaOCols <- grep('o__',colnames(inDF2))
  taxaO <- inDF2[,taxaOCols]
  if (length(taxaO) > 0) {inDF2 <- inDF2[,-taxaOCols]}
  # pick classes (C)
  taxaCCols <- grep('c__',colnames(inDF2))
  taxaC <- inDF2[,taxaCCols]
  if (length(taxaC) > 0) {inDF2 <- inDF2[,-taxaCCols]}
  # pick phyla (P)
  taxaPCols <- grep('p__',colnames(inDF2))
  taxaPColsKeep <- grep('p__',colnames(inDF2),invert = T)
  taxaPColsKeepN <- colnames(inDF2)[grep('p__',colnames(inDF2),invert = T)]
  taxaP <- inDF2[,taxaPCols]
  if (length(taxaP) > 0) {inDF2 <- as.data.frame(inDF2[,taxaPColsKeep])}
  colnames(inDF2) <- taxaPColsKeepN
  taxaK <- inDF2
  # pick 
  oDF <- inDFnonTaxa
  if (verbose) {print ('Keeping following taxonomic levels:'); print(keepLevels)}
  if ('K' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaK),'Kingdoms'))}
    if (rescaleTaxa) {taxaK <- taxaK/rowSums(taxaK)}
    oDF <- cbind(oDF,taxaK)}
  if ('P' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaP),'Phyla'))} 
    if (rescaleTaxa) {taxaP <- taxaP/rowSums(taxaP)}
    oDF <- cbind(oDF,taxaP)}
  if ('C' %in% keepLevels) {
    if (verbose) {print(paste(' -> kept',ncol(taxaC),'Classes'))} 
    if (rescaleTaxa) {taxaC <- taxaC/rowSums(taxaC)}
    oDF <- cbind(oDF,taxaC)}
  if ('O' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaO),'Orders'))}
    if (rescaleTaxa) {taxaO <- taxaO/rowSums(taxaO)}
    oDF <- cbind(oDF,taxaO)}
  if ('F' %in% keepLevels) {
    if (verbose){print(paste(' -> kept',ncol(taxaF),'Families'))}
    if (rescaleTaxa) {taxaF <- taxaF/rowSums(taxaF)}
    oDF <- cbind(oDF,taxaF)}
  if ('G' %in% keepLevels) {if (verbose){ print(paste(' -> kept',ncol(taxaG),'Genera'))}
    if (rescaleTaxa) {taxaG <- taxaG/rowSums(taxaG)}
    oDF <- cbind(oDF,taxaG)}
  if ('S' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaS),'Species'))}
    if (rescaleTaxa) {taxaS <- taxaS/rowSums(taxaS)}
    oDF <- cbind(oDF,taxaS)}
  if ('T' %in% keepLevels) {if (verbose){print(paste(' -> kept',ncol(taxaT),'Strains'))}
    if (rescaleTaxa) {taxaT <- taxaT/rowSums(taxaT)}
    oDF <- cbind(oDF,taxaT)}
  
  if (verbose) {print ('data processing done, returning Dataframe')}
  oDF
}
# ================================================================================


# ================================================================================
# function that extracts list of taxa from metaphlan file, returns dataframe
# ================================================================================
metaPhlanExtractTaxa <- function(inFileA,taxaList) {
  taxaA=read.table(inFileA, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaA=as.data.frame(t(taxaA))
  taxaA[,taxaList]
}

# ================================================================================
# function that finds species in metaphlan file
# ================================================================================
metaPhlanDFGetTaxa <- function(inDF, Ttresh=0.001,ret="S",FUN=mean) {
  taxaA=as.data.frame(inDF)
  # remove strains
  s_speciesA=taxaA[,grep("t__",colnames(taxaA),invert = TRUE)]
  # remove genera
  s_generaA = s_speciesA[,grep('s__',colnames(s_speciesA),invert = TRUE)]
  # remove families
  s_familyA = s_generaA[,grep('g__',colnames(s_generaA),invert = TRUE)]
  # remove orders
  s_orderA = s_familyA[,grep('f__',colnames(s_familyA),invert = TRUE)]
  # remove classes
  s_classA = s_orderA[,grep('o__',colnames(s_orderA),invert = TRUE)]
  #remove phyla
  s_phylaA = s_classA[,grep('c__',colnames(s_classA),invert = TRUE)]
  # grab species and drop 0s
  s_speciesA=s_speciesA[,grep("s__",colnames(s_speciesA))]
  s_speciesAFN=colnames(s_speciesA[,apply(s_speciesA,2,FUN) >= Ttresh])
  # grab genera and drop 0s
  s_generaA=s_generaA[,grep("g__",colnames(s_generaA))]
  s_generaAFN=colnames(s_generaA[,apply(s_generaA,2,FUN) >= Ttresh])
  # grab families and drop 0s
  s_familyA=s_familyA[,grep("f__",colnames(s_familyA))]
  s_familyAFN=colnames(s_familyA[,apply(s_familyA,2,FUN) >= Ttresh])
  # grab orders and drop 0s
  s_orderA=s_orderA[,grep("o__",colnames(s_orderA))]
  s_orderAFN=colnames(s_orderA[,apply(s_orderA,2,FUN) >= Ttresh])
  # grab classes and drop 0s
  s_classA=s_classA[,grep("c__",colnames(s_classA))]
  s_classAFN=colnames(s_classA[,apply(s_classA,2,FUN) >= Ttresh])
  # grab phyla and drop 0s
  s_phylaA=s_phylaA[,grep("p__",colnames(s_phylaA))]
  s_phylaAFN=colnames(s_phylaA[,apply(s_phylaA,2,FUN) >= Ttresh])
  
  if (ret=="P") s_phylaAFN
  else if (ret=="C") s_classAFN
  else if (ret=="O") s_orderAFN
  else if (ret=="F") s_familyAFN
  else if (ret=="G") s_generaAFN
  else if (ret=="S") s_speciesAFN
  else s_speciesAFN
  
}

# ================================================================================
# find overlapping species between two metaphlan files
# ================================================================================
metaPhlanGetOverlapSp <- function(inFileA,inFileB,Ttresh = 0.001,taxa="S") {
  # -- process 1st input (A)
  # ----------------------------
  taxaA=read.table(inFileA, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaA=as.data.frame(t(taxaA))
  # remove strains
  s_speciesA=taxaA[,grep("t__",colnames(taxaA),invert = TRUE)]
  # remove genera
  s_generaA = s_speciesA[,grep('s__',colnames(s_speciesA),invert = TRUE)]
  # remove families
  s_familyA = s_generaA[,grep('g__',colnames(s_generaA),invert = TRUE)]
  # remove orders
  s_orderA = s_familyA[,grep('f__',colnames(s_familyA),invert = TRUE)]
  # remove classes
  s_classA = s_orderA[,grep('o__',colnames(s_orderA),invert = TRUE)]
  #remove phyla
  s_phylaA = s_classA[,grep('c__',colnames(s_classA),invert = TRUE)]
  # grab species and drop 0s
  s_speciesA=s_speciesA[,grep("s__",colnames(s_speciesA))]
  s_speciesAFN=colnames(s_speciesA[,apply(s_speciesA,2,mean) > Ttresh])
  # grab genera and drop 0s
  s_generaA=s_generaA[,grep("g__",colnames(s_generaA))]
  s_generaAFN=colnames(s_generaA[,apply(s_generaA,2,mean) > Ttresh])
  # grab families and drop 0s
  s_familyA=s_familyA[,grep("f__",colnames(s_familyA))]
  s_familyAFN=colnames(s_familyA[,apply(s_familyA,2,mean) > Ttresh])
  # grab orders and drop 0s
  s_orderA=s_orderA[,grep("o__",colnames(s_orderA))]
  # grab classes and drop 0s
  s_classA=s_classA[,grep("c__",colnames(s_classA))]
  # grab phyla and drop 0s
  s_phylaA=s_phylaA[,grep("p__",colnames(s_phylaA))]
  
  # -- process 2nd input 
  # ----------------------------
  taxaB=read.table(inFileB, sep = "\t", row.names = 1, header = T, check.names = F,blank.lines.skip = T)
  taxaB=as.data.frame(t(taxaB))
  # remove strains
  s_speciesB=taxaB[,grep("t__",colnames(taxaB),invert = TRUE)]
  # remove genera
  s_generaB = s_speciesB[,grep('s__',colnames(s_speciesB),invert = TRUE)]
  # remove families
  s_familyB = s_generaB[,grep('g__',colnames(s_generaB),invert = TRUE)]
  # remove orders
  s_orderB = s_familyB[,grep('f__',colnames(s_familyB),invert = TRUE)]
  # remove classes
  s_classB = s_orderB[,grep('o__',colnames(s_orderB),invert = TRUE)]
  #remove phyla
  s_phylaB = s_classB[,grep('c__',colnames(s_classB),invert = TRUE)]
  # grab species and drop 0s
  s_speciesB=s_speciesB[,grep("s__",colnames(s_speciesB))]
  s_speciesBFN=colnames(s_speciesB[,apply(s_speciesB,2,mean) > Ttresh])
  # grab genera and drop 0s
  s_generaB=s_generaB[,grep("g__",colnames(s_generaB))]
  s_generaBFN=colnames(s_generaB[,apply(s_generaB,2,mean) > Ttresh])
  # grab families and drop 0s
  s_familyB=s_familyB[,grep("f__",colnames(s_familyB))]
  s_familyBFN=colnames(s_familyB[,apply(s_familyB,2,mean) > Ttresh])
  # grab orders and drop 0s
  s_orderB=s_orderB[,grep("o__",colnames(s_orderB))]
  # grab classes and drop 0s
  s_classB=s_classB[,grep("c__",colnames(s_classB))]
  # grab phyla and drop 0s
  s_phylaB=s_phylaB[,grep("p__",colnames(s_phylaB))]  
  
  venn(list(inFileA=s_speciesAFN,inFileB=s_speciesBFN),names = c(inFileA,inFileB))
}

# ================================================================================
# function that preps metaphlan dataframe
# ================================================================================
# - grab MetaPhlan DF, returns dataframe with
# numbers of taxonomical units and shannon diversity
metaPhlanDFgetStats <- function(taxa_df) {
  # remove strains
  s_species=taxa_df[,grep("t__",colnames(taxa_df),invert = TRUE)]
  # remove genera
  s_genera = s_species[,grep('s__',colnames(s_species),invert = TRUE)]
  # remove families
  s_family = s_genera[,grep('g__',colnames(s_genera),invert = TRUE)]
  # remove orders
  s_order = s_family[,grep('f__',colnames(s_family),invert = TRUE)]
  # remove classes
  s_class = s_order[,grep('o__',colnames(s_order),invert = TRUE)]
  #remove phyla
  s_phyla = s_class[,grep('c__',colnames(s_class),invert = TRUE)]
  
  # grab species
  s_species=s_species[,grep("s__",colnames(s_species))]
  # grab genera
  s_genera=s_genera[,grep("g__",colnames(s_genera))]
  # grab families
  s_family=s_family[,grep("f__",colnames(s_family))]
  # grab orders
  s_order=s_order[,grep("o__",colnames(s_order))]
  # grab classes
  s_class=s_class[,grep("c__",colnames(s_class))]
  # grab phyla
  s_phyla=s_phyla[,grep("p__",colnames(s_phyla))]
  # count taxonomy stuff per sample
  # ---------------------------------------------------
  s_taxnr <- data.frame(cbind(row.names(s_species), apply(X = s_species,MARGIN = 1,FUN = nsum), 
                              apply(X = s_genera,MARGIN = 1,FUN = nsum),
                              apply(X = s_family,MARGIN = 1,FUN = nsum),
                              apply(X = s_order,MARGIN = 1,FUN = nsum),
                              apply(X = s_class,MARGIN = 1,FUN = nsum),
                              apply(X = s_phyla,MARGIN = 1,FUN = nsum)))
  
  # get rid of rows with no species detected
  s_taxnr <- s_taxnr[apply(s_species,MARGIN = 1,sum) != 0,]
  s_species <- s_species[apply(s_species,MARGIN = 1,sum) != 0,]
  
  # add diversity measures
  s_taxnr <- cbind(s_taxnr,diversity(s_species,index="shannon"))
  colnames(s_taxnr) <- c("sample","NrSpecies","NrGenera","NrFamilies","NrOrders","NrClasses","NrPhyla","Shannon")
  s_taxnr$NrSpecies <- as.numeric(as.character(s_taxnr$NrSpecies))
  s_taxnr$NrGenera <- as.numeric(as.character(s_taxnr$NrGenera))
  s_taxnr$NrFamilies <- as.numeric(as.character(s_taxnr$NrFamilies))
  s_taxnr$NrClasses <- as.numeric(as.character(s_taxnr$NrClasses))
  s_taxnr$NrOrders <- as.numeric(as.character(s_taxnr$NrOrders))
  s_taxnr$NrPhyla <- as.numeric(as.character(s_taxnr$NrPhyla))
  s_taxnr$sample <- rownames(s_taxnr)
  
  s_taxnr
}


# ================================================================================
# metaphlan clean load
# -> loads metaphlan file, removes weird stuff (like random # and blank rows),
# removes useless stuff from sample names and converts fields to numeric
# if purge = 1, purges taxa, leaving only species!
# if purge = 2, also leaves strains
# ================================================================================
metaPhlanCleanLoad <- function(inFile,outFile,purge=0) {
  mfRaw <- read.csv(inFile, sep = "\t",check.names = F, blank.lines.skip = T, row.names = 1, comment.char = '')
  colnames(mfRaw) <- gsub('_taxonomic_profile','',colnames(mfRaw))
  colnames(mfRaw) <- gsub('_metaphlan','',colnames(mfRaw))
  mfRaw <- mfRaw[!grepl('#',row.names(mfRaw)),]
  rnames <- row.names(mfRaw)
  cnames <- sapply(colnames(mfRaw), FUN=function(x) gsub('\\.','',x))
  cnames <- sapply(cnames, FUN=function(x) gsub('-','',x))
  mfRaw <- as.data.frame(lapply(mfRaw,function(x) as.numeric(as.character(x)) ))
  colnames(mfRaw) <- cnames
  row.names(mfRaw) <- rnames
  mfRaw = as.data.frame(t(mfRaw))
  
  # purge 1: remove non-bacteria, extract species only
  if (purge == 1) {
    # grab only bacteria
    mfRaw <- mfRaw[,grep("k__Bact",colnames(mfRaw),invert=F)]
    # extract species
    mfRaw <- mfRaw[,grep("s__",colnames(mfRaw),invert=F)]
    # get rid of strains
    mfRaw <- mfRaw[,grep("t__",colnames(mfRaw),invert=T)]
  }
  # purge 1: remove non-bacteria, extract specie & strains
  else if (purge == 2) {
    # grab only bacteria
    mfRaw <- mfRaw[,grep("k__Bact",colnames(mfRaw),invert=F)]
    # extract species
    mfRaw <- mfRaw[,grep("s__",colnames(mfRaw),invert=F)]
  }
  
  mfRaw
}

# =================================================
# -> HELPER FUNCTIONS
# =================================================
nzmean <- function(a){
  mean(a[a!=0])
}

##Function to calculate nº of 0
zsum <- function(a){
  sum (a==0)
}

##Function to calculate nº of non-0
nsum <- function(a){
  sum (a!=0)
}

# ====================================================================================================
# function that tests abundancy of pathways or taxa in one dataset with two-four classes of diagnosis, 
# cohorts or equivalent 
# - does t-test
# - calculates p nominal, p adjusted (holm), FDR, and 
# - plots and outputs to saveFolder
# => stTest should be "wilcox" or "permutation"
# ====================================================================================================
testAbundances <- function(dataIn,saveFolder,pathway=T,doSave=T,display="P",onlyShowSig=T, na.rem = T,alpha = 0.25,
                           responseVar="Diagnosis",stTest="wilcox",doPlots = T,nrTests=-1,xLab=NA,sX=9,sY=6)
{
  if (is.na(xLab)) {
    xLab = responseVar
  }
  if (pathway) {
    pws <- grep('PWY',colnames(dataIn))
    ftr="pathway"
  } else {
    pws <- grep('__',colnames(dataIn))
    ftr="microbial_taxa"
  }
  
  ret = F
  tmpData <- dataIn
  for (i in pws) {
    dataIn <- tmpData
    if (na.rem) {
      dataIn <- dataIn[!is.na(dataIn[[i]]),]
      dataIn <- dataIn[!is.na(responseVar),]
    }
    row1sig = F
    row2sig = F
    row3sig = F
    row4sig = F
    rowCnt = 0
    sigP = ""
    sigN = 0
    pw <- dataIn[,i]
    classes <- unique(dataIn[[responseVar]])
    # extract classes
    pwC <- list()
    for (j in classes) {
      pwC[[j]] <- dataIn[dataIn[[responseVar]] == j,i]
    }
    # do pairwise combination tests
    #      calculate combinations
    combos <- as.data.frame(t(combn(as.character(classes),2)))
    combos$V1 <- as.factor(combos$V1)
    combos$V2 <- as.factor(combos$V2)
    nCombos <- nrow(combos)
    #      do tests
    tTests <- list()
    tTestsP <- c()
    tTestsFDR <- c()
    tTestsPCorr <- c()
    if (nrTests == -1) {
      nTests <- nrow(combos)*length(pws)
    } else {
      nTests = nrTests
    }
    for (r in seq(1,nrow(combos))) {
      #print(paste("testing",combos$V1[r],"VS",combos$V2[r]))
      if (stTest != "permutation") {
        test <- wilcox.test(pwC[[as.character(combos$V1[r]) ]],pwC[[as.character(combos$V2[r]) ]])
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,test$p.value)
      } else {
        testA <- pwC[[as.character(combos$V1[r]) ]]
        testB <- pwC[[as.character(combos$V2[r]) ]]
        test <- independence_test(value ~ group, rbind.data.frame(data.frame(value=testA,group="A"),data.frame(value=testB,group="B")) )
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,pvalue(test))
      }
      #print(paste(" -> p-value:"))
    }
    # add more data to combos
    combos$pValue <- tTestsP
    combos$pValue[is.na(combos$pValue)] <- 1.0
    combos$FDR <- p.adjust(tTestsP,method="fdr",n=nTests)
    combos$FDR[is.na(combos$FDR)] <- 1.0
    combos$pCorr <-  p.adjust(tTestsP,method="holm",n=nTests)
    combos$pCorr[is.na(combos$pCorr)] <- 1.0
    
    # set positioning for plot
    dodge = position_dodge(width=0.8)
    jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
    # NR rows
    nRows <- 1
    if (nCombos == 3) {nRows <- 2}
    if (nCombos == 6) {nRows <- 4}
    if (nCombos > 6) {nRows <- -1}
    classes <- sort(classes)
    
    if (doPlots) {
      g <- ggplot(data=dataIn,aes(x=Diagnosis, y=dataIn[[i]],col=Diagnosis)) +
        #scale_color_manual(values=c("darkgreen","red")) +
        geom_jitter(alpha=alpha,position=jdodge) +
        geom_violin(alpha=alpha,width=0.9,position=dodge,size=0.75) +
        geom_boxplot(alpha=alpha,width=0.3,position=dodge,size=0.75,outlier.colour = NA) +
        ggtitle(colnames(dataIn)[i]) + ylab("Normalised relative abundancy") + xlab(xLab)
      # add annotations
      lineSize = 0.75
      # first row:
      if (nRows > 0) {
        for (f in seq(1,length(classes)-1)) {
          v1 <- as.character(classes[f])
          v2 <- as.character(classes[f+1])
          
          if (display == "FDR") {
            fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
            lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")
          } else if (display == "P") {
            fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
            lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")
          } else if (display == "Pn") {
            fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
            lbl <- paste("Pn = ",formatC(fdr, format = "e", digits = 1),sep="")
          }
          doShow = F
          if (!onlyShowSig) {doShow = T; rowCnt <- 1;ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
          if (fdr < 0.05) {rowCnt <- 1; ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
          if (fdr < 0.05) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
          if (fdr < 0.0005) {lbl <- paste(lbl,"*",sep="");sigP = "s"; sigN = sigN + 1; doShow = T}
          if (doShow) {
            g <- g + annotate("text", x=f+0.5,y=ypos, label=lbl,vjust=0 )
            g <- g + annotate("segment", x = f+0.01, y = yposLine, xend = f+1-0.01,lineend = "round",
                              yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
              annotate("segment", x = f+1-0.01, y = yposLine, xend = f+0.01,lineend = "round", 
                       yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
          }
        }
      }
      # second row: 1 vs 3
      if (nRows > 1) {
        v1 <- as.character(classes[1])
        v2 <- as.character(classes[3])
        #ypos <- max(dataIn[[i]])*(1+0.25*2-0.15)
        #yposLine <- max(dataIn[[i]])*(1+0.25*2-0.2)
        if (display == "FDR") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
          lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "P") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
          lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "Pn") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
          lbl <- paste("Pn = ",formatC(fdr, format = "e", digits = 1),sep="")
        }
        doShow = F
        if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
        if (fdr < 0.05) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
        if (fdr < 0.05) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (fdr < 0.0005) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (doShow) {
          g <- g + annotate("text", x=1.5+0.5,y=ypos, label=lbl,vjust=0 )
          g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 3-0.01,lineend = "round",
                            yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
            annotate("segment", x = 3-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                     yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
        }
      }
      # third row: 2 vs 4
      if (nRows > 2) {
        v1 <- as.character(classes[2])
        v2 <- as.character(classes[4])
        #ypos <- max(dataIn[[i]])*(1+0.25*3-0.15)
        #yposLine <- max(dataIn[[i]])*(1+0.25*3-0.2)
        if (display == "FDR") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
          lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "P") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
          lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "Pn") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
          lbl <- paste("Pn = ",formatC(fdr, format = "e", digits = 1),sep="")
        }
        doShow = F
        if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}  
        if (fdr < 0.05) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
        if (fdr < 0.05) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (fdr < 0.0005) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (doShow) {
          g <- g + annotate("text", x=2.5+0.5,y=ypos, label=lbl,vjust=0 )
          g <- g + annotate("segment", x = 2+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                            yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
            annotate("segment", x = 4-0.01, y = yposLine, xend = 2.01,lineend = "round", 
                     yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
        }
      }
      # fourth row: 1 vs 4
      if (nRows > 3) {
        v1 <- as.character(classes[1])
        v2 <- as.character(classes[4])
        #ypos <- max(dataIn[[i]])*(1+0.25*4-0.15)
        #yposLine <- max(dataIn[[i]])*(1+0.25*4-0.2)
        if (display == "FDR") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
          lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "P") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
          lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")
        } else if (display == "Pn") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
          lbl <- paste("Pn = ",formatC(fdr, format = "e", digits = 1),sep="")
        }
        doShow = F
        if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
        if (fdr < 0.05) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.15); yposLine <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt-0.2)}
        if (fdr < 0.05) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (fdr < 0.0005) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (doShow) {
          g <- g + annotate("text", x=2.0,y=ypos, label=lbl,vjust=0 )
          g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                            yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
            annotate("segment", x = 4-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                     yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
        }
        minY <- min(dataIn[[i]],na.rm = T)
        maxY <- max(dataIn[[i]],na.rm = T)   
        maxYadj <- max(dataIn[[i]])
        if (nRows != -1) {maxYadj <- max(dataIn[[i]],na.rm = T)*(1+0.25*rowCnt)  }
        g <- g + ylim(minY,maxYadj)
      }
      if (doSave) {
        ggsave(g,filename = paste(saveFolder,'_',sigP,sigN,'_',i,'.png',sep=''),width = sX,height = sY)    
      }
    }
    
    
    if (class(ret) != "data.frame") {
      ret <- combos
      ret$feature <- ftr
      ret$name <- colnames(dataIn)[i]
      ret$mcorr <- nTests
    } else {
      rett <- combos
      rett$feature <- ftr
      rett$name <- colnames(dataIn)[i]
      rett$mcorr <- nTests
      ret <- rbind(ret,rett)
    }
  }
  ret <- ret[order(ret$FDR),]
  ret
}


# ====================================================================================
# PREVALENCE / PROPORTION TESTING FUNCTION FOR ONE FEATURE
#
# - does chi-squared test (or fisher's exact test)
# - calculates p nominal, p adjusted (holm) and FDR (with input number of tests)
# - plots and outputs to saveFolder
# ====================================================================================
testOneFeaturePrevalence <- function(dataIn,saveFolder,feature=NA,doSave=T,display="P",onlyShowSig=T,yLab=NA,title=NA,na.rem=T,
                             responseVar="Diagnosis",stTest="chisquared",doPlots = T,nrTests=-1,xLab=NA,retPlot=F,ylim=NA,addNRs=F,
                             cutoff=0.05)
{
  ret = NULL
  # check input
  if (is.na(feature)) {
    print ('needs a feature as input!')
    break
  }
  if (is.na(xLab)) {
    xLab = responseVar
  }
  if (is.na(yLab)) {
    yLab = "Prevalence"
  }
  if (feature %in% colnames(dataIn)) {
    i <- grep(paste0("^",feature,"$"),colnames(dataIn))[1]
    if (class(dataIn[[i]]) != "factor") {
      print ('feature MUST BE factor with 2 classes')
      return()
    } else if (length(unique(dataIn[[i]])) != 2 ) { 
      print(paste0('ERROR: ',feature,' has ',length(unique(dataIn[[i]])),' classes!'))
      return()
      }
  } else {
    print (paste("ERROR:",feature,"does not exist in input dataset!"))
    return()
  }
  
  # do tests
  ftr=feature
  ymax = 1000000
  if (!is.na(ylim)) {ymax = max(ylim)}
  ret = F
  row1sig = F
  row2sig = F
  row3sig = F
  row4sig = F
  rowCnt = 0
  sigP = ""
  sigN = 0
  pw <- dataIn[,i]
  if (na.rem) {
    dataIn <- dataIn[!is.na(dataIn[[responseVar]]),]
  }
  if (class(dataIn[[responseVar]]) == "factor" ) {
    classes <- levels(dataIn[[responseVar]])
  } else {
    classes <- sort(na.exclude( unique(dataIn[[responseVar]]) ))
  }
  # extract classes
  pwC <- list()
  for (j in classes) {
    pwC[[j]] <- dataIn[dataIn[[responseVar]] == j,i]
  }
  # do pairwise combination tests
  #      calculate combinations
  combos <- as.data.frame(t(combn(as.character(classes),2)))
  combos$V1 <- as.character(combos$V1)
  combos$V2 <- as.character(combos$V2)
  nCombos <- nrow(combos)
  #      do tests
  tTests <- list()
  tTestsP <- c()
  tTestsFDR <- c()
  tTestsPCorr <- c()
  if (nrTests == -1) {
    nTests <- nrow(combos)
  } else {
    nTests = nrTests
  }
  V1Npos = c()
  V1Nneg = c()
  V2Npos = c()
  V2Nneg = c()
  for (r in seq(1,nrow(combos))) {
    print (paste0('testing',r))
    cTableRdy <- dataIn[dataIn[[responseVar]] %in% c(combos$V1[r],combos$V2[r]), c(feature,responseVar)]
    cTableRdy[[responseVar]] <- as.factor(as.character(cTableRdy[[responseVar]]))
    cTableRdy[[feature]] <- as.factor(as.character(cTableRdy[[feature]]))
    cTable <- table(cTableRdy)
    V1Npos = c(V1Npos,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V1[r]] == "1",na.rm = T) )
    V1Nneg = c(V1Nneg,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V1[r]] == "0",na.rm = T))
    V2Npos = c(V2Npos,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V2[r]] == "1",na.rm = T))
    V2Nneg = c(V2Nneg,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V2[r]] == "0",na.rm = T))
    if (nrow(cTable) != 2) {
      print(paste0(' > WARNING: ',feature,'has != 2 levels (should always be 2); tests will be p-value = 1.0'))
      tTestsP <- c(tTestsP,1.0)
    } else {
      #/debug print(cTable)
      if (stTest == "chisquared") {
        #print('chi-squared')
        test <- chisq.test(cTable)
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,test$p.value)
      } else if (stTest == "fisher") {
        #print('fisher')
        test <- fisher.test(cTable)
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,test$p.value)
      } else {
        #print('ind test')
        test <- independence_test(cTable,distribution=approximate(nresample = 10000))
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,pvalue(test))
      }
    }
    #/debug print(test)
  }
  # add more data to combos
  combos$pValue <- tTestsP
  combos$FDR <- p.adjust(tTestsP,method="fdr",n=nTests)
  combos$pCorr <-  p.adjust(tTestsP,method="holm",n=nTests)
  combos$V1Npos <- V1Npos
  combos$V1Nneg <- V1Nneg
  combos$V1prop <- combos$V1Npos/(combos$V1Npos+combos$V1Nneg)
  combos$V2Npos <- V2Npos
  combos$V2Nneg <- V2Nneg
  combos$V2prop <- combos$V2Npos/(combos$V2Npos+combos$V2Nneg)
  
  # set positioning for plot
  dodge = position_dodge(width=0.8)
  jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
  # NR rows
  nRows <- 1
  if (nCombos == 3) {nRows <- 2}
  if (nCombos == 6) {nRows <- 4}
  if (nCombos > 6) {nRows <- -1}
  classes <- sort(classes)
  
  # make plot
  if (doPlots) {
    dfForPlot <- data.frame()
    for (u in unique(dataIn[[responseVar]])) {
      dfForPlot <- rbind.data.frame(dfForPlot,data.frame(Feature=as.character(u),
                                                         Proportion=sum(dataIn[[feature]] == "1" & dataIn[[responseVar]] == u)/sum(dataIn[[responseVar]] == u) 
                                                         ))
    }
    dfForPlot$Feature = factor(as.character(dfForPlot$Feature),levels = levels(dataIn[[responseVar]]))
    g <- ggplot(data=dfForPlot,aes_string(x="Feature", y="Proportion",col="Feature",fill="Feature")) + geom_col() + 
      ggtitle(colnames(dataIn)[i]) + ylab(yLab) + xlab(xLab)
    if (!is.na(ylim)) {
      print (ylim)
      g <- g + ylim(ylim)
    }
    # add annotations
    yLim <- max(dfForPlot$Proportion)
    lineSize = 0.66
    rowMove = 0.08
    rowLineAdj = 0.05
    rowTextAdj = 0.07
    yAdj = 0.02
    ymax = 1000000
    rowCnt <- 0
    ftDF <- combos
    if (!is.na(yLim)) {ymax = yLim}
    # first row:
    if (nRows > 0) {
      for (f in seq(1,length(classes)-1)) {
        v1 <- as.character(classes[f])
        v2 <- as.character(classes[f+1])
        if (display == "FDR") {
          fdr <- as.numeric( combos[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$FDR)
          if (fdr > 0.00095) {
            lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
        } else if (display == "P") {
          fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pCorr)
          if (fdr > 0.00095) {
            lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
        } else if (display == "Pn") {
          fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pValue)
          if (fdr > 0.00095) {
            lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
        }
        doShow = F
        if (!onlyShowSig) {doShow = T}
        if (fdr <= cutoff) {doShow = T}
        if (fdr <= cutoff) {lbl <- paste(lbl,"*");doShow = T}
        if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep="");doShow = T}
        if (doShow) {
          rowCnt <- 1; 
          ypos <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowLineAdj+yAdj); 
          yposLine <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowTextAdj+yAdj)        
          g <- g + annotate("text", x=f+0.5,y=ypos, label=lbl,vjust=0 )
          g <- g + annotate("segment", x = f+0.01, y = yposLine, xend = f+1-0.01,lineend = "round",
                            yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
            annotate("segment", x = f+1-0.01, y = yposLine, xend = f+0.01,lineend = "round", 
                     yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
        }
      }
    }
    # second row: 1 vs 3
    if (nRows > 1) {
      v1 <- as.character(classes[1])
      v2 <- as.character(classes[3])
      if (display == "FDR") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      doShow = F
      if (!onlyShowSig) {doShow = T}
      if (fdr <= cutoff) {doShow = T}
      if (fdr <= cutoff) {lbl <- paste(lbl,"*"); doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep=""); doShow = T}
      if (doShow) {
        rowCnt <- rowCnt+1; 
        ypos <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowLineAdj+yAdj); 
        yposLine <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowTextAdj+yAdj)        
        g <- g + annotate("text", x=1.5+0.5,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 3-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 3-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }
    
    # third row: 2 vs 4
    if (nRows > 2) {
      v1 <- as.character(classes[2])
      v2 <- as.character(classes[4])
      #ypos <- max(dataIn[[i]])*(1+0.25*3-0.15)
      #yposLine <- max(dataIn[[i]])*(1+0.25*3-0.2)
      if (display == "FDR") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      doShow = F
      if (!onlyShowSig) {doShow = T}
      if (fdr <= cutoff) {doShow = T}
      if (fdr <= cutoff) {lbl <- paste(lbl,"*"); doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep=""); doShow = T}
      if (doShow) {
        rowCnt <- rowCnt+1; 
        ypos <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowLineAdj+yAdj); 
        yposLine <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowTextAdj+yAdj)        
        g <- g + annotate("text", x=2.5+0.5,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 2+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 4-0.01, y = yposLine, xend = 2.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }
    # fourth row: 1 vs 4
    if (nRows > 3) {
      v1 <- as.character(classes[1])
      v2 <- as.character(classes[4])
      #ypos <- max(dataIn[[i]])*(1+0.25*4-0.15)
      #yposLine <- max(dataIn[[i]])*(1+0.25*4-0.2)
      if (display == "FDR") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( ftDF[(ftDF$V1 == v1 & ftDF$V2 == v2) | (ftDF$V1 == v2 & ftDF$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      doShow = F
      if (!onlyShowSig) {doShow = T}
      if (fdr <= cutoff) {doShow = T}
      if (fdr <= cutoff) {lbl <- paste0(lbl," *"); doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste0(lbl,"*",sep=""); doShow = T}
      if (doShow) {
        rowCnt <- rowCnt+1; 
        ypos <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowLineAdj+yAdj); 
        yposLine <- min(ymax,max(dfForPlot[[2]],na.rm = T)*100)*(1+rowMove*rowCnt-rowTextAdj+yAdj)        
        g <- g + annotate("text", x=2.0,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 4-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }
    if (addNRs) {
      namesWnrs <- data.frame()
      
      V1Npos = c(V1Npos,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V1[r]] == "1",na.rm = T) )
      V1Nneg = c(V1Nneg,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V1[r]] == "0",na.rm = T))
      V2Npos = c(V2Npos,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V2[r]] == "1",na.rm = T))
      V2Nneg = c(V2Nneg,sum(dataIn[[feature]][dataIn[[responseVar]] == combos$V2[r]] == "0",na.rm = T))
      
      lbls <- c()
      for (lvl in levels(dfForPlot$Feature) ) {
        lbls <- c(lbls,paste0(lvl,'\n','pos=',
                              sum(dataIn[[feature]][dataIn[[responseVar]] == lvl] == "1",na.rm=T),
                              '\n','neg=',sum(dataIn[[feature]][dataIn[[responseVar]] == lvl] == "0",na.rm=T) ) )
      }
      g <- g + scale_x_discrete(labels = lbls)
    }
    
    if (doSave) {
      fn <- paste(saveFolder,'.png',sep='')
      print (fn)
      ggsave(g,filename = fn)    
    }
  }
  
  if (class(ret) != "data.frame") {
    ret <- combos
    ret$feature <- ftr
    #ret$name <- colnames(dataIn)[i]
    #ret$mcorr <- nTests
  } else {
    rett <- combos
    rett$feature <- ftr
    #rett$name <- colnames(dataIn)[i]
    #rett$mcorr <- nTests
    ret <- rbind(ret,rett)
  }
  ret$responseVar <- responseVar
  ret$V1toV2 = ret$V1prop/ret$V2prop
  ret$FDR <- NULL
  ret$pCorr <- NULL
  ret$mcorr <- NULL
  #ret <- ret[order(ret$FDR),]
  if (retPlot) {
    ret <- list(ret,g)
  }
  ret 
}


# ====================================================================================================
# function that tests one feature of samples with 2 to 4 diagnosis, cohorts or equivalent 
# - does t-test
# - calculates p nominal, p adjusted (holm), FDR, and 
# - plots and outputs to saveFolder
# => stTest should be "wilcox" or "permutation"
# -> if controlClass is not NA, does only one-way testing (all classes VS control class)
# ====================================================================================================
testOneFeature <- function(dataIn,saveFolder,feature=NA,doSave=T,display="P",onlyShowSig=T,yLab=NA, doViolin=T,title=NA,na.rem=T,
                           responseVar="Diagnosis",stTest="wilcox",doPlots = T,nrTests=-1,xLab=NA,retPlot=F,ylim=NA,alpha=0.25,
                           discardZeros=F,cutoff=0.05,controlClass=NA,doAnnotation=T,xRot=0)
{
  ret = NULL
  if (is.na(feature)) {
    stop('needs a feature as input!')
  }
  if (is.na(xLab)) {
    xLab = responseVar
  }
  if (is.na(yLab)) {
    yLab = "Normalised relative abundancy"
  }
  if (feature %in% colnames(dataIn)) {
    i <- grep(paste0("^",feature,"$"),colnames(dataIn))[1]
  } else {
    stop("ERROR:",feature,"does not exist in input dataset!")
  }
  ftr=feature
  ymax = 1000000
  if (!is.na(ylim)) {ymax = max(ylim)}
  ret = F
  row1sig = F
  row2sig = F
  row3sig = F
  row4sig = F
  rowCnt = 0
  sigP = ""
  sigN = 0
  pw <- dataIn[,i]
  if (na.rem) {
    dataIn <- dataIn[!is.na(dataIn[[responseVar]]),]
  }
  if (discardZeros) {
    dataIn <- dataIn[!dataIn[[i]] == 0,]
  }
  if (class(dataIn[[responseVar]]) == "factor" ) {
    classes <- levels(dataIn[[responseVar]])
  } else {
    classes <- sort(na.exclude( unique(dataIn[[responseVar]]) ))
  }
  # check if any response-var factors are all-0, if so, break
  if (length(classes) <= 1) {
    stop(paste0('ERROR: ',responseVar,' has < 2 classes for ',feature,'! quitting!'))
  }
  # extract classes
  pwC <- list()
  for (j in classes) {
    pwC[[j]] <- dataIn[dataIn[[responseVar]] == j,i]
  }
  
  # do pairwise combination tests
  #      calculate combinations
  if (is.na(controlClass)) {
    combos <- as.data.frame(t(combn(as.character(classes),2)))
    combos$V1 <- as.character(combos$V1)
    combos$V2 <- as.character(combos$V2)
  } else {
    combos <- cbind.data.frame(as.character(classes),rep(controlClass,length(classes) ))
    colnames(combos) <- c("V1","V2")
    combos <- combos[combos$V1 != controlClass,]
  }
  nCombos <- nrow(combos)
  #      do tests
  tTests <- list()
  tTestsP <- c()
  tTestsFDR <- c()
  tTestsPCorr <- c()
  if (nrTests == -1) {
    nTests <- nrow(combos)
  } else {
    nTests = nrTests
  }
  
  V1zeros = c()
  V1nonZeros = c()
  V1means = c()
  V1sds = c()
  V1medians = c()
  V2zeros = c()
  V2nonZeros = c()
  V2means = c()
  V2sds = c()
  V2medians = c()
  
  for (r in seq(1,nrow(combos))) {
    #print(paste("testing",combos$V1[r],"VS",combos$V2[r]))
    if (length(pwC[[as.character(combos$V1[r]) ]]) < 3 | length(pwC[[as.character(combos$V2[r]) ]]) < 3) {
      print (paste0(' >> WARNING: not enough data-point to do test for ',feature,' returning pV = 1'))
      tTestsP <- c(tTestsP,1.0)
    }
    else {
      if (stTest != "permutation") {
        test <- wilcox.test(pwC[[as.character(combos$V1[r]) ]],pwC[[as.character(combos$V2[r]) ]])
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,test$p.value)
      } else {
        testA <- pwC[[as.character(combos$V1[r]) ]]
        testB <- pwC[[as.character(combos$V2[r]) ]]
        test <- independence_test(value ~ group, 
                                  rbind.data.frame(data.frame(value=testA,group="A"),
                                                   data.frame(value=testB,group="B")),
                                  distribution=approximate(nresample = 10000))
        tTests[[r]] <- test
        tTestsP <- c(tTestsP,pvalue(test))
      }
    }
    #print(paste(" -> p-value:"))
    V1zeros =    c(V1zeros, sum(pwC[[as.character(combos$V1[r]) ]] == 0,na.rm = T) )
    V1nonZeros = c(V1nonZeros,sum(pwC[[as.character(combos$V1[r]) ]] != 0,na.rm = T) )
    V1means =    c(V1means,mean(pwC[[as.character(combos$V1[r]) ]],na.rm=T))
    V1sds =      c(V1sds,sd(pwC[[as.character(combos$V1[r]) ]],na.rm = T))
    V1medians =  c(V1medians,median(pwC[[as.character(combos$V1[r]) ]],na.rm = T))
    V2zeros =    c(V2zeros, sum(pwC[[as.character(combos$V2[r]) ]] == 0,na.rm=T) )
    V2nonZeros = c(V2nonZeros,sum(pwC[[as.character(combos$V2[r]) ]] != 0,na.rm=T) )
    V2means =    c(V2means,mean(pwC[[as.character(combos$V2[r]) ]],na.rm=T) )
    V2medians =  c(V2medians,median(pwC[[as.character(combos$V2[r]) ]],na.rm=T) )
    V2sds =      c(V2sds,sd(pwC[[as.character(combos$V2[r]) ]],na.rm = T) )
  }
  # add more info to combos
  combos$pValue <- tTestsP
  combos$FDR <- p.adjust(tTestsP,method="fdr",n=nTests)
  combos$pCorr <-  p.adjust(tTestsP,method="holm",n=nTests)
  combos$V1zeros <- V1zeros
  combos$V1nonZeros <- V1nonZeros
  combos$V1mean <- V1means
  combos$V1sd <- V1sds
  combos$V1median <- V1medians
  combos$V2zeros <- V2zeros
  combos$V2nonZeros <- V2nonZeros
  combos$V2mean <- V2means
  combos$V2sd <- V2sds
  combos$V2median <- V2medians
  combos$V1minusV2 <- combos$V1mean-combos$V2mean
  
  # set positioning for plot
  dodge = position_dodge(width=0.8)
  jdodge = position_jitterdodge(jitter.width = 0.35, jitter.height = 0, dodge.width = 0.8)
  # NR rows
  nRows <- 1
  if (nCombos == 3) {nRows <- 2}
  if (nCombos == 6) {nRows <- 4}
  if (nCombos > 6) {nRows <- -1}
  #classes <- sort(classes)
  
  if (doPlots) {
    if (!doViolin) {
      g <- ggplot(data=dataIn,aes_string(x=responseVar, y=colnames(dataIn)[i],col=responseVar)) +
        #scale_color_manual(values=c("darkgreen","red")) +
        geom_jitter(alpha=alpha,position=jdodge) +
        geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
        ggtitle(colnames(dataIn)[i]) + ylab(yLab) + xlab(xLab)
      
      if (!is.na(ylim)) {
        print (ylim)
        g <- g + ylim(ylim)
      }
    }
    if (doViolin) {
      g <- ggplot(data=dataIn,aes_string(x=responseVar, y=colnames(dataIn)[i],col=responseVar)) +
        #scale_color_manual(values=c("darkgreen","red")) +
        geom_jitter(alpha=alpha,position=jdodge) +
        geom_boxplot(alpha=0.2,width=0.45,position=dodge,size=0.75,outlier.colour = NA) +
        geom_violin(alpha=0.2,width=0.9,position=dodge,size=0.75) +
        ggtitle(colnames(dataIn)[i]) + ylab(yLab) + xlab(xLab)
      if (!is.na(ylim)) {
        print (ylim)
        g <- g + ylim(ylim)
      }
    }
    if (!is.na(title)) {g <- g + ggtitle(title)}
    
    # add annotations
    lineSize = 0.66
    rowMove = 0.10
    rowTextAdj = 0.17
    
    combos$FDR[is.na(combos$FDR)] <- 1.0
    combos$pCorr[is.na(combos$pCorr)] <- 1.0
    combos$pValue[is.na(combos$pValue)] <- 1.0
    
    # add styling
    if (xRot != 0) {
      g <- g + theme(axis.text.x = element_text(angle = xRot,hjust = 1))
    }
    
    # first row:
    if (doAnnotation) {
    if (nRows > 0) {
      for (f in seq(1,length(classes)-1)) {
        v1 <- as.character(classes[f])
        v2 <- as.character(classes[f+1])
        
        if (display == "FDR") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
          if (fdr > 0.00095) {
            lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
        } else if (display == "P") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
          if (fdr > 0.00095) {
            lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
        } else if (display == "Pn") {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
          if (fdr > 0.00095) {
            lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
          } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
        }
        if (display=='') {
          fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
          lbl <- ''
        }
        doShow = F
        if (!onlyShowSig) {doShow = T; rowCnt <- 1;ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
        yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
        if (fdr <= cutoff) {rowCnt <- 1; ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
        yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
        if (fdr <= cutoff) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
        if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep="");sigP = "s"; sigN = sigN + 1; doShow = T}
        if (doShow) {
          g <- g + annotate("text", x=f+0.5,y=ypos, label=lbl,vjust=0 )
          g <- g + annotate("segment", x = f+0.01, y = yposLine, xend = f+1-0.01,lineend = "round",
                            yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
            annotate("segment", x = f+1-0.01, y = yposLine, xend = f+0.01,lineend = "round", 
                     yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
        }
      }
    }
    # second row: 1 vs 3
    if (nRows > 1) {
      v1 <- as.character(classes[1])
      v2 <- as.character(classes[3])
      #ypos <- max(dataIn[[i]])*(1+0.25*2-0.15)
      #yposLine <- max(dataIn[[i]])*(1+0.25*2-0.2)
      if (display == "FDR") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      if (display=='') {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        lbl <- ''
      }
      doShow = F
      if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
      if (fdr <= cutoff) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
      if (fdr <= cutoff) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (doShow) {
        g <- g + annotate("text", x=1.5+0.5,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 3-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 3-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }
    # third row: 2 vs 4
    if (nRows > 2) {
      v1 <- as.character(classes[2])
      v2 <- as.character(classes[4])
      #ypos <- max(dataIn[[i]])*(1+0.25*3-0.15)
      #yposLine <- max(dataIn[[i]])*(1+0.25*3-0.2)
      if (display == "FDR") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      if (display=='') {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        lbl <- ''
      }
      doShow = F
      if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <-min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}  
      if (fdr <= cutoff) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
      if (fdr <= cutoff) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (doShow) {
        g <- g + annotate("text", x=2.5+0.5,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 2+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 4-0.01, y = yposLine, xend = 2.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }

    # fourth row: 1 vs 4
    if (nRows > 3) {
      v1 <- as.character(classes[1])
      v2 <- as.character(classes[4])
      #ypos <- max(dataIn[[i]])*(1+0.25*4-0.15)
      #yposLine <- max(dataIn[[i]])*(1+0.25*4-0.2)
      if (display == "FDR") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        if (fdr > 0.00095) {
          lbl <- paste("FDR = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("FDR = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "P") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pCorr)
        if (fdr > 0.00095) {
          lbl <- paste("Pc = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("Pc = ",formatC(fdr, format = "e", digits = 1),sep="")}
      } else if (display == "Pn") {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$pValue)
        if (fdr > 0.00095) {
          lbl <- paste("P = ",formatC(fdr, digits = 1),sep="")
        } else {lbl <- paste("P = ",formatC(fdr, format = "e", digits = 1),sep="")}
      }
      if (display=='') {
        fdr <- as.numeric( combos[(combos$V1 == v1 & combos$V2 == v2) | (combos$V1 == v2 & combos$V2 == v1),]$FDR)
        lbl <- ''
      }
      doShow = F
      
      if (!onlyShowSig) {doShow = T; rowCnt <- rowCnt + 1;ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
      if (fdr <= cutoff) {if (onlyShowSig) {rowCnt <- rowCnt + 1}; ypos <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-0.15); 
      yposLine <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*rowCnt-rowTextAdj)}
      if (fdr <= cutoff) {lbl <- paste(lbl,"*"); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (fdr <= cutoff/100) {lbl <- paste(lbl,"*",sep=""); sigP = "s"; sigN = sigN + 1; doShow = T}
      if (doShow) {
        g <- g + annotate("text", x=2.0,y=ypos, label=lbl,vjust=0 )
        g <- g + annotate("segment", x = 1+0.01, y = yposLine, xend = 4-0.01,lineend = "round",
                          yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc"))) + 
          annotate("segment", x = 4-0.01, y = yposLine, xend = 1.01,lineend = "round", 
                   yend = yposLine,size=lineSize,colour="black",arrow = arrow(length = unit(0.02, "npc")))
      }
    }
    }
    # final adjustment for axis
    minY <- min(dataIn[[i]],na.rm = T)
    maxY <- min(ymax,max(dataIn[[i]],na.rm = T))  
    maxYadj <- max(dataIn[[i]],na.rm = T)
    if (nRows != -1 & is.na(ylim)) {
      maxYadj <- max(dataIn[[i]],na.rm = T)*(1+rowMove*(rowCnt-1)-0.00)  
      g <- g + ylim(minY,maxYadj)}
    if (nRows != -1 & !is.na(ylim)) {
      maxYadj <- min(ymax,max(dataIn[[i]],na.rm = T))*(1+rowMove*(rowCnt-1)-0.10); #max(dataIn[[i]],na.rm = T)*(1+rowMove*(rowCnt-0.5))  
      g <- g + ylim(ylim[1],maxYadj*1.1)
    }
    # save plot
    if (doSave) {
      fn <- paste(saveFolder,'.png',sep='')
      print (fn)
      ggsave(g,filename = fn)    
    }
  }
  if (class(ret) != "data.frame") {
    ret <- combos
    ret$feature <- ftr
    #ret$name <- colnames(dataIn)[i]
    #ret$mcorr <- nTests
  } else {
    rett <- combos
    rett$feature <- ftr
    #rett$name <- colnames(dataIn)[i]
    #rett$mcorr <- nTests
    ret <- rbind(ret,rett)
  }
  #ret <- ret[order(ret$FDR),]
  if (retPlot) {
    ret <- list(ret,g)
  }
  ret$FDR <- NULL
  ret$pCorr <- NULL
  ret$mcorr <- NULL
  ret$name <- NULL
  ret
}

# ==================================================================================
# goes through abundancy tests [generated by testAbundances], picks "best" candidates
# calculates score for each feature as follows:
# -> for each test in feature, score = 1 if below sigCut[1], 2 if below sigCut[2]
# -> totScore = sum(testScores); norScore = totScore/nr tests
# if norScore >= minScore, picks the feature
# finally sorts by score and picks best maxFeatures features (based on score)
# ==================================================================================
pickTopFeatures <- function(tests,sigMetric="FDR",sigCut=c(0.05,0.005,0.0005,0.00005,0.000005),sigValue=c(1,1.25,1.5,1.75,2),minScore=1.0,maxScore=1000,maxFeatures=1000) {
  cModel <- data.frame(name=as.character(),norScore=as.numeric(),stringsAsFactors = F)
  un <- unique(tests$name)
  for (u in un) {
    uTest <- tests[tests$name == u,]
    totScore = 0
    for (t in seq(1,nrow(uTest))) {
      for (r in seq(1,length(sigCut))) {
        if (uTest[t,][[sigMetric]] <= rev(sigCut)[r]) {
          totScore <- totScore + rev(sigValue[r])
          break
        }
      }
    }
    norScore <- totScore/nrow(uTest)
    if (norScore >= minScore & norScore <= maxScore) {
      cModel <- rbind.data.frame(cModel,data.frame(name=as.character(u),norScore=as.numeric(norScore),stringsAsFactors = F),stringsAsFactors = F)
    }
  }
  cModel <- cModel[order(cModel$norScore,decreasing = T),]
  if (nrow(cModel) > maxFeatures) {cModel <- cModel[1:maxFeatures,]}
  cModel
}

# ==================================================================================
# keeps "toKeep" pathways & taxa, also keeps phenotype data
# ==================================================================================
filterMGPWYmodel <- function(inData,toKeep) {
  for (c in colnames(inData)) {
    if (c %in% toKeep | (length(grep("PWY", c)) == 0 & length(grep("__", c)) == 0) ) {
    } else {
      inData[[c]] <- NULL
    }
  }
  inData
}