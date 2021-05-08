**Packages & Functions**
```{r}

setwd('/Users/ibd/Desktop/MD:PhD/Chapter5-virome/CrAssPhage_project/Scripts/')

## Functions
library(ggplot2)
library(reshape2)
library(forcats)
library(dplyr)
library(RColorBrewer)
library(ggridges) # install.packages("ggridges")
library(VennDiagram) # install.packages("VennDiagram") 



# 1. Prevalence by cohort ......

PrevalenceByCohort <- function(x, y = 3, byVariable){
  my_matrix<-  matrix(nrow = length(y:ncol(x)), ncol = 11)
  colnames(my_matrix) <- c('CrAss_Phage', 'Prev_All', 'Prev_LLD', 'Prev_LLD2','Prev_1000IBD','Prev_300OB','n_All', 'n_LLD','n_LLD2','n_1000IBD','n_300OB')
  a=1
  groups <- levels(x[,byVariable])
  for (k in y:ncol(x)) {
    my_matrix[a,1] <- colnames(x)[k] 
    my_matrix[a,2] <- colSums(x[k] != 0, na.rm = T) / nrow(x) *100
    my_matrix[a,3] <- colSums(subset(x, x[c(byVariable)] == groups[1], select = k) != 0, na.rm = T) /sum(x[byVariable] == groups[1]) *100
    my_matrix[a,4] <- colSums(subset(x, x[c(byVariable)] == groups[2], select = k) != 0, na.rm = T) /sum(x[byVariable] == groups[2]) *100
    my_matrix[a,5] <- colSums(subset(x, x[c(byVariable)] == groups[3], select = k) != 0, na.rm = T) /sum(x[byVariable] == groups[3]) *100
    my_matrix[a,6] <- colSums(subset(x, x[c(byVariable)] == groups[4], select = k) != 0, na.rm = T) /sum(x[byVariable] == groups[4]) *100
    my_matrix[a,7] <- nrow(x)
    my_matrix[a,8] <- sum(x[byVariable] == groups[1])
    my_matrix[a,9] <- sum(x[byVariable] == groups[2])
    my_matrix[a,10] <- sum(x[byVariable] == groups[3])
    my_matrix[a,11] <-sum(x[byVariable] == groups[4])
    a=a+1
  }
  return(as.data.frame(my_matrix))
}


# 2. Prevalence by diagnosis ......

PrevalenceByDiagnosis <- function(x, y = 3, level){
  my_matrix<-  matrix(nrow = length(y:ncol(x)), ncol = 13)
  colnames(my_matrix) <- c('CrAss_Phage', 'Prev_All', 'Prev_LLD', 'Prev_LLD2','Prev_CD','Prev_UC','Prev_300OB','n_All', 'n_LLD','n_LLD2','n_CD','n_CD','n_300OB')
  a=1
  groups <- levels(x[,level])
  for (k in y:ncol(x)) {
    my_matrix[a,1] <- colnames(x)[k] 
    my_matrix[a,2] <- colSums(x[k] != 0, na.rm = T) / nrow(x) *100
    my_matrix[a,3] <- colSums(subset(x, x[c(level)] == groups[1], select = k) != 0, na.rm = T) /sum(x[level] == groups[1]) *100
    my_matrix[a,4] <- colSums(subset(x, x[c(level)] == groups[2], select = k) != 0, na.rm = T) /sum(x[level] == groups[2]) *100
    my_matrix[a,5] <- colSums(subset(x, x[c(level)] == groups[3], select = k) != 0, na.rm = T) /sum(x[level] == groups[3]) *100
    my_matrix[a,6] <- colSums(subset(x, x[c(level)] == groups[4], select = k) != 0, na.rm = T) /sum(x[level] == groups[4]) *100
    my_matrix[a,7] <- colSums(subset(x, x[c(level)] == groups[5], select = k) != 0, na.rm = T) /sum(x[level] == groups[5]) *100
    my_matrix[a,8] <- nrow(x)
    my_matrix[a,9] <- sum(x[level] == groups[1])
    my_matrix[a,10] <- sum(x[level] == groups[2])
    my_matrix[a,11] <- sum(x[level] == groups[3])
    my_matrix[a,12] <-sum(x[level] == groups[4])
    my_matrix[a,13] <-sum(x[level] == groups[5])
    a=a+1
  }
  return(as.data.frame(my_matrix))
}


# 3. Prevelance Filter....

Filter_func <- function(df, start, thres = 5){
  remove_cols <- vector()
  for (i in start:ncol(df)) {
    cname <- colnames(df)[i]
    if(colSums(df[i] != 0, na.rm = T) / nrow(df) *100 < thres){
      remove_cols <- c(remove_cols, cname) 
    }
  }
  df_1 <- df %>% select(-remove_cols)
  print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
  return(df_1)
}

#4. Transformation & filter function...

transform_and_filter_taxa=function(x, samples_row=T, method="asin", missing_filter=0){
  x[x=="NA"]=0
  x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
    x=as.data.frame(t(x))
  }
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100")
  }
  x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  if (method=="asin"){
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))
  } else if (method=="log"){
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0]/2)
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)
  }else if (method=="clr"){
    x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0]/2)
    x_filt=x_filt+my_min
    #clr transformation (adapted from microbiome package function)
    d <- apply(x_filt, 2, function(b) {
      log(b) - mean(log(b))
    })
    x_filt=d
  }
  return(as.data.frame(x_filt))
}


#5. Pairwise wilcoxon test function...

wilcox_pairwise_function <- function(df, start = 2, pheno_col = 1, ncombos){
  
  my_matrix<-  matrix(nrow = ncombos * ncol(df[c(start:ncol(df))]), ncol = 14)
  colnames(my_matrix) <- c('CrAss_Phage', 'Variable', 'Cohort1', 'Cohort2', 'n_Cohorts', 'n_SamplesCohort1', 'n_SamplesCohort2', 'P-value', 'BonferroniAdjust', 'FDRadjust',  'MeanAbundance_Cohort1', 'MeanAbundance_Cohort2',  'MedianAbundance_Cohort1', 'MedianAbundance_Cohort2')
  
  #__Add correlation analysis data to the matrix 
  a=1 
  #Iterate through taxonomy columns
  for (j in start:ncol(df)) {
    #iterate through phenotype columns
    for (k in pheno_col) {
      #following conditions applied to all factor variables
      combos <- combn(levels(df[,k]),2)
      #pairwise.wilcox.test 
      wilcox <- (pairwise.wilcox.test(df[,j], df[,k], p.adjust.method = 'none'))$p.value
      #condition applied to all factor variables with more than 2 groups
      for (z in 1:ncol(combos)) {
        #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
        my_matrix[a,1] <- colnames(df)[j]
        #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
        my_matrix[a,2] <- colnames(df)[k]
        #adding the group combinations to my_matrix
        my_matrix[a,3] <- combos[1,z]
        my_matrix[a,4] <- combos[2,z]
        my_matrix[a,5] <- nlevels(df[,k])
        my_matrix[a,6] <- sum(df[,k]==combos[1,z], na.rm = T) #Number of samples in group1
        my_matrix[a,7] <- sum(df[,k]==combos[2,z], na.rm = T) #Number of samples in group2
        my_matrix[a,8] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
        my_matrix[a,11] <- mean(subset(df, df[c(k)] == combos[1,z], select = j)[,1],na.rm = T)
        my_matrix[a,12] <- mean(subset(df, df[c(k)] == combos[2,z], select = j)[,1],na.rm = T)
        my_matrix[a,13] <- median(subset(df, df[c(k)] == combos[1,z], select = j)[,1],na.rm = T)
        my_matrix[a,14] <- median(subset(df, df[c(k)] == combos[2,z], select = j)[,1],na.rm = T)
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,8] <- as.numeric(as.character(my_matrix[,8]))
  my_matrix[,9] <- p.adjust(c(my_matrix[,8]), method = 'bonferroni')
  my_matrix[,10] <- p.adjust(c(my_matrix[,8]), method = 'fdr')
  return(my_matrix)
}

# 6. testOneFeaturePrevalence (Ranko's code)
function(dataIn,saveFolder,feature=NA,doSave=T,display="P",onlyShowSig=T,yLab=NA,title=NA,na.rem=T,
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
    g <- ggplot(data=dfForPlot,aes_string(x="Feature", y="Proportion",fill="Feature")) + geom_col(colour="black") + ylab(yLab) + xlab(xLab)
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


# 7. Logistic Regression funtion (presence/absence)

LogisticRegression.function <- function(DF, nPheno, colTaxa, my_preds){
  Nlevels <- vector()
  NonFactors <- vector()
  for (i in 1:nPheno) {
    if(is.factor(DF[,i])) {
      Nlevels <- c(Nlevels, (nlevels(DF[,i])-1))
    }
    else {
      NonFactors <- c(NonFactors, colnames(DF)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Nlevels)) * ncol(DF[c(colTaxa:ncol(DF))]), ncol = 9)
    colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Estimate', 'Std.Error', 'z_value','Pr(>|z|)','BonferroniAdjust','FDRAdjust', 'Prevalence')
  }
  
  a=1
  for (j in colTaxa:ncol(DF)) {
    my_f <- as.formula(paste(colnames(DF)[j], paste(my_preds, collapse = " + "), sep = " ~ "))
    my_lm <- glm(my_f, family = binomial(link="logit"), data = DF)
    for(k in 2:nrow(summary(my_lm)$coefficients)){
      my_matrix[a,1] <- colnames(DF)[j]
      my_matrix[a,2] <- rownames(summary(my_lm)$coefficients)[k]
      my_matrix[a,3] <- summary(my_lm)$coefficients[k,1]
      my_matrix[a,4] <- summary(my_lm)$coefficients[k,2]
      my_matrix[a,5] <- summary(my_lm)$coefficients[k,3]
      my_matrix[a,6] <- summary(my_lm)$coefficients[k,4]
      my_matrix[a,9] <- sum(DF[,j])/nrow(DF) *100
      a=a+1
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,6] <- as.numeric(as.character(my_matrix[,6]))
  my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
  my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')
  return(my_matrix)
}



# 8. Descriptive statistics for logisitic regression 

LgRgr.DescriptiveStats <- function(DF, nPheno, colTaxa){
  Nlevels <- vector()
  NonFactors <- vector()
  for (i in 1:nPheno) {
    if(is.factor(DF[,i])) {
      Nlevels <- c(Nlevels, (nlevels(DF[,i])-1))
    }
    else {
      NonFactors <- c(NonFactors, colnames(DF)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Nlevels)) * ncol(DF[c(colTaxa:ncol(DF))]), ncol = 8)
    colnames(my_matrix) <- c('Taxonomy', 'Phenotype','Group', 'nGroup', 'nGroupOther', 'PrevalenceWithinGroup', 'PrevalenceWithinGroupOther', 'TotalPrevalence')
  }
  
  a=1
  for (j in colTaxa:ncol(DF)) {
    for (i in 1:nPheno) {
      if(is.factor(DF[,i])){
        Levels <- levels(DF[,i])[2:nlevels(DF[,i])]
          for (m in 1:length(Levels)) {
          my_matrix[a,1] <- colnames(DF)[j]
          my_matrix[a,2] <- colnames(DF)[i]
          my_matrix[a,3] <- Levels[m]
          my_matrix[a,4] <- sum(DF[,i] == Levels[m], na.rm = T)
          n_samples <- sum(DF[,i] == Levels[m], na.rm = T)
          my_matrix[a,5] <- sum(DF[,i] != Levels[m], na.rm = T)
          n_samplesother <- sum(DF[,i] != Levels[m], na.rm = T)
          my_matrix[a,6] <- colSums(subset(DF, DF[c(i)] == Levels[m], select = j) == 1)/n_samples *100 #bacteria prevalence within the group/level 'm'
          my_matrix[a,7] <- colSums(subset(DF, DF[c(i)] != Levels[m], select = j) == 1)/n_samplesother *100 #bacteria prevalence within the groups/levels other than m
          my_matrix[a,8] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria within whole cohort
          a=a+1
          }
      }
      else{
        my_matrix[a,1] <- colnames(DF)[j]
        my_matrix[a,2] <- colnames(DF)[i]
        my_matrix[a,3] <- colnames(DF)[i]
        my_matrix[a,4] <- nrow(DF)
        my_matrix[a,5] <- nrow(DF)
        my_matrix[a,6] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria among the non NA individuals per numerical phenotype 
        my_matrix[a,7] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria among the non NA individuals per numerical phenotype 
        my_matrix[a,8] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria within whole cohort
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  return(my_matrix)
}

```



**Data Import & process**
```{r}
options(scipen = 999)

# ---Import---
CrAss_abundances <- read.table(file = "../from_Peregrine/LLD_LLD2_300OB_IBD_crAss_abundance_table_metaphlan_style_RR.txt", header = TRUE, sep = "\t", row.names = 1)# CrAss-like phage abundance data 

Samples_cohort <- read.table(file = "../from_Peregrine/sample_cohort_RR.txt", header = TRUE, sep = "\t",row.names = 1) # The cohort of each sample 

new_ids <- read.table(file = "../from_Peregrine/rename_IBD.txt", header = TRUE, sep = "\t") # The cohort of each sample 

Metadata <- read.table(file = "../from_Peregrine/LLD_IBD_meta_201020.txt", sep = "\t", header = T,  quote = "\\")

Associations <- read.table(file = "../from_Peregrine/LLD_crass_pheno_associations_noSpecies.txt", header = TRUE, sep = "\t")
Associations_1 <- Associations[(Associations$FDR.quant<0.05),] # for covariates

BiomarkerVals <- read.csv("../from_Peregrine/LLDIBDMIBS_biomarkers_raw.csv", header = TRUE)

LLDeepIDs <- read.table(file = "../from_Peregrine/LLDeepIDs.txt", header = TRUE, sep = "\t")

OB_key <- read.table(file = "../from_Peregrine/key_300OB.txt", header = TRUE, sep = "\t")
OB_meta <- read.table(file = "../from_Peregrine/300OB_AllData_CVON study_def_correction.txt", header = TRUE, sep = "\t")
OB_MetabolSyn <- read.table(file = "../from_Peregrine/300OB_metabolicSyndrome.txt", header = TRUE, sep = "\t")

LLDfup_key <- read.table(file = "../from_Peregrine/key_LLD_baseline_fup_338sample_pairs.txt", header = TRUE, sep = "\t")
LLDfup_meta <- read.table(file = "../from_Peregrine/data_pheno_LLD_base_fup_338pairs_45pheno_log_imput_min_10participants.txt", header = TRUE, sep = "\t")

LLD_reads <- read.table(file = "../from_Peregrine/LLD_raw_clean_reads_number.txt", sep = "\t", header = T)
IBD_reads <- read.table(file = "../from_Peregrine/IBD_raw_clean_reads_number.txt", sep = "\t", header = T)
OB_reads <- read.table(file = "../from_Peregrine/300OB_raw_clean_reads_number.txt", sep = "\t", header = T)



# ---Process---
  #Phage data
colnames(CrAss_abundances) <- gsub("X","",colnames(CrAss_abundances))
CrAss_abundances <- as.data.frame(t(CrAss_abundances))
colnames(CrAss_abundances) <- gsub("\\|","_",colnames(CrAss_abundances))
CrAss_abundances <- CrAss_abundances[-c(1775,1876,2042),] # remove duplicate samples (they are the same as samples IBDFEC0745 and IBDFEC0471.2)


CrAss_abundances_new <- merge(new_ids,CrAss_abundances,by.y = "row.names", by.x = "old", all.y = T)
CrAss_abundances_new$Classic <- as.character(CrAss_abundances_new$Classic)
CrAss_abundances_new$old <- as.character(CrAss_abundances_new$old)
CrAss_abundances_new[c(102:2288),2] <- CrAss_abundances_new[c(102:2288),1]
rownames(CrAss_abundances_new) <- CrAss_abundances_new$Classic

  #metadata
Samples_cohort$IDs <- rownames(Samples_cohort)
Samples_cohort <- Samples_cohort[-c(1775,1876,2042),]  # remove duplicate samples (they are the same as samples IBDFEC0745 and IBDFEC0471.2)
Samples_cohort_new <- merge(Samples_cohort,new_ids,by.x = "IDs", by.y = "old", all.x = T)
Samples_cohort_new$Classic <- as.character(Samples_cohort_new$Classic)
Samples_cohort_new[c(1:1771,1873:2288),3] <- Samples_cohort_new[c(1:1771,1873:2288),1]
Samples_cohort_new_1 <- subset(Samples_cohort_new, cohort != "LLD2")
Samples_cohort_new_1$cohort <- factor(Samples_cohort_new_1$cohort, levels = c("LLD", "IBD", "300OB"))


Metadata$ID <- as.character(Metadata$ID)
Metadata[460,1] <- "IBDFEC0745"   #Change sample ID from IBDFEC9745 to IBDFEC0745

# Remove Stoma and Pouch samples
IBD_BowelPheno <- Metadata[c(1,137)]
IDs <- subset(IBD_BowelPheno, CurrentStomaOrPouchType == "ileostoma" |CurrentStomaOrPouchType == "colostoma" | CurrentStomaOrPouchType == "pouch", select = "ID")
CrAss_abundances_new_1 <- subset(CrAss_abundances_new, !(Classic %in% IDs$ID)) 


#Subset 
#Genus level only
CrAss_G <- CrAss_abundances_new_1[,grepl('g__',colnames(CrAss_abundances_new_1),ignore.case = T)]
CrAss_G<- CrAss_G[,!grepl('s__',colnames(CrAss_G),ignore.case = T)]
colnames(CrAss_G) <- gsub('.*g__','',colnames(CrAss_G)) #Keep taxonomy level name only 

#Family level only
CrAss_F <- CrAss_abundances_new_1[,grepl('f__',colnames(CrAss_abundances_new_1),ignore.case = T)]
CrAss_F<- CrAss_F[,!grepl('g__',colnames(CrAss_F),ignore.case = T)]
colnames(CrAss_F) <- gsub('.*f__','',colnames(CrAss_F)) #Keep taxonomy level name only 

#Order level only
CrAss_O <- CrAss_abundances_new_1[,grepl('o__',colnames(CrAss_abundances_new_1),ignore.case = T)]
CrAss_O<- CrAss_O[!grepl('f__',colnames(CrAss_O),ignore.case = T)]
colnames(CrAss_O) <- gsub('.*o__','',colnames(CrAss_O)) #Keep taxonomy level name only 


#Merge
  #Cohort stratification data

CrAss_G_PerCohort <- merge(Samples_cohort_new_1, CrAss_G, by.x = "Classic", by.y = "row.names")
rownames(CrAss_G_PerCohort) <- CrAss_G_PerCohort$Classic
CrAss_G_PerCohort$Classic <- NULL

CrAss_F_PerCohort <- merge(Samples_cohort_new_1, CrAss_F, by.x = "Classic", by.y = "row.names")
rownames(CrAss_F_PerCohort) <- CrAss_F_PerCohort$Classic
CrAss_F_PerCohort$Classic <- NULL

CrAss_O_PerCohort <- merge(Samples_cohort_new_1, CrAss_O, by.x = "Classic", by.y = "row.names")
rownames(CrAss_O_PerCohort) <- CrAss_O_PerCohort$Classic
CrAss_O_PerCohort$Classic <- NULL


  #merge with IBD diagnosis data
IBD_diagnosis <- Metadata[c(1,35)]
IBD_diagnosis$DiagnosisCurrent <- as.character(IBD_diagnosis$DiagnosisCurrent)
IBD_diagnosis$DiagnosisCurrent <- ifelse(IBD_diagnosis$DiagnosisCurrent == "IBDI" |IBD_diagnosis$DiagnosisCurrent == "IBDU"|IBD_diagnosis$DiagnosisCurrent == "MicroscopicColitis"| IBD_diagnosis$DiagnosisCurrent == "ReconsideringDiagnosis", "UC",IBD_diagnosis$DiagnosisCurrent)  #Put IBDU, IBDI etc into UC group
IBD_diagnosis$DiagnosisCurrent <- as.factor(IBD_diagnosis$DiagnosisCurrent)  

#Genus
CrAss_G_Meta <- merge(IBD_diagnosis,CrAss_G_PerCohort, by.x = "ID", by.y = "row.names", all.y = T)
rownames(CrAss_G_Meta) <- CrAss_G_Meta$ID
CrAss_G_Meta$ID <- NULL

#Add additional levels to 'DiagnosisCurrent' 
levels <- levels(CrAss_G_Meta$DiagnosisCurrent)
levels <- c("LLD", levels,"300OB")

# refactor DiagnosisCurrent to include additional levels
CrAss_G_Meta$DiagnosisCurrent <- factor(CrAss_G_Meta$DiagnosisCurrent, levels = levels)
CrAss_G_Meta[1:1433,1] <- CrAss_G_Meta[1:1433,3] # which(rownames(CrAss_G_Meta) == "G105118")
CrAss_G_Meta$DiagnosisCurrent <- fct_drop(CrAss_G_Meta$DiagnosisCurrent)


#Family
CrAss_F_Meta <- merge(IBD_diagnosis,CrAss_F_PerCohort, by.x = "ID", by.y = "row.names", all.y = T)
rownames(CrAss_F_Meta) <- CrAss_F_Meta$ID
CrAss_F_Meta$ID <- NULL

#Add additional levels to 'DiagnosisCurrent' 
levels <- c("LLD","CD","UC","300OB")

# refactor DiagnosisCurrent to include additional levels
CrAss_F_Meta$DiagnosisCurrent <- factor(CrAss_F_Meta$DiagnosisCurrent, levels = levels)
CrAss_F_Meta[1:1433,1] <- CrAss_F_Meta[1:1433,3]


#Order
CrAss_O_Meta <- merge(IBD_diagnosis,CrAss_O_PerCohort, by.x = "ID", by.y = "row.names", all.y = T)
rownames(CrAss_O_Meta) <- CrAss_O_Meta$ID
CrAss_O_Meta$ID <- NULL

#Add additional levels to 'DiagnosisCurrent' 
levels <- c("LLD","CD","UC","300OB")

# refactor DiagnosisCurrent to include additional levels
CrAss_O_Meta$DiagnosisCurrent <- factor(CrAss_O_Meta$DiagnosisCurrent, levels = levels)
CrAss_O_Meta[1:1433,1] <- CrAss_O_Meta[1:1433,3]


```


**Data Processing - Filtering**
```{r}

#Filter out taxa present in <5% LLD samples  

FilteredTaxa <- Filter_func(subset(CrAss_G_PerCohort, cohort == "LLD"), 3)
Names_FilteredTaxa <- colnames(FilteredTaxa)[3:ncol(FilteredTaxa)]


```



**Logistic Regression - presence/absence analysis**
```{r}

#=============
# ORDER LEVEL
#=============


#cofactors <- Metadata[,c(1,3,5,9,34,565,433)]
cofactors <- Metadata[,c(1,3,9,433,565)]

BiomarkerVals <- read.csv("../from_Peregrine/LLDIBDMIBS_biomarkers_raw.csv", header = TRUE)

CrAss_O_Meta$newID <- rownames(CrAss_O_Meta)
CrAss_O_Meta_1 <- merge(BiomarkerVals,CrAss_O_Meta,by.x = "ID",by.y = "IDs",all.y = T)
CrAss_O_Meta_1 <- merge(LLDeepIDs, CrAss_O_Meta_1,by.x = "LLDeepID2", by.y = "ID", all.y = T)
CrAss_O_Meta_1$LLDeepID1 <- as.character(CrAss_O_Meta_1$LLDeepID1)
CrAss_O_Meta_1$newID <- as.character(CrAss_O_Meta_1$newID)
CrAss_O_Meta_1[1136:1888,2] <- CrAss_O_Meta_1[1136:1888,11] 

CrAss_O_Meta_1 <- merge(cofactors, CrAss_O_Meta_1, by.x = "ID", by.y = "LLDeepID1", all.y = T)


# LLD vs IBD
#===========

CrAss_O_Meta_IBDLLD <- subset(CrAss_O_Meta_1, cohort != "300OB", select = c(1:5,9,12:15))
CrAss_O_Meta_IBDLLD$DiagnosisCurrent <- fct_drop(CrAss_O_Meta_IBDLLD$DiagnosisCurrent)
CrAss_O_Meta_IBDLLD$cohort <- fct_drop(CrAss_O_Meta_IBDLLD$cohort)
rownames(CrAss_O_Meta_IBDLLD) <- CrAss_O_Meta_IBDLLD$newID
CrAss_O_Meta_IBDLLD$newID <- NULL

CrAss_O_Meta_IBDLLD[c(9)] <- as.data.frame(ifelse(CrAss_O_Meta_IBDLLD[c(9)] > 0, 1, 0))

LLD.IBD_lmO <- CrAss_O_Meta_IBDLLD[c(2:6,8:ncol(CrAss_O_Meta_IBDLLD))]

my_preds=c(colnames(LLD.IBD_lmO)[1:6])

LgR <- LogisticRegression.function(LLD.IBD_lmO,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.IBD_lmO, 6, 7)

LogisticRegrLLDIBD_O <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLDIBD_O, "CrAssOrder_LogisticRegression_LLDvsIBD.txt", sep = "\t", col.names = T)



# LLD vs 300OB
#=============

OB_keymeta <- merge(OB_key[c(4:6)],OB_meta,by.x = "ID", by.y = "C_ID")
CrAss_O_Meta_300OBLLD <- merge(OB_keymeta, CrAss_O_Meta_1,by.x = "G_id", by.y = "ID", all.y = T)
CrAss_O_Meta_300OBLLD_1 <- CrAss_O_Meta_300OBLLD[c(1,2,5,6,69,70,79:82)]
CrAss_O_Meta_300OBLLD_1$sex <- ifelse(CrAss_O_Meta_300OBLLD_1$sex == 0, "female", "male") 

CrAss_O_Meta_300OBLLD_1[1:298,5] <- CrAss_O_Meta_300OBLLD_1[1:298,3] 
CrAss_O_Meta_300OBLLD_1[1:298,6] <- CrAss_O_Meta_300OBLLD_1[1:298,4] 

CrAss_O_Meta_300OBLLD_2 <- subset(CrAss_O_Meta_300OBLLD_1, cohort != "IBD", select = c(5,6,8:10))
CrAss_O_Meta_300OBLLD_2$cohort <- fct_drop(CrAss_O_Meta_300OBLLD_2$cohort)

rownames(CrAss_O_Meta_300OBLLD_2) <- CrAss_O_Meta_300OBLLD_2$newID
CrAss_O_Meta_300OBLLD_2$newID <- NULL

CrAss_O_Meta_300OBLLD_2[c(4)] <- as.data.frame(ifelse(CrAss_O_Meta_300OBLLD_2[c(4)] > 0, 1, 0))

LLD.300OB_lmO <- CrAss_O_Meta_300OBLLD_2

my_preds=c(colnames(LLD.300OB_lmO)[1:3])

LgR <- LogisticRegression.function(LLD.300OB_lmO,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.300OB_lmO, 3,4)

LogisticRegrLLD300OB_o <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLD300OB_o, "CrAssOrder_LogisticRegression_LLDvs300OB.txt", sep = "\t", col.names = T)


#=============
# FAMILY LEVEL
#=============


CrAss_F_Meta$newID <- rownames(CrAss_F_Meta)
CrAss_F_Meta_1 <- merge(BiomarkerVals,CrAss_F_Meta,by.x = "ID",by.y = "IDs",all.y = T)
CrAss_F_Meta_1 <- merge(LLDeepIDs, CrAss_F_Meta_1,by.x = "LLDeepID2", by.y = "ID", all.y = T)
CrAss_F_Meta_1$LLDeepID1 <- as.character(CrAss_F_Meta_1$LLDeepID1)
CrAss_F_Meta_1$newID <- as.character(CrAss_F_Meta_1$newID)
CrAss_F_Meta_1[1136:1888,2] <- CrAss_F_Meta_1[1136:1888,15] 

CrAss_F_Meta_1 <- merge(cofactors, CrAss_F_Meta_1, by.x = "ID", by.y = "LLDeepID1", all.y = T)


# LLD vs IBD
#===========

CrAss_F_Meta_IBDLLD <- subset(CrAss_F_Meta_1, cohort != "300OB", select = c(1:5,9,12:19))
CrAss_F_Meta_IBDLLD$DiagnosisCurrent <- fct_drop(CrAss_F_Meta_IBDLLD$DiagnosisCurrent)
CrAss_F_Meta_IBDLLD$cohort <- fct_drop(CrAss_F_Meta_IBDLLD$cohort)
rownames(CrAss_F_Meta_IBDLLD) <- CrAss_F_Meta_IBDLLD$newID
CrAss_F_Meta_IBDLLD$newID <- NULL

CrAss_F_Meta_IBDLLD[c(9:13)] <- as.data.frame(ifelse(CrAss_F_Meta_IBDLLD[c(9:13)] > 0, 1, 0))

LLD.IBD_lm <- CrAss_F_Meta_IBDLLD[c(2:6,8:ncol(CrAss_F_Meta_IBDLLD))]

my_preds=c(colnames(LLD.IBD_lm)[1:6])

LgR <- LogisticRegression.function(LLD.IBD_lm,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.IBD_lm, 6, 7)

LogisticRegrLLDIBD <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLDIBD, "CrAssFamily_LogisticRegression_LLDvsIBD.txt", sep = "\t", col.names = T)



# LLD vs 300OB
#=============

CrAss_F_Meta_300OBLLD <- merge(OB_keymeta, CrAss_F_Meta_1,by.x = "G_id", by.y = "ID", all.y = T)
CrAss_F_Meta_300OBLLD_1 <- CrAss_F_Meta_300OBLLD[c(1,2,5,6,69,70,79:86)]
CrAss_F_Meta_300OBLLD_1$sex <- ifelse(CrAss_F_Meta_300OBLLD_1$sex == 0, "female", "male") 

CrAss_F_Meta_300OBLLD_1[1:298,5] <- CrAss_F_Meta_300OBLLD_1[1:298,3] 
CrAss_F_Meta_300OBLLD_1[1:298,6] <- CrAss_F_Meta_300OBLLD_1[1:298,4] 

CrAss_F_Meta_300OBLLD_2 <- subset(CrAss_F_Meta_300OBLLD_1, cohort != "IBD", select = c(5,6,8:14))
CrAss_F_Meta_300OBLLD_2$cohort <- fct_drop(CrAss_F_Meta_300OBLLD_2$cohort)

rownames(CrAss_F_Meta_300OBLLD_2) <- CrAss_F_Meta_300OBLLD_2$newID
CrAss_F_Meta_300OBLLD_2$newID <- NULL


CrAss_F_Meta_300OBLLD_2[c(4:8)] <- as.data.frame(ifelse(CrAss_F_Meta_300OBLLD_2[c(4:8)] > 0, 1, 0))

LLD.300OB_lm <- CrAss_F_Meta_300OBLLD_2

my_preds=c(colnames(LLD.300OB_lm)[1:3])

LgR <- LogisticRegression.function(LLD.300OB_lm,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.300OB_lm, 3,4)

LogisticRegrLLD300OB <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLD300OB, "CrAssFamily_LogisticRegression_LLDvs300OB.txt", sep = "\t", col.names = T)




#============
# GENUS LEVEL
#============


CrAss_G_Meta$newID <- rownames(CrAss_G_Meta)
CrAss_G_Meta_1 <- merge(BiomarkerVals,CrAss_G_Meta,by.x = "ID",by.y = "IDs",all.y = T)
CrAss_G_Meta_1 <- merge(LLDeepIDs, CrAss_G_Meta_1,by.x = "LLDeepID2", by.y = "ID", all.y = T)
CrAss_G_Meta_1$LLDeepID1 <- as.character(CrAss_G_Meta_1$LLDeepID1)
CrAss_G_Meta_1$newID <- as.character(CrAss_G_Meta_1$newID)
CrAss_G_Meta_1[1136:1888,2] <- CrAss_G_Meta_1[1136:1888,142] 

CrAss_G_Meta_1 <- merge(cofactors, CrAss_G_Meta_1, by.x = "ID", by.y = "LLDeepID1", all.y = T)


# LLD vs IBD
#===========

CrAss_G_Meta_IBDLLD <- subset(CrAss_G_Meta_1, cohort != "300OB", select = c(1:5,9,12:ncol(CrAss_G_Meta_1)))
CrAss_G_Meta_IBDLLD$DiagnosisCurrent <- fct_drop(CrAss_G_Meta_IBDLLD$DiagnosisCurrent)
CrAss_G_Meta_IBDLLD$cohort <- fct_drop(CrAss_G_Meta_IBDLLD$cohort)
rownames(CrAss_G_Meta_IBDLLD) <- CrAss_G_Meta_IBDLLD$newID
CrAss_G_Meta_IBDLLD$newID <- NULL


CrAss_G_Meta_IBDLLD[c(9:ncol(CrAss_G_Meta_IBDLLD))] <- as.data.frame(ifelse(CrAss_G_Meta_IBDLLD[c(9:ncol(CrAss_G_Meta_IBDLLD))] > 0, 1, 0))

#Filter taxa present in less than 5% of LLD samples (see previous section)
Names_FilteredTaxa_1 <- c(colnames(CrAss_G_Meta_IBDLLD)[1:8],Names_FilteredTaxa)
CrAss_G_Meta_IBDLLD_fil <- CrAss_G_Meta_IBDLLD[, Names_FilteredTaxa_1]

# Logistic regression
LLD.IBD_lmG <- CrAss_G_Meta_IBDLLD_fil[c(2:6,8:ncol(CrAss_G_Meta_IBDLLD_fil))]

my_preds=c(colnames(LLD.IBD_lmG)[1:6])

LgR <- LogisticRegression.function(LLD.IBD_lmG,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.IBD_lmG, 6, 7)

LogisticRegrLLDIBD_g <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLDIBD_g, "CrAssGenera_LogisticRegression_LLDvsIBD.txt", sep = "\t", col.names = T)



# LLD vs 300OB
#=============

CrAss_G_Meta_300OBLLD <- merge(OB_keymeta, CrAss_G_Meta_1,by.x = "G_id", by.y = "ID", all.y = T)
CrAss_G_Meta_300OBLLD_1 <- CrAss_G_Meta_300OBLLD[c(1,2,5,6,69,70,79:ncol(CrAss_G_Meta_300OBLLD))]
CrAss_G_Meta_300OBLLD_1$sex <- ifelse(CrAss_G_Meta_300OBLLD_1$sex == 0, "female", "male") 

CrAss_G_Meta_300OBLLD_1[1:298,5] <- CrAss_G_Meta_300OBLLD_1[1:298,3] 
CrAss_G_Meta_300OBLLD_1[1:298,6] <- CrAss_G_Meta_300OBLLD_1[1:298,4] 

CrAss_G_Meta_300OBLLD_2 <- subset(CrAss_G_Meta_300OBLLD_1, cohort != "IBD", select = c(5,6,8:ncol(CrAss_G_Meta_300OBLLD_1)))
CrAss_G_Meta_300OBLLD_2$cohort <- fct_drop(CrAss_G_Meta_300OBLLD_2$cohort)

rownames(CrAss_G_Meta_300OBLLD_2) <- CrAss_G_Meta_300OBLLD_2$newID
CrAss_G_Meta_300OBLLD_2$newID <- NULL

CrAss_G_Meta_300OBLLD_2[c(4:ncol(CrAss_G_Meta_300OBLLD_2))] <- as.data.frame(ifelse(CrAss_G_Meta_300OBLLD_2[c(4:ncol(CrAss_G_Meta_300OBLLD_2))] > 0, 1, 0))

#Filter taxa present in less than 5% of LLD samples (see previous section)
Names_FilteredTaxa_1 <- c(colnames(CrAss_G_Meta_300OBLLD_2)[1:3],Names_FilteredTaxa)
CrAss_G_Meta_300OBLLD_2_fil <- CrAss_G_Meta_300OBLLD_2[, Names_FilteredTaxa_1]

#Logistic regression
LLD.300OB_lmG <- CrAss_G_Meta_300OBLLD_2_fil

my_preds=c(colnames(LLD.300OB_lmG)[1:3])

LgR <- LogisticRegression.function(LLD.300OB_lmG,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(LLD.300OB_lmG, 3,4)

LogisticRegrLLD300OB_g <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrLLD300OB_g, "CrAssGenera_LogisticRegression_LLDvs300OB.txt", sep = "\t", col.names = T)


```



**Metabolic Syndrome_OLD**
```{r}

#=============
# WITHIN 300OB
#=============


# ------------
# FAMILY LEVEL
# ------------

CrAss_F_300OB_MetabolSyn <- subset(CrAss_F_Meta_300OBLLD, cohort == "300OB", select = c(1,5,6,27,81:ncol(CrAss_F_Meta_300OBLLD)))
rownames(CrAss_F_300OB_MetabolSyn) <- CrAss_F_300OB_MetabolSyn$newID
CrAss_F_300OB_MetabolSyn$newID <- NULL
CrAss_F_300OB_MetabolSyn$sex <- ifelse(CrAss_F_300OB_MetabolSyn$sex == 0, "female", "male") 
CrAss_F_300OB_MetabolSyn$sex <- factor(CrAss_F_300OB_MetabolSyn$sex, levels = c("female", "male"))

CrAss_F_300OB_MetabolSyn[c(4)] <- as.data.frame(ifelse(CrAss_F_300OB_MetabolSyn[c(4)] > 0, "yes", "no"))

colnames(CrAss_F_300OB_MetabolSyn)[4] <- "PlaquesYesNo"

CrAss_F_300OB_MetabolSyn[c(5:9)] <- as.data.frame(ifelse(CrAss_F_300OB_MetabolSyn[c(5:9)] > 0, 1, 0))

within300OB_lm_f <- CrAss_F_300OB_MetabolSyn[-c(1)]

my_preds=c(colnames(within300OB_lm_f)[1:3])

LgR <- LogisticRegression.function(within300OB_lm_f,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(within300OB_lm_f, 3,4)

LogisticRegr300OB_f <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegr300OB_f, "CrAssFamily_LogisticRegression_within300OB.txt", sep = "\t", col.names = T)


# -----------
# GENUS LEVEL
# -----------

CrAss_G_300OB_MetabolSyn <- subset(CrAss_G_Meta_300OBLLD, cohort == "300OB", select = c(1,5,6,27,81:ncol(CrAss_G_Meta_300OBLLD)))
rownames(CrAss_G_300OB_MetabolSyn) <- CrAss_G_300OB_MetabolSyn$newID
CrAss_G_300OB_MetabolSyn$newID <- NULL
CrAss_G_300OB_MetabolSyn$sex <- ifelse(CrAss_G_300OB_MetabolSyn$sex == 0, "female", "male") 
CrAss_G_300OB_MetabolSyn$sex <- factor(CrAss_G_300OB_MetabolSyn$sex, levels = c("female", "male"))

CrAss_G_300OB_MetabolSyn[c(4)] <- as.data.frame(ifelse(CrAss_G_300OB_MetabolSyn[c(4)] > 0, "yes", "no"))


colnames(CrAss_G_300OB_MetabolSyn)[4] <- "PlaquesYesNo"


CrAss_G_300OB_MetabolSyn[c(5:136)] <- as.data.frame(ifelse(CrAss_G_300OB_MetabolSyn[c(5:136)] > 0, 1, 0))

#Filter taxa present in less than 5% of LLD samples (see previous section)
Names_FilteredTaxa_1 <- c(colnames(CrAss_G_300OB_MetabolSyn)[1:4],Names_FilteredTaxa)
CrAss_G_300OB_MetabolSyn_fil <- CrAss_G_300OB_MetabolSyn[, Names_FilteredTaxa_1]

within300OB_lm <- CrAss_G_300OB_MetabolSyn_fil[-c(1)]

my_preds=c(colnames(within300OB_lm)[1:3])

LgR <- LogisticRegression.function(within300OB_lm,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(within300OB_lm, 3,4)

LogisticRegr300OB <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegr300OB, "CrAssGenera_LogisticRegression_within300OB.txt", sep = "\t", col.names = T)


```

**Metabolic Syndrome_NEW**
```{r}

#=============
# WITHIN 300OB
#=============

# ------------
# ORDER LEVEL
# ------------

CrAss_O_300OB_MetabolSyn <- subset(CrAss_O_Meta_300OBLLD_1, cohort == "300OB", select = c(1,2,5,6,9:ncol(CrAss_O_Meta_300OBLLD_1)))
CrAss_O_300OB_MetabolSyn <- merge(OB_MetabolSyn,CrAss_O_300OB_MetabolSyn,by = "ID")
rownames(CrAss_O_300OB_MetabolSyn) <- CrAss_O_300OB_MetabolSyn$newID
CrAss_O_300OB_MetabolSyn$newID <- NULL
CrAss_O_300OB_MetabolSyn$Sex <- factor(CrAss_O_300OB_MetabolSyn$Sex, levels = c("female", "male"))

CrAss_O_300OB_MetabolSyn[c(2)] <- as.data.frame(ifelse(CrAss_O_300OB_MetabolSyn[c(2)] >= 3, "yes", "no"))

colnames(CrAss_O_300OB_MetabolSyn)[2] <- "MetabolSyn_YesNo"

CrAss_O_300OB_MetabolSyn[c(6)] <- as.data.frame(ifelse(CrAss_O_300OB_MetabolSyn[c(6)] > 0, 1, 0))

within300OB_lm_o <- CrAss_O_300OB_MetabolSyn[-c(1,3)]

my_preds=c(colnames(within300OB_lm_o)[1:3])

LgR <- LogisticRegression.function(within300OB_lm_o,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(within300OB_lm_o, 3,4)

LogisticRegr300OB_o <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegr300OB_o, "CrAssOrder_LogisticRegression_within300OB_New.txt", sep = "\t", col.names = T)


# ------------
# FAMILY LEVEL
# ------------

CrAss_F_300OB_MetabolSyn <- subset(CrAss_F_Meta_300OBLLD_1, cohort == "300OB", select = c(1,2,5,6,9:ncol(CrAss_F_Meta_300OBLLD_1)))
CrAss_F_300OB_MetabolSyn <- merge(OB_MetabolSyn,CrAss_F_300OB_MetabolSyn,by = "ID")
rownames(CrAss_F_300OB_MetabolSyn) <- CrAss_F_300OB_MetabolSyn$newID
CrAss_F_300OB_MetabolSyn$newID <- NULL
CrAss_F_300OB_MetabolSyn$Sex <- factor(CrAss_F_300OB_MetabolSyn$Sex, levels = c("female", "male"))

CrAss_F_300OB_MetabolSyn[c(2)] <- as.data.frame(ifelse(CrAss_F_300OB_MetabolSyn[c(2)] >= 3, "yes", "no"))

colnames(CrAss_F_300OB_MetabolSyn)[2] <- "MetabolSyn_YesNo"

CrAss_F_300OB_MetabolSyn[c(6:10)] <- as.data.frame(ifelse(CrAss_F_300OB_MetabolSyn[c(6:10)] > 0, 1, 0))

within300OB_lm_f <- CrAss_F_300OB_MetabolSyn[-c(1,3)]

my_preds=c(colnames(within300OB_lm_f)[1:3])

LgR <- LogisticRegression.function(within300OB_lm_f,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(within300OB_lm_f, 3,4)

LogisticRegr300OB_f <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegr300OB_f, "CrAssFamily_LogisticRegression_within300OB_New.txt", sep = "\t", col.names = T)


# -----------
# GENUS LEVEL
# -----------

CrAss_G_300OB_MetabolSyn <- subset(CrAss_G_Meta_300OBLLD_1, cohort == "300OB", select = c(1,2,5,6,9:ncol(CrAss_G_Meta_300OBLLD_1)))
CrAss_G_300OB_MetabolSyn <- merge(OB_MetabolSyn,CrAss_G_300OB_MetabolSyn,by = "ID")
rownames(CrAss_G_300OB_MetabolSyn) <- CrAss_G_300OB_MetabolSyn$newID
CrAss_G_300OB_MetabolSyn$newID <- NULL
CrAss_G_300OB_MetabolSyn$Sex <- factor(CrAss_G_300OB_MetabolSyn$Sex, levels = c("female", "male"))

CrAss_G_300OB_MetabolSyn[c(2)] <- as.data.frame(ifelse(CrAss_G_300OB_MetabolSyn[c(2)] >= 3, "yes", "no"))

colnames(CrAss_G_300OB_MetabolSyn)[2] <- "MetabolSyn_YesNo"

CrAss_G_300OB_MetabolSyn[c(6:137)] <- as.data.frame(ifelse(CrAss_G_300OB_MetabolSyn[c(6:137)] > 0, 1, 0))

within300OB_lm_g <- CrAss_G_300OB_MetabolSyn[-c(1,3)]

my_preds=c(colnames(within300OB_lm_g)[1:3])

LgR <- LogisticRegression.function(within300OB_lm_g,3,4,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(within300OB_lm_g, 3,4)

LogisticRegr300OB_g <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegr300OB_g, "CrAssGenus_LogisticRegression_within300OB_New.txt", sep = "\t", col.names = T)


```


**Within IBD: CD vs UC**
```{r}

#=============
# ORDER LEVEL
#=============

CrAss_O_Meta_IBD <- subset(CrAss_O_Meta_1, cohort == "IBD", select = c(1:5,9,12:15))
CrAss_O_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_O_Meta_IBD$DiagnosisCurrent)
CrAss_O_Meta_IBD$cohort <- fct_drop(CrAss_O_Meta_IBD$cohort)
rownames(CrAss_O_Meta_IBD) <- CrAss_O_Meta_IBD$newID
CrAss_O_Meta_IBD$newID <- NULL

CrAss_O_Meta_IBD[c(9)] <- as.data.frame(ifelse(CrAss_O_Meta_IBD[c(9)] > 0, 1, 0))

IBD_lm_o <- CrAss_O_Meta_IBD[c(2:7,9:ncol(CrAss_O_Meta_IBD))]

my_preds=c(colnames(IBD_lm_o)[1:6])

LgR <- LogisticRegression.function(IBD_lm_o,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm_o, 6, 7)

LogisticRegrCDvUC_o <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrCDvUC_o, "CrAssOrder_LogisticRegression_CDvsUC.txt", sep = "\t", col.names = T)

#=============
# FAMILY LEVEL
#=============

CrAss_F_Meta_IBD <- subset(CrAss_F_Meta_1, cohort == "IBD", select = c(1:5,9,12:19))
CrAss_F_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_F_Meta_IBD$DiagnosisCurrent)
CrAss_F_Meta_IBD$cohort <- fct_drop(CrAss_F_Meta_IBD$cohort)
rownames(CrAss_F_Meta_IBD) <- CrAss_F_Meta_IBD$newID
CrAss_F_Meta_IBD$newID <- NULL

CrAss_F_Meta_IBD[c(9:13)] <- as.data.frame(ifelse(CrAss_F_Meta_IBD[c(9:13)] > 0, 1, 0))

IBD_lm <- CrAss_F_Meta_IBD[c(2:7,9:ncol(CrAss_F_Meta_IBD))]

my_preds=c(colnames(IBD_lm)[1:6])

LgR <- LogisticRegression.function(IBD_lm,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm, 6, 7)

LogisticRegrCDvUC <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrCDvUC, "CrAssFamily_LogisticRegression_CDvsUC.txt", sep = "\t", col.names = T)


#============
# GENUS LEVEL
#============


CrAss_G_Meta_IBD <- subset(CrAss_G_Meta_1, cohort == "IBD", select = c(1:5,9,12:ncol(CrAss_G_Meta_1)))
CrAss_G_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_G_Meta_IBD$DiagnosisCurrent)
CrAss_G_Meta_IBD$cohort <- fct_drop(CrAss_G_Meta_IBD$cohort)
rownames(CrAss_G_Meta_IBD) <- CrAss_G_Meta_IBD$newID
CrAss_G_Meta_IBD$newID <- NULL

CrAss_G_Meta_IBD[c(9:ncol(CrAss_G_Meta_IBD))] <- as.data.frame(ifelse(CrAss_G_Meta_IBD[c(9:ncol(CrAss_G_Meta_IBD))] > 0, 1, 0))

#Filter taxa present in less than 5% of LLD samples (see previous section)
Names_FilteredTaxa_1 <- c(colnames(CrAss_G_Meta_IBD)[1:8],Names_FilteredTaxa)
CrAss_G_Meta_IBD_fil <- CrAss_G_Meta_IBD[, Names_FilteredTaxa_1]

IBD_lm_g <- CrAss_G_Meta_IBD_fil[c(2:7,9:ncol(CrAss_G_Meta_IBD_fil))]

my_preds=c(colnames(IBD_lm_g)[1:6])

LgR <- LogisticRegression.function(IBD_lm_g,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm_g, 6, 7)

LogisticRegrCDvUC_g <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrCDvUC_g, "CrAssGenera_LogisticRegression_CDvsUC.txt", sep = "\t", col.names = T)



```



**Within IBD: Disease Location**
```{r}

cofactors <- Metadata[,c(1,38)]
cofactors$DiseaseLocation <- ifelse(cofactors$DiseaseLocation == "both", "ileum",as.character(cofactors$DiseaseLocation)) 
cofactors$DiseaseLocation <- factor(cofactors$DiseaseLocation, levels = c("colon","ileum"))
cofactors$DiseaseLocation <- fct_drop(cofactors$DiseaseLocation)

#=============
# ORDER LEVEL
#=============

CrAss_O_Meta_IBD <- subset(CrAss_O_Meta_1, cohort == "IBD", select = c(1:5,9,12:15))
CrAss_O_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_O_Meta_IBD$DiagnosisCurrent)
CrAss_O_Meta_IBD$cohort <- fct_drop(CrAss_O_Meta_IBD$cohort)
rownames(CrAss_O_Meta_IBD) <- CrAss_O_Meta_IBD$newID
CrAss_O_Meta_IBD$newID <- NULL

CrAss_O_Meta_IBD[c(9)] <- as.data.frame(ifelse(CrAss_O_Meta_IBD[c(9)] > 0, 1, 0))

CrAss_O_Meta_IBD_1 <- merge(cofactors, CrAss_O_Meta_IBD, by = "ID")

IBD_lm_o <- CrAss_O_Meta_IBD_1[c(1:7,10:ncol(CrAss_O_Meta_IBD_1))]
rownames(IBD_lm_o) <- IBD_lm_o$ID
IBD_lm_o$ID <- NULL

my_preds=c(colnames(IBD_lm_o)[1:6])

LgR <- LogisticRegression.function(IBD_lm_o,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm_o, 6, 7)

LogisticRegrDiseaseLoc_o <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrDiseaseLoc_o, "CrAssOrder_LogisticRegression_DiseaseLocation.txt", sep = "\t", col.names = T)


#=============
# FAMILY LEVEL
#=============

CrAss_F_Meta_IBD <- subset(CrAss_F_Meta_1, cohort == "IBD", select = c(1:5,9,12:19))
CrAss_F_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_F_Meta_IBD$DiagnosisCurrent)
CrAss_F_Meta_IBD$cohort <- fct_drop(CrAss_F_Meta_IBD$cohort)
rownames(CrAss_F_Meta_IBD) <- CrAss_F_Meta_IBD$newID
CrAss_F_Meta_IBD$newID <- NULL

CrAss_F_Meta_IBD[c(9:13)] <- as.data.frame(ifelse(CrAss_F_Meta_IBD[c(9:13)] > 0, 1, 0))

CrAss_F_Meta_IBD_1 <- merge(cofactors, CrAss_F_Meta_IBD, by = "ID")

IBD_lm <- CrAss_F_Meta_IBD_1[c(1:7,10:ncol(CrAss_F_Meta_IBD_1))]
rownames(IBD_lm) <- IBD_lm$ID
IBD_lm$ID <- NULL

my_preds=c(colnames(IBD_lm)[1:6])

LgR <- LogisticRegression.function(IBD_lm,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm, 6, 7)

LogisticRegrDiseaseLoc <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrDiseaseLoc, "CrAssFamily_LogisticRegression_DiseaseLocation.txt", sep = "\t", col.names = T)


#============
# GENUS LEVEL
#============

CrAss_G_Meta_IBD <- subset(CrAss_G_Meta_1, cohort == "IBD", select = c(1:5,9,12:ncol(CrAss_G_Meta_1)))
CrAss_G_Meta_IBD$DiagnosisCurrent <- fct_drop(CrAss_G_Meta_IBD$DiagnosisCurrent)
CrAss_G_Meta_IBD$cohort <- fct_drop(CrAss_G_Meta_IBD$cohort)
rownames(CrAss_G_Meta_IBD) <- CrAss_G_Meta_IBD$newID
CrAss_G_Meta_IBD$newID <- NULL

CrAss_G_Meta_IBD[c(9:ncol(CrAss_G_Meta_IBD))] <- as.data.frame(ifelse(CrAss_G_Meta_IBD[c(9:ncol(CrAss_G_Meta_IBD))] > 0, 1, 0))

CrAss_G_Meta_IBD_1 <- merge(cofactors, CrAss_G_Meta_IBD, by = "ID")

#Filter taxa present in less than 5% of LLD samples (see previous section)
Names_FilteredTaxa_1 <- c(colnames(CrAss_G_Meta_IBD_1)[1:9],Names_FilteredTaxa)
CrAss_G_Meta_IBD_1_fil <- CrAss_G_Meta_IBD_1[, Names_FilteredTaxa_1]

#Logistic Regression
IBD_lm_g <- CrAss_G_Meta_IBD_1_fil[c(1:7,10:ncol(CrAss_G_Meta_IBD_1_fil))]
rownames(IBD_lm_g) <- IBD_lm_g$ID
IBD_lm_g$ID <- NULL

my_preds=c(colnames(IBD_lm_g)[1:6])

LgR <- LogisticRegression.function(IBD_lm_g,6,7,my_preds)
LgR.decrp <- LgRgr.DescriptiveStats(IBD_lm_g, 6, 7)

LogisticRegrDiseaseLoc_g <- cbind(LgR.decrp, LgR[c(3:8)])

write.table(LogisticRegrDiseaseLoc_g, "CrAssGenera_LogisticRegression_DiseaseLocation.txt", sep = "\t", col.names = T)



```


