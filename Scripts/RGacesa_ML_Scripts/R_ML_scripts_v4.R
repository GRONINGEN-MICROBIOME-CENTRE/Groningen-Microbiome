# Machine Learning scripts for microbiome analysis
#
# by: Ranko Gacesa,
#     Weersma Group, 
#     UMCG (2019)
# ==============================================

library(caret)
library(plotROC)
library(pROC)

# helper function for exctracting best tuning prediction data from Caret model
# ===========================================================================
extractBestTunePred <- function(inCaret) {
  bt <- inCaret$bestTune
  if (bt[1] == "none") {
    res <- inCaret$pred
  } else {
    res <- inCaret$pred
    for (v in colnames(bt)) {
      print(v)
      res <- res[res[[v]] == bt[[v]][1], ]
    }
  }
  res
}

# analysis of RFE profiles
# ======================================
rfeMakePlots <- function(rfeProfile) {
  g1 <- ggplot(rfeProfile[[3]],aes(x=Variables,y=Kappa)) + geom_line(size=1.15) + geom_point(col="blue",size=3) +
    geom_errorbar(aes(ymin=Kappa-KappaSD,ymax=Kappa+KappaSD),width=0.1)
  
  g2 <- ggplot(rfeProfile[[3]],aes(x=Variables,y=Accuracy)) + geom_line(size=1.15) + geom_point(col="blue",size=3) +
    geom_errorbar(aes(ymin=Accuracy-AccuracySD,ymax=Accuracy+AccuracySD),width=0.1)
  
  list(g1,g2)
}

# do train & testing, make confusion matrices
# ======================================================================
doTrainTest <- function(dataI,mdl,trainPerc=0.7) {
  inTrain <- createDataPartition(y=dataI$Diagnosis,p=trainPerc,list=F)
  trainSet <-dataI[inTrain,]
  testSet <- dataI[-inTrain,]
  t <- train(Diagnosis ~ ., data=trainSet, method=mdl)
  cTest <- confusionMatrix(predict(t,testSet),testSet$Diagnosis)
  cTrain <- confusionMatrix(predict(t,trainSet),trainSet$Diagnosis)
  return (list(t,cTrain,cTest))
}


# ============================
#
# plot all learning curves
#
# ===========================
plotAllLearningCurves <- function(dataI,dataName="",dataNameShort="",
                                  mdls=c("glm","gbm","rpart2","avNNet","nnet","pcaNNet","hdrda","svmRadial","svmLinear3",
                                         "rda","rrlda","regLogistic","rf","RRFglobal"),
                                  bts=3,sStep=2,sMin=50,outTable=T,outPlots=T,outFolder='plots',trainRepNR = 2,trainBootNR = 10,
                                  doROC = T,posClass="IBS",allowParallel=T,excludeVars=NULL) {
  # outtable
  outTbl = data.frame(Model=character(),Metric=numeric(),Value=numeric(),SD=numeric())
  # go over all models
  for (md in mdls) {
    startTime <- Sys.time()
    print (paste(">>>>>> PREPPING LEARNING CURVE FOR ",md," <<<<<<"))
    #registerDoSNOW(makeCluster(4, type = "SOCK"))
    lCurve <- prepLearningCurve(dataIn = dataI,mdl = md, minSam = sMin,lcBoots = bts,samStep = sStep,trNumber = trainRepNR,trBoot = trainBootNR,
                                saveVarImp=paste(outFolder,'/lCurve_',dataNameShort,'_',md,'_varImp.png',sep=''),
                                saveVarImpTit=paste('Covariate importance (',md,')',sep=''),
                                posClass=posClass,ROCtitle=paste('ROC: ',dataName,' [',md,']',sep=''),
                                allowParallel=allowParallel,excludeVars=excludeVars)
    g <- plotLearningCurves(lCurve,tit=paste("L-Curve",dataName,md))
    ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'_',md,'.png',sep=''),width = 9,height = 9)
    endTime <- Sys.time()
    fStats <- as.data.frame(lCurve[nrow(lCurve),3:ncol(lCurve)])
    fStats$TIME <- round(as.numeric(difftime(time1 = endTime, time2= startTime,units = "sec")))
    fStats$TIME.SD <- 0.0
    
    fStats$Model <- md
    fStatsSD <- fStats[,grep("\\.SD",colnames(fStats))]
    fStatsSD$Model <- md
    fStats <- fStats[,grep("\\.SD",colnames(fStats),invert = T)]
    fStatsLong <- gather(fStats,key = "Metric",value = "Value",Acc:TIME)
    fStatsSDLong <- gather(fStatsSD,key = "Metric",value = "SD",Acc.SD:TIME.SD)
    t <- cbind.data.frame(fStatsLong,fStatsSDLong$SD)
    t[is.na(t)] <- 0.0
    colnames(t) <- c("Model","Metric","Value","SD")
    
    outTbl <- rbind.data.frame(outTbl,t)
    
    #outTbl$Model <- as.character(outTbl$Model)
    #print(outTbl)
  }
  #for (i in c(2:9)) {outTbl[[i]] <- as.numeric(outTbl[[i]])}
  
  g <- ggplot(data=outTbl[outTbl$Metric!="TIME",],aes(x=reorder(Model,Value),y=Value,col=reorder(Model,Value))) + geom_point(size=3) + facet_grid(Metric ~ .) + 
    geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), width=0.15) + ylim(0,1) + ggtitle(paste(dataName," Metrics",sep="")) + 
    xlab("ML algorithm") + ylab("Prediction metric value") + guides(col=guide_legend(title="ML algorithm")) + theme(axis.text.x = element_text(size = rel(0.75), angle = 00))
  
  ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'__COMP.png',sep=''),width = 12,height = 10)
  
  g <- ggplot(data=outTbl[outTbl$Metric=="TIME",],aes(x=reorder(Model,Value),y=Value,fill=reorder(Model,Value))) + geom_col() + 
    ggtitle(paste(dataName," Training Time",sep="")) + xlab("ML algorithm") + ylab("Training time (seconds)") + 
    guides(fill=guide_legend(title="ML algorithm")) + theme(axis.text.x = element_text(size = rel(0.75), angle = 00))
  ggsave(g,filename = paste(outFolder,'/lCurve_',dataNameShort,'__TIME.png',sep=''),width = 12,height = 8)
  
  
  if (outTable) {
    write.csv(outTbl,file=paste(outFolder,'/lCurve_',dataNameShort,'__COMP_TABLE.csv',sep=''))
  }

  if (outTable) {
    return(outTbl)
  }
}

# ===============================================================================
# generate ROC based on RF trained on given model (70% of data)
#
# ===============================================================================
generateROC <- function(dataModel,trSet=0.7,klasses=c("IBS","HC"),tit="") {
  inTrain <- createDataPartition(y=dataModel$Diagnosis,p=trSet,list=F)
  trainSet <-dataModel[inTrain,]
  testSet <- dataModel[-inTrain,]
  fittedModel <- train(Diagnosis ~ ., data=trainSet,method="rf",metric="ROC",
                trControl=trainControl(savePrediction=T,classProbs=T,summaryFunction=twoClassSummary))
  sI <- fittedModel$pred$mtry == fittedModel$bestTune$mtry
  
  grug <- fittedModel$pred[sI,grep(klasses[1],colnames(fittedModel$pred))]
  g <- ggplot(fittedModel$pred[sI,], aes(m=grug, d=factor(obs, levels = klasses))) + 
    style_roc(xlab = "1 - Specificity / False Positive rate",
              ylab = "Sensitivity / True Positive rate") + 
    geom_roc(n.cuts = 0,color="blue",size=1.5) + 
    coord_equal() 
  auc <- fittedModel$results[fittedModel$results$mtry==fittedModel$bestTune$mtry,]$ROC
  g <- g + annotate("text", x=0.70, y=0.3, label=paste("AUC =", round(auc, 3)))
  g <- g + ggtitle(tit)
  g
}

# ===============================================================================
# generate ROC for listed model (mdl) and testSet (testSet)
# positive class is inputed one (posClass), other class is considered to be negative
# ===============================================================================
generateMdlROC <- function(mdl,testSet,posClass="IBS",tit="") {
  predP <- predict(mdl,testSet,type="prob")
  r <- roc(testSet$Diagnosis,predP[[posClass]],auc=T,percent=F)
  #cir <- ci(r,of="se", sp = seq(0, 100, 5), boot.n=100)
  #plot(r,print.auc=T,col="darkgreen")
  
  g <- pROC::ggroc(r,col="darkblue",size=1.5,alpha=0.75) +
  annotate("text", x = .75, y = .25, label = paste("AUC =",round(r$auc[1],2) )) +
    ggtitle(tit)

}
# ===============================================================================
# generate ROC for listed model (mdl) and testSet (testSet) + trainSet (trainSet)
# positive class is inputed one (posClass), other class is considered to be negative
# ===============================================================================
generateMdlROCTrainTest <- function(mdl,trainSet,testSet,posClass="IBS",tit="",doSmooth=F) {
  
  predTr <- predict(mdl,trainSet,type="prob")
  rTr <- roc(trainSet$Diagnosis,predTr[[posClass]],auc=T,percent=F)
  
  predTst <- predict(mdl,testSet,type="prob")
  rTst <- roc(testSet$Diagnosis,predTst[[posClass]],auc=T,percent=F)
  
  #g <- pROC::ggroc(r,col="darkblue",size=1.5,alpha=0.75) +
  #  annotate("text", x = .75, y = .25, label = paste("AUC =",round(r$auc[1],2) )) +
  #  ggtitle(tit)  
  if (doSmooth) {
    g <- pROC::ggroc(list(Train=rTr,Test=rTst), size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(list(Train=rTr,Test=rTst), size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    annotate("text", x = .45, y = .23, label = paste("Training Set AUC =",round(rTr$auc[1],2)),hjust = 0) +
    annotate("text", x = .45, y = .30, label = paste("Test Set AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    labs(color='Dataset')  
  g
}
# ===============================================================================
# generator of comparative ROC for MULTIPLE MODELS and ONE DATASET (test set)
# ===============================================================================
generateMdlROCComparative <- function(mdls,modelNames,testSet,posClass="IBS",tit="",doSmooth=F) {
  rocs = list()
  c = 0
  for (m in mdls) {
    c = c + 1
    namen = modelNames[c]
    pr <- predict(m,testSet,type="prob")
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
  }
  
  if (doSmooth) {
  g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  g
}
# ===============================================================================
# generator of comparative ROC for MULTIPLE DATASETS and ONE MODEL
# ===============================================================================
generateMdlROCCompareData <- function(mdl,testSetNames,testSets,posClass="IBS",tit="",doSmooth=F) {
  rocs = list()
  c = 0
  m <- mdl
  for (testSet in testSets) {
    c = c + 1
    namen = testSetNames[c]
    pr <- predict(m,testSet,type="prob")
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
  }
  
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(tit) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    geom_abline(slope = 1,intercept = 1,col="grey") +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(testSetNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(testSetNames[i],"AUC =",round(rocs[[testSetNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  g
}

# ===============================================================================
# compares model performance on different datasets
# ===============================================================================
compareMdlDatasets <- function(mdl,dataSets,dataSetNames,posClass="IBS",roc.smooth=F,tit="",
                               response="Diagnosis",roc.conf = T,roc.conf.boot = 100,target = F) {
  trainedModels <- list()
  c = 0
  # list of ROC curves
  rocs = list()
  # data frame of predictions
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  # generate results for each dataset
  for (ds in dataSets) {
    c = c + 1
    testSet <- ds
    fittedModel <- mdl
    trainedModels[[c]] <- fittedModel
    pr <- predict(fittedModel,newdata = testSet,type="prob")
    namen = dataSetNames[c]
    if (roc.smooth) {
      rocs[[namen]] <- smooth(roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F))
    } else {
      rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    }
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet[[response]])
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataSet"=namen,"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # extract coordinates
  if (!roc.conf) {
    rocsDFprecalc = NULL
    rocStep = 0.01
    for (i in seq(1,length(rocs))) {
      if (class(rocsDFprecalc) != "data.frame") {
        if (!roc.smooth) {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        } else {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        }
        rocsDFprecalc$Dataset <- names(rocs)[i]
        rocsDFprecalc$sz = 1
      } else {
        if (!roc.smooth) {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        } else {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x=c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep)),ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        }
        rocsDFprecalctmp$Dataset <- names(rocs)[i]
        rocsDFprecalctmp$sz = 1
        rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
        rownames(rocsDFprecalc) <- NULL
      }
    }
  }
  
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = roc.conf.boot
    meanConf <- F
    for (i in seq(1,length(rocs))) {
      if (roc.smooth) {
        sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities =c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots,progress = "none")
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots,progress="none")
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- as.numeric(rownames(rocConf))
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf$Dataset <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }
      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
      if (!roc.smooth) {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Dataset=names(rocs)[i],sz=1))
      }
    }
    rocsDFprecalc <- rocDFConf
  }
  
  if (roc.conf) {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$specificity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes(y=specificity,x=sensitivity,col=Dataset)) +
      geom_ribbon(data=rocsDFprecalc,aes(ymin=sp.low,ymax=sp.high,fill=Dataset),alpha=0.2,colour=NA) + geom_line(size=1.25)
  } else {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$sensitivity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes(y=specificity,x=sensitivity,col=Dataset)) + geom_line(size=1.25) #pROC::ggroc(rocs) + geom_line(size=1.25)
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1) +
    ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +  
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(dataSetNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(dataSetNames[i],"AUC =",round(rocs[[dataSetNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF$DataSet )) {
      if (u != "multi.class") {
        if (class(tarDF) != "data.frame") {
          tarDF <- data.frame(DataSet=u,
                              Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                              Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                              F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                              B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value
          )
        } else {
          tarDF <- rbind.data.frame(tarDF,data.frame(DataSet=u,
                                                     Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                                                     Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                                                     F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                                                     B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value))
        }
      }
    }
    g <- g + geom_point(data=tarDF,aes(x=Specificity,y=Sensitivity,col=DataSet),size=4,shape=4,stroke=2)
  }
  
  # do metrics plot
  unik <- length(dataSets)
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataSet,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(dataSetNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}


# ===============================================================================
# compares training CV of multiple trained models
# -> dataMdls should be list( dataframes(Diagnosis ~ covariate) )
# -> note: models should predict Diagnosis, with positive class = posClass
# ===============================================================================
compareModelsTrainingCV <- function(fittedMdls,modelNames,mtd="glm",posClass="IBS",roc.smooth=F,tit="",
                                    annotateAUConly=F,roc.conf = T,conf.lvl = 0.95,roc.conf.boot=10,
                                    annotS = 5,removeLegend=T,textSize=15) {
  c = 0
  rocs = list()
  for (m in fittedMdls) {
    c = c + 1
    namen = modelNames[c]
    bestPred <- extractBestTunePred(m)
    rocs[[namen]] <- roc(predictor = bestPred[[posClass]],response = bestPred$obs,auc=T,percent=F)
    
    # rocs[[namen]] <- roc(predictor = fittedMdls[[c]]$pred[[posClass]], response = fittedMdls[[c]]$pred$obs,auc=T,percent = F,
    #                      smooth = roc.smooth)
  }
  
  # extract coordinates
  if (!roc.conf) {
    rocsDFprecalc = NULL
    rocStep = 0.01
    senses <- c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep))
    for (i in seq(1,length(rocs))) {
      if (class(rocsDFprecalc) != "data.frame") {
        if (!roc.smooth) {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        } else {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x=senses,ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        }
        rocsDFprecalc$name <- names(rocs)[i]
        rocsDFprecalc$sz = 1
      } else {
        if (!roc.smooth) {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        } else {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x=senses,ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        }
        rocsDFprecalctmp$name <- names(rocs)[i]
        rocsDFprecalctmp$sz = 1
        rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
        rownames(rocsDFprecalc) <- NULL
      }
    }
  }
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = roc.conf.boot
    meanConf <- F
    for (i in seq(1,length(rocs))) {
      senses <- c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep))
      if (roc.smooth) {
        sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities =senses,
                         conf.level=conf.lvl,boot.n = boots,progress = "none")
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=senses,
                         conf.level=conf.lvl,boot.n = boots,progress = "none")
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- senses
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf$name <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }
      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
      rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,name=names(rocs)[i],sz=1))
    }
    rocsDFprecalc <- rocDFConf
  }
  
  if (roc.conf) {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$specificity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col="name")) +
      geom_ribbon(data=rocsDFprecalc,aes_string(ymin="sp.low",ymax="sp.high",fill="name"),alpha=0.2,colour=NA) + geom_line(size=1.25)
  } else {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$sensitivity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col="name")) + geom_line(size=1.25) #pROC::ggroc(rocs) + geom_line(size=1.25)
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1) +
    ggtitle(paste(tit,"ROC <xval>") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +  
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color="Model",fill="Model") + theme(text = element_text(size = textSize))
  
  # annotate w auc
  c = -1
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    if (!annotateAUConly) {
      if (length(modelNames) == 1) {
        g <- g + annotate("text", x = .6, y = 0.00+0.06*c, label = paste(" AUC = ",round(rocs[[modelNames[i] ]]$auc[1],2),
                                                                       "; Kappa = ",round(getTrainPerf(fittedMdls[[i]])[[2]],2),
                                                                       "; ACC = ",round(getTrainPerf(fittedMdls[[i]])[[1]],2),sep=""),hjust = 0, size=annotS)
      } else {
        g <- g + annotate("text", x = .6, y = 0.00+0.06*c, label = paste(modelNames[i]," AUC = ",round(rocs[[modelNames[i] ]]$auc[1],2),
                                                                         "; Kappa = ",round(getTrainPerf(fittedMdls[[i]])[[2]],2),
                                                                         "; ACC = ",round(getTrainPerf(fittedMdls[[i]])[[1]],2),sep=""),hjust = 0, size=annotS)
      }
      
    } else {
      if (roc.conf) {
        rocAUCconf <- ci.auc(rocs[[ modelNames[i] ]],conf.level=conf.lvl,boot.n = boots,progress = "none")
        rocAUCm <- rocAUCconf[[2]]
        rocAUCdelta <- rocAUCconf[[2]] - rocAUCconf[[1]]
        if (length(modelNames) == 1) {
          g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste("AUC =",round(rocAUCm,2),"\u00b1",round(rocAUCdelta,2)),hjust = 0, size=annotS)
        } else {
          g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocAUCm,2),"\u00b1",round(rocAUCdelta,2)),hjust = 0, size=annotS)  
        }
      } else {
        if (length(modelNames) == 1) {
          g <- g + annotate("text", x = .6, y = 0.00+0.06*c, label = paste("AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0, size=annotS)
        } else {
          g <- g + annotate("text", x = .6, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0, size=annotS)
        }
      }
    }
  }
  if (removeLegend) {
    g <- g + theme(legend.position = "none")
  }
  #return plot
  g
}

# ===============================================================================
# compares multiple data models - each one is built from its own dataset
# does the model training
# -> dataMdls should be list( dataframes(Diagnosis ~ covariate) )
# -> note: models should predict Diagnosis, with positive class = posClass
# ===============================================================================
trainCompareModels <- function(dataMdls,modelNames,mtd="glm",posClass="IBS",doSmooth=F,trSet=0.7,tit="") {
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  #trainedModels = list()
  for (m in dataMdls) {
    c = c + 1
    inTrain <- createDataPartition(y=m$Diagnosis,p=trSet,list=F)
    trainSet <- m[inTrain,]
    testSet <- m[-inTrain,]
    fittedModel <- train(Diagnosis ~ ., data=trainSet,method=mtd,metric="Kappa")
    pr <- predict(fittedModel,testSet,type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"ROC (",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Data Model') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + xlab("Data Model") +
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics (",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}

# ===============================================================================
# compares multiple models on multiple test datasets
# (each model with paired with appropriate test-set in order they are entered)
# does not do model training
# -> fittedMdls should be list( fitted models ) generated by caret train
# ===============================================================================
compareModelsOnDatasets <- function(fittedMdls,testSets,modelNames,mtd="glm",posClass="IBS",doSmooth=F,tit="",annotateAUConly=F) {
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  for (m in fittedMdls) {
    c = c + 1
    #trainedModels[[c]] <- m
    pr <- predict(m,testSets[[c]],type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSets[[c]]$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(m,testSets[[c]])
    conf <- confusionMatrix(pr2,testSets[[c]]$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"(",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = -1.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    if (annotateAUConly) {
      g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    } else {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(modelNames[i]," AUC = ",round(rocs[[modelNames[i] ]]$auc[1],2),
                                                                     "; Kappa = ",round(predDF[predDF$DataModel==modelNames[i] & predDF$Metric=="Kappa",]$Value,2),
                                                                     "; ACC = ",round(predDF[predDF$DataModel==modelNames[i] & predDF$Metric=="Acc",]$Value,1),sep=""),hjust = 0)
    }
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"(",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}



# ===============================================================================
# compares multiple models on one dataset
# does not do model training
# -> fittedMdls should be list( fitted models ) generated by caret train
# ===============================================================================
compareModelsOnDataset <- function(fittedMdls,testSet,modelNames,mtd="glm",posClass="IBS",doSmooth=F,tit="") {
  #trainedModels <- list()
  c = 0
  rocs = list()
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  for (m in fittedMdls) {
    c = c + 1
    #trainedModels[[c]] <- m
    pr <- predict(m,testSet,type="prob")
    namen = modelNames[c]
    rocs[[namen]] <- roc(testSet$Diagnosis,pr[[posClass]],auc=T,percent=F)
    pr2 <- predict(m,testSet)
    conf <- confusionMatrix(pr2,testSet$Diagnosis)
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataModel"=modelNames[c],"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Det.Rate","Value"=conf$byClass[["Detection Rate"]],stringsAsFactors=F))
    #predDF <- rbind.data.frame(predDF,
    #                           data.frame("DataModel"=modelNames[c],"Metric"="Prevalence","Value"=conf$byClass[["Prevalence"]],stringsAsFactors=F))
  }
  
  # do combined ROC plot
  if (doSmooth) {
    g <- pROC::ggroc(rocs, size=1.25,alpha=0.0) + geom_smooth(size=1.5,alpha=0.08)
  } else { 
    g <- pROC::ggroc(rocs, size=1.25,alpha=1.0) 
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1.25,linetype="longdash") +
    ggtitle(paste(tit,"ROC (",mtd,")") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color='Dataset') 
  # annotate w auc
  c = 0.0
  for (i in order(modelNames,decreasing = T)) {
    c = c + 1.0
    g <- g + annotate("text", x = .35, y = 0.00+0.06*c, label = paste(modelNames[i],"AUC =",round(rocs[[modelNames[i] ]]$auc[1],2)),hjust = 0)
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # do metrics plot
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataModel,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics (",mtd,")") ) +
    geom_vline(xintercept =  c(1:length(modelNames)) ,size=40,alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}

# ===============================================================================
# find correlated covariates
# note: input data should not have NAs
# ===============================================================================
findCorrelatedCovariates <- function(inData,cutOff) 
{
  numa <- sapply(inData, is.numeric)
  #print (numa)
  if (sum(numa) > 1) {
    correlationMatrix <- cor(inData[,sapply(inData, is.numeric)])
    cs <- colnames(inData)[sapply(inData, is.numeric)]
    # calculate correlation matrix
    #print(correlationMatrix)
    correlationMatrix[is.na(correlationMatrix)] <- 0.0
    # find attributes that are highly corrected (ideally >0.75)
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=cutOff,verbose = TRUE)
    # print indexes of highly correlated attributes
    print(cs[highlyCorrelated])
  } else {
    print ('less then 2 numeric covariates!')
  }
}

findCorrelatedCovariates2 <- function(inData,cutOff) 
{
  correlationMatrix <- cor(inData[,sapply(inData, is.numeric)])
  cs <- colnames(inData)[sapply(inData, is.numeric)]
  for (i in 1:length(cs)) {
    cs[i] <- strsplit(cs[i],'\\.')[[1]][length(strsplit(cs[i],'\\.')[[1]])]
    #print (cs[i])
  }
  # calculate correlation matrix
  #print(correlationMatrix)
  correlationMatrix[is.na(correlationMatrix)] <- 0.0
  cntr = 0
  for (i in 1:nrow(correlationMatrix)){
    for (j in 1:i) {
      if ( !is.na(correlationMatrix[i,j]) & !(i == j)){
        if ( correlationMatrix[i,j] > cutOff) {
          cntr = cntr + 1
          print(paste(cntr,cs[i], "-" , cs[j], ": ", correlationMatrix[i,j]))
        }
      }
    }
  }
}


#' Prepares learning curves (but does not draw them)
#'
#'
#'
prepLearningCurve <- function(dataIn,mdl,lcBoots=1,trSet=0.7,samStep=1.5,minSam=25,scaleCenter=T,trBoot = 25,trNumber = 5,trainType="boot",
                              pProc = c("center","scale"), saveVarImp="", saveVarImpTit="",responseVar="Diagnosis",posClass="IBS",ROCtitle='',
                              allowParallel=T,optimiseMetric="Kappa",excludeVars=NULL)
{
  inTrain <- createDataPartition(y=dataIn[[responseVar]],p=trSet,list=F)
  trainSet <- dataIn[inTrain,]
  testSet <- dataIn[-inTrain,]
  
  p <- c()
  pp <- minSam / nrow(trainSet)
  if (pp <= 0) {pp = 0.01}
  while (pp <= 1.5) {
    p <- c(p,pp)
    pp <- pp * samStep
    if (pp > 1.0) {pp <- 1.0; p <- c(p,pp); break}
  }
  results = data.frame()
  #NR, Name, "Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD"
  results <- rbind.data.frame(results,c(0.0,"Train",1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0,1.0,0.0),stringsAsFactors = F)
#  results <- as.data.frame(apply(results,MARGIN = 2,FUN= function(x) as.character(x)))
  results <- rbind.data.frame(results,c(0.0,"Test",0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0),stringsAsFactors = F)
#  results <- as.data.frame(apply(results,MARGIN = 2,FUN= function(x) as.character(x)))
  
  cnt = 0
  for (n in p) { # numbers loop
    #print(n)
    # vectors for accuracy measures
    cnt = cnt + 1
    acc = c()
    acctrain = c()
    kappa = c()
    kappatrain = c()
    sens <- c()
    spec <- c()
    f1 <- c()
    ppv <- c()
    npv <- c()
    senstrain <- c()
    spectrain <- c()
    f1train <- c()
    ppvtrain <- c()
    npvtrain <- c()
    auctrain <- c()
    auc <- c()
    
    for (b in seq(1:lcBoots)) { # bootstraps loop
      # another partition
      inTrain <- createDataPartition(y=dataIn[[responseVar]],p=trSet,list=F)
      trainSet <- dataIn[inTrain,]
      testSet <- dataIn[-inTrain,]
      # subselect
      #smplTr <- createDataPartition(y=trainSet[[responseVar]],p=trSet,list=F)
      #smpl <- trainSet[sample(nrow(trainSet), n), ]
      smplTake <- createDataPartition(y=trainSet[[responseVar]],p=n,list=F)
      smpl <- trainSet[smplTake,]
      
      print(paste("Test",cnt,"(",nrow(smpl),"cases);","bootstrap",b))
      # prep traincontrol
      if (optimiseMetric=="ROC") {
        tC = trainControl(method=trainType,repeats = trNumber,number = trBoot,classProbs=T,savePredictions = T,allowParallel = allowParallel,summaryFunction = twoClassSummary)
      } else {
        tC = trainControl(method=trainType,repeats = trNumber,number = trBoot,classProbs=T,savePredictions = T,allowParallel = allowParallel)
      }
      #print(n)
      
      result <- tryCatch({
        print (paste('train Set size =',nrow(smpl)))
        if (mdl %in% c("gbm")) {
          mFit <- train(reformulate(response=responseVar,termlabels = '.'), data=smpl,method=mdl,metric="Kappa",verbose=F, trControl=tC)
        } else {
          mFit <- train(reformulate(response=responseVar,termlabels = '.'), data=smpl,method=mdl,metric="Kappa", trControl=tC)
        }
        
        print (paste('test Set size =',nrow(testSet)))
        
        # ROC (test)
        # ==========================
        # > exclusion for testing 
        if (!is.null(excludeVars)) {
          testSet2 <- testSet[,!colnames(testSet) %in% excludeVars]
        } else {
          testSet2 <- testSet
        }
        
        predP <- predict(mFit,newdata=testSet2,type="prob")
        # resolve multi-class AUC
        #    2-class prediction
        if (length(levels(testSet2[[responseVar]])) == 2) {
          cauc <- roc(testSet2[[responseVar]],predP[[posClass]],auc=T)$auc+0.0
        } else {
        #    3+ class prediction: we use multiclass.roc
          cauc <- as.numeric(multiclass.roc(testSet2[[responseVar]],predP)$auc+0.0)
        }
        # confusion (test)
        pred <- predict(mFit,newdata=testSet2)
        conf <- confusionMatrix(pred,testSet2[[responseVar]])
        cacc <- conf$overall[['Accuracy']]
        ckappa <- conf$overall[['Kappa']]
        # > 2-class
        if (length(levels(testSet2[[responseVar]])) == 2) {
          csens <- conf$byClass[['Sensitivity']]
          cspec <- conf$byClass[['Specificity']]
          cppv <- conf$byClass[['Pos Pred Value']]
          cnpv <- conf$byClass[['Neg Pred Value']]
          cf1 <- conf$byClass[['Balanced Accuracy']]
        } else {
          # multi-class prediction: get averages of Sensitivity/Specificity/PPV/NPV/...
          csens <- mean(as.data.frame(conf$byClass)[['Sensitivity']],na.rm = T)
          cspec <-  mean(as.data.frame(conf$byClass)[['Specificity']],na.rm = T)
          cppv <-  mean(as.data.frame(conf$byClass)[['Pos Pred Value']],na.rm = T)
          cnpv <-  mean(as.data.frame(conf$byClass)[['Neg Pred Value']],na.rm = T)
          cf1 <-  mean(as.data.frame(conf$byClass)[['Balanced Accuracy']],na.rm = T)
        }
        
        
        # ROC (train)
        # > exclusion for testing 
        if (!is.null(excludeVars)) {
          smpl <- smpl[,!colnames(smpl) %in% excludeVars]
        }
        predtP <- predict(mFit,newdata=smpl,type="prob")
        
        if (length(levels(testSet2[[responseVar]])) == 2) {
          cauct <- as.numeric(roc(smpl[[responseVar]],predtP[[posClass]],auc=T)$auc)+0.0
        } else {
          #    3+ class prediction: we use multiclass.roc
          cauct <- as.numeric(multiclass.roc(smpl[[responseVar]],predtP)$auc)+0.0
        }
        
        # confusion (train)
        predt <- predict(mFit,newdata=smpl)
        conft <- confusionMatrix(predt,smpl[[responseVar]])
        cacct <- conft$overall[['Accuracy']]
        ckappat <- conft$overall[['Kappa']]
        if (length(levels(testSet2[[responseVar]])) == 2) {
          csenst <- conft$byClass[['Sensitivity']]
          cspect <- conft$byClass[['Specificity']]
          cppvt <- conft$byClass[['Pos Pred Value']]
          cnpvt <- conft$byClass[['Neg Pred Value']]
          cf1t <- conft$byClass[['Balanced Accuracy']]
        } else {
          csenst <- mean(as.data.frame(conft$byClass)[['Sensitivity']],na.rm = T)
          cspect <- mean(as.data.frame(conft$byClass)[['Specificity']],na.rm = T)
          cppvt <- mean(as.data.frame(conft$byClass)[['Pos Pred Value']],na.rm = T)
          cnpvt <- mean(as.data.frame(conft$byClass)[['Neg Pred Value']],na.rm = T)
          cf1t <- mean(as.data.frame(conft$byClass)[['Balanced Accuracy']],na.rm = T)
        }
        
        # debug:
        #print (paste("TRAIN: SPEC",cspec))
        #print (paste("TRAIN: SENS",csens))       
        #print (paste("TEST: SPEC",cspect))
        #print (paste("TEST: SENS",csenst))       
        
        print ('Success: training OK!')
        #print(paste("train: acc:",cacct,"kappa",ckappat,"test: acc:",cacc,'kappa:',ckappa))        
        #return(prepLearningCurveTrainTest(smpl,mdl,testSet))
        c(cacc,ckappa,cacct,ckappat,csens,cspec,cf1,cppv,cnpv,csenst,cspect,cf1t,cppvt,cnpvt,cauc,cauct)
      }, error = function(err) {
        print (paste('Error: training failed:',err))
        return(c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0))
      })
      #print(result)
      acc <- c(acc,result[1])
      kappa <- c(kappa,result[2])
      acctrain <- c(acctrain,result[3])
      kappatrain <-c(kappatrain,result[4])
      sens <- c(sens,result[5])
      spec <- c(spec,result[6])
      f1 <- c(f1,result[7])
      ppv <- c(ppv,result[8])
      npv <- c(npv,result[9])
      senstrain <- c(senstrain,result[10])
      spectrain <- c(spectrain,result[11])
      f1train <- c(f1train,result[12])
      ppvtrain <- c(ppvtrain,result[13])
      npvtrain <- c(npvtrain,result[14])
      auc <- c(auc,result[15])
      auctrain <- c(auctrain,result[16])
      
      #print (acc)
      #print (acctrain)
      #print (kappa)
      #print (kappatrain)
      #print (paste(" -> acc:",acc,";kappa:",kappa))
      #results <- rbind(results,c(nrow(smpl),acc,kappa))
    }
    accM = mean(acc)
    accSD = sd(acc)
    kappaM = mean(kappa)
    kappaSD = sd(kappa)
    accTrainM = mean(acctrain)
    accTrainSD = sd(acctrain)
    kappaTrainM = mean(kappatrain)
    kappaTrainSD = sd(kappatrain)
    
    sensM <- mean(sens)
    sensSD <- sd(sens)
    specM <- mean(spec)
    specSD <- sd(spec)
    f1M <- mean(f1)
    f1SD <- sd(f1)
    ppvM <- mean(ppv)
    ppvSD <- sd(ppv)
    npvM <- mean(npv)
    npvSD <- sd(npv)
    sensTrainM <- mean(senstrain)
    sensTrainSD <- sd(senstrain)
    specTrainM <- mean(spectrain)
    specTrainSD <- sd(spectrain)
    f1TrainM <- mean(f1train)
    f1TrainSD <- sd(f1train)
    ppvTrainM <- mean(ppvtrain)
    ppvTrainSD <- sd(ppvtrain)
    npvTrainM <- mean(npvtrain)
    npvTrainSD <- sd(npvtrain)
    
    aucTrainM <- mean(auctrain)
    aucTrainSD <- sd(auctrain)
    aucTestM <- mean(auc)
    aucTestSD <- sd(auc)
    
    #print(paste(accM,accSD,kappaM,kappaSD))
    results <- rbind.data.frame(results,c(nrow(smpl),"Train",accTrainM,accTrainSD,kappaTrainM,kappaTrainSD,
                               sensTrainM,sensTrainSD,specTrainM,specTrainSD,f1TrainM,f1TrainSD,
                               ppvTrainM,ppvTrainSD,npvTrainM,npvTrainSD,aucTrainM,aucTrainSD))
    results <- rbind.data.frame(results,c(nrow(smpl),"Test",accM,accSD,kappaM,kappaSD,
                               sensM,sensSD,specM,specSD,f1M,f1SD,
                               ppvM,ppvSD,npvM,npvSD,aucTestM,aucTestSD))

    #print (paste(" -> Test ",cnt,'done;','Mean Acc:',accM,'; Mean Kappa:',kappaM))
    colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
    #print(results)
    for (i in c(1,3:ncol(results))) {
      #print(i)
      #print(results[[i]])
      results[[i]] = as.numeric(as.character(results[[i]]))
    }
    results[[2]] <- as.character(results[[2]])    
    colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
  }
  results <- as.data.frame(results)
  for (i in c(1,3:ncol(results))) {
    results[[i]] <- as.numeric(as.character(results[[i]]))
  }
  colnames(results) <- c("N","Dataset","Acc","Acc.SD","Kappa","Kappa.SD","SENS","SENS.SD","SPEC","SPEC.SD","BACC","BACC.SD","PPV","PPV.SD","NPV","NPV.SD","AUC","AUC.SD")
  results[is.na(results)] <- 0.0
  #results[is.nan(results)] <- 0.0
  results
}

plotLearningCurves <- function(lVector,tit,metrics=c("ACC","B.ACC","Kappa","NPV","PPV","Sensitivity","Specificity")) {
  vec1 <- lVector[,c(1,2,3,4)]
  vec1$Metric <- "ACC"
  colnames(vec1) <- c("N","Dataset","Value","SD","Metric")
  vec2 <- lVector[,c(1,2,5,6)]
  vec2$Metric <- "Kappa"
  colnames(vec2) <- c("N","Dataset","Value","SD","Metric")
  vec3 <- lVector[,c(1,2,7,8)]
  vec3$Metric <- "Sensitivity"
  colnames(vec3) <- c("N","Dataset","Value","SD","Metric")
  vec4 <- lVector[,c(1,2,9,10)]
  vec4$Metric <- "Specificity"
  colnames(vec4) <- c("N","Dataset","Value","SD","Metric")
  vec5 <- lVector[,c(1,2,11,12)]
  vec5$Metric <- "B.ACC"
  colnames(vec5) <- c("N","Dataset","Value","SD","Metric")
  vec6 <- lVector[,c(1,2,13,14)]
  vec6$Metric <- "PPV"
  colnames(vec6) <- c("N","Dataset","Value","SD","Metric")
  vec7 <- lVector[,c(1,2,15,16)]
  vec7$Metric <- "NPV"
  colnames(vec7) <- c("N","Dataset","Value","SD","Metric")
  vec8 <- lVector[,c(1,2,17,18)]
  vec8$Metric <- "AUC"
  colnames(vec8) <- c("N","Dataset","Value","SD","Metric")
  
  vecLong <- rbind.data.frame(vec1,vec2)
  vecLong <- rbind.data.frame(vecLong,vec3)
  vecLong <- rbind.data.frame(vecLong,vec4)
  vecLong <- rbind.data.frame(vecLong,vec5)
  vecLong <- rbind.data.frame(vecLong,vec6) #NPV
  vecLong <- rbind.data.frame(vecLong,vec7) #PPV
  vecLong <- rbind.data.frame(vecLong,vec8) #AUC
  vecLong$Metric <- as.factor(vecLong$Metric)
  vecLong$N <- as.numeric(vecLong$N)
  vecLong$Value <- as.numeric(vecLong$Value)
  vecLong$SD <- as.numeric(vecLong$SD)
  vecLong <- vecLong[vecLong$Metric %in% metrics,]
  
  #pd <- position_dodge(0.25)
  ggplot(data=vecLong,aes(x=N,y=Value,col=Dataset)) + 
    geom_line(linetype="dashed") +
    geom_point(size=3)+
    geom_errorbar(aes(ymin=Value-SD, ymax=Value+SD), colour="black", width=3) + 
    ylim(0,1) +
    xlab("Number of cases in Training set") +
    ylab("Prediction success Metric") +
    facet_grid(Metric ~ .) +
    ggtitle(tit)
}

#' Calculates prediction metrics for one model and test set
#' 
#' \code{calcPredictionMetrics} returns vector of prediction metrics
#' includes Accuracy, Sensitivity, Specificity, F1 and Kappa
#' 
#' @param testSet : test set (must be compatible with training set)
#' @param modelNames : vector of strings, names of used models, must match modles
#' @param ... : fitted models (produced by caret train)
calcPredictionMetrics <- function(mdl,testSet)
{
  pred <- predict(mdl,testSet)
  conf <- confusionMatrix(pred,testSet$Diagnosis)
  acc <- conf$overall[["Accuracy"]]
  sens <- conf$byClass[["Sensitivity"]]
  spec <- conf$byClass[["Specificity"]]
  f1 <- conf$byClass[["F1"]]
  kappa <- conf$overall[["Kappa"]]
  res <- c("Acc"=round(acc,3) ,"Kappa"=round(kappa,3),"Sens"=round(sens,3),"Spec"=round(spec,3),"F1"=round(f1,3))
  res
}

#' Calculates prediction metrics for set of models and test set
#' 
#' \code{calcPredictionMetricsTable} returns data frame of calculated ML prediction metrics
#' 
#' @param testSet : test set (must be compatible with training set)
#' @param modelNames : vector of strings, names of used models, must match modles
#' @param ... : fitted models (produced by caret train)
calcPredictionMetricsTable <- function(testSet,modelNames,...)
{
  #print (list(...))
  c = 0
  for (m in list(...)) {
    c = c + 1
    r <- calcPredictionMetrics(m,testSet)
    r <- c("Model"=modelNames[c],r)
    if (c == 1) {
      res <- r
    } else {
      res <- rbind(res,r)
    }
  }
  res <- as.data.frame(res)
  rownames(res) <- c(1:nrow(res))
  res
}

# examines list of trained models (made by Caret)
# and outputs predictive power metrics for each of them
# on internal cross-validation set
# ==========================================================
getModelsTrainingCVstats <- function(fittedMdls,modelNames,posClass)
{
  cnt = 0
  resultsCV <- NULL
  for (fittedMdl in fittedMdls) {
    cnt = cnt + 1
    md <- modelNames[cnt]
    print(paste0(' >> analyzing x-validation performance for ',md))
    predBest <- extractBestTunePred(fittedMdl)
    
    confM <- confusionMatrix(predBest$pred,predBest$obs)
    # acc, kappa
    acc <- confM[[3]][1]
    kappa <- confM[[3]][2]
    # sens spec ppv npv prec recall F1 prev bacc
    sens <- confM[[4]][1]
    spec <- confM[[4]][2]
    ppv <- confM[[4]][3]
    npv <- confM[[4]][4]
    f1 <- confM[[4]][7]
    bacc <- confM[[4]][11]
    # get auc
    r <- roc(predictor = predBest[[posClass]], response = predBest$obs, auc=T,percent=T,smooth=F)
    #r <- roc(predictor = fittedMdl$pred[[posClass]], response = fittedMdl$pred$obs, auc=T,percent=T,smooth=F)
    auc <- as.numeric(r$auc)/100
    nrVars <- length(fittedMdl$coefnames)
    oneRow <- data.frame(Model=md,Method=fittedMdl$method,NrVars=nrVars,ACC=acc,Kappa=kappa,Sensitivity=sens,Specificity=spec,PPV=ppv,NPV=npv,B.ACC=bacc,F1=f1,AUC=auc)
    resultsCV <- rbind.data.frame(resultsCV,oneRow)
    row.names(resultsCV) <- NULL
  }
  resultsCV
}

# MULTI CLASS PREDICTION FUNCTIONS
# ===========================================

#' Calculates prediction metrics and ROC curves for fitted ML/Caret multi-class prediction model, 
#' 
#' \code{generateMultiClassROC} returns list of ROC gg2plot, metrics gg2plot and metrics dataframe
#' 
#' @param fittedML : fitted Caret ML model, must enable multi-class prediction and be able to return prediction likelihood value
#' @param testSet : test set (must be compatible with training set)
#' @param tit : name of model, used for chart titles
#' @param respName : COLUMN name of response variable (to be predicted) in test set
#' @param multiROC : if TRUE, will generate averaged/multiROC curve(s)
#' @param oneAllROC : if TRUE, will generate one-vs-all ROC curves
generateMultiClassROC <- function(fittedML,testSet,tit,respName,roc.conf=T,multiROC=T,oneAllROC = T, 
                                  multiROCPriority = F, smooth=F,dottedNonPriority=T,confBoots=250,
                                  mcROCmetric = T,mcKappa = T, mcAcc = T,target = T,tarAcc = T) {
  
  # predict on test set
  pred <- predict(fittedML,testSet,type="prob")
  predAg <- as.numeric(predict(fittedML,testSet,type="raw"))
  mcAUC <- auc(multiclass.roc(testSet[[respName]], predAg))  
  
  # MAKE ROC CURVE(S):
  rocs <- list()
  rocsDFprecalc = NULL
  for (i in seq(1,ncol(pred))) {
    rocs[[i]] <- roc(testSet[[respName]]==colnames(pred)[i],pred[[i]],percent=F)
  }
  names(rocs) <- colnames(pred)
  for (i in seq(1,ncol(pred))) {
    if (class(rocsDFprecalc) != "data.frame") {
      rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"),transpose = TRUE)))
      rocsDFprecalc$outcome <- names(rocs)[i]
      rocsDFprecalc$sz = 1
    } else {
      rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"),transpose = TRUE )))
      rocsDFprecalctmp$outcome <- names(rocs)[i]
      rocsDFprecalctmp$sz = 1
      rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
    }
  }
  
  # make data frame with results
  # ==========================================
  # conf matrix
  conf <- confusionMatrix(testSet[[respName]],predict(fittedML,testSet,type="raw"))
  # correct names of classes
  rownames(conf$byClass) <- gsub("Class: ","",rownames(conf$byClass))
  #rownames <- names(rocs)
  # prep results DF
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  # extract overall statistics
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
  predDF <- rbind.data.frame(predDF,
                             data.frame("DataSet"="multi.class","Metric"="ROC.AUC","Value"=round(mcAUC,2),stringsAsFactors=F))
  
  # extract per class statistics
  for (i in seq(1,nrow(conf$byClass))) {
    cc <- conf$byClass[i,]
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="Sensitivity","Value"=cc[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="Specificity","Value"=cc[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="PPV","Value"=cc[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="NPV","Value"=cc[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="B.Acc","Value"=cc[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="F1","Value"=cc[["F1"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("DataSet"=rownames(conf$byClass)[i],"Metric"="ROC.AUC","Value"=round(auc(rocs[[rownames(conf$byClass)[i]]]),2),stringsAsFactors=F))
  }
  dataSetNames <- rownames(conf$byClass)
  unik <- length(conf$byClass)
  gg2 <- ggplot(data=predDF,aes(col=Metric,y=Value,x=DataSet,shape=Metric)) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Multi-class prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(dataSetNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  
  # aggregate curve to calculate mean curve
  # note: we still use non-mean curves for plot of non-mean curves
  rocCoords <- F
  rocDF <- F
  rocStep = 0.01
  for (i in seq(1,ncol(pred))) {
    if (smooth) {
      rtmp <- as.data.frame(t(coords(smooth(rocs[[i]]),x=seq(0,1,rocStep),input="sensitivity", ret = c("specificity","sensitivity"),transpose = T )))
    } else {
      rtmp <- as.data.frame(t(coords((rocs[[i]]),x=seq(0,1,rocStep),input="sensitivity", ret = c("specificity","sensitivity"),transpose = T )))
    }
    if (class(rocCoords) != "data.frame") {
      rocDF <- rtmp
      rocCoords <- rtmp
      rocDF$outcome <- names(rocs)[i]
    } else {
      rocCoords <- cbind.data.frame(rocCoords,rtmp)
      rocDFtmp <- rtmp
      rocDFtmp$outcome <- names(rocs)[i]
      rocDF <- rbind.data.frame(rocDF,rocDFtmp)
    }
  }
  rocDF$sz = 1
  rocMeanSens <- apply(rocCoords[,grep("sensitivity",colnames(rocCoords))],MARGIN = 1,mean)
  rocMeanSpec <- apply(rocCoords[,grep("specificity",colnames(rocCoords))],MARGIN = 1,mean)
  rocMDF <- data.frame(specificity=rocMeanSpec,sensitivity=rocMeanSens,outcome="mean")
  rocMDF$outcome <- as.character(rocMDF$outcome)
  rocMDF$sz = 1.1
  
  # prep dataframe for plotting
  rocDFAll <- rocsDFprecalc
  if (smooth) {
    rocDFAll <- rocDF
  }
  if (multiROC & oneAllROC) {
    rocDFAll <- rbind.data.frame(rocDFAll,rocMDF)
  } else if (multiROC & !oneAllROC) {
    rocDFAll <- rocMDF
  } else if (!multiROC & !oneAllROC) {
    print ("multiROC and oneAllROC must not all be FALSE!")
    stop()
  }
  # not doing confidence intervals
  if (!roc.conf) {
    rocDFAll <- rocDFAll[order(rocDFAll$sensitivity),]
    rocDFPlot <- rocDFAll
    rocDFPlot$ltype = rocDFPlot$sz
    if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot$ltype[rocDFPlot$sz == 1] <- 1.25
      rocDFPlot$ltype[rocDFPlot$sz == 1.25] <- 1
    } 
    if (!dottedNonPriority) {
      rocDFPlot$ltype <- 1
    }
    gg <- ggplot(data=rocDFAll,aes(x=sensitivity,y=specificity,col=outcome)) +
      scale_size_continuous(range = c(1.25, 2)) + 
      style_roc() +
      scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
      scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
      geom_abline(slope = 1,intercept = 1,col="grey") +
      geom_line(aes(size=as.numeric(sz),linetype=as.factor(sz))) +
      guides(size = FALSE,linetype=FALSE) + #theme(legend.position="bottom") +
      ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
      labs(color='Dataset') + guides(colour = guide_legend(override.aes = list(size=3)))
  }
  # doing confidence intervals
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = confBoots
    meanConf <- F
    for (i in seq(1,ncol(pred))) {
      if (smooth) {
      sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities = c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                       conf.level=0.66,boot.n = boots,progress = "none")
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep)),
                         conf.level=0.66,boot.n = boots,progress = "none")
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- as.numeric(rownames(rocConf))
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf$outcome <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }

      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
    }
    meanConfSpLow <- apply(meanConf[,grep('sp.low',colnames(meanConf))],MARGIN = 1,mean)
    meanConfSpHigh <- apply(meanConf[,grep('sp.high',colnames(meanConf))],MARGIN = 1,mean)
    meanConfMedian <- apply(meanConf[,grep('specificity',colnames(meanConf))],MARGIN = 1,mean)
    # prep for plotting
    if (oneAllROC & multiROC & !multiROCPriority) {
      rocMDF$sp.low = 0.0
      rocMDF$sp.high = 0.0
      rocDFPlot <- rbind.data.frame(rocDFConf,rocMDF)
    } else if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot <- rocDFConf
      rocDFPlot$sp.low = 0.0
      rocDFPlot$sp.high = 0.0
      rocDFPlotTmp <- cbind.data.frame(meanConfSpLow,meanConfMedian,meanConfSpHigh,rocConf$sensitivity)
      colnames(rocDFPlotTmp) <- c("sp.low","specificity","sp.high","sensitivity")
      rocDFPlotTmp$outcome = "multi.class"
      rocDFPlotTmp$sz = 1.25
      rocDFPlot <- rbind(rocDFPlotTmp,rocDFPlot)
    } else if (oneAllROC & !multiROC) {
      rocDFPlot <- rocDFConf
    } else if (!oneAllROC & multiROC) {
        rocDFPlot <- cbind.data.frame(meanConfSpLow,meanConfMedian,meanConfSpHigh,rocConf$sensitivity)
        colnames(rocDFPlot) <- c("sp.low","specificity","sp.high","sensitivity")
        rocDFPlot$outcome = "multi.class"
        rocDFPlot$sz = 1.25
    }
    rocDFPlot$ltype = rocDFPlot$sz
    if (oneAllROC & multiROC & multiROCPriority) {
      rocDFPlot$ltype[rocDFPlot$sz == 1] <- 1.25
      rocDFPlot$ltype[rocDFPlot$sz == 1.25] <- 1
    } 
    if (!dottedNonPriority) {
      rocDFPlot$ltype <- 1
    }
    # add 0s
    for (oc in unique(rocDFPlot$outcome)) {
      szt = rocDFPlot[rocDFPlot$outcome==oc,]$sz[1]
      ltt = rocDFPlot[rocDFPlot$outcome==oc,]$ltype[1]
      rocDFPlot <- rbind.data.frame(rocDFPlot,
                                    data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,outcome=oc,sz=szt,ltype=ltt) )
    }
    rocDFPlot <- rocDFPlot[order(rocDFPlot$specificity),]
    gg <- ggplot(rocDFPlot,aes(x=sensitivity,y=specificity,col=outcome)) +
      geom_ribbon(data=rocDFPlot,aes(ymin=sp.low,ymax=sp.high,fill=outcome),alpha=0.2,colour=NA) + 
      scale_size_continuous(range = c(1.25, 2)) + 
      style_roc() + geom_abline(slope = 1,intercept = 1,col="grey") + 
      scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
      scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
      guides(size = FALSE,linetype=FALSE) + #theme(legend.position="bottom") +
      ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
      labs(color='Dataset') + guides (fill = FALSE) + guides(colour = guide_legend(override.aes = list(size=3)))
    gg <- gg + geom_line(aes(size=as.numeric(sz),linetype=as.factor(ltype)))
  }
  
  # annotation
  # ====================
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF$DataSet )) {
      if (u != "multi.class") {
        if (class(tarDF) != "data.frame") {
          tarDF <- data.frame(DataSet=u,
                              Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                              Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                              F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                              B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value
                              )
        } else {
          tarDF <- rbind.data.frame(tarDF,data.frame(DataSet=u,
                                                     Sensitivity=predDF[predDF$DataSet==u & predDF$Metric=="Sensitivity",]$Value,
                                                     Specificity=predDF[predDF$DataSet==u & predDF$Metric=="Specificity",]$Value,
                                                     F1=predDF[predDF$DataSet==u & predDF$Metric=="F1",]$Value,
                                                     B.Acc=predDF[predDF$DataSet==u & predDF$Metric=="B.Acc",]$Value))
        }
      }
    }
    gg <- gg + geom_point(data=tarDF,aes(x=Specificity,y=Sensitivity,col=DataSet),size=4,shape=4,stroke=2)
  }
  
  # separate AUCses
  mvDown = 0
  if (mcROCmetric) {
    # multi-class AUC
    gg <- gg + annotate("text", x = .5, y = .11, label = paste("multi.class AUC = ",round(mcAUC,2),sep=""),hjust = 0)
    mvDown = mvDown + 1
  }
  if (mcAcc) {
    gg <- gg + annotate("text", x = .5, y = .11-0.07*mvDown, label = paste("multi.class Acc = ",round(conf$overall[["Accuracy"]],2),sep=""),hjust = 0)
    mvDown = mvDown + 1
  }
  if (mcKappa) {
    gg <- gg + annotate("text", x = .5, y = .11-0.07*mvDown, label = paste("multi.class Kappa = ",round(conf$overall[["Kappa"]],2),sep=""),hjust = 0)
  }
  if (oneAllROC) {
    posit = 0
    for (i in order(names(rocs),decreasing = T)) {
      posit = posit + 1
      cauc <- ci.auc(rocs[[i]],conf.level=0.66,progress = "none")
      caucsd <-  round(cauc[2]-cauc[1],2)
      if (roc.conf) {
        #lbl = paste(names(rocs)[[i]]," AUC = ",round(auc(rocs[[i]]),2),", sd = ",caucsd,sep="")
        lbl = paste(names(rocs)[[i]],": AUC = ",round(auc(rocs[[i]]),2),sep="")
      } else {
        lbl = paste(names(rocs)[[i]],": AUC = ",round(auc(rocs[[i]]),2),sep="")
      }
      if (target) {
        lbl = paste(lbl,
                    "; B.Acc = ",round(tarDF[tarDF$DataSet==names(rocs)[i], ]$B.Acc,2),
                    "; F1 = ",round(tarDF[tarDF$DataSet==names(rocs)[i], ]$F1,2),sep="")
      }
      gg <- gg + annotate("text", x = .5, y = .11+0.07*posit, label = lbl,hjust = 0)
    }
  }
  
  list(gg,gg2,predDF)
}


# normalisation function (asinsqrt, for use on metagenomes and/or pathways)
asinSqrtNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  mN[,rowsToNor] <- asin(sqrt(mN[,rowsToNor]))
  mN
}
# normalisation function (divide by max & center, for use on metagenomes and/or pathways)
divideByMaxCenterNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  for (r in rowsToNor) {
    mN[[r]] <- mN[[r]]-mean(mN[[r]])
    mN[[r]] <- mN[[r]]/max(abs(mN[[r]]))
  }
  mN
}
# normalisation function (divide by stdev & center, for use on metagenomes and/or pathways)
divideBySdevCenterNormalise <- function(mN,norTaxa=T,norPWY=T) {
  rowsToNor <- c()
  if (norTaxa) {rowsToNor <- c(rowsToNor,grep('__',colnames(mN)))}
  if (norPWY) {rowsToNor <- c(rowsToNor,grep('PWY',colnames(mN)))}
  for (r in rowsToNor) {
    mN[[r]] <- mN[[r]]-mean(mN[[r]])
    mN[[r]] <- mN[[r]]/sd(abs(mN[[r]]))
  }
  mN
}

# do linear correction for phenotypes
linearCorrectMGPwy <- function (iData,corrNames,corrMG=TRUE,corrPWY=TRUE,corrVFDB=T,corrCARD=T,
                                correctZeros=F,removeCorCol=T,negIsZero=F,verbose=T,debug=F) {
  
  testNames <- make.names(colnames(iData))
  origNames <- colnames(iData)
  if (sum(testNames == colnames(iData)) != ncol(iData) ) {
    print ('WARNING: columns names are not R-happy, consider using make.names(colnames(input)) !!')
    print ('WARNING:  >> will attempt fix, code WILL DIE if corrNames are not R-compatible! ')
    colnames(iData) <- testNames
  }
  
  # determine which columns to correct
  #print(sum(iData==0))
  tt <- iData
  corrCols <- c()
  if (corrMG) {
    corrCols <- c(corrCols,grep("^[dgtspcfko]__",colnames(iData)))
    if (verbose) {
      print(paste0(' > correcting ',length(grep("^[dgtspcfko]__",colnames(iData))),' Taxa!' ))
    }
  }
  if (corrPWY) {
    corrCols <- c(corrCols,grep("PWY",colnames(iData)))
    if (verbose) {
      print(paste0(' > correcting ',length(grep("PWY",colnames(iData))),' PWYs!' ))
    }
  }
  if (corrVFDB) {
    corrCols <- c(corrCols,grep("^VF",colnames(iData)))
    if (verbose) {
      print(paste0(' > correcting ',length(grep("^VF",colnames(iData))),' VFs!' ))
    }
  }
  if (corrCARD) {
    corrCols <- c(corrCols,grep("CARD_",colnames(iData)))
    if (verbose) {
      print(paste0(' > correcting ',length(grep("CARD_",colnames(iData))),' CARDs!' ))
    }
  }
  # iterate over columns, correct them for selected variables (entered as corrNames)
  for (c in corrCols) {
    # build formula
    if (debug) {
      print(paste0('  >> correcting ',colnames(iData)[c]))
    }
    if (!correctZeros) {
      corrRows <- rownames(iData)[iData[[c]] != 0.0]
    } else {
      corrRows <- rownames(iData)
    }
    if (length(corrRows) == 0) {
      print(paste0(' > WARNING: ',colnames(iData)[[c]],' is all zero, will not be corrected!'))
    } else {
      toCor <- iData[corrRows,c(colnames(iData)[c], corrNames)]
      # check for factors having only 1 level:
      corCorrNames = c()
      for (cn in corrNames) {
        if (is.factor(toCor[[cn]])) {
          #print (cn)
          doDrop = F
          zer = c()
          for (lvl in levels(toCor[[cn]])) {
            if (sum(toCor[[cn]]==lvl) == 0) {
              zer = c(zer,F)
            } else {zer = c(zer,T)}
          }
          #print (zer)
          if (sum(zer) >= 2) {
            corCorrNames = c(corCorrNames,cn)
          } else {
            print(paste('WARNING:',colnames(iData)[c],': variable',cn,'has less then 2 levels! Dropping it from linear model!'))
          }
        } else {
          corCorrNames <- c(corCorrNames,cn)
        }
      }
      #print (corCorrNames)
      if (is.null(corCorrNames)) {
        print(paste('WARNING:',colnames(iData)[c],': Will not be corrected - all variables dropped from model'))
      } else {
        frm <- reformulate(termlabels = corCorrNames,response=colnames(iData)[c])
        #print (paste0(' > ',frm))
        # do linear model
        m <- lm(data = toCor,frm)
        #print(m)
        # correct for linear model
        corTo <- m$coefficients[["(Intercept)"]] + resid(m)
        if (negIsZero) {
          corTo[corTo < 0] <- 0
        }
        iData[corrRows,colnames(iData)[c]] <- corTo
      }
    }
  }
  # get rid of columns used in model
  colnames(iData) <- origNames
  if (removeCorCol) {
    for (c in corrNames) {
      iData[[c]] <- NULL
    } 
  }
  iData
}

# ===============================================
# prepares dataset by refactoring as necessary
# and taking only input covariates
# ===============================================
prepData <- function(inData,covars) {
  # grab covariates
  mdl <- inData[,covars]
  # omit NAs
  mdl <- na.omit(mdl)
  # refactor factors
  fak <- sapply(mdl, is.factor)
  c = 0
  for (i in fak) {
    c = c + 1
    if (i) {
      mdl[,c] <- as.factor(as.character(mdl[,c]))
    }
  }
  fak <- sapply(mdl, is.character)
  c = 0
  for (i in fak) {
    c = c + 1
    if (i) {
      mdl[,c] <- as.factor(as.character(mdl[,c]))
    }
  }
  mdl
}

# >>> MODEL PREP
prepModel <- function(dFrame,mdl) {
  if (mdl=="MG_PWY") {
    dFrame <- purgeMGNames(dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )])
  } else if (mdl=="P") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="P_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("P"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="F") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("F"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="F_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("F"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="G") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="G_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="S") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="S_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="GS") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S","G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="GS_PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('k__',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(filterMetaGenomeDF(dFrame,keepLevels=c("S","G"),presPerc = -1,minMRelAb = -1,minMedRelAb = -1))
  } else if (mdl=="PWY") {
    dFrame <- dFrame[,c(grep('Diagnosis',colnames(dFrame)),grep('PWY',colnames(dFrame)) )]
    dFrame <- purgeMGNames(dFrame)
  } 
  colnames(dFrame) <- make.names(colnames(dFrame),unique=T)
  dFrame
}

#' =================================================================================
#' dataModelOptimiseRFE function
#' ===============================================================================
#' 
#'\code{dataModelOptimiseRFE} optimises data model using Recursive Feature Elimination algorithm
#'
#'@param responseVar : what is to be predicted [def: Diagnosis]
#'@param trP : percentage of data to put into in training sets [def: 0.75]
#'@param positive : what value is considered "positive" [def: 'IBD']
#'@param rfeMethod : which method to use for optimisation [def: glm]
#'@param testMethod : which method to use for testing [def: glm]
#'@param xvnr : # of x-validations [def=5]
#'@param xvreps : how many times to repeat x-validation [def=1]
#'@param parallel : if parallel processing allowed [def=T]
#'@param verb : should it prints out debug info [def=T]
#'@param optimiseMetric : what to optimise (Kappa/ROC) [def='Kappa']
#'@return list of 1) max accuracy vars, 2) 2% tolerance vars, 3) 5% tolerance vars, 4) test set ROCs, 5) x-validation ROCs, 6) RFE plot
#'}
# =================================================================================
dataModelOptimiseRFE <- function(dModel,dModelName="",responseVar="Diagnosis",trP=0.75,positive="IBD",rfeMethod="glm",
                                 testMethod="glm",xvnr=5,xvreps=1,rfeReps=3,
                                 parallel=T,verb=T,szs=F,trainMethod="repeatedcv",optimiseMetric="Kappa",
                                 saveRFE=T,saveRFEpath='rfe_model.RDS',doComparisons=T) {
  
  #trC <- trainControl(method=rfeMethod,number=tBut,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = parallel)
  inTrain <- createDataPartition(dModel[[responseVar]],p=trP,list=F)
  dModel[[responseVar]] <- as.factor(as.character(dModel[[responseVar]]))
  trSet <- dModel[inTrain,]
  testSet <- dModel[-inTrain,]
  
  if (optimiseMetric == "ROC") {
    caretFuncs$summary <- twoClassSummary
    trC <- trainControl(method="repeatedcv",number=xvnr,
                        repeats = xvreps,savePredictions = T,
                        classProbs = T,allowParallel = parallel,
                        summaryFunction = twoClassSummary)
    rfeCtrl <- rfeControl(method = "repeatedcv", repeats = rfeReps,functions = caretFuncs,
                          number=xvnr, verbose = T,allowParallel = parallel,returnResamp="final")
  } else {
    trC <- trainControl(method="repeatedcv",number=xvnr,
                        repeats = xvreps,savePredictions = T,
                        classProbs = T,allowParallel = parallel)
    
    # NOTE: caretFuncs have BUG and randomly crash !!!
    # rfeCtrl <- rfeControl(functions = caretFuncs, method = "repeatedcv", repeats = rfeReps,
    #                       number=xvnr, verbose = T,allowParallel = parallel)
    
    rfeCtrl <- rfeControl(functions = rfFuncs, method = "repeatedcv", repeats = rfeReps,
                          number=xvnr, verbose = T,allowParallel = parallel)
  }
  if (!szs) {
    if (ncol(trSet) <= 50) {szs=seq(1,ncol(trSet)-1,1)
    } else if (ncol(trSet) <= 105) {szs=c(seq(1,50,1),seq(50,ncol(trSet)-1,2))  
    } else if (ncol(trSet) <= 210) {szs=c(seq(1,50,1),seq(50,100,2),seq(105,ncol(trSet)-1,5))
    } else if (ncol(trSet) > 210) {szs=c(seq(1,50,1),seq(50,100,2),seq(105,200,5), seq(200,ncol(trSet)-1,10))  }
  }
  szs <- unique(szs)
  
  #debug
  # trSet2 <- trSet
  # nrRes <- grep(paste0('^',responseVar,'$'),colnames(trSet))
  # colnames(trSet2) <- paste0("VAR",seq(1,ncol(trSet2)))
  # colnames(trSet2)[nrRes] <- "RESPONSE"
  # rfeProfile <- rfe(x=trSet2[,-grep('^RESPONSE$',colnames(trSet2))], y=trSet2$RESPONSE, sizes=szs, rfeControl = rfeCtrl,
  #                   metric=optimiseMetric,method=rfeMethod,maximize = T,trControl = trC)
    
  #print(szs)
  print ('  >> doing RFE profile')
  rfeProfile <- rfe(x=trSet[,-grep(paste0('^',responseVar,'$'),colnames(trSet))], y=trSet[[responseVar]], sizes=szs, rfeControl = rfeCtrl,
                    metric=optimiseMetric,method=rfeMethod,maximize = T,trControl = trC)
  print ('    >>> DONE!')
  # save output
  varImpPlot <- cbind.data.frame(rownames(varImp(rfeProfile)),varImp(rfeProfile))
  colnames(varImpPlot) <- c("Variable","Importance")
  varImpPlot$Importance <- varImpPlot$Importance/max(varImpPlot$Importance)
  varImpPlot$Variable <- factor(varImpPlot$Variable, levels = varImpPlot$Variable[order(varImpPlot$Importance,decreasing = T)])
  rownames(varImpPlot) <- NULL
  
  # select best NR of vars within tol% margin of error, while keeping number as low as possible
  varNrMax <- pickSizeTolerance(rfeProfile$results, metric=optimiseMetric,maximize = T,tol=0.1)
  if (is.na(varNrMax)) {varNrMax = ncol(trSet)-1}
  varSelMax <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNrMax)
  varNr2 <- pickSizeTolerance(rfeProfile$results, metric=optimiseMetric,maximize = T,tol=2)
  if (is.na(varNr2)) {varNr2 = ncol(trSet)-1}
  varSel2 <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNr2)
  varNr5 <- pickSizeTolerance(rfeProfile$results, metric=optimiseMetric,maximize = T,tol=6)
  if (is.na(varNr5)) {varNr5 = ncol(trSet)-1}
  varSel5 <- caretFuncs$selectVar(y = rfeProfile$variables,size=varNr5)
  # build new model with these only, compare to original model
  if (optimiseMetric == "ROC") {
    trC <- trainControl(method=trainMethod,number=xvnr,repeats = xvreps,savePredictions = T,classProbs = T,allowParallel = parallel,
                        summaryFunction = twoClassSummary)
  } else {
    trC <- trainControl(method=trainMethod,number=xvnr,repeats = xvreps,savePredictions = T,classProbs = T,allowParallel = parallel)
  }
  trSetMax <- trSet[,c(responseVar,varSelMax)]
  trSetSel2 <- trSet[,c(responseVar,varSel2)]
  trSetSel5 <- trSet[,c(responseVar,varSel5)]
  
  trForm <- reformulate(response=responseVar,termlabels = '.')
  
  fitAll <- train(trForm, trSet,trControl = trC,method = testMethod,metric=optimiseMetric)
  fitMax <- train(trForm, trSetMax,trControl = trC,method = testMethod,metric=optimiseMetric)
  fitRfe2 <- train(trForm, trSetSel2,trControl = trC,method = testMethod,metric=optimiseMetric)
  fitRfe5 <- train(trForm, trSetSel5,trControl = trC,method = testMethod,metric=optimiseMetric)
  # compare
  mNames <- c(paste("All:",ncol(testSet)-1,sep=''),paste("M.max:",varNrMax,sep=''),
              paste("M.2:",varNr2,sep=''),paste("M.5:",varNr5,sep=''))
  
  if (length(unique(trSet[[responseVar]])) == 2 & doComparisons) {
    testComp <- compareMdlsDatasets(mdls = list(fitAll,fitMax,fitRfe2,fitRfe5),mdNames = mNames,dataSets = list(testSet),
                                    posClass = positive, tit = "",response = responseVar,annotS = 4,textSize = 15)
    
    cvComp <- compareModelsTrainingCV(fittedMdls = list(fitAll,fitMax,fitRfe2,fitRfe5),mtd = paste('x-val',testMethod), 
                                      modelNames = mNames,posClass = positive)
  }
  rfeplot <- ggplot(rfeProfile) + geom_vline(xintercept = varNrMax,linetype='longdash') + geom_vline(xintercept = varNr2,linetype='longdash')+
    geom_vline(xintercept = varNr5,linetype='longdash')+ggtitle(paste("RFE plot (",dModelName,' / ',rfeMethod,')',sep=""))
  
  
  if (saveRFE) {
    saveRDS(rfeProfile,file = saveRFEpath)
  }
  if (length(unique(trSet[[responseVar]])) == 2 & doComparisons) {
    return(list(varSelMax,varSel2,varSel5,testComp,cvComp,rfeplot,varImpPlot,c("Vars:Max","Vars:2%","Vars:5%","plot:testset","plot:X-val","RFE-plot","VarImp-DF")))
  } else {
    return(list(varSelMax,varSel2,varSel5,rfeplot,varImpPlot,c("Vars:Max","Vars:2%","Vars:5%","RFE-plot","VarImp-DF")))
  }
}

# ===============================================================================
# compares model(s) performance on one or more datasets
# ===============================================================================
compareMdlsDatasets <- function(mdls,dataSets,mdNames,posClass="IBS",roc.smooth=F,tit="",specSensAnnot=T,
                                response="Diagnosis",roc.conf = T,roc.conf.boot = 10,target = F,conf.lvl=0.66,
                                textSize=15,annotS=4,removeLegend=F) {
  
  nrModels <- length(mdls)
  nrDataSets <- length(dataSets)
  #print (nrModels)
  #print (nrDataSets)
  # list of ROC curves
  rocs = list()
  # data frame of predictions
  predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
  if (nrModels == 1 & nrDataSets == 1) {
    # ===================================================================================
    # ================== 1 MODEL, 1 DATASET =====================================
    # ===================================================================================
    leg = "Model"
    print (' -> 1 model, 1 dataset')
    trainedModels <- list()
    # generate results for each dataset
    testSet <- dataSets[[1]]
    fittedModel <- mdls[[1]]
    trainedModels[[1]] <- fittedModel
    pr <- predict(fittedModel,newdata = testSet,type="prob")
    namen = mdNames[1]
    if (roc.smooth) {
      rocs[[namen]] <- smooth(roc(testSet[[response]],pr[[posClass]],auc=T,percent=F))
    } else {
      rocs[[namen]] <- roc(testSet[[response]],pr[[posClass]],auc=T,percent=F)
    }
    pr2 <- predict(fittedModel,testSet)
    conf <- confusionMatrix(pr2,testSet[[response]])
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
    predDF <- rbind.data.frame(predDF,
                               data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
  } else if (nrModels == 1 & nrDataSets > 1) {
    # ===================================================================================
    # ================== 1 MODEL, MULTIPLE DATASETS =====================================
    # ===================================================================================
    print (' -> 1 model, multiple datasets')
    trainedModels <- list()
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    leg = "Dataset"
    for (ds in dataSets) {
      c = c + 1
      testSet <- dataSets[[c]]
      fittedModel <- mdls[[1]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet[[response]],pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet[[response]],pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
    
  } else if (nrModels > 1 & nrDataSets == 1) {
    # ===================================================================================
    # ================== MULTIPLE MODELs, 1 DATASET =====================================
    # ===================================================================================
    leg = "Model"
    print (' -> multiple models, 1 dataset')
    trainedModels <- list()
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("DataModel"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    for (fittedModel in mdls) {
      c = c + 1
      testSet <- dataSets[[1]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet[[response]],pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet[[response]],pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
    
  } else if (nrModels == nrDataSets & nrModels > 1) {
    print (' -> multiple dataset - model pairs')
    trainedModels <- list()
    leg = "Data.Model"
    c = 0
    # list of ROC curves
    rocs = list()
    # data frame of predictions
    predDF <- data.frame("Data.Model"=as.character(),"Metric"=as.character(),"Value"=as.numeric(),stringsAsFactors=F)
    # generate results for each dataset
    for (fittedModel in mdls) {
      c = c + 1
      testSet <- dataSets[[c]]
      trainedModels[[c]] <- fittedModel
      pr <- predict(fittedModel,newdata = testSet,type="prob")
      namen = mdNames[c]
      if (roc.smooth) {
        rocs[[namen]] <- smooth(roc(testSet[[response]],pr[[posClass]],auc=T,percent=F))
      } else {
        rocs[[namen]] <- roc(testSet[[response]],pr[[posClass]],auc=T,percent=F)
      }
      pr2 <- predict(fittedModel,testSet)
      conf <- confusionMatrix(pr2,testSet[[response]])
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Acc","Value"=conf$overall[["Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Kappa","Value"=conf$overall[["Kappa"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Sensitivity","Value"=conf$byClass[["Sensitivity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="Specificity","Value"=conf$byClass[["Specificity"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="PPV","Value"=conf$byClass[["Pos Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="NPV","Value"=conf$byClass[["Neg Pred Value"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="B.Acc","Value"=conf$byClass[["Balanced Accuracy"]],stringsAsFactors=F))
      predDF <- rbind.data.frame(predDF,
                                 data.frame("Data.Model"=namen,"Metric"="F1","Value"=conf$byClass[["F1"]],stringsAsFactors=F))
    }
  }
  
  # extract coordinates
  if (!roc.conf) {
    rocStep = 0.01
    senses <- c(seq(0,0.1,rocStep/10),seq(0.1,1,rocStep))
    rocsDFprecalc = NULL
    for (i in seq(1,length(rocs))) {
      if (class(rocsDFprecalc) != "data.frame") {
        if (!roc.smooth) {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        } else {
          rocsDFprecalc <- as.data.frame(t(coords(rocs[[i]],x=senses,ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalc) <- NULL
        }
        rocsDFprecalc[[leg]] <- names(rocs)[i]
        rocsDFprecalc$sz = 1
      } else {
        if (!roc.smooth) {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x="all",ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        } else {
          rocsDFprecalctmp <- as.data.frame(t(coords(rocs[[i]],x=senses,ret= c("specificity","sensitivity"))))
          rownames(rocsDFprecalctmp) <- NULL
        }
        rocsDFprecalctmp[[leg]] <- names(rocs)[i]
        rocsDFprecalctmp$sz = 1
        rocsDFprecalc <- rbind(rocsDFprecalc, rocsDFprecalctmp)
        rownames(rocsDFprecalc) <- NULL
      }
    }
  }
  
  if (roc.conf) {
    rocDFConf <- F
    rocStep = 0.01
    boots = roc.conf.boot
    meanConf <- F
    senses <- c(seq(0,0.1,rocStep/10),seq(0.1, 1, rocStep))
    for (i in seq(1,length(rocs))) {
      if (roc.smooth) {
        sens.ci <- ci.sp(smooth(rocs[[i]]), sensitivities = senses,
                         conf.level=conf.lvl,boot.n = boots,progress = "none")
      } else {
        sens.ci <- ci.sp((rocs[[i]]), sensitivities=senses,
                         conf.level=conf.lvl,boot.n = boots,progress = "none")
      }
      rocConf <- as.data.frame(sens.ci)
      rocConf$sensitivity <- senses
      colnames(rocConf) <- c("sp.low","specificity","sp.high","sensitivity")
      rocConf[[leg]] <- names(rocs)[i]
      rocConf$sz = 1.0
      rownames(rocConf) <- NULL
      if (class(meanConf) != "data.frame") {
        meanConf <- rocConf[,c(4,1,2,3)]
      } else {
        meanConf <- cbind(meanConf,rocConf[,c(1,2,3)])
      }
      if (class(rocDFConf) != "data.frame") {
        rocDFConf <- rocConf
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,rocConf)
      }
      if (leg=="Model") {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Model=names(rocs)[i],sz=1))
      } else if (leg=="Dataset") {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Dataset=names(rocs)[i],sz=1))
      } else {
        rocDFConf <- rbind.data.frame(rocDFConf,data.frame(sp.low=0,specificity=0,sp.high=0,sensitivity=1,Data.Model=names(rocs)[i],sz=1))
      }
    }
    rocsDFprecalc <- rocDFConf
  }
  
  if (roc.conf) {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$specificity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col=leg)) +
      geom_ribbon(data=rocsDFprecalc,aes_string(ymin="sp.low",ymax="sp.high",fill=leg),alpha=0.2,colour=NA) + geom_line(size=1.25)
  } else {
    rocsDFprecalc <- rocsDFprecalc[order(rocsDFprecalc$sensitivity,decreasing = F),]
    g <- ggplot(rocsDFprecalc,aes_string(y="specificity",x="sensitivity",col=leg)) + geom_line(size=1.25) #pROC::ggroc(rocs) + geom_line(size=1.25)
  }
  g <- g + 
    annotate("segment", x = 1, xend = 0, y = 0, yend = 1,colour = "gray",size=1) +
    ggtitle(paste(tit,"ROC") ) + xlab("Specificity") + ylab("Sensitivity") +
    style_roc() +  
    scale_y_continuous(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0, 1, 0.1)) +
    scale_x_reverse(minor_breaks = c(seq(0,0.8,0.1),seq(0.75,0.9,0.05),seq(0.91,1,0.01)), breaks = seq(0,1,0.1)) +
    labs(color=leg) + theme(text = element_text(size = textSize))
  # annotate w auc
  c = 0.0
  for (i in order(mdNames,decreasing = T)) {
    c = c + 1.0
    if (roc.conf) {
      rocAUCconf <- ci.auc(rocs[[ mdNames[i] ]],conf.level=conf.lvl,boot.n = boots,progress="none")
      rocAUCm <- rocAUCconf[[2]]
      rocAUCdelta <- rocAUCconf[[2]] - rocAUCconf[[1]]
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(mdNames[i],"AUC =",round(rocAUCm,2),"\u00b1",round(rocAUCdelta,2)),hjust = 0, size=annotS)
    } else {
      g <- g + annotate("text", x = .5, y = 0.00+0.06*c, label = paste(mdNames[i],"AUC =",round(rocs[[mdNames[i] ]]$auc[1],2)),hjust = 0,size=annotS)
    }
    #annotate("text", x = .5, y = .30, label = paste("AUC =",round(rTst$auc[1],2) ),hjust = 0) +
    #
  }
  # targets
  tarDF <- F  
  if (target) {
    for (u in unique(predDF[[leg]] )) {
      if (class(tarDF) != "data.frame") {
        tarDF <- data.frame(DataSet=u,
                            Sensitivity=predDF[predDF[[leg]]==u & predDF$Metric=="Sensitivity",]$Value,
                            Specificity=predDF[predDF[[leg]]==u & predDF$Metric=="Specificity",]$Value,
                            F1=predDF[predDF[[leg]]==u & predDF$Metric=="F1",]$Value,
                            B.Acc=predDF[predDF[[leg]]==u & predDF$Metric=="B.Acc",]$Value)
      }
      g <- g + geom_point(data=tarDF,aes_string(x="Specificity",y="Sensitivity",col="DataSet"),size=4,shape=4,stroke=2)
    }
  }
  # legend removal
  if (removeLegend) {
    g <- g + theme(legend.position = "none")
  }

  
  # do metrics plot
  unik <- max(length(dataSets),length(mdls))
  gg2 <- ggplot(data=predDF,aes_string(col="Metric",y="Value",x=leg,shape="Metric")) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4)) + theme_minimal() + ylim(0.0,1) + 
    scale_shape_manual(values = c(0,1,2,3,4,5,6,7,9,10,12,13)) +  ggtitle(paste(tit,"Prediction metrics") ) +
    geom_vline(xintercept =  c(1:length(mdNames)) ,size=max(35,130-unik*20),alpha=0.05) +
    geom_point(stroke=2,size=3,position=position_dodge(width=0.4))
  list(g,gg2,predDF)
}


# =====================================================================================
# Predictive Modelling wrapper script
# - Does:
#   > separation into training & test sets
#      > saves them
#   > data pre-processing 
#   > learning curves
#   > feature optimisation (RFE)
#   > optimised learning curves
#   > raw model training
#   > optimised model training
#   > comparison of raw vs optimised models
# =====================================================================================
doMLModelling <- function(outFolderName,
                        allDM,         # list of data models
                        allDMn,        # list of names of data models 
                        responseVar,   # what is response variable
                        posC,          # positive class
                        inTrainP=0.75, # how much goes into training set
                        keepTrSplit=T, # keep training / test split the same between datasets
                        doDataPrep=T,     # do data preparation (scaling, centering...)
                        doDataPrepPreselection = F, # do pre-selection?
                        doDataPrepPreselectP = 0.2, # pre-selection p-value of independance_test 
                        doDataPrepPreselectMaxFeatures = 50, #pre-selection max # of features
                        doDataPrepPreselectMinFeatures = 5,
                        dataPrepMtd = c("scale","center"),  # how to pre-process data
                        doDataSep=T,   
                        doLC = T,      # do learning curves
                        doRFE = T,     # do recursive feature optimisation
                        doRFEpreselection = T, # do pre-selection?
                        doRFEpreselectP = 0.1, # pre-selection p-value of independance_test 
                        doRFEpreselectMaxFeatures = 50, #pre-selection max # of features
                        doRFEpreselectMinFeatures = 5,
                        doRfeComparisons = T,
                        doOptLC = T,   # do optimised learning curves
                        doRawMdls = T,  # do raw models
                        doOptMdls = T,  # do optimised models
                        doMdlComparison = T, # do model comparison (raw vs optimised)
                        doOptMax = T,  # do MAX optimisation
                        doOpt2 = T,    # do V2 optimisation
                        tButLC=10,     # bootstraps of learning curves
                        tButXV=10,       # training xvalidations (or bootstraps) for all steps
                        tRep=3,        # training repeats (for all steps)
                        #lcB=5,         # learning curves bootstraps
                        rfeR=5,        # repeats of RFE
                        saveMdls = T,   # save fitted models
                        smoothROCs = F, # if T, make smooth ROC curves <warning: buggy!>
                        mtds = c("glm","svmRadial"),
                        optimiseMetric="Kappa",  # what to optimise (usually Kappa, TODO: ROC)
                        allowParallel=T,
                        excludeVars = NULL # which variables to exclude FROM TESTING (they are still used for TRAINING!)
) 
{
  # root output folder, create/clean as necessary
  # ====================================================
  if (!dir.exists(outFolderName)) {dir.create(outFolderName)}
  
  # separate data into training and testets 
  # ====================================================
  if (doDataSep) {
    print(' > SEPARATING DATA INTO TRAINING & TEST')
    if (!dir.exists(paste0(outFolderName,'/inputData'))) {dir.create(paste0(outFolderName,'/inputData'))}
    if (keepTrSplit) {
      inTrain <- createDataPartition(y=allDM[[1]][[responseVar]],p=inTrainP,list=F)
    }
    for (c in c(1:length(allDM))) {
      if (!keepTrSplit) {
        inTrain <- createDataPartition(y=allDM[[c]][[responseVar]],p=inTrainP,list=F)
      }
      trainSet <-allDM[[c]][inTrain,]
      testSet <- allDM[[c]][-inTrain,]
      mN <- allDMn[[c]]
      write.table(trainSet,file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'), sep=',',row.names = F)
      if (inTrainP != 1) {
        write.table(testSet, file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),sep=',',row.names = F)
      }
    }
  }

  # preprocess data
  # ====================================================
  allDMtr <- list()
  if (doDataPrep) {
    print(' > PRE-PROCESSING DATA')
    for (c in c(1:length(allDM))) {
      mN <- allDMn[[c]]
      allDMtr[[c]] <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')
      prepModel <- preProcess(allDMtr[[c]],method=dataPrepMtd)
      #print(prepModel)
      preppedData <- predict(prepModel, allDMtr[[c]])
      write.table(preppedData, file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),sep=',',row.names = F)
      saveRDS(prepModel,file=paste0(outFolderName,'/inputData/',mN,'_trSet_prepMdl.RData'))
    }
  } else {
    for (c in c(1:length(allDM))) {
      mN <- allDMn[[c]]
      allDMtr[[c]] <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')
      write.table(allDMtr[[c]], file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),sep=',',row.names = F)
    }
  }
  # do pre-selection?
  if (doDataPrepPreselection) {
    print(' > PRE-SELECTING FEATURES')
    if (doDataPrep) {
      inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')  
    } else {
      inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')  
    }
    if (ncol(inData) > 2) {
      print (paste("    >> DOING Pre-processsing PRE-SELECTION"))
      toTest <- colnames(inData)[colnames(inData) != responseVar]
      testRes <- data.frame()
      for (t in toTest) {
        rVC <- unique(inData[[responseVar]])
        if ( class(inData[[t]]) == "factor" | class(inData[[t]]) == "character" )
        {
          tst <- chisq.test(table(inData[[t]],inData[[responseVar]]))
        } else {
          t1 <- inData[[t]][inData[[responseVar]]==rVC[1]]
          t2 <- inData[[t]][inData[[responseVar]]==rVC[2]]
          tst <- wilcox.test(t1,t2)
        }
        testRes <- rbind.data.frame(testRes,data.frame(feature=t,pv=tst$p.value))
      }
      testRes$pv[is.na(testRes$pv) | is.nan(testRes$pv)] <- 1.0
      testRes <- testRes[order(testRes$pv),]
      touse <- as.character(testRes$feature)
      touse <- touse[testRes$pv <= doDataPrepPreselectP]
      if (length(touse) < doDataPrepPreselectMinFeatures) {
        touse <- as.character(testRes$feature[1:min(ncol(testRes),doDataPrepPreselectMinFeatures)])
      }
      if (length(touse) > doDataPrepPreselectMaxFeatures) {
        touse <- as.character(touse[1:doDataPrepPreselectMaxFeatures])
      }
      inData <- inData[,colnames(inData) %in% c(responseVar,touse)]
      print (paste('   >> PRE-SELECTION DONE, selected',length(touse),'FEATURES!'))
    }
    write.table(c(responseVar,touse), file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep_features.csv'),sep=',',row.names = F)
    write.table(inData, file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),sep=',',row.names = F)
  }
    
  # make learning curves
  # =============================================
  if (doLC) {
    print (" > MAKING LEARNING CURVES:")
    if (!dir.exists(paste(outFolderName,'/lcurves',sep=''))) {dir.create(paste(outFolderName,'/lcurves',sep=''))}
    #trB = tButLC; trRep = tRep # training parameters for each step of learning curve
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        print (paste('prepping learning curve for',mN,'; method:',mtd))
        if (doDataPrep) {
          inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')  
        } else {
          inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')  
        }
        lCurve <- prepLearningCurve(dataIn = inData,mdl = mtd,lcBoots = tButLC,trBoot = 3,trNumber = 3,trainType="cv",trSet = 0.75,
                                    saveVarImp = F,responseVar = responseVar,posClass = posC,minSam = 100,samStep = 2,
                                    pProc = T,excludeVars=excludeVars)
        write.table(lCurve,file=paste0(outFolderName,"/lcurves/",mN,"_lc_",mtd,".csv"),sep=',',row.names = F)
        lp <- plotLearningCurves(lCurve,tit=paste('Learning curve (',mN,' [',ncol(inData)-1,'f], ',mtd,')',sep=''),metrics = c("Kappa"))
        ggsave(plot = lp,file=paste0(outFolderName,"/lcurves/",mN,"_lc_",mtd,"_k.png"),width = 8,height = 6)
        lp <- plotLearningCurves(lCurve,tit=paste('Learning curve (',mN,' [',ncol(inData)-1,'f], ',mtd,')',sep=''),metrics = c("Sensitivity","Specificity","ACC"))
        ggsave(plot = lp,file=paste0(outFolderName,"/lcurves/",mN,"_lc_",mtd,"_acc.png"),width = 8,height = 6)
        lp <- plotLearningCurves(lCurve,tit=paste('Learning curve (',mN,' [',ncol(inData)-1,'f], ',mtd,')',sep=''),metrics = c("NPV","PPV","AUC"))
        ggsave(plot = lp,file=paste0(outFolderName,"/lcurves/",mN,"_lc_",mtd,"_npv.png"),width = 8,height = 6)
      }
    }
  }
  
  # run RFE to determine most relevant genera, species & pathways
  # ===================================================================
  # general training parameters
  # -> optimise model using RFE and compare to full model & save results
  # ===================================================================
  if (doRFE) {
    print (" > DOING RFE RUNS:")
    if (!dir.exists(paste(outFolderName,'/opt_RFE',sep=''))) {dir.create(paste(outFolderName,'/opt_RFE',sep=''))}
    # RFE parameters
    xvn = rfeR; xvr = tButXV
    # do RFE
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        if (doDataPrep) {
          inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')
        } else {
          inData <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')
        }
        if (ncol(inData) > 2) {
          print (paste("  >> OPTIMISING (RFE)",allDMn[[c]],"with",mtd))
          if (doRFEpreselection) { # do pre-selection?
            print (paste("    >> DOING PRE-SELECTION"))
            toTest <- colnames(inData)[colnames(inData) != responseVar]
            testRes <- data.frame()
            rVC <- unique(inData[[responseVar]])
            combos <- combn(rVC,m = 2)
            # do all pair-wise combos, get p-values for them
            for (cmb in c(1:ncol(combos))) {
              v1 <- combos[,cmb][1]
              v2 <- combos[,cmb][2]
              #print(paste(v1,',',v2))
              for (t in toTest) {
                rVC <- unique(inData[[responseVar]])
                if ( class(inData[[t]]) == "factor" | class(inData[[t]]) == "character" )
                {
                  tst <- chisq.test(table(inData[[v1]],inData[[v2]]))
                } else {
                  t1 <- inData[[t]][inData[[responseVar]]==v1]
                  t2 <- inData[[t]][inData[[responseVar]]==v2]
                  tst <- wilcox.test(t1,t2)
                }
                testRes <- rbind.data.frame(testRes,data.frame(feature=t,v1=v1,v2=v2,pv=tst$p.value))
              }
            }
            testRes$pv[is.na(testRes$pv) | is.nan(testRes$pv)] <- 1.0
            testRes <- testRes[order(testRes$pv),]
            #testRes$touse <- 1
            #testRes <- 0
            touse <- as.character(testRes$feature)
            touse <- touse[testRes$pv <= doRFEpreselectP]
            if (length(touse) < doRFEpreselectMinFeatures) {
              touse <- as.character(testRes$feature[1:min(ncol(testRes),doRFEpreselectMinFeatures)])
            }
            if (length(touse) > doRFEpreselectMaxFeatures) {
              touse <- as.character(touse[1:doRFEpreselectMaxFeatures])
            }
            inData <- inData[,colnames(inData) %in% c(responseVar,touse)]
            print (paste('   >> PRE-SELECTION DONE, selected',length(touse),'FEATURES!'))
            write.table(testRes,paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_preselection.csv"),row.names = F,col.names = F)
          }
          
          # debug
          #dModel=allDM[[c]];xvnr=xvn;xvreps = xvr; rfeMethod = mtd; parallel = T;testMethod = mtd; dModelName=allDMn[[c]];positive = posC
          rfeRes <- dataModelOptimiseRFE(dModel=inData,xvnr=xvn,xvreps = xvr,rfeReps = rfeR,rfeMethod = mtd,optimiseMetric = optimiseMetric,
                                         testMethod = mtd, dModelName=allDMn[[c]],positive = posC,responseVar = responseVar,
                                         parallel = allowParallel,saveRFE = T,saveRFEpath = paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_RFE.RDS"),
                                         doComparisons = doRfeComparisons)
          print('RFE done!')
          # return plots (2 classes & doRfeComparisons = T)
          if (length(unique(inData[[responseVar]])) == 2 &  doRfeComparisons ) {
            write.table(rfeRes[[1]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_vmax.csv"),row.names = F,col.names = F)
            write.table(rfeRes[[2]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v2.csv"),row.names = F,col.names = F)
            write.table(rfeRes[[3]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v5.csv"),row.names = F,col.names = F)
            write.table(rfeRes[[4]][[3]],  paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_opt_results.csv"),row.names = F)
            ggsave(plot = rfeRes[[4]][[1]],paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_roc_test.png"),width = 8,height = 6)
            ggsave(plot = rfeRes[[5]]     ,paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_roc_xval.png"),width = 8,height = 6)
            ggsave(plot = rfeRes[[6]]     ,paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_opt.png"),width = 8,height = 6)
            write.table(rfeRes[[7]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_opt_varImp.csv"),row.names = F)
            vip <- ggplot(rfeRes[[7]], aes(x=Variable,y=Importance,col=Variable,fill=Variable)) + geom_col() + ylab('Relative Importance') +
              theme(axis.text.x=element_blank())
            ggsave(plot = vip     ,paste(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_varImp.png",sep = ''),width = 8,height = 6)
          } else {
          # return tables
            write.table(rfeRes[[1]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_vmax.csv"),row.names = F,col.names = F)
            write.table(rfeRes[[2]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v2.csv"),row.names = F,col.names = F)
            write.table(rfeRes[[3]],       paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v5.csv"),row.names = F,col.names = F)
            ggsave(plot = rfeRes[[4]]     ,paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_opt.png"),width = 8,height = 6)
            
          }
        } else {
          write.table(colnames(inData)[colnames(inData) != responseVar],paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_vmax.csv"),row.names = F,col.names = F)
          write.table(colnames(inData)[colnames(inData) != responseVar],paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v2.csv"),row.names = F,col.names = F)
          write.table(colnames(inData)[colnames(inData) != responseVar],paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",allDMn[[c]],"_vars_v5.csv"),row.names = F,col.names = F)
        }
      }
    }
    print ("RFE RUNS DONE!")
  }
  
  # remake models using RFE-identified genera, species & pathways
  # make new learning curves with those
  # ===================================================================
  if (doOptLC) {
    if (!dir.exists(paste(outFolderName,'/opt_lcurves',sep=''))) {dir.create(paste(outFolderName,'/opt_lcurves',sep=''))}
    defMethod = "glm"
    #trB = tButLC; trRep = tRep # training parameters for each step of learning curve
    print ("BUILDING OPTIMISED LEARNING CURVES")
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        print(paste0(' > dataset = ',mN,", method = ",mtd))
        # optimisation level (2% unless if it taxa & pathways, then 5%)
        if (doOpt2) {tarOpt = "v2"} else {tarOpt="vmax"}
        if (doDataPrep) {
          dm <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')
        } else {
          dm <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')
        }
        inFileOpt = paste(outFolderName,"/opt_RFE/rfe_",mtd,"_",mN,"_vars_",tarOpt,".csv",sep = '')
        #print (inFileOpt)
        inFileDef = paste(outFolderName,"/opt_RFE/rfe_",defMethod,"_",mN,"_vars_",tarOpt,".csv",sep = '')
        #print (inFileDef)
        # select optimal RFE, if it doesn't exist select default one; if that ones doesn't exist, don't do it    
        #print (paste('loading RFE data for',allDMn[[c]],'; method:',mtd))
        inFile = inFileOpt
        if (file.exists(inFileOpt)) {
          print (paste(' -> using',inFileOpt,'for optimisation!'))
        } else if (file.exists(inFileDef)) {
          print (paste(' -> using',inFileDef,'for optimisation!'))
          inFile = inFileDef
        }
        if (file.exists(inFile)) {
          optVars <- read.table(inFile,header= F,stringsAsFactors = F)
          optVarsV <- optVars$V1
          #print (paste('building optimised model for',allDMn[[c]],'; method:',mtd))
          dmOpt <- dm[,c(responseVar,optVarsV)]
          print (paste('prepping learning curve for',mN,'; method:',mtd))
          lCurve <- prepLearningCurve(dataIn = dmOpt,mdl = mtd,lcBoots = tButLC,trBoot = tButXV,trNumber = tRep,trainType="cv",trSet = 0.75,
                                      saveVarImp = F,responseVar = responseVar,posClass = posC,minSam = 10,samStep = 1.5,allowParallel=allowParallel)
          lp <- plotLearningCurves(lCurve,tit=paste('Opt learning curve (',mN,' [',ncol(dmOpt)-1,'f], ',mtd,')',sep=''),metrics = c("Kappa"))
          ggsave(plot = lp,paste(outFolderName,"/opt_lcurves/lc_",mtd,"_",mN,'_k.png',sep=''),width = 8,height = 6)
          lp <- plotLearningCurves(lCurve,tit=paste('Opt learning curve (',mN,' [',ncol(dmOpt)-1,'f], ',mtd,')',sep=''),metrics = c("Sensitivity","Specificity","ACC"))
          ggsave(plot = lp,paste(outFolderName,"/opt_lcurves/lc_",mtd,"_",mN,'_ssa.png',sep=''),width = 8,height = 6)
          lp <- plotLearningCurves(lCurve,tit=paste('Opt learning curve (',mN,' [',ncol(dmOpt)-1,'f], ',mtd,')',sep=''),metrics = c("NPV","PPV","AUC"))
          ggsave(plot = lp,paste(outFolderName,"/opt_lcurves/lc_",mtd,"_",mN,'_npa.png',sep = ''),width = 8,height = 6)
        } else {
          print (' WARNING: optimisation file does not exist, skipping optimised learning curves')
        }
      }
    }
    print ("OPTIMISED LEARNING CURVES DONE!")
  }

  # ======================================================================
  # MODEL TRAINING - prep raw models and test them
  # ======================================================================
  if (doRawMdls) {
    if (!dir.exists(paste(outFolderName,'/raw_fittedmdls',sep=''))) {dir.create(paste(outFolderName,'/raw_fittedmdls',sep=''))}
    if (!dir.exists(paste(outFolderName,'/raw_mdlResults',sep=''))) {dir.create(paste(outFolderName,'/raw_mdlResults',sep=''))}
    # general training parameters
    if (optimiseMetric=="ROC") {
      trC <- trainControl(method="repeatedcv",number=tButXV,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = allowParallel,summaryFunction = twoClassSummary)
    } else {
      trC <- trainControl(method="repeatedcv",number=tButXV,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = allowParallel)
    }
    # go through methods and data models
    fitRaw <- list()
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        # training set (pre-processed)
        if (doDataPrep) {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')  
        } else {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')  
        }
        # pre-process test set
        if (doDataPrep & inTrainP != 1) {
          tstSetRaw <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
          prepper <- readRDS(paste0(outFolderName,'/inputData/',mN,'_trSet_prepMdl.RData'))
          tstSet <- predict(prepper,tstSetRaw)
        } else if (inTrainP != 1) {
          tstSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
        }
        # train non-optimised model
        print (paste(" >> Training raw model for",allDMn[[c]],'/',mtd))
        # do training
        frm <- reformulate(response = responseVar, termlabels = ".")
        if (mtd == "glm") {
          fitRaw[[c]] <- train(frm,data = trSet,trControl = trC,method=mtd,metric=optimiseMetric,family="binomial") # logistic reg.
        } else {
          fitRaw[[c]] <- train(frm,data = trSet,trControl = trC,method=mtd,metric=optimiseMetric)
        }
        # save models
        if (saveMdls) {
          saveRDS(fitRaw[[c]], paste(outFolderName,"/raw_fittedmdls/mdl_raw_",mtd,"_",allDMn[[c]],'.rds',sep = ''))
        }
        # test set
        if (inTrainP != 1) {
          # 2-class
          if (length(unique(trSet[[responseVar]])) == 2) {
            compTest <- compareMdlsDatasets(mdls=list(fitRaw[[c]]),dataSets=list(tstSet),mdNames=c("test"),
                                          posClass=posC,tit = paste0(posC," Mdl <test set> (",allDMn[[c]],", ",mtd,")"),
                                          roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=T,response = responseVar)
          # save it
            ggsave(plot = compTest[[1]],paste(outFolderName,"/raw_mdlResults/roc_raw_",mtd,"_",allDMn[[c]],'_test.png',sep = ''),width = 8,height = 6)
            ggsave(plot = compTest[[2]],paste(outFolderName,"/raw_mdlResults/roc_raw_",mtd,"_",allDMn[[c]],'_test_val.png',sep = ''),width = 8,height = 6)
            write.table(compTest[[3]],paste(outFolderName,"/raw_mdlResults/roc_raw_",mtd,"_",allDMn[[c]],'_test.csv',sep = ''),row.names = F,sep=',')
          } else {
            # multi-class
            
          }
        }
        # x-validation
        # 2-class
        if (length(unique(trSet[[responseVar]])) == 2) {
          compXV <- compareModelsTrainingCV(fittedMdls = list(fitRaw[[c]]),modelNames=c("Raw-XV"),mtd=mtd,posClass=posC,
                                          tit = paste0(posC," Mdl (",allDMn[[c]],", ",mtd,")"),
                                          roc.smooth = F,roc.conf = T,annotateAUConly = T,removeLegend=T)
        # save it
          ggsave(plot = compXV,paste(outFolderName,"/raw_mdlResults/roc_raw_",mtd,"_",allDMn[[c]],'_xv.png',sep = ''),width = 8,height = 6)
        } else {
          # multi-class
          
        }
      }
    }
  }
    
  # ======================================================================
  # MODEL TRAINING - prep final (optimised) models and test them
  # ======================================================================
  if (doOptMdls) {
    if (!dir.exists(paste(outFolderName,'/opt_fittedmdls',sep=''))) {dir.create(paste(outFolderName,'/opt_fittedmdls',sep=''))}
    if (!dir.exists(paste(outFolderName,'/opt_mdlResults',sep=''))) {dir.create(paste(outFolderName,'/opt_mdlResults',sep=''))}
    defMethod = "glm"
    # general training parameters
    if (optimiseMetric=="ROC") {
      trC <- trainControl(method="repeatedcv",number=tButXV,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = allowParallel,summaryFunction = twoClassSummary)
    } else {
      trC <- trainControl(method="repeatedcv",number=tButXV,repeats = tRep,savePredictions = T,classProbs = T,allowParallel = allowParallel)
    }
    print ("BUILDING FINAL OPTIMISED MODELS:")
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        # training set (pre-processed)
        if (doDataPrep) {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')  
        } else {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')  
        }
        # pre-process test set
        if (doDataPrep & inTrainP != 1) {
          tstSetRaw <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
          prepper <- readRDS(paste0(outFolderName,'/inputData/',mN,'_trSet_prepMdl.RData'))
          tstSet <- predict(prepper,tstSetRaw)
        } else if (inTrainP != 1) {
          tstSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
        }
        # train optimised model (V2)
        if (doOpt2) {
          print (paste(" >> Training optimised model for",mN,'/',mtd,' (V2 optimisation) '))
          # get file with optimised variables
          tarOpt = "v2"
          inFileOpt = paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",mN,"_vars_",tarOpt,".csv")
          inFileDef = paste0(outFolderName,"/opt_RFE/rfe_",defMethod,"_",mN,"_vars_",tarOpt,".csv")
          # select optimal RFE, if it doesn't exist select default one; if that ones doesn't exist, don't do it    
          #print (paste('loading RFE data for',allDMn[[c]],'; method:',mtd))
          inFile = inFileOpt
          if (file.exists(inFileOpt)) {
            print (paste(' -> using',inFileOpt,'for optimisation!'))
            optVars <- read.table(inFile,header= F,stringsAsFactors = F)[[1]]
          } else if (file.exists(inFileDef)) {
            print (paste(' -> using',inFileDef,'for optimisation!'))
            optVars <- read.table(inFileDef,header= F,stringsAsFactors = F)[[1]]
          } else {
            optVars <- colnames(trSet)[!colnames(trSet) %in% responseVar]
          }
          print (paste("   -> using ",length(optVars)," / ",ncol(trSet)-1," features"))
          trSetT <- trSet[,c(responseVar,optVars)]
          # do training
          frm <- reformulate(response = responseVar, termlabels = ".")
          if (mtd == "glm") {
            fit <- train(frm,data = trSetT,trControl = trC,method=mtd,metric="Kappa",family="binomial") # logistic reg.
          } else {
            fit <- train(frm,data = trSetT,trControl = trC,method=mtd,metric="Kappa")
          }
          # save models
          if (saveMdls) {
            saveRDS(fit, paste(outFolderName,"/opt_fittedmdls/mdl_opt_",mtd,"_",mN,'_v2.rds',sep = ''))
          }
          # compare on test set
          # > 2 class prediction
          if (inTrainP != 1) {
            if (length(unique(trSet[[responseVar]])) == 2) {
              compTest <- compareMdlsDatasets(mdls=list(fit),dataSets=list(tstSet),mdNames=c("test"),
                                              posClass=posC,tit = paste0(posC," v2 opt mdl [",length(optVars),"f] <test set> (",mN,", ",mtd,")"),
                                              roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=T,response = responseVar)
              # save it
              ggsave(plot = compTest[[1]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_v2_test.png',sep = ''),width = 8,height = 6)
              ggsave(plot = compTest[[2]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_v2_test_val.png',sep = ''),width = 8,height = 6)
              write.table(compTest[[3]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_v2_test.csv',sep = ''),row.names = F,sep=',')
            } else {
              # > multi-class prediction
            }
          } 
          # x-validation
          if (length(unique(trSet[[responseVar]])) == 2) {
            compXV <- compareModelsTrainingCV(fittedMdls = list(fit),modelNames=c("Raw-XV"),mtd=mtd,posClass=posC,
                                              tit = paste0(posC," v2 opt mdl [",length(optVars),"f] <x-val> (",mN,", ",mtd,")"),
                                              roc.smooth = smoothROCs,roc.conf = T,annotateAUConly = T,removeLegend=T,roc.conf.boot = 10)
            # save it
            ggsave(plot = compXV,paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_v2_xv.png',sep = ''),width = 8,height = 6)
          } else {
            
          }
        }
        if (doOptMax) {
          print (paste(" >> Training optimised model for",mN,'/',mtd,' (max optimisation) '))
          # get file with optimised variables
          tarOpt = "vmax"
          inFileOpt = paste0(outFolderName,"/opt_RFE/rfe_",mtd,"_",mN,"_vars_",tarOpt,".csv")
          inFileDef = paste0(outFolderName,"/opt_RFE/rfe_",defMethod,"_",mN,"_vars_",tarOpt,".csv")
          # select optimal RFE, if it doesn't exist select default one; if that ones doesn't exist, don't do it    
          #print (paste('loading RFE data for',allDMn[[c]],'; method:',mtd))
          inFile = inFileOpt
          if (file.exists(inFileOpt)) {
            print (paste(' -> using',inFileOpt,'for optimisation!'))
            optVars <- read.table(inFile,header= F,stringsAsFactors = F)[[1]]
          } else if (file.exists(inFileDef)) {
            print (paste(' -> using',inFileDef,'for optimisation!'))
            optVars <- read.table(inFileDef,header= F,stringsAsFactors = F)[[1]]
          } else {
            optVars <- colnames(trSet)[!colnames(trSet) %in% responseVar]
          }
          print (paste("   -> using ",length(optVars)," / ",ncol(trSet)-1," features"))
          trSetT <- trSet[,c(responseVar,optVars)]
          # do training
          frm <- reformulate(response = responseVar, termlabels = ".")
          if (mtd == "glm") {
            fit <- train(frm,data = trSetT,trControl = trC,method=mtd,metric="Kappa",family="binomial") # logistic reg.
          } else {
            fit <- train(frm,data = trSetT,trControl = trC,method=mtd,metric="Kappa")
          }
          # save models
          if (saveMdls) {
            saveRDS(fit, paste(outFolderName,"/opt_fittedmdls/mdl_opt_",mtd,"_",mN,'_vMax.rds',sep = ''))
          }
          # compare on test set
          if (inTrainP != 1) {
            if (length(unique(trSet[[responseVar]])) == 2) {
              compTest <- compareMdlsDatasets(mdls=list(fit),dataSets=list(tstSet),mdNames=c("test"),
                                              posClass=posC,tit = paste0(posC," max opt mdl [",length(optVars),"f] <test set> (",mN,", ",mtd,")"),
                                              roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=T,response = responseVar)
              # save it
              ggsave(plot = compTest[[1]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_vmax_test.png',sep = ''),width = 8,height = 6)
              ggsave(plot = compTest[[2]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_vmax_test_val.png',sep = ''),width = 8,height = 6)
              write.table(compTest[[3]],paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_vmax_test.csv',sep = ''),row.names = F,sep=',')
            }
          }
          # x-validation
          # > 2 class prediction
          if (length(unique(trSet[[responseVar]])) == 2) {
            compXV <- compareModelsTrainingCV(fittedMdls = list(fit),modelNames=c("Raw-XV"),mtd=mtd,posClass=posC,
                                              tit = paste0(posC," max opt mdl [",length(optVars),"f] <x-val> (",mN,", ",mtd,")"),
                                              roc.smooth = smoothROCs,roc.conf = T,annotateAUConly = T,removeLegend=T,roc.conf.boot = 10)
            # save it
            ggsave(plot = compXV,paste(outFolderName,"/opt_mdlResults/roc_opt_",mtd,"_",mN,'_vmax_xv.png',sep = ''),width = 8,height = 6)
          } else {
          # > multi-class prediction
          }
        }
      }
    }
  }
  # ======================================================================
  # MODEL COMPARISON
  # ======================================================================
  #print(doMdlComparison)
  if (doMdlComparison) {
    print(" >> Comparing Optimised vs Raw Models")
    if (!dir.exists(paste(outFolderName,'/opt_vs_raw',sep=''))) {dir.create(paste(outFolderName,'/opt_vs_raw',sep=''))}
    for (mtd in mtds) {
      for (c in c(1:length(allDM))) {
        mN <- allDMn[[c]]
        print(paste0('   -> ',mN,' / ',mtd))
        # training set (pre-processed)
        if (doDataPrep) {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet_prep.csv'),header = T,sep=',')  
        } else {
          trSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_trSet.csv'),header = T,sep=',')  
        }
        # pre-process test set
        if (inTrainP != 1) {
          if (doDataPrep) {
            tstSetRaw <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
            prepper <- readRDS(paste0(outFolderName,'/inputData/',mN,'_trSet_prepMdl.RData'))
            tstSet <- predict(prepper,tstSetRaw)
          } else {
            tstSet <- read.table(file = paste0(outFolderName,'/inputData/',mN,'_tstSet.csv'),header = T,sep=',')  
          }
          # load raw and optimised models
          fitRaw <- readRDS(paste(outFolderName,"/raw_fittedmdls/mdl_raw_",mtd,"_",mN,'.rds',sep = ''))
          if (doOptMax) {
            fitOptM <- readRDS(paste(outFolderName,"/opt_fittedmdls/mdl_opt_",mtd,"_",mN,'_vMax.rds',sep = ''))
          } else {fitOptM <- NULL}
          if (doOpt2) {
            fitOpt2 <- readRDS(paste(outFolderName,"/opt_fittedmdls/mdl_opt_",mtd,"_",mN,'_v2.rds',sep = ''))
          } else {fitOpt2 <- NULL}
          # compare on test set
          # > 2 class prediction
          if (length(unique(trSet[[responseVar]])) == 2) {
            if (!is.null(fitOpt2) & !is.null(fitOptM)) {
              compTest <- compareMdlsDatasets(mdls=list(fitRaw,fitOptM,fitOpt2),dataSets=list(tstSet),mdNames=c("Raw","Opt[2]","Opt[M]"),
                                              posClass=posC,tit = paste0(posC,", Raw[",length(predictors(fitRaw)),"] vs Opt[",length(predictors(fitOpt2)),",",length(predictors(fitOptM)),"] <test set> (",mN,", ",mtd,")"),
                                              roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=F,response = responseVar)
            } else if (is.null(fitOpt2) & !is.null(fitOptM)) {
              compTest <- compareMdlsDatasets(mdls=list(fitRaw,fitOptM),dataSets=list(tstSet),mdNames=c("Raw","Opt[M]"),
                                              posClass=posC,tit = paste0(posC,", Raw[",length(predictors(fitRaw)),"] vs Opt[",length(predictors(fitOptM)),"] <test set> (",mN,", ",mtd,")"),
                                              roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=F,response = responseVar)
            } else if (!is.null(fitOpt2) & is.null(fitOptM)) {
              compTest <- compareMdlsDatasets(mdls=list(fitRaw,fitOpt2),dataSets=list(tstSet),mdNames=c("Raw","Opt[2]"),
                                              posClass=posC,tit = paste0(posC,", Raw[",length(predictors(fitRaw)),"] vs Opt[",length(predictors(fitOpt2)),"] <test set> (",mN,", ",mtd,")"),
                                              roc.conf.boot = 100,roc.smooth = smoothROCs,removeLegend=F,response = responseVar)
            }
            # save it
            ggsave(plot = compTest[[1]],paste(outFolderName,"/opt_vs_raw/roc_",mtd,"_",mN,'_test.png',sep = ''),width = 8,height = 6)
            ggsave(plot = compTest[[2]],paste(outFolderName,"/opt_vs_raw/roc_",mtd,"_",mN,'_test_val.png',sep = ''),width = 8,height = 6)
            write.table(compTest[[3]],paste(outFolderName,"/opt_vs_raw/roc_",mtd,"_",mN,'_test.csv',sep = ''),row.names = F,sep=',')
          } else {
            # > multi-class prediction
          }
          # x-validation
          # > 2 class prediction
          if (length(unique(trSet[[responseVar]])) == 2) {
            if (!is.null(fitOpt2) & !is.null(fitOptM)) {
              compXV <- compareModelsTrainingCV(fittedMdls = list(fitRaw,fitOptM,fitOpt2),modelNames=c("Raw-XV","Opt[M]","Opt[2]"),mtd=mtd,posClass=posC,
                                                tit = paste0(posC,", Raw[",length(predictors(fitRaw)),"] vs Opt[",length(predictors(fitOpt2)),",",length(predictors(fitOptM)),"] (",mN,", ",mtd,")"),
                                                roc.smooth = T,roc.conf = T,annotateAUConly = T,removeLegend=F,roc.conf.boot = 10)
            } else if (is.null(fitOpt2) & !is.null(fitOptM)) {
              compXV <- compareModelsTrainingCV(fittedMdls = list(fitRaw,fitOptM),modelNames=c("Raw-XV","Opt[M]"),mtd=mtd,posClass=posC,
                                                tit = paste0(posC," Opt mdl [",length(optVars),"f] <x-val> (",mN,", ",mtd,")"),
                                                roc.smooth = T,roc.conf = T,annotateAUConly = T,removeLegend=F,roc.conf.boot = 10)
            } else if (!is.null(fitOpt2) & is.null(fitOptM)) {
              compXV <- compareModelsTrainingCV(fittedMdls = list(fitRaw,fitOpt2),modelNames=c("Raw-XV","Opt[2]"),mtd=mtd,posClass=posC,
                                                tit = paste0(posC," Opt mdl [",length(optVars),"f] <x-val> (",mN,", ",mtd,")"),
                                                roc.smooth = T,roc.conf = T,annotateAUConly = T,removeLegend=F,roc.conf.boot = 10)
            }
            # save it
            ggsave(plot = compXV,paste(outFolderName,"/opt_vs_raw/roc_",mtd,"_",mN,'_xv.png',sep = ''),width = 8,height = 6)
          }
        } else {
          # multi-class prediction
          
        }
      }
    }
  }
}
 