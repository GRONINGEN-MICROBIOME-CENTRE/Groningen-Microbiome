# load libs & stuff
setwd('C:/Users/ranko/Documents/UMCG/DAG3_stats/')
source('leafletPatch.R')
source('../myLibs/R_Microbiome_scripts.R')
source('../myLibs/R_ML_scripts_v3b.R')
source('../myLibs/R_Misc.R')

library(vegan)

dag3 <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filt_rsc_asnor_corrtech.csv',sep=',',header = T)  
dag3 <- dag3[dag3$PROVINCE %in% c("Drenthe","Friesland","Groningen"),]

dag3raw <- read.table('data_processed/profiling/DAG3_taxa_pwys_CARD_VFDB_filtered_rescaled.csv',sep=',',header = T)  
dag3raw <- purgeMGNames(dag3raw)
dag3raw <- dag3raw[dag3raw$DAG3_sampleID %in% dag3$DAG3_sampleID,]

# fix factors
dag3$PROVINCE <- as.factor(as.character(dag3$PROVINCE))
# purge names
dag3 <- purgeMGNames(dag3)
toTest <- colnames(subsetMicrobiomeDF(dag3,getPWYs = T,getVFs = T,getCARDs = T,getTaxa = T,getPhenos = F))
toTestSpec <- colnames(dag3)[grep('^s__',colnames(dag3))]
phenos <- colnames(subsetMicrobiomeDF(dag3,getPWYs = F,getVFs = F,getCARDs = F,getTaxa = F,getPhenos = T))
allCARDS <- colnames(subsetMicrobiomeDF(dag3,getPWYs = F,getVFs = F,getCARDs = T,getTaxa = F,getPhenos = F))
allVFs <- colnames(subsetMicrobiomeDF(dag3,getPWYs = F,getVFs = T,getCARDs = F,getTaxa = F,getPhenos = F))

dag3$AllCARDS <- rowSums(dag3[,colnames(dag3) %in% allCARDS])
dag3$AllVFs <- rowSums(dag3[,colnames(dag3) %in% allVFs])

dag3$URBAN.INDEX.F <- as.factor(dag3$URBAN.INDEX)
dag3$URBAN.INDEX.4 <- dag3$URBAN.INDEX
dag3$URBAN.INDEX.4[dag3$URBAN.INDEX==5] <- 4
dag3$URBAN.INDEX.4.F <- as.factor(dag3$URBAN.INDEX.4)

dag3$div.shannon <- diversity(dag3raw[,toTestSpec])
dag3 <- dag3[dag3$div.shannon > mean(dag3$div.shannon) - 3 * sd(dag3$div.shannon),]


# male VS female microbiome
# ========================================
# look for biases (age, BMI)
ggplot(dag3,aes(x=age,col=gender)) + geom_density(size=2)
ggplot(dag3,aes(x=BMI,col=gender)) + geom_density(size=2)
ggplot(dag3,aes(y=div.shannon,x=gender,col=gender)) + geom_boxplot()
# correct for Age + BMI
dag3rdy <- linearCorrectMGPwy(dag3,corrNames = c("age","BMI"),corrMG = T,corrPWY = T,corrVFDB = T,corrCARD = T,
                              correctZeros = F,removeCorCol = T,negIsZero = F)
ggplot(dag3rdy,aes(x=age,col=gender)) + geom_density(size=2)
ggplot(dag3rdy,aes(y=div.shannon,x=gender,col=gender)) + geom_boxplot()

testsMFcorr <- data.frame()
toTestAll <- c("AllCARDS","AllVFs","div.shannon",toTest)
for (f in toTestAll) {
  print (paste(' > testing ',f,' VS gender'))
  testsMFcorr <- rbind.data.frame(testsMFcorr,testOneFeature(dataIn = dag3rdy,feature = f,responseVar = "gender",doSave = F,doPlots = F))
}
testsMFcorr$FDR <- p.adjust(testsMFcorr$pValue)
testsMFcorr <- testsMFcorr[order(testsMFcorr$FDR),]
c = 0
for (t in testsMF$feature[testsMF$FDR < 1.0e-10] ) {
  c = c + 1
  p <- testOneFeature(dag3rdy,saveFolder = paste0(""),feature = t,responseVar = "gender",
                      doPlots = T,discardZeros = F,retPlot = T,doSave = F)
  print(p[[2]])
}

# AGE associations
# ========================================
# correct for Gender + BMI
dag3rdy <- linearCorrectMGPwy(dag3,corrNames = c("BMI","gender"),corrMG = T,corrPWY = T,corrVFDB = T,corrCARD = T,
                              correctZeros = F,removeCorCol = T,negIsZero = F)
toTestAll <- c("AllCARDS","AllVFs","div.shannon",toTest)
testsAgeLMwZ <- testsLinearFit(inDF = dag3rdy,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = F)
testsAgeLMnZ <- testsLinearFit(inDF = dag3rdy,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = T)

plotTestsLM(inDF = dag3rdy,tests = testsAgeLMwZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 1.0e-2)
plotTestsLM(inDF = dag3rdy,tests = testsAgeLMnZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 1.0e-2)

# AGE associations (< 20) 
# ========================================
dag3rdyA20 <- dag3rdy[dag3rdy$age < 20,]
toTestAll <- c("AllCARDS","AllVFs","div.shannon",toTest)
testsAgeB20LMwZ <- testsLinearFit(inDF = dag3rdyA20,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = F)
testsAgeB20LMnZ <- testsLinearFit(inDF = dag3rdyA20,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = T)

plotTestsLM(inDF = dag3rdyA20,tests = testsAgeB20LMwZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 0.1)


# AGE associations ( 20 +) 
# ========================================
dag3rdy20plus <- dag3rdy[dag3rdy$age >= 20,]
toTestAll <- c("AllCARDS","AllVFs","div.shannon",toTest)
testsAge20PluswZ <- testsLinearFit(inDF = dag3rdy20plus,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = F)
testsAge20PlusnZ <- testsLinearFit(inDF = dag3rdy20plus,responseVar = "age",varsToTest = toTestSpec,correctors = c(),nonZero = T)

plotTestsLM(inDF = dag3rdy20plus,tests = testsAge20PluswZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 1.0e-2)
plotTestsLM(inDF = dag3rdy20plus,tests = testsAge20PlusnZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 1.0e-2)

# AGE associations ( 20 + prevalence only)
# ========================================
#dag3rdy20plus$age.F <- as.factor(dag3rdy20plus$age)
# ... <- testProportionLM(dag3rdy20plus,responseVar = "age",varsToTest = toTestSpec)

# BMI associations <AGE 18+>
# ========================================
dag3rdy <- linearCorrectMGPwy(dag3,corrNames = c("age","gender"),corrMG = T,corrPWY = T,corrVFDB = T,corrCARD = T,
                              correctZeros = F,removeCorCol = F,negIsZero = F)
dag3rdy <- dag3rdy[dag3rdy$BMI <= 45 & dag3rdy$age >= 18,]
toTestAll <- c("AllCARDS","AllVFs","div.shannon",toTest)
testsBMIwZ <- testsLinearFit(inDF = dag3rdy,responseVar = "BMI",varsToTest = toTestAll,correctors = c(),nonZero = F)
testsBMInZ <- testsLinearFit(inDF = dag3rdy,responseVar = "BMI",varsToTest = toTestAll,correctors = c(),nonZero = T)

# all
plotTestsLM(inDF = dag3rdy,tests = testsBMIwZ,responseVarLM = "BMI",responseVarX = "BMI",pointColVar = "BMI",FDRcutoff = 1.0e-5)
# non-zero
plotTestsLM(inDF = dag3rdy,tests = testsAge20PlusnZ,responseVarLM = "age",responseVarX = "age",pointColVar = "age",FDRcutoff = 1.0e-2)
# proportions

# >> BMI split into obesity levels
# underweight: BMI < 18.5; normal 18.5 - 25; pre-obese 25-30; obese: 30+; morbid-obese: 40+ (35+ and bad health ??)
dag3rdy$BMI.F <- "N"
dag3rdy$BMI.F[dag3rdy$BMI < 18.5] <- "UWt"
dag3rdy$BMI.F[dag3rdy$BMI > 25] <- "OWt"
dag3rdy$BMI.F[dag3rdy$BMI > 30] <- "OB"

dag3rdyNOB <- dag3rdy[dag3rdy$BMI.F %in% c("N","OB"),]
dag3rdyNOB$BMI.F <- as.factor(dag3rdyNOB$BMI.F)
testsNOBcorr <- data.frame()
for (f in toTestAll) {
  print (paste(' > testing ',f,' VS gender'))
  testsNOBcorr <- rbind.data.frame(testsNOBcorr,testOneFeature(dataIn = dag3rdyNOB, feature = f,responseVar = "BMI.F",doSave = F,doPlots = F))
}
testsNOBcorr$FDR <- p.adjust(testsNOBcorr$pValue)
testsNOBcorr <- testsNOBcorr[order(testsNOBcorr$pValue),]
#testsBMIwZ <- testsLinearFit(inDF = dag3rdy,responseVar = "BMI.F",varsToTest = toTestSpec,correctors = c(),nonZero = F)
c = 0
for (t in testsNOBcorr$feature[testsNOBcorr$FDR < 1.0e-10] ) {
  c = c + 1
  p <- testOneFeature(dag3rdyNOB,saveFolder = paste0(""),feature = t,responseVar = "BMI.F",
                      doPlots = T,discardZeros = F,retPlot = T,doSave = F)
  print(p[[2]])
}

# Test << URBAN.INDEX >>
# ================================================
# >> ALL DATA
responseVar = "URBAN.INDEX"
#correctors = c("age","gender","BMI","postclean.reads","vol.ul","conc.ng.ul")

correctors = c()
testsZ <- testsLinearFit(inDF = dag3,responseVar = "URBAN.INDEX",varsToTest = toTest,correctors = c(),nonZero = F)
testsNZ <- testsLinearFit(inDF = dag3,responseVar = "URBAN.INDEX",varsToTest = toTest,correctors = c(),nonZero = T)

testsZ <- testsLinearFit(inDF = dag3,responseVar = "LOG.POP.KM2",varsToTest = toTestSpec,correctors = c(),nonZero = T)
plotTestsLM(inDF = dag3,tests = testsZ,responseVarLM = "LOG.POP.KM2",responseVarX = "LOG.POP.KM2",pointColVar = "URBAN.INDEX.F",FDRcutoff = 1)
testsNZ <- testsLinearFit(inDF = dag3,responseVar = "LOG.POP.KM2",varsToTest = toTest,correctors = c(),nonZero = F)
plotTestsLM(inDF = dag3,tests = testsNZ,responseVarLM = "LOG.POP.KM2",responseVarX = "LOG.POP.KM2",pointColVar = "URBAN.INDEX.F")

dag3$POP.DENSITY.CUT <- cut(dag3$LOG.POP.KM2,breaks=5)

testsNZ <- testsLinearFit(inDF = dag3,responseVar = "LOG.POP.KM2",varsToTest = toTest,correctors = c(),nonZero = T)
plotTestsLM(inDF = dag3,tests = testsNZ,responseVarLM = "LOG.POP.KM2",responseVarX = "LOG.POP.KM2",pointColVar = "POP.DENSITY.CUT")

testsT <- testsLinearFit(inDF = dag3,responseVar = "LOG.POP.KM2",varsToTest = c("AllCARDS","AllVFs","div.shannon"),correctors = c(),nonZero = T)
plotTestsLM(inDF = dag3,tests = testsT,responseVarLM = "LOG.POP.KM2",responseVarX = "LOG.POP.KM2",pointColVar = "POP.DENSITY.CUT",FDRcutoff = 1,pVcutoff = 1)

# merge with new metadata table
demoTable <- read.table('data_processed/dutch_demo_table_2015d.csv',sep=',',header=T,quote = '"')
dag3wdemo <- merge(dag3,demoTable,by.x="municipaly",by.y="GAMENTEE")
dag3wdemo <- dag3wdemo[dag3wdemo$div.shannon > 1.5,]

phenosdemo <- colnames(demoTable)

# mass tester
res <- data.frame()
dag3wdemo[is.na(dag3wdemo)] <- 1.0
dag3wdemo[dag3wdemo <= 0] <- 1.0
for (geometa in colnames(demoTable)[2:110] ) {
  demoTable[[geometa]][is.na(demoTable[[geometa]])] <- 1.0
  demoTable[[geometa]][ demoTable[[geometa]] <= 0] <- 1.0
  print(paste('> testing',geometa))
  frm <- reformulate(response = "div.shannon",termlabels = paste0("log(",geometa,")"))
  mdl <- lm(data = dag3wdemo,formula = frm)
  rs <- summary(mdl)$r.squared
  pv <- summary(mdl)$coefficients[paste0("log(",geometa,")"),4]
  coeff <- summary(mdl)$coefficients[paste0("log(",geometa,")") ,1]
  res <- rbind.data.frame(res,data.frame(feature=geometa,rsquared=rs,coef=coeff,pvalue=pv ))
}
res <- res[order(res$pvalue),]
res$FDR <- p.adjust(res$pvalue)

for (r in res$feature[res$pvalue < 0.05]) {
  g <- ggplot(dag3wdemo,aes_string(x=paste0("log(",r,")"),y="div.shannon")) + geom_point() + geom_smooth(method="lm")
  print(g)
}

ggplot(dag3wdemo,aes(x=municipaly,y=div.shannon)) + geom_boxplot()
dfMUNI <- data.frame()
for (gam in unique(dag3wdemo$municipaly)) {
  print(gam)
  shan <- dag3wdemo$div.shannon[dag3wdemo$municipaly==gam]
  dutchness <- dag3wdemo$POP.ETHNIC.DUTCH.PERC[dag3wdemo$municipaly==gam]
  dfMUNI <- rbind.data.frame(dfMUNI,data.frame(gameteen=gam,nr=length(shan),shannon.avg=mean(shan),
                                               shannon.sd=sd(shan),shannon.med=median(shan),dutchness=mean(dutchness)))
}
dfMUNI <- dfMUNI[order(dfMUNI$shannon.med,decreasing = T),]


# DO PROPORTION TESTING 
# ==========================
# > data -> 0/1
# >> ALL DATA
testProportionLM <- function(inDF,responseVar,varsToTest) {
  testsProp <- data.frame()
  inDF[,varsToTest][inDF[,varsToTest] > 0] <- 1
  for (f in varsToTest) {
    propDF <- data.frame()
    print(paste(' >> testing',f))
    for (u in sort(unique(inDF[[responseVar]]))) {
      propF <- sum(inDF[inDF[[responseVar]]==u,f])/length(inDF[inDF[[responseVar]]==u,f])
      propDF <- rbind.data.frame(propDF,data.frame(RVAR=u,prop=propF))
    }
    propDF$RVARNR <- as.numeric(row.names(propDF))
    mdl <- lm(data = propDF, formula =  prop ~ RVARNR)
    rs <- summary(mdl)$r.squared
    if (length(summary(mdl)$coefficients[,4]) > 1) {
      coeff <- summary(mdl)$coefficients[,1][[2]]
      pv <- summary(mdl)$coefficients[,4][[2]]
      testsProp <- rbind.data.frame(testsProp,data.frame(feature=f,coef=coeff,rsquared=rs,pvalue=pv))
    }
  }
  testsProp <- testsProp[order(testsProp$pvalue),]
  testsProp$FDR <- p.adjust(testsProp$pvalue)
  testsProp
}

testsP <- testProportionLM(inDF = dag3,responseVar = "POP.DENSITY.CUT",varsToTest = toTestSpec)

# select 'good looking' subset for plotting
teststest <- testsP[testsP$coef > 0.005 & testsP$rsquared > 0.8,]
for (f in teststest$feature) {
  print(f)
  propDF <- data.frame()
  for (u in sort(unique(inDF[[responseVar]]))) {
    propF <- sum(inDF[inDF[[responseVar]]==u,f])/length(inDF[inDF[[responseVar]]==u,f])
    propDF <- rbind.data.frame(propDF,data.frame(ResponseVar=u,prop=propF))
  }
  #propDF$URBAN.INDEX.F <- as.factor(propDF[[responseVar]])
  g <- ggplot(propDF,aes_string(x="ResponseVar",y="prop",fill="ResponseVar")) + geom_col() + 
    ylab(paste0("Prevalence of ",f)) + xlab(responseVar)
  print(g)
}
