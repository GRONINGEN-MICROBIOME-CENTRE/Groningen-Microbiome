# ===============================================
# ===============================================
# > DMP code for building heritability models
# ===============================================
# ===============================================
# by: R.Gacesa (UMCG, 2021)
#
# implementation of pathway heritability analysis on mock data

# =============================================
# function for permutation analysis of model
# =============================================
# input: 
#  - model data (dataframe)
#  - model formula
#  - kinship matrix
# return data frame with row = proportions of variance explained by random variables
# =======================================================================================
doOnePerm <- function(mdlData,mdlFormula,kinMatrix,permuteVars=T,permuteKinship=F,
                      toPermute=c("ID"),seed=NA) {
  inDFrdyPerm <- mdlData
  # permute vars
  # > if seed != NA: we use fixed seed to ensure same permutation structure
  if (!is.na(seed)) {
    set.seed(seed)
  }
  if (permuteVars) {
    for (v in toPermute) {
      inDFrdyPerm[[v]] <- sample(inDFrdyPerm[[v]],replace = F)
      #debug: test if fixing seeds worked [it does]
      #print(inDFrdyPerm[[v]][1:5])
    }
  }
  kinMatPerm <- kinMatrix
  # permute kinship matrix
  if (permuteKinship) {
    pID <- sample(rownames(kinMatPerm),replace = F)
    rownames(kinMatPerm) <- pID
    colnames(kinMatPerm) <- pID
  }
  mdlPerm <- relmatLmer(mdlFormula, data = inDFrdyPerm, relmat = list(ID = kinMatPerm),REML=F)
  #    get variance explained by random effects
  varPropMdl <- lme4qtl::VarProp(mdlPerm)
  result <- varPropMdl[,c("grp","prop")]
  rownames(result) <- result$grp
  resultt <- as.data.frame(t(result))
  resultt <- resultt[2,]
  rownames(resultt) <- NULL
  for (cc in colnames(resultt)) {resultt[[cc]] <- as.numeric(as.character(resultt[[cc]]))}
  resultt
}

# install packages:
# install.packages(c("compositions","coxme","devtools","lmerTest"))
# devtools::install_github("variani/lme4qtl")

# LOAD LIBS
library(compositions)
library(coxme)
library(lme4qtl)
library(RLRsim)
library(lmerTest)
library(foreach)
library(plyr)

# === INITIALIZATION ===
# =================================
# NOTE: inFld should be changed to appropriate path on the hard drive
# example: #inFld <- "D:/Vbox/shared/dag/git_14_05/DMP/heritability_analysis_v2/"
inFld <- "."

setwd(inFld)

# === OPTIONS ===
# should we calculate confidence intervals?
# NOTE: it is slow and somewhat unstable but usually works
doConfInt <- T
# should we run permutation test?
doPermutations <- T
#  > number of permutations
nPerm <- 10
# do we keep random seeds fixed (for reproducibility)
useFixedSeeds <- T
# MODEL SELECTION: if T, builds selection of models and takes lowest AIC one
doModelSelection <- T
# keepOnlyRelated: if T: removes zero-sum columns and rows from kinship matrix
keepOnlyRelated = F
# dropSingletFamsd: if T: removes all families with only 1 individual
dropSingletFams = T

# Load Data
# =======================
# >> phenotypes
inPhenos3 <- read.table('mock_data/DMP_heritability_phenos_mock.csv', sep=',',header=T, stringsAsFactors = F)
# shuffle IDs for merging
rownames(inPhenos3) <- inPhenos3$PSEUDOIDEXT
inPhenos3$PSEUDOIDEXT <- NULL; inPhenos3$ID <- NULL; inPhenos3$ID <- rownames(inPhenos3)

# >> microbiome (untransformed)
inMBmfraw <- read.table('mock_data/DMP_mock_microbiome_pwys_filtered.csv',sep=',',header=T,row.names = 1)

# >> microbiome (CLR transformed)
inMBmfclr <- read.table('mock_data/DMP_mock_microbiome_pwys_CLR_filtered.csv',sep=',',header=T,row.names = 1)

# >> kinship matrix
# NOTE: x 2 multiplication is to turn it into actual genetics shared matrix from kinship implementation
kinMat <- readRDS('mock_data/DMP_heritability_kinmatrix_mock.RDS')*2

# >> load seeds (if using fixed seeds)
if (useFixedSeeds) {
  fixedSeeds <- read.table('mock_data/seeds.txt',header=T)
}
# >> merge everything
# =======================
inDFrdyCLR <- merge(inPhenos3,inMBmfclr,by="row.names")

# >> iterate over taxa
# ============================
#for (taxNR in c(10:12)) {# c(1:ncol(inMBmfclr))) {
for (taxNR in c(173:ncol(inMBmfclr))) {
  tax <- colnames(inMBmfclr)[[taxNR]]
  # set permutation variable
  varsToPermute <- tax
  # Lme4qtl implementation of model fitting
  # =================================================
  # > define formula for model
  # primary model
  modelParams1 <- c("Age",
                    "Age^2",
                    "Sex",
                    "META.POOP.Freq","META.DNA.postclean.reads",
                    "(1|ID)",
                    "(1|famID)",
                    "(1|COHOUSING.ID_DMP)")
  if (doModelSelection) {
    # reduced model: no family grouping
    modelParams2 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|ID)",
                      "(1|COHOUSING.ID_DMP)")
    # reduced model: no cohousing grouping
    modelParams3 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|ID)",
                      "(1|famID)")
    # reduced model: ID only
    modelParams4 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|ID)")
    # reduced model: no ID
    modelParams5 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|famID)",
                      "(1|COHOUSING.ID_DMP)")
    # reduced model: family only
    modelParams6 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|famID)")
    # reduced model: cohousing only
    modelParams7 <- c("Age",
                      "Age^2",
                      "Sex",
                      "META.POOP.Freq","META.DNA.postclean.reads",
                      "(1|COHOUSING.ID_DMP)")
  }
  frm1 <- reformulate(response = tax,termlabels = modelParams1)
  if (doModelSelection) {
    frm2 <- reformulate(response = tax,termlabels = modelParams2)
    frm3 <- reformulate(response = tax,termlabels = modelParams3)
    frm4 <- reformulate(response = tax,termlabels = modelParams4)
    frm5 <- reformulate(response = tax,termlabels = modelParams5)
    frm6 <- reformulate(response = tax,termlabels = modelParams6)
    frm7 <- reformulate(response = tax,termlabels = modelParams7)
  }
  
  # > final data check & prep:
  # - NAs have to be purged from META.POOP.Freq
  # - META.DNA.postclean.reads has to be transformed, otherwise model refuses to converge
  # - we have to get rid of NAs in random effect(s)
  inDFrdyR <- inDFrdyCLR
  # impute poop freq to median
  if ("META.POOP.Freq" %in% modelParams1) {
    inDFrdyR$META.POOP.Freq[is.na(inDFrdyR$META.POOP.Freq)] <- median(inDFrdyR$META.POOP.Freq,na.rm = T)
  }
  if ("META.DNA.postclean.reads" %in% modelParams1) {
    # log-transform read depth
    inDFrdyR$META.DNA.postclean.reads <- log(inDFrdyR$META.DNA.postclean.reads)
  }
  # Remove samples without housing information
  if ("(1|COHOUSING.ID_DMP)" %in% modelParams1) {
    inDFrdyR <- inDFrdyR[!is.na(inDFrdyR$COHOUSING.ID_DMP),]
  }
  # prep for model building
  res <- NULL # results go here
  resList <- list()
  
  # MODEL BUILDING
  print(paste0(' >> FITTING HERITABILITY MODELS FOR ',tax))
  
  # drop NA rows from data 
  inDFrdyR2 <- inDFrdyR[!is.na(inMBmfclr[[tax]]),]
  # double check for NAs
  # clean bad metadata samples
  inDFrdyR2 <- inDFrdyR2[!is.na(inDFrdyR2$Age),]
  inDFrdyR2 <- inDFrdyR2[!is.na(inDFrdyR2$Sex),]
  inDFrdyR2 <- inDFrdyR2[!is.na(inDFrdyR2$BMI),]
  
  # if keepOnlyRelated:
  #   keep only samples for which kinship matrix is avaliable and non-0
  if (keepOnlyRelated) {
    inDFrdyR2 <- inDFrdyR2[inDFrdyR2$ID %in% row.names(kinMat),]
    rr <- rowSums(kinMat)
    rr2 <- rr[rr != 1]
    inDFrdyR2 <- inDFrdyR2[inDFrdyR2$ID %in% names(rr2),]
    inMBmfraw2 <- inMBmfraw[inMBmfraw$ID %in% inDFrdyR2$ID]
  } else {
    inMBmfraw2 <- inMBmfraw
  }
  
  # if dropSingletFams:
  #   drop singlet families
  if (dropSingletFams) {
    print ('dropping singlet families!')
    ft <- as.data.frame(table(inDFrdyR2$famID))
    ft <- ft[ft$Freq > 1,]
    inDFrdyR2 <- inDFrdyR2[inDFrdyR2$famID %in% ft$Var1,]
    print (paste0(' > retained ',nrow(inDFrdyR2),' samples'))
    
  }
  
  # fit the models
  # =========================
  # > primary (all random effects)
  fittedMdl1 <- relmatLmer(frm1, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
  
  if (doModelSelection) {
    # > reduced model 1 
    fittedMdlr2 <- relmatLmer(frm2, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
    # > reduced model 2
    fittedMdlr3 <- relmatLmer(frm3, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
    # > reduced model 3
    fittedMdlr4 <- relmatLmer(frm4, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
    # > reduced model 4
    fittedMdlr5 <- relmatLmer(frm5, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
    # > reduced model 5
    fittedMdlr6 <- relmatLmer(frm6, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
    # > reduced model 6
    fittedMdlr7 <- relmatLmer(frm7, data = inDFrdyR2, relmat = list(ID = kinMat),REML=F)
  }
  
  if (doModelSelection) {
    # compare models
    # =========================
    aic1 <- extractAIC(fittedMdl1)[2]
    aic2 <- extractAIC(fittedMdlr2)[2]
    aic3 <- extractAIC(fittedMdlr3)[2]
    aic4 <- extractAIC(fittedMdlr4)[2]
    aic5 <- extractAIC(fittedMdlr5)[2]
    aic6 <- extractAIC(fittedMdlr6)[2]
    aic7 <- extractAIC(fittedMdlr7)[2]
    # debug
    print(paste0('Full model      AIC = ',aic1))
    print(paste0(' no family      AIC = ',aic2))
    print(paste0(' no cohousing   AIC = ',aic3))
    print(paste0(' herit only.    AIC = ',aic4))
    print(paste0(' no herit.      AIC = ',aic5))
    print(paste0(' family only .  AIC = ',aic6))
    print(paste0(' housing only.  AIC = ',aic7))
    if (aic1 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 1 (all random effects)')
      fittedMdl <- fittedMdl1
      mName <- "full model"
    } else if (aic2 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 2 [no family ID]')
      fittedMdl <- fittedMdlr2
      mName <- "no-family"
    } else if (aic3 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 3 [no cohousing ID]')
      fittedMdl <- fittedMdlr3
      mName <- "no-cohousing"
    } else if (aic4 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 4 [heritability only]')
      fittedMdl <- fittedMdlr4
      mName <- "heritability.only"
    } else if (aic5 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 5 [no heritability]]')
      fittedMdl <- fittedMdlr5
      mName <- "no heritability"
    } else if (aic6 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 6 family only]')
      fittedMdl <- fittedMdlr6
      mName <- "family only"
    } else if (aic7 == min(aic1,aic2,aic3,aic4,aic5,aic6,aic7)) {
      print('selecting model 7 cohousing only]')
      fittedMdl <- fittedMdlr7
      mName <- "cohousing only"
    }
  } else {
    mName <- "full model"
    fittedMdl <- fittedMdl0
  }
  
  # Extract random effects stats
  # =============================
  # > get variance explained by random effects
  # NOTE: variance explained by ID is H2
  print ('   >> extracting random effects & calculating statistics')
  rndEffProp <- lme4qtl::VarProp(fittedMdl)
  # clean tables
  rownames(rndEffProp) <- rndEffProp$grp
  rndEffProp$var1 <- NULL; rndEffProp$var2 <- NULL; rndEffProp$grp <- NULL
  colnames(rndEffProp) <- c("Variance","Variance.SD","VarExp.Prop")
  rndEffProp$Variance <- NULL; rndEffProp$Variance.SD <- NULL
  
  # RND EFFECTS: profile confidence intervals
  # ============================================
  # NOTE: try-catch block makes sure that loop survives
  #       if profiling crashes (it tends to with low prevalence taxa)
  if (doConfInt) {
    rndEffProp2 = tryCatch({
      print('   >> calculating confidence intervals')
      # > profile results
      print('      >> profiling')
      prof <- profile(fittedMdl)
      #  # > convert to proportions
      # >> drop duplicate values to prevent some instances of spline error:
      prof2 <- prof
      for (cc in colnames(prof2)[!colnames(prof2) %in% ".par"]) {
        #        print(cc)  
        prof2 <- prof2[!duplicated(prof2[[cc]]),]
      }
      print('      >> getting proportions')
      profProp <- lme4qtl::varpropProf(prof2,ranef=T) # convert to proportions
      # # > get confidence intervals
      print('      >> calculating CI')
      confInt <- as.data.frame(confint(profProp,level = 0.95))
      confInt2 <- confInt[c(1:nrow(rndEffProp)),]
      rownames(confInt2) <- rownames(rndEffProp)
      confInt2$CI <- paste0(round(confInt2[[1]],2),'-',round(confInt2[[2]],2))
      colnames(confInt2) <- c("CI.l","CI.h","CI")
      # add these to random effect variance explained
      rndEffProp2 <- merge(rndEffProp,confInt2,by="row.names")
      rndEffProp2$CI.h <- NULL; rndEffProp2$CI.l <- NULL
      row.names(rndEffProp2) <- rndEffProp2$Row.names; rndEffProp2$Row.names <- NULL        
      rndEffProp2
    }, error = function(e) {
      print(e)
      rndEffProp2 <- rndEffProp
      rndEffProp2$CI <- NA
      rndEffProp2
    })
  } else {
    rndEffProp2 <- rndEffProp
  }
  # RND EFFECTS: do permutations:
  # =================================
  if (doPermutations) {
    # NOTE: try-catch block makes sure that loop survives
    #       if one permutation crashes (it is rare, but happens)
    #       in that case, that permutation is dropped from the permutations table
    print('   >> running permutation analysis')
    # >> NOTE: only permute taxon -> ID relationship
    #   >> then run all analyses for the permutation,
    #     >> then repeat
    #       >> then calculate H2s
    # >> NOTE: only permute taxa, we want to keep predictor structure!!
    permRndEff <- foreach(i=c(1:nPerm),.packages=c("lme4qtl"),.combine=rbind.data.frame) %do% {
      tmpOut = tryCatch({
        if (useFixedSeeds) {
          seedN <- useFixedSeeds[i]
        } else {
          seedN <- NA
        }
        oneP <- doOnePerm(mdlData = inDFrdyR2,mdlFormula=frm1,kinMatrix=kinMat,
                          permuteVars = T,toPermute = varsToPermute,
                          permuteKinship = F,seed = seedN)
        oneP
      }, error = function(e) {
        print(e)
        oneP <- NULL
        oneP
      })
    }
    # get empirical FDR for one permutation
    empPVtbl <- data.frame()
    print('     >> calculating empirical FDR')
    for (cc in c(colnames(permRndEff))) {
      if (cc != "Residual") {
        empPV <- sum(rndEffProp2[rownames(rndEffProp2)==cc,]$VarExp.Prop <= permRndEff[[cc]])/nPerm
      } else {
        empPV <- sum(rndEffProp2[rownames(rndEffProp2)==cc,]$VarExp.Prop >= permRndEff[[cc]])/nPerm
      }
      empPVtbl <- rbind.data.frame(empPVtbl,data.frame(Var=cc,empPV=empPV))
    }
    # merge with main table
    print('     >> saving permutation run')
    rownames(empPVtbl) <- empPVtbl$Var
    rndEffProp3 <- merge(rndEffProp2,empPVtbl,by="row.names")
    rownames(rndEffProp3) <- rndEffProp3$Row.names
    rndEffProp3$Var <- NULL
    rndEffProp3$Row.names <- NULL
    # save permutation run
    permRun <- list()
    permRun[["tax"]] <- tax
    permRun[["nPerm"]] <- nPerm
    permRun[["varsToPermute"]] <- varsToPermute
    permRun[["realFittedMdl"]] <- fittedMdl
    permRun[["rndEffProp"]] <- rndEffProp2
    permRun[["permutations"]] <- permRndEff
    saveRDS(permRun,file=paste0('mock_data_permutation_runs/','permrun_pwy_',taxNR,'.RDS'))
  } else {
    rndEffProp3 <- rndEffProp2
  }
  
  # RND EFFECTS: get rest of stats for random effects
  # =========================================================
  #print(paste0('H2 = ',round(h2,2),', sd = ',round(h2sd,2) ))
  print('   >> calculating log-likelihood tests for random effects')
  rndRanova <- as.data.frame(lmerTest::ranova(fittedMdl))
  # make names compatible for merging
  row.names(rndRanova) <- gsub('^<none>','Residual',row.names(rndRanova))
  row.names(rndRanova) <- gsub('^\\(1 \\| ','',row.names(rndRanova))
  row.names(rndRanova) <- gsub(')$','',row.names(rndRanova))
  #rownames(rndRanova)[rownames(rndRanova)=="ID"] <- "Kinship"
  rndTable <- merge(rndEffProp3,rndRanova,by="row.names")
  print(rndTable)
  # drop unnecessary data
  rndTable$AIC <- NULL; rndTable$Df <- NULL
  rndTable$npar <- NULL; rndTable$LRT <- NULL
  if (doConfInt & doPermutations) {
    colnames(rndTable) <- c("VAR","VE","CI","empPV","LogL","Pv")
  } else if (doConfInt & !doPermutations) {
    colnames(rndTable) <- c("VAR","VE","CI","LogL","Pv")
  } else {
    colnames(rndTable) <- c("VAR","VE","LogL","Pv")
  }
  resList[[paste0(tax,'_rnd')]] <- rndTable
  # make wide
  rndTable$ID <- "X"
  # make wide
  rndEffwide <- reshape(rndTable, idvar="ID", timevar="VAR", direction="wide", sep="_")  
  rndEffwide$ID <- NULL; rndEffwide$Pv_Residual <- NULL
  
  # Extract fixed effect stats, quantify p-values
  # =================================================
  # Analysis of Deviance (Type II Wald chisquare tests)
  print('   >> extracting fixed effects & calculating statistics')
  fixedEffectsPv <- car::Anova(mod = fittedMdl)
  # > grab rest of numbers for fixed effects
  fixedEffT2 <- summary(fittedMdl)[[10]]
  # > fix names to allow merging
  rownames(fixedEffT2) <- gsub('SexM','Sex',row.names(fixedEffT2))
  rownames(fixedEffT2) <- gsub('AddressChangedY','AddressChanged',row.names(fixedEffT2))
  # > merge
  fixedEffTable <- merge(fixedEffT2,fixedEffectsPv,by='row.names')
  # > clean
  fixedEffTable$Df <- NULL ; fixedEffTable$Chisq <- NULL; fixedEffTable$`t value` <- NULL
  rownames(fixedEffTable) <- fixedEffTable$Row.names; fixedEffTable$Row.names <- NULL
  colnames(fixedEffTable) <- c("Est","SE","P.Chisq")
  print(fixedEffTable)
  # save to list
  resList[[paste0(tax,'_fixed')]] <- fixedEffTable
  # dummy vars for transformation
  fixedEffTable$dummy <- row.names(fixedEffTable)
  fixedEffTable$ID <- "X"
  # make wide
  fixedEffwide <- reshape(fixedEffTable, idvar="ID", timevar="dummy", direction="wide", sep="_")
  # beautify
  fixedEffwide$ID <- NULL
  
  # merge with random effect table and fixed effect table
  # =============================================
  print('   >> preparing final table for output')
  tblRdy <- cbind.data.frame(rndEffwide,fixedEffwide)
  
  # other stuff to add to table
  # =============================
  # trait name and short name (for taxa)
  tblRdy$Trait <- tax
  ss <- strsplit(tblRdy$Trait,'\\.')[[1]]
  tblRdy$Trait.short <- ss[length(ss)]
  
  # number statistics (prevalence of taxon)
  nZ <- sum(inMBmfraw2[[tax]] == 0)
  nnZ <- sum(inMBmfraw2[[tax]] != 0)
  prev <- nnZ/(nZ+nnZ)
  tblRdy$N <- length(inDFrdyR2[[tax]])
  tblRdy$Nz <- nZ
  tblRdy$Prev <- round(prev*100,2)
  tblRdy$Model <- mName    
  # add to results
  res <- rbind.fill(res,tblRdy)
  print(paste0(' >> DONE with ',tax))
  print('========================================')
  
  #saveRDS(resList,paste0('mock_data_heritability_results/results_pwys_',taxNR,'.RDS'))
  write.table(res,paste0('mock_data_heritability_results/results_pwys_',taxNR,'.csv'),sep=',',row.names = F)
}
