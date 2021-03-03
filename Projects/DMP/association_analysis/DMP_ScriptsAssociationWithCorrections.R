### Script to correct the phenotype-species associations for mutiple covariates###
## Author: A.Kur, 12 th August, 2020
## Updates: T.Sinha, 16th August, 2020
##          R.Gacesa, 22/09/2020
## =================================================
library(data.table)
library(foreach)
library(doParallel)

## CLR transformation function 
do_clr_externalWeighting = function(interest_matrix, core_matrix) {
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  #estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

# ========================================
# ========================================
#   ASSOCIATION ANALYSIS FOR DAG3 
# ========================================
# ========================================
dag3AssociationsWithCorrections <- function(dataType='taxa',
                                            covariates=c(),
                                            idsToKeep = NULL,
                                            dataPath='/groups/umcg-lifelines/tmp01/users/umcg-akurilshchikov') {
  # define parameters
  dataType = dataType
  #phenoOut = outFile
  pheno2covar = covariates
  # error check
  if (!(dataType %in% c("pathways","taxa","CARD","VFDB") )) {
    stop("ERROR: dataType must be one of [taxa,pathways,CARD,VFDB]")
  }
  ## LOAD DATA
  # ========================================
  print(' > Loading & Preparing Data !')
  workPWD = dataPath
  print ('  >> Taxa')
  # > taxa
  taxa = read.table(paste0(workPWD,"/Association.analysis.input/taxa.txt"))
  if (!is.null(idsToKeep)) {
    taxa <- taxa[rownames(taxa) %in% idsToKeep,]
    if (nrow(taxa) == 0) {
      stop("ERROR: No samples kept after subsetting by idsToKeep, quitting!")
    }
  }
  #   > do CLR
  taxa_transformed = do_clr_externalWeighting(taxa,taxa[,grep("[.]s__",colnames(taxa))]) # Performing transformation 
  taxa_transformed = taxa_transformed[,colSums(taxa>0)>nrow(taxa) * 0.05] # Filtering 
  print ('    >> Done!')
  
  # > pathways
  print ('  >> PWYs')
  pathways = read.table(paste0(workPWD,"/Association.analysis.input/pathways_metacyc.txt"))
  if (!is.null(idsToKeep)) {
    pathways <- pathways[rownames(pathways) %in% idsToKeep,]
    if (nrow(pathways) == 0) {
      stop("ERROR: No samples kept after subsetting by idsToKeep, quitting!")
    }
  }
  #   > do CLR
  pathways_transformed = do_clr_externalWeighting(pathways,pathways) # Performing transformation 
  pathways_transformed = pathways_transformed[,colSums(pathways>0)>nrow(pathways)*0.05]# Filtering 
  print ('    >> Done!')
  
  # > CARD
  print ('  >> CARDs')
  CARD = read.table(paste0(workPWD,"/Association.analysis.input/DAG3_CARD_nofiltering.txt"),header=T,sep='\t',stringsAsFactors = F,row.names = 1)
  CARD = CARD[intersect(rownames(CARD),rownames(taxa)),]
  if (!is.null(idsToKeep)) {
    CARD <- CARD[rownames(CARD) %in% idsToKeep,]
    if (nrow(CARD) == 0) {
      stop("ERROR: No samples kept after subsetting by idsToKeep, quitting!")
    }
  }
  #  > do CLR
  CARD_transformed = do_clr_externalWeighting(CARD,CARD)
  #  > filter by prevalence
  CARD_transformed = CARD_transformed[,colSums(CARD>0)>0.05 * nrow(CARD)]
  print ('    >> Done!')
  
  # VFs
  print ('  >> VFs')
  VFDB = read.table(paste0(workPWD,"/Association.analysis.input/DAG3_VFDB_VFs_nofiltering.txt"),header=T,sep='\t',stringsAsFactors = F,row.names = 1)
  VFDB = VFDB[intersect(rownames(taxa),rownames(VFDB)),]
  if (!is.null(idsToKeep)) {
    VFDB <- VFDB[rownames(VFDB) %in% idsToKeep,]
    if (nrow(VFDB) == 0) {
      stop("ERROR: No samples kept after subsetting by idsToKeep, quitting!")
    }
  }
  #  > do CLR
  VFDB_transformed = do_clr_externalWeighting(VFDB,VFDB)
  #  > filter by prevalence
  VFDB_transformed = VFDB_transformed[,colSums(VFDB>0)>0.05 * nrow(VFDB)]
  print ('    >> Done!')
  
  # > load phenotypes
  print ('  >> Phenotypes')
  load(paste0(workPWD,"/Association.analysis.input/pheno_release27.RData")) 
  pheno = pheno27[rownames(taxa_transformed),]
  if (!is.null(idsToKeep)) {
    pheno <- pheno[rownames(pheno) %in% idsToKeep,]
    if (nrow(pheno) == 0) {
      stop("ERROR: No samples kept after subsetting by idsToKeep, quitting!")
    }
  }
  print ('    >> Done!')
  
  # ===============================================================
  # prepare covariate data frame
  # ===============================================================
  # > check for missing covariates
  missingCovar <- pheno2covar[!(pheno2covar %in% colnames(pheno))]
  if (length(missingCovar) > 0) {
    print(paste0('WARNING: Covariates < ',paste(missingCovar,collapse=', '),' > missing from phenotype file!'))
  }
  # > technical covariates (these should always be included)
  covar = read.table(paste0(workPWD,"/Association.analysis.input/covariates.txt"),sep='\t') 
  covar$GeoMean <- NULL #Remove geometric mean
  # > TODO: allow retention of these
  # remove covariates from phenotypes
  pheno$META.BATCH <- NULL
  pheno$META.DNA.conc.ng.ul <- NULL
  pheno$META.POOP.COLLECTION_SEASON <- NULL
  pheno$ANTHRO.AGE <- NULL
  pheno$ANTHRO.Sex <- NULL
  
  # clean covariate matrix & add extra covariates
  covar = covar[rownames(taxa_transformed),]
  covar = data.frame(covar,pheno[,which(colnames(pheno) %in% pheno2covar),drop = F])
  
  # drop extra covariates
  pheno <- pheno[,!(colnames(pheno) %in% pheno2covar)]
  pheno <- pheno[,!(colnames(pheno) %in% colnames(covar) )]
  
  # ===========================================
  # RUN ANALYSIS
  # ===========================================
  print (paste0(" > Building associations models for ",dataType,""))
  # ===========================================
  
  # ERROR TESTS
  # ==========================================
  if (dataType=="pathways") {
    ftrs_transformed = pathways_transformed
  } else if (dataType=="CARD"){
    ftrs_transformed = CARD_transformed
  } else if (dataType == "VFDB") {
    ftrs_transformed = VFDB_transformed
  } else if (dataType == "taxa") {
    ftrs_transformed = taxa_transformed
  } else {
    quit ("ERROR:wrong data type requested!")
  }
  
  # ====================================================
  # Run multivariate models, multi-thread implementation
  # ======================================================
  # prep parallelization
  #registerDoSEQ() # debug mode = single threaded
  registerDoParallel(makeCluster(8))
  # debug: timer
  # t1 <- Sys.time()
  # loop over all phenotypes
  result_ftrs = foreach(i = 1:ncol(pheno),.combine = rbind) %:%
    
    #result_ftrs = foreach(i = 1:50,.combine = rbind) %:%     # debug/test
    
    # loop over all features of requested type (ftrs, VFs, PWYs or CARDs)
    # ==========================================================================
    foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %do% { # single-threaded implementation for debug
    
    #foreach(j = 1:ncol(ftrs_transformed),.combine = rbind) %dopar% {  # parallel implementation
      #debug output/mode
      #foreach(j = 1:50,.combine = rbind) %dopar% {  #debug/test
      #print(i)
      #print(j)
      predictors = data.frame(covar[!is.na(pheno[,i]),],
                              model.matrix(
                                as.formula(paste0("~ ",colnames(pheno)[i])),data = pheno)[,-1,drop = F])
      
      cleaned_data = predictors[complete.cases(predictors),]
      rn <- rownames(cleaned_data)
      rn <- rn[rn %in% rownames(ftrs_transformed)]
      ftrs.cleaned = ftrs_transformed[rn,]
      cleaned_data = cleaned_data[rn,]
      if (nrow(cleaned_data) > 3) {
        # debug: print model
        #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data))))
        
        # make model
        s1 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)))),
          data = cleaned_data
        )
        # debug: print model
        #print(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)])))
        
        # make model with extra covariates
        s0 = lm(
          as.formula(paste("ftrs.cleaned[,j] ~ ",paste(collapse=" + ",colnames(cleaned_data)[1:ncol(covar)]))),
          data = cleaned_data
        )
        
        # compare models
        an1 = anova(s1,s0)
        output = data.frame(
          phenotype = colnames(pheno)[i],
          taxon = colnames(ftrs.cleaned)[j],
          Nsamples = nrow(cleaned_data),
          levels = if(class(pheno[,i]) == "factor") paste(collapse=":",levels(pheno[,i])) else "Not Applicable",
          levels_SampleSize = 
            if(class(pheno[,i]) == "factor" | length(table(pheno[,i]))==2) paste(collapse= ":",table(pheno[,i])) else "Not Applicable",
          effect.size = if(class(pheno[,i]) == "factor") {paste(collapse = ":",c(0,round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))])))}
          else round(digits = 5, s1$coef[grep(colnames(pheno)[i],names(s1$coef))]) ,
          
          R2 = summary(s1)$r.squared- summary(s0)$r.squared,
          F.stat = an1[2,5],
          Pvalue = an1[2,6]
        )
        #add covariates
        output
      }
    }
  # debug
  #t2 <-  Sys.time()
  print("ftrs_done")
  # debug
  #print(t2-t1)
  rownames(result_ftrs) <- NULL
  result_ftrs$FDR <- p.adjust(result_ftrs$Pvalue)
  #write.table(,file = phenoOut,row.names=F,quote = F,sep="\t")
  
  #close parallelization
  registerDoSEQ()
  # return results
  result_ftrs
}


