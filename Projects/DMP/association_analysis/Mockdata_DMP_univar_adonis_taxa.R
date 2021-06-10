# ============================================================
# By: Weersma Group, UMCG (2021)
#
# DMP Variance explained analysis,
# implementation on MOCK DATA
#
# script performs univariate ADONIS analysis for set of phenotypes
#
# Notes: 
#  - working directory must be set to appropriate path (DMP folder)
#  - this is implementation on MOCK DATA, it does not reproduce results from the DMP study
#    and is intended for code testing only
#    reproducing the study requires phenotype data that can be requested from Lifelines Biobank
#    (please see manuscript for details)
#   - this is slightly different implementation from DMP_do_univar_adonis_taxa.R code 
#     which performs only adonis run for one phenotype and is intended for use on HPC cluster    
# ==========================================================================================

# set WD: please set the path to appropriate location (DMP folder of github repo)
# ============================================
#setwd('D:/Vbox/shared/dag/git_10_06/DMP/')
setwd('.')
# ============================================

# load libraries
library(vegan)

# ==========================================

# ==========================================
# >> MAIN <<
# ==========================================
# set number of permutations (here 100 for quick demonstration run, 
# 1,000+ is appropriate for "real" run, DMP study used 20,000 
permNR <- 100

# select how many samples to use in analysis 
# (2000 is selected here for quick demonstration run)
# if set to -1, uses everything (appropriate for "real" run)
nrSamples <- 2000 

# output folder
outFolder <- "microbiome_variance/mockdata_output_adonis"
if (!dir.exists(outFolder)) {dir.create(outFolder)}

# ===================================

# load helper scripts
source('r_scripts_library/R_Microbiome_scripts.R')

# load phenotypes
print(' >> LOADING PHENOTYPES')
inPhenos <- read.table('Mock_data/MOCK_phenotypes_Fig1c.csv',sep=',',header=T)
inPhenos$DAG3_sampleID <- paste0('Sample',c(1:nrow(inPhenos)))
rownames(inPhenos) <- inPhenos$DAG3_sampleID
if (nrSamples == -1) {nrSamples <- nrow(inPhenos)}
inPhenos <- inPhenos[c(1:nrSamples),]

# load data (microbiome)
print(' >> LOADING MICROBIOME')
inDF <- read.table('heritability_analysis_v2/mock_data/DMP_mock_microbiome_taxa_filtered.csv',sep=',',header=T)

# grab species only
inMB <- filterMetaGenomeDF(inDF,presPerc = -1,minMRelAb = -1,minMedRelAb = -1,
                           rescaleTaxa = T,verbose = T,keepLevels = c("S"),keepDomains = c("Archaea","Bacteria"))
inMB$DAG3_sampleID <- paste0('Sample',c(1:nrow(inMB)))
inMB <- inMB[inMB$DAG3_sampleID %in% inPhenos$DAG3_sampleID,]

# load phenotypes
# ==========================================
# > remove IDs and other (irrelevant) phenotypes
adonisVarsTouse <- colnames(inPhenos)
toExclude = c("DAG3_sampleID","ID")
adonisVarsTouse <- adonisVarsTouse[!adonisVarsTouse %in% toExclude]

# =================================================================
# PART I: do UNIVARIATE ADONIS for all variables
# =================================================================
print (' >> STARTING ADONIS CALCULATIONS')
print (paste0('   >> NR of premut ', permNR))
print (' ============================================ ')
adonisResults <- NULL
for (i in adonisVarsTouse) {
  print (paste(' >>> ANALYSING VARIABLE <',i,'>    <<<'))
  print ('  >> collecting complete cases')
  inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(i,"DAG3_sampleID")]
  allDF <- merge(x=inPhenosOneVarID,by.x="DAG3_sampleID",y=inMB,by.y="DAG3_sampleID")
  rownames(allDF) <- allDF$DAG3_sampleID
  allDF$DAG3_sampleID <- NULL
  allDF <- allDF[complete.cases(allDF),]
  av <- allDF[[i]]
  allDF[[i]] <- NULL
  nrRows <- length(av)
  print(paste0(' Numbers: Total: ',nrRows, ', NAs: ',sum(is.na(av))))
  print ('  >> calculating B/C distance')
  inBC <- vegdist(allDF,method = "bray",parallel=4)
  print ('  >> doing adonis')
  if (length(av) < 3 | length(unique(av)) < 2) {
    print(paste0(' >> WARNING: ',i,' has no useful data!!'))
  } else {
    ad <- adonis(inBC ~ av,permutations=permNR)
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
    print (paste0('--- ',i,' DONE! ---'))
  }
}
print (' >> DONE WITH UNIVARIATE ADONIS CALCULATIONS')

adonisResults$FDR.BH=p.adjust(adonisResults$pval, method = "BH")
adonisResults$Significant="No"
adonisResults$Significant[adonisResults$FDR.BH<0.05]="Yes"
adonisResults <- adonisResults[order(adonisResults$pval),]

write.table(adonisResults,paste0(outFolder,"/mockdata_adonis_Taxa_results_table.csv"),sep=",",row.names=F)

# ==================================================================
#    ========= PART II: Multivariate ADONIS ===================
# ==================================================================
# > loads results of part I (univariate analysis)
# > selects significant variables
# > groups them by phenotype groups (diseases, anthropometrics...)
# > does multivariate ADONIS per group

# helper functions:

# maxFreqClass function 
# ================
# returns most common class
# in vector or matrix of class factor
# ex: maxFreqClass(c("A","A","B","C","C","C")) = "C" 
#     maxFreqClass(matrix(c("A","A","B","A","B","C"),2)) = "A"
maxFreqClass <- function(InVec, mult = FALSE) {
  if (!is.factor(InVec)) InVec <- factor(InVec)
  A <- tabulate(InVec)
  if (isTRUE(mult)) {
    levels(InVec)[A == max(A)]
  }
  else levels(InVec)[which.max(A)]
}

# ============ SETTINGS  ============
permNR <- 100
# type of adonis to perform: seq = sequential; sum = summary; marg = marginal
# summary adonis is adequate for group-level analysis
analysisType <- "sum"

# univariate ADONIS significance cutoff, should be set to 0.05 for real analysis
fdrCutoff <- 1 
# ====================================

# === LOAD STEP 1 RESULTS ===
adonisResults1 <- read.table(paste0(outFolder,"/mockdata_adonis_Taxa_results_table.csv"),sep=",",header=T,stringsAsFactors = F)

# === LOAD PHENOTYPE GROUPS ===
phenoGroups <- read.table('association_analysis/data/phenotype_groups.csv',sep=',',header=T)

# select significant phenotypes per group
# ===========================================
adonisRes1grp <- merge(adonisResults1,phenoGroups,by.x="Var",by.y="Phenotype")

# results
resAdonisMV <- NULL

# iterate over groups
# ===========================================
for (ug in unique(adonisRes1grp$Group.smaller)) {
  print(paste0(' >> PERFORMING MULTIVARIATE ADONIS FOR Phenotype group: ',ug))
  sigPhenos <- adonisRes1grp$Var[adonisRes1grp$Group.smaller==ug & adonisRes1grp$FDR.BH < fdrCutoff]
  print(paste0('   >> NR of significant phenotypes to include: ',length(sigPhenos)))
  if (length(sigPhenos) ==0) {
    print(' WARNING: no phenotypes are significant in Univar ADONIS, skipping the group!')
  }
  else if (length(sigPhenos) >= 1) {
    # impute phenotypes
    for (cn in colnames(inPhenos)) {
      if (cn != "DAG3_sampleID") {
        if (class(inPhenos[[cn]]) == "numeric" | class(inPhenos[[cn]])=="integer") {
          inPhenos[[cn]][is.na(inPhenos[[cn]])] <- median(inPhenos[[cn]],na.rm=T)
        } else if (class(inPhenos[[cn]]) == "factor") {
          inPhenos[[cn]][is.na(inPhenos[[cn]])] <- maxFreqClass(inPhenos[[cn]][!is.na(inPhenos[[cn]])])
        } else {
          print(paste0(cn,' = ',class(inPhenos[[cn]])))
        }
      }
    }
    
    # =================================================================
    # do adonis
    # =================================================================
    phenosList <- sigPhenos
    print (' >> STARTING ADONIS CALCULATIONS')
    print (paste0('   >> phenotypes: '))
    print (phenosList)
    print (paste0('   >> NR of premut ', permNR))
    print (' ============================================ ')
    adonisResults <- NULL
    inPhenosOneVarID <- inPhenos[,colnames(inPhenos) %in% c(phenosList,"DAG3_sampleID")]
    colnames(inMB)[colnames(inMB)=="ID"] <- "DAG3_sampleID"
    allDF <- merge(x=inPhenosOneVarID,by.x="DAG3_sampleID",y=inMB,by.y="DAG3_sampleID")
    rownames(allDF) <- allDF$DAG3_sampleID
    allDF$DAG3_sampleID <- NULL
    
    # === VARIABLE CHECK ===
    # drop vars with < X cases
    xCases <- 0
    toDrop <- c()
    for (cc in phenosList) {
      nrNA <- sum(is.na(allDF[[cc]]))
      nrTot <- length(allDF[[cc]])
      nrNonNA <- nrTot - nrNA
      if (nrNonNA < xCases) {
        print(paste0(' > WARNING: ',cc,' is ',nrNA,' NA; ',nrNonNA,' non-NA; ',nrTot,' Total samples, dropping it!'))
        toDrop <- c(toDrop,cc)
      }    
    }
    phenosList <- phenosList[!(phenosList %in% toDrop)]
    allDF <- allDF[,!(colnames(allDF) %in% toDrop)]
    # use complete cases
    allDF <- allDF[complete.cases(allDF),]
    # =============================
    av <- as.data.frame(allDF[,colnames(allDF) %in% phenosList])
    if (ncol(av) == 1) {colnames(av) <- phenosList}
    # grab taxa for Bray-Curtis calculation
    allDF <- allDF[,!(colnames(allDF) %in% phenosList)]
    print ('  >> calculating B/C distance')
    inBC <- vegdist(allDF,method = "bray")
    print ('  >> doing multivariate adonis')
    if (length(av[[1]]) < 3 | length(unique(av[[1]])) < 2) {
      print(paste0(' >> WARNING: ',i,' has no useful data!!'))
    } else {
      frm <- reformulate(termlabels=colnames(av),response='inBC')
      if (analysisType == 'seq') {
        print ('  >> running sequential adonis')
        ad <- adonis(formula=frm,data=av,permutations=permNR)
        print(ad)
      } else if (analysisType == 'sum') {
        print ('  >> running summary adonis')
        ad <- adonis2(formula=frm,data=av,permutations=permNR,by=NULL)
        print(ad)
      } else if (analysisType == 'marg') {
        print ('  >> running marginal adonis')
        ad <- adonis2(formula=frm,data=av,permutations=permNR,by="margin")
        print(ad)
      }
      
      inR <- ad
      R2 <- ad$R2[1]
      pV <- ad$`Pr(>F)`[1]
      resAdonisMV <- rbind.data.frame(resAdonisMV,
                                      data.frame(Group=ug,R2=R2,pV=pV))
      print (paste0('--- DONE with group ',ug))
    }
  }
}
print(' >> DONE WITH ALL GROUPS << ')
resAdonisMV

# make plot
inDFss <- resAdonisMV[,c("Group","R2")]
inDFss$Group <- factor(as.character(inDFss$Group),levels = inDFss$Group[order(inDFss$R2,decreasing = F)])
inDFss <- inDFss[order(inDFss$R2, decreasing = T),]
inDFss$csumR2 <- inDFss$R2[1]
for (rn in c(2:nrow(inDFss))) {inDFss$csumR2[rn] <- inDFss$csumR2[rn-1] + inDFss$R2[rn]}
inDFss$ypos <- inDFss$csumR2
for (rn in c(1:nrow(inDFss))) {inDFss$ypos[rn] <- inDFss$ypos[rn] - inDFss$R2[rn]/2}
inDFss$R2lbl <- paste0(format(inDFss$R2*100,digits = 1),'%')
inDFss$Data <- "Taxa"
rownames(inDFss) <- inDFss$Group

g <- ggplot(inDFss,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col() + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + ylim(0,0.015) + theme_classic()
g


