# ======================================================================================================
# ======================================================================================================
# By: Weersma Group, UMCG, 2020

# Runs association analysis, including age, sex and BMI, Bistol stool scale of the faecal sample,
# and technical factors (DNA concentration, sequencing read depth, sequencing batch and sampling season)
# into the model
# ======================================================================================================
# ======================================================================================================
setwd('#placeholder#')
# association analysis is done by the following code:
source('DMP_ScriptsAssociationWithCorrections.R')

# define covariates
# NOTE: Age, Sex, Sequencing Batch, Season of sample collection and DNA concentration are 
#       included by default [see dag3ScriptsAssociationWithCorrections.R]

covs <- c("ANTHRO.BMI","META.POOP.BristolMean","META.DNA.postclean.reads")

# loop over data layer, run associations for each:
for (dataType in c("pathways","taxa","CARD","VFDB")) {
  phenoOut <- paste0("associations_",dataType,"_corrected_BMI_Bristol.txt")
  res <- dag3AssociationsWithCorrections(dataType = dataType,
                                         dataPath = '#placeholder#',
                                         covariates = covs)
  write.table(res,paste0('Association.analysis.output/',phenoOut),sep=',',row.names = F)
}
