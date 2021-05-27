# ===================================
# ==================================
# BY: Weersma Group (UMCG)
# DMP association plots (Fig 4/a; Fig 3/b, Fig S7-S10)
# ==================================
# ===================================
# function for trimming phenotype names
trimPhenos <- function(inPh) {
  inPh <- gsub('META\\.','',inPh)
  inPh <- gsub('ANTHRO\\.','',inPh)
  inPh <- gsub('SOCIOEC\\.','',inPh)
  inPh <- gsub('MED\\.','',inPh)
  inPh <- gsub('MEDS\\.','',inPh)
  inPh <- gsub('SOCIOECONOMIC\\.','',inPh)
  inPh <- gsub('INCOME\\.','',inPh)
  inPh <- gsub('WORK\\.','',inPh)
  inPh <- gsub('DISEASES\\.','',inPh)
  inPh <- gsub('Cardiovascular\\.','',inPh)
  inPh <- gsub('BLOOD\\.','',inPh)
  inPh <- gsub('EXP\\.','',inPh)
  inPh <- gsub('EARLYLIFE\\.','',inPh)
  inPh <- gsub('DIET\\.','',inPh)
  inPh <- gsub('Scores\\.','',inPh)
  inPh <- gsub('BIRTH\\.','',inPh)
  inPh <- gsub('DEMOGRAPHICS\\.','',inPh)
  inPh <- gsub('Scores\\.','',inPh)
  inPh <- gsub('HEALTH\\.','',inPh)
  inPh <- gsub('RAND\\.','',inPh)
  inPh <- gsub('Autoimmune\\.','',inPh)
  inPh
}
# function for linking phenotypes to phenotype groups
remapPhenos <- function(inPh) {
  out <- c()
  cc = 0
  for (inP in inPh) {
    cc = cc + 1
    if (cc %% 5000 == 0) {
      print(cc)
    }
    mm <- grep(paste0('^',inP,'$'),inRename$Phenotype)
    if (length(mm) == 1) {
      inPh2 <- inRename$Group.renamed[mm]
    } else {
      inPh2 <- inP 
    }
    out <- c(out,inPh2)
  }
  out
}

# =========== MAIN ============
# > set wd (note: make sure this is set to appropriate path)
setwd('D:/Vbox/shared/dag/git_14_05/DMP')

#  > plotting scripts live in R_Misc.R
source('r_scripts_library/R_Misc.R')
#  > other microbiome-related scripts live in R_Microbiome_scripts.R
source('r_scripts_library/R_Microbiome_scripts.R')
# regular R libraries
library(plyr)

# what to plot? (plotters include options for "Pvalue",	"PadjBH","PadjBonferroni","F.stat","R2","effect.size","z.score")
statToPlot <- "z.score"

# file with phenotype groups
inRename <- read.table('association_analysis/data/phenotype_groups.csv',sep=',',header=T,stringsAsFactors = F)

# =========================================================================
# =================== HEALTHY MICROBIOME PLOTS (Fig 4/a) ======================
# =========================================================================
print(' >> loading data ...')
inData <- read.table('association_analysis/data/DAG3_v27d_supplementary_tables_S3B.txt',sep='\t',header=T,stringsAsFactors = F)
inDF <- inData
inDF$effect.size.asInteger <- inDF$effect.size
inDF$PadjBH <- inDF$FDR
inDF$PadjBonferroni <- inDF$FDR
print(' >> grouping and cleaning phenotypenames (takes a bit) ...')
inDF$phenotype <- remapPhenos(inDF$phenotype)

phenosToPlot <- c(
                  "MED.HEALTH.RAND.Health.General.asInteger",
                  "MED.DISEASES.Endocrine.DiabetesT2",
                  "MED.DISEASES.Gastrointestinal.Rome3_IBS.Any",
                  "MED.DISEASES.Gastrointestinal.Discomfort.Pain",
                  "MED.DISEASES.Gastrointestinal.Discomfort.Unwell",
                  "MED.DISEASES.Pulmonary.COPD",
                  "MED.DISEASES.Neurological.Mental.Fibromyalgia",
                  "MED.DISEASES.Gastrointestinal.Discomfort.Nausea",
                  "MED.DISEASES.Cardiovascular.Heartrate.complains",
                  "MED.DISEASES.Mental.Burn.Out",
                  "MED.DISEASES.Other.Chronic.Muscle.Weakness",
                  "MED.DISEASES.Hepatologic.Gallstones",
                  "MED.DISEASES.Skin.Autoimmune.Atopic.dermatitis",
                  "MED.DISEASES.Cardiovascular.Colesterol.high",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.UC",
                  "MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD",
                  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD",
                  "MED.DISEASES.Pulmonary.Autoimmune.Asthma",
                  "MED.DISEASES.Other.Chronic.Inflammation.Throatnose",
                  "MED.DISEASES.Cardiovascular.Arrythmia.MedDiagnosed",
                  "MED.DISEASES.Gastrointestinal.Discomfort.Burping",
                  "MED.DISEASES.Cancer.AnyNonBasal",
                  "MED.DISEASES.Other.Autoimmune.Rheumatoid.Artritis",
                  "MED.DISEASES.Other.Osteoarthritis",
                  "MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM",
                  "MED.DISEASES.Gastrointestinal.Discomfort.Flatulence",
                  "MED.DISEASES.Cardiovascular.Hypertension",
                  "MED.DISEASES.Mental.Depression",
                  "MED.DISEASES.Gastrointestinal.Stomach.Ulcer",
                  "MED.DISEASES.Mental.Bipolar",
                  "MED.DISEASES.Mental.Eating.disorder",
                  "MED.DISEASES.Neurological.Migraine",
                  "MED.DISEASES.Other.Chronic.cystitis",
                  "MED.DISEASES.Blood.Thrombosis",
                  "MED.DISEASES.Gastrointestinal.FGID.Bloating"
)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)

# SPECIES plot (sorted) [Fig 4/a]
# ======================================
# make the plot
nrFeatures = 40
phenosToPlot <- phenosToPlot[order(phenosToPlot)]
print(paste0(' >> PLOTTING Healthy microbiome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = statToPlot,
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = F, # do not cluster phenotypes
                                     sortPhenos = T,
                                     retData = T) 
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig_4a_healthy_microbiome_species_',statToPlot,'.png'),width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig_4a_healthy_microbiome_species_',statToPlot,'_datatable_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig_4a_healthy_microbiome_species_',statToPlot,'_datatable_directions.csv'),sep=',',row.names = F)

# ======================================
# JOINED EXPOSOME PLOTS, Fig 3/b
# ======================================
phenosToPlot <- c(
  # ======================
  # diet
  # ======================
  "EXP.DIET.Amounts.GlycLoadTotal",
  "EXP.DIET.Amounts.KcalTotal",
  "EXP.DIET.ECorrected.AlchAll",
  "EXP.DIET.ECorrected.CarbsAll",
  "EXP.DIET.ECorrected.FatAll",
  "EXP.DIET.ECorrected.ProtAnimal",
  "EXP.DIET.ECorrected.ProtPlant",
  "EXP.DIET.Probiotics",
  "EXP.DIET.Scores.LLDietScore",
  # ==========================
  # baby phenotypes
  # ==========================
  "EXP.EARLYLIFE.Preg.Mother.Smoking",
  "EXP.EARLYLIFE.PretermBorn",
  "EXP.EARLYLIFE.Birth.Mode",
  "EXP.BIRTH.Breastfed.6.M.Plus",
  # ==========================
  # childhood phenotypes
  # ==========================
  "EXP.EARLYLIFE.LIVINGPLACE.child.1.4",
  "EXP.EARLYLIFE.Parent.smoker.childhood",
  "EXP.EARLYLIFE.PETS.child.0.15.Any",
  # ==========================
  # current exposome
  # ==========================
  "EXP.GREENSPACE.NDVI.100m",
  "EXP.PETS.now.Any",
  "EXP.POLLUTION.NO2",
  "EXP.POLLUTION.PM2.5",
  "SOCIOECONOMIC.BUURT.Urbanicity",
  "EXP.SMOKING.Passive",
  "EXP.SMOKING.Smoker.Now",
  "EXP.SMOKING.Smoker.OneYear.Ever",
  "EXP.SMOKING.Smoker.stopped",
  # ==============================
  # socioeconomics
  # ==============================
  "SOCIOEC.INCOME.Income.Month",
  "SOCIOEC.STATUS.Retired",
  "SOCIOEC.WORK.Housewife.husband",
  "SOCIOEC.WORK.Paidwork.any",
  "SOCIOECONOMIC.BUURT.Highincome.Pop.Prop",
  # ========= health for comparison =======
  "MED.HEALTH.RAND.Health.General.asInteger"
)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)
#phenosToPlot <- trimPhenos(phenosToPlot)

# ==== Species Plot, sorted ===
nrFeatures = 40
#phenosToPlot <- phenosToPlot[order(phenosToPlot)]
print(paste0(' >> PLOTTING Exposome vs Species'))
#phenosToPlot <- as.character(phenosToPlot)
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = statToPlot,
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = F,
                                     clusterFeatures = T,
                                     sortPhenos = F,
                                     retData = T,
                                     keepOrder = T,
                                     flipCoords = F) # do not cluster phenotypes
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig3b_Exposome_vs_Species_',statToPlot,'.png'),height = 1400+70*nrPhenos,width = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig3b_Exposome_vs_Species_',statToPlot,'_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig3b_Exposome_vs_Species_',statToPlot,'_directions.csv'),sep=',',row.names = F)

# =========================================================================================
# =================== HEALTHY MICROBIOME + DRUGS PLOTS [Fig S/7] ======================
# =========================================================================================

phenosToPlot <- c(
  "MED.HEALTH.RAND.Health.General.asInteger",
  "MED.DISEASES.Endocrine.DiabetesT2",
  "MED.DISEASES.Gastrointestinal.Rome3_IBS.Any",
  "MED.DISEASES.Gastrointestinal.Discomfort.Pain",
  "MED.DISEASES.Gastrointestinal.Discomfort.Unwell",
  "MED.DISEASES.Pulmonary.COPD",
  "MED.DISEASES.Neurological.Mental.Fibromyalgia",
  "MED.DISEASES.Gastrointestinal.Discomfort.Nausea",
  "MED.DISEASES.Cardiovascular.Heartrate.complains",
  "MED.DISEASES.Mental.Burn.Out",
  "MED.DISEASES.Other.Chronic.Muscle.Weakness",
  "MED.DISEASES.Hepatologic.Gallstones",
  "MED.DISEASES.Skin.Autoimmune.Atopic.dermatitis",
  "MED.DISEASES.Cardiovascular.Colesterol.high",
  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.UC",
  "MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeD",
  "MED.DISEASES.Gastrointestinal.Autoimmune.IBD.CD",
  "MED.DISEASES.Pulmonary.Autoimmune.Asthma",
  "MED.DISEASES.Other.Chronic.Inflammation.Throatnose",
  "MED.DISEASES.Cardiovascular.Arrythmia.MedDiagnosed",
  "MED.DISEASES.Gastrointestinal.Discomfort.Burping",
  "MED.DISEASES.Cancer.AnyNonBasal",
  "MED.DISEASES.Other.Autoimmune.Rheumatoid.Artritis",
  "MED.DISEASES.Other.Osteoarthritis",
  "MED.DISEASES.Gastrointestinal.Rome3_IBS.TypeM",
  "MED.DISEASES.Gastrointestinal.Discomfort.Flatulence",
  "MED.DISEASES.Cardiovascular.Hypertension",
  "MED.DISEASES.Mental.Depression",
  "MED.DISEASES.Gastrointestinal.Stomach.Ulcer",
  "MED.DISEASES.Mental.Bipolar",
  "MED.DISEASES.Mental.Eating.disorder",
  "MED.DISEASES.Neurological.Migraine",
  "MED.DISEASES.Other.Chronic.cystitis",
  "MED.DISEASES.Blood.Thrombosis",
  "MED.DISEASES.Gastrointestinal.FGID.Bloating"
)

phDrugs <- c(
  "PPIs_ATC_A02BC",
  "Biguanides_ATC_A10BA",
  "Laxatives_osmotic_ATC_A06AD",
  "Antibacterials_ATC_J01",
  "Beta_blokkers_ATC_C07AB",
  "Sulfonylureumderivaties_ATC_A10BB",
  "Vitamine_B_ATC_A11DA",
  "Corticosteroids_Eyedrops_ATC_S03A",
  "Salycilates_ATC_B01AC",
  "Intestinal_antiinflammatory_agents_ATC_A07E",
  "Laxatives_volume_increasing_ATC_A06AC",
  "Seratonine_uptake_inhibitors_selective_ATC_N06AB",
  "Statines_ATC_C10AA",
  "Other_Alimentary_Metabolism_ATC_A",
  "Tricyclic_antidepressants_ATC_N06AA",
  "Antipsychotics_atypical_ATC_N05AH",
  "ACE_inhibitors_ATC_C09AA",
  "Fluorochinolons_ATC_C08DA",
  "Calcium_antagonists_ATC_C08CA",
  "Corticosteroids_Oral_ATC_H02"
)

phenosToPlot <- c(phenosToPlot,phDrugs)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)
inDF$phenotype <- trimPhenos(inDF$phenotype)

# === SPECIES PLOT [Fig S7] ===
# save features table
# make actual plot
nrFeatures = 50
print(paste0(' >> PLOTTING Healthy microbiome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = statToPlot,
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = F,
                                     nonLogLimit = 10,
                                     retData = T) # do not cluster phenotypes
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig_S7_microbiome_drug_diseases_species_',statToPlot,'.png'),width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig_S7_microbiome_drug_diseases_species_',statToPlot,'_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig_S7_microbiome_drug_diseases_species_',statToPlot,'_directions.csv'),sep=',',row.names = F)

# =========================================================================
# =================== Early-life exposome & health (Supplement S8) ==
# =========================================================================
phenosToPlot <- c(
  "EXP.EARLYLIFE.LIVINGPLACE.child.1.4",
  "EXP.EARLYLIFE.Parent.smoker.childhood",
  "EXP.BIRTH.Breastfed.6.M.Plus",
  "EXP.EARLYLIFE.Birth.wt.g",
  "EXP.EARLYLIFE.Birth.lt",
  "EXP.BIRTH.Breastfed.1.3.M",
  "EXP.BIRTH.Breastfed.0.4.W",
  "EXP.EARLYLIFE.Birth.Mode",
  "EXP.BIRTH.Breastfed.3.6.M",
  "EXP.EARLYLIFE.PretermBorn",
  "EXP.EARLYLIFE.Preg.Mother.Smoking",
  "EXP.EARLYLIFE.PETS.child.0.15.Cat",
  "EXP.EARLYLIFE.PETS.child.0.15.Any",
  "EXP.EARLYLIFE.PETS.child.0.15.Dog"
#  "MED.HEALTH.RAND.Health.General.asInteger"
)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)
#inDF$phenotype <- trimPhenos(inDF$phenotype)

# make plot
nrFeatures = 50
print(paste0(' >> PLOTTING Current Early-life exposome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = statToPlot,
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = T,
                                     retData = T) # do not cluster phenotypes
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig_S8_Earlylife_exposome_microbiome_',statToPlot,'.png'),width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig_S8_Earlylife_exposome_microbiome_',statToPlot,'_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig_S8_Earlylife_exposome_microbiome_',statToPlot,'_directions.csv'),sep=',',row.names = F)


# =========================================================================
# =================== Current exposome & health (Supplement S9) ==
# =========================================================================
phenosToPlot <- c(
  "EXP.GREENSPACE.NDVI.100m",
  "EXP.PETS.now.Any",
  "EXP.PETS.now.Cat",
  "EXP.PETS.now.Dog",
  "EXP.POLLUTION.NO2",
  "EXP.POLLUTION.PM2.5",
  "EXP.SMOKING.Passive",
  "EXP.SMOKING.Smoker.Now",
  "EXP.SMOKING.Smoker.OneYear.Ever",
  "EXP.SMOKING.Smoker.stopped",
  "SOCIOECONOMIC.BUURT.Urbanicity",
  "MED.HEALTH.RAND.Health.General.asInteger"
)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)
#inDF$phenotype <- trimPhenos(inDF$phenotype)

# make plot
nrFeatures = 50
print(paste0(' >> PLOTTING Current Exposome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = statToPlot,
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = T,
                                     retData = T) # do not cluster phenotypes
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig_S9_Exposome_microbiome_',statToPlot,'.png'),width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig_S9_Exposome_microbiome_',statToPlot,'_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig_S9_Exposome_microbiome_',statToPlot,'_directions.csv'),sep=',',row.names = F)


# =========================================================================
# =================== DIET & health (Supplement S10) ==
# =========================================================================
phenosToPlot <- c(
  "EXP.DIET.Amounts.GlycLoadTotal",
  "EXP.DIET.Amounts.KcalTotal",
  "EXP.DIET.ECorrected.AlchAll",
  "EXP.DIET.ECorrected.CarbsAll",
  "EXP.DIET.ECorrected.CarbsFructose",
  "EXP.DIET.ECorrected.CarbsGlucose",
  "EXP.DIET.ECorrected.CarbsLactose",
  "EXP.DIET.ECorrected.CarbsMaltose",
  "EXP.DIET.ECorrected.CarbsPoly",
  "EXP.DIET.ECorrected.CarbsSuchrose",
  "EXP.DIET.ECorrected.FatAll",
  "EXP.DIET.ECorrected.ProtAll",
  "EXP.DIET.ECorrected.ProtAnimal",
  "EXP.DIET.ECorrected.ProtPlant",
  "EXP.DIET.Probiotics",
  "EXP.DIET.Ratios.CarbToFat",
  "EXP.DIET.Ratios.CarbToProtein",
  "EXP.DIET.Scores.LLDietScore",
  "EXP.DIET.Scores.LLProtein",
  "EXP.DIET.Scores.LowCarb",
  "EXP.DIET.Scores.ProtPlantScore",
  # ========= health for comparison =======
  "MED.HEALTH.RAND.Health.General.asInteger"
)
nrPhenos <- length(phenosToPlot)
phenosToPlot <- remapPhenos(phenosToPlot)
#inDF$phenotype <- trimPhenos(inDF$phenotype)

# make plot
nrFeatures = 50
print(paste0(' >> PLOTTING Current Exposome vs Species'))
phm <- plotAssociationsDag3HeatmapV2(inData = inDF, # data
                                     phenosToPlot = phenosToPlot, # what is plotted
                                     statToPlot = "z.score", # plot z-score (rather then p-value)
                                     featuresToPlot = "S", # plot genera (All plots everthing, Taxa plots all taxa, S/O/F/... plots appropriate levels)
                                     nrFeaturesToPlot = nrFeatures, # how many top shared features to plot
                                     nrColors = 51, # increase number of colors for gradient
                                     clusterPhenos = T,
                                     retData = T) # do not cluster phenotypes
save_pheatmap_png(phm[[3]],filename = paste0('association_analysis/plots/Fig_S10_Diet_microbiome_',statToPlot,'.png'),width = 1400+70*nrPhenos,height = 70*nrFeatures+1800,res = 300)
# save data table for plot
write.table(phm[[1]],file = paste0('association_analysis/plots/Fig_S10_Diet_microbiome_',statToPlot,'_numbers.csv'),sep=',',row.names = F)
write.table(phm[[2]],file = paste0('association_analysis/plots/Fig_S10_Diet_microbiome_',statToPlot,'_directions.csv'),sep=',',row.names = F)