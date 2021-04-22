# ==================================================
# By: Weersma Group, UMCG (2020)
# DMP project, training and testing of
# models for prediction of health and diseases
# ==================================================


# this script generates some random data for testing disease prediction models. 
# parameters for generating are taken from real data.  
# script requires real "taxa" and "pathways" table which can be obtained on EGA

ids = paste("Sample",1:2000,sep = "")

covariates_prediction = data.frame(
  row.names=ids,
  Age = round(rnorm(2000,mean = 48,sd =14.77)),
  Sex = rbinom(n = 2000,size = 1,prob = 0.4204435),
  BMI = round(digits = 1,rnorm(2000,mean = 25.556,sd = 4.406)))

write.table(covariates_prediction,file = "mock_data/covariates.txt",sep="\t",quote = F)

dis2_prediction = data.frame(
                                         row.names = ids,
                         MED.DISEASES.Blood.Anemia = rbinom(n = 2000,size = 1,prob = 0.01279240),
                     MED.DISEASES.Blood.Thrombosis = rbinom(n = 2000,size = 1,prob = 0.14485867),
                           MED.DISEASES.Cancer.Any = rbinom(n = 2000,size = 1,prob = 0.05653021),
MED.DISEASES.Cardiovascular.Arrythmia.MedDiagnosed = rbinom(n = 2000,size = 1,prob = 0.07894737),
       MED.DISEASES.Cardiovascular.Colesterol.high = rbinom(n = 2000,size = 1,prob = 0.12524366),
          MED.DISEASES.Cardiovascular.Heart.Attack = rbinom(n = 2000,size = 1,prob = 0.01327973),
MED.DISEASES.Cardiovascular.Heart.Failure.Disorder = rbinom(n = 2000,size = 1,prob = 0.01522904),
   MED.DISEASES.Cardiovascular.Heartrate.complains = rbinom(n = 2000,size = 1,prob = 0.22161306),
          MED.DISEASES.Cardiovascular.Hypertension = rbinom(n = 2000,size = 1,prob = 0.20845517),
                 MED.DISEASES.Endocrine.DiabetesT2 = rbinom(n = 2000,size = 1,prob = 0.02205166),
       MED.DISEASES.Gastrointestinal.Stomach.Ulcer = rbinom(n = 2000,size = 1,prob = 0.03691520),
               MED.DISEASES.Hepatologic.Gallstones = rbinom(n = 2000,size = 1,prob = 0.04300682),
                           MED.DISEASES.Mental.Any = rbinom(n = 2000,size = 1,prob = 0.13949805),
                      MED.DISEASES.Mental.Burn.Out = rbinom(n = 2000,size = 1,prob = 0.05250975),
                    MED.DISEASES.Mental.Depression = rbinom(n = 2000,size = 1,prob = 0.04191033),
                 MED.DISEASES.Mental.Other.anxiety = rbinom(n = 2000,size = 1,prob = 0.02631579),
                MED.DISEASES.Mental.Panic.disorder = rbinom(n = 2000,size = 1,prob = 0.02546296),
       MED.DISEASES.Neurological.Dizziness.Falling = rbinom(n = 2000,size = 1,prob = 0.01461988),
     MED.DISEASES.Neurological.Mental.Fibromyalgia = rbinom(n = 2000,size = 1,prob = 0.03252924),
                MED.DISEASES.Neurological.Migraine = rbinom(n = 2000,size = 1,prob = 0.17202729),
 MED.DISEASES.Other.Autoimmune.Rheumatoid.Artritis = rbinom(n = 2000,size = 1,prob = 0.01985867),
               MED.DISEASES.Other.Chronic.cystitis = rbinom(n = 2000,size = 1,prob = 0.01632554),
MED.DISEASES.Other.Chronic.Inflammation.Throatnose = rbinom(n = 2000,size = 1,prob = 0.04958577),
        MED.DISEASES.Other.Chronic.Muscle.Weakness = rbinom(n = 2000,size = 1,prob = 0.05275341),
                MED.DISEASES.Other.Fractures.Other = rbinom(n = 2000,size = 1,prob = 0.02083333),
                   MED.DISEASES.Other.Incontinence = rbinom(n = 2000,size = 1,prob = 0.02107700),
                  MED.DISEASES.Other.Kidney.Stones = rbinom(n = 2000,size = 1,prob = 0.01230507),
                 MED.DISEASES.Other.Osteoarthritis = rbinom(n = 2000,size = 1,prob = 0.12426901),
                   MED.DISEASES.Other.Osteoporosis = rbinom(n = 2000,size = 1,prob = 0.01900585),
                            MED.DISEASES.Other.RSI = rbinom(n = 2000,size = 1,prob = 0.01803119),
          MED.DISEASES.Pulmonary.Autoimmune.Asthma = rbinom(n = 2000,size = 1,prob = 0.04653996),
                       MED.DISEASES.Pulmonary.COPD = rbinom(n = 2000,size = 1,prob = 0.03204191),
    MED.DISEASES.Skin.Autoimmune.Atopic.dermatitis = rbinom(n = 2000,size = 1,prob = 0.14741715),
            MED.DISEASES.Skin.Autoimmune.Psoriasis = rbinom(n = 2000,size = 1,prob = 0.02655945),
          MED.DISEASES.Skin.Autoimmune.Severe.acne = rbinom(n = 2000,size = 1,prob = 0.02570663),
       MED.DISEASES.Gastrointestinal.Rome3_IBS.Any = rbinom(n = 2000,size = 1,prob = 0.10804521),
                     MED.DISEASES.None.No.Diseases = rbinom(n = 2000,size = 1,prob = 0.35782164)
)
write.table(dis2_prediction, file = "Mock_data/diseases.txt",sep="\t",quote = F)

#taxa & pathways
taxa_random = taxa
taxa_random = taxa_random[1:2000,]
rownames(taxa_random) = ids
for(i in 1:ncol(taxa_random)){
  taxa_random[,i] = rbinom(2000, size = 1000000,prob = mean(taxa[,i]))/1000000
}
pathways_random = pathways
pathways_random = pathways_random[1:2000,]
rownames(pathways_random) = ids
for(i in 1:ncol(pathways_random)){
  pathways_random[,i] = rbinom(2000, size = 1000000,prob = mean(taxa[,i]))/1000000
}

write.table(taxa_random, file = "Mock_data/taxa.txt",sep="\t",quote = F)
write.table(pathways_random, file = "Mock_data/pathways.txt",sep="\t",quote = F)
