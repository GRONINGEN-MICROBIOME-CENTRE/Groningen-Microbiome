library(tidyverse)

#All Health
Health_info = read_tsv("TABLES/Data_v2.tsv")
print("Samples with health information")
length(unique(Health_info$ID))
# IBD baseline, IBD devm, healthy
print("From those:")
print("Healthy 164666")
print("UC patients BL 813")
print("Develped UC 238")
print("CD patients BL 465")
print("Developed CD 94")
print("Both UC and CD BL 36")
print("Developed both UC and CD 25")
#Dataset of samples to test
To_test = read_tsv("TABLES/Dataset_IBD_test.tsv")
print("Keeping Only UC/CD developers and Healthy")
To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())
print("Availability of Cov")
To_test  %>% gather(Covariate, Value,  Age_n:Smoking_n) -> Covariates
Covariates %>% group_by(Covariate, is.na(Value)) %>% summarise(n())
#All diet
Diet_info = read_tsv("TABLES/Table_participant_GDAG.tsv_script_eva_input")
print("Samples with Diet information")
dim(Diet_info)[1]



#Filtered info
read_tsv("TABLES/Samples_BMR_filtered.tsv") -> To_keep
print("Filtering samples based on cal intake")
dim(To_keep)[1]
print("Applied to health data")
To_test %>% filter(! ID %in% To_keep$Participant) %>% group_by(Status_UC, Status_CD) %>% summarise(n())
print("Remaining after filtering")
To_test %>% filter(ID %in% To_keep$Participant) -> To_test
To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())
##How many with PCA
PCA = read_tsv("TABLES/PCs.tsv")
PCA_robust = read_tsv("TABLES/PCs_robust.tsv")
print("Participants in PCA")
PCA %>% drop_na() %>% dim()
print("Paticipants in robust PCA")
PCA_robust %>% drop_na() %>% dim()
#Diet scores
print("Different scores")
Scores = read_tsv("TABLES/DIET_SCORES.tsv")
Scores %>% gather(Score_t, Value, LL:mMED) -> Scores
Scores %>% group_by(Score_t, is.na(Value)) %>% summarise(n())

print("Covariates after caloric filtering")
To_test  %>% gather(Covariate, Value,  Age_n:Smoking_n) -> Covariates
Covariates %>% group_by(Covariate, is.na(Value)) %>% summarise(n())

#Regressions
print("Data available for regressions - Removed if missing cov or regressors")
print("UC regression")
To_test %>% select(-Status_CD) %>% drop_na() %>% group_by(Status_UC) %>% summarise(n())
print("CD regression")
To_test %>% select(-Status_UC) %>% drop_na() %>% group_by(Status_CD) %>% summarise(n())


