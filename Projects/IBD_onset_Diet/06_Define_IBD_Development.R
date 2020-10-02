library(tidyverse)
Data = read_tsv("TABLES/Data_v2.tsv")
Data %>% mutate(Condition = ifelse(Disease=="HEALTH100E1", "UC_d", ifelse(Disease=="HEALTH72C2" , "UC", ifelse(Disease=="HEALTH100E2", "CD_d", ifelse(Disease == "HEALTH72C3", "CD", "" ))))) -> Data


Data %>% filter(Condition == "UC") %>% filter(`Baseline assessment, questionnaire 1` == "Yes") -> UC_patients
Data %>% filter(Condition == "CD") %>% filter(`Baseline assessment, questionnaire 1` == "Yes") -> CD_patients
#UC
Data %>% filter(Condition == "UC_d") %>% filter(! ID %in% c(UC_patients$ID, CD_patients$ID) ) %>% filter(`Second assessment, questionnaire 1` == "Yes" | `Second assessment, questionnaire 2` == "Yes" | `Follow-up questionnaire (1B)` == "Yes" | `Follow-up questionnaire (1C)` =="Yes") -> UC_developed
#CD
Data %>% filter(Condition == "CD_d") %>% filter(! ID %in% c(UC_patients$ID, CD_patients$ID)) %>% filter(`Second assessment, questionnaire 1` == "Yes" | `Second assessment, questionnaire 2` == "Yes" | `Follow-up questionnaire (1B)` == "Yes" | `Follow-up questionnaire (1C)` =="Yes") -> CD_developed
#Healthy
REMOVE = c("Disease", "Condition", "Baseline assessment, questionnaire 1",  "Second assessment, questionnaire 1",  "Second assessment, questionnaire 2","Follow-up questionnaire (1B)","Follow-up questionnaire (1C)" )
Data %>% filter(! ID %in%  c(UC_patients$ID, UC_developed$ID, CD_patients$ID, CD_developed$ID)) %>% spread(Disease, c(`Follow-up questionnaire (1B)`)) %>% select(-one_of(REMOVE)) -> Healthy
#Unknown
intersect(UC_patients$ID, CD_patients$ID) -> Unknown_patients_ID
Data %>% filter(Condition == "UC") %>% filter(ID %in% Unknown_patients_ID) %>% select(-one_of(REMOVE)) -> Unknown_patients
#Remove unknown from the other
UC_patients %>% filter(!ID %in% Unknown_patients_ID) %>% select(-one_of(REMOVE)) -> UC_patients
CD_patients %>% filter(!ID %in% Unknown_patients_ID) %>% select(-one_of(REMOVE)) -> CD_patients
#Unknown developmnet
intersect(UC_developed$ID, CD_developed$ID) -> Unknown_patients_d_ID
Data %>% filter(Condition == "UC_d") %>% filter(ID %in% Unknown_patients_d_ID) %>% select(-one_of(REMOVE)) -> Unknown_patients_developed
#Remove unknown from the other
UC_developed %>% filter(!ID %in% Unknown_patients_d_ID) %>% select(-one_of(REMOVE)) -> UC_developed
CD_developed %>% filter(!ID %in% Unknown_patients_d_ID) %>% select(-one_of(REMOVE)) -> CD_developed



###Numbers
print("Healthy")
dim(Healthy)[1]
print("UC patients BL")
dim(UC_patients)[1]
print("Develped UC")
dim(UC_developed)[1]
print("CD patients BL")
dim(CD_patients)[1]
print("Developed CD")
dim(CD_developed)[1]
print("Both UC and CD BL")
dim(Unknown_patients)[1]
print("Developed both UC and CD")
dim(Unknown_patients_developed)[1]

print("Keeping only Healthy and people who developed CD or UC (not both)")

###TEST CD and UC
Healthy %>% mutate(Status_CD = 0 , Status_UC= 0) -> Healthy
CD_developed %>% mutate(Status_CD = 1, Status_UC =NA) -> Dev_CD
UC_developed %>%  mutate(Status_UC = 1, Status_CD = NA) -> Dev_UC


Healthy %>% select(colnames(Dev_CD)) %>% unique() -> Healthy
rbind(Healthy, Dev_CD) -> To_test
rbind(To_test, Dev_UC) -> To_test

print("Dataset")
To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())

print("Adding Scores")

#Now need to add Diet Quality info as the regressor
read_tsv("TABLES/DIET_SCORES.tsv") -> Diet_scores
colnames(Diet_scores)[1]  = "ID"
left_join(To_test, Diet_scores, by="ID") -> To_test
To_test = To_test[!duplicated(To_test$ID),]

print("Adding PCs")
#Add PCs
read_tsv("TABLES/PCs.tsv") -> PCs
left_join(To_test, PCs, by="ID") -> To_test

print("Adding PCs or robust PCA")
read_tsv("TABLES/PCs_robust.tsv") -> PCs_r
left_join(To_test, PCs_r, by="ID") -> To_test

#Add Food items
read_tsv("TABLES/Table_participant_GDAG.tsv_PCA_input") -> Food_groups
col_remove = str_split(str="SumOfkcal       SumOfkJ SumOfeiwittot   SumOfeiwitplant SumOfeiwitdier  SumOfvettot     SumOfkhtot      SumOfmodis      SumOfpolys     SumOfalcohol     SumOffree_sug   SumOfadd_sug    SumOfglucose    SumOffructose   SumOflactose    SumOfmaltose    SumOfsucrose    SumOfGlyc_ind   SumOfGlyc_load  leeftijd  GESLACHT        GEWICHT LENGTE  BMI     BMR     EI_BMR",pattern=" ")
col_remove = unique(col_remove[[1]][col_remove[[1]] != ""])
Food_groups %>% select(-col_remove) %>% mutate(ID = Participant) %>% select(-Participant)  -> Food_groups

left_join(To_test, Food_groups) -> To_test


print("Saving complete dataset in TABLES/Dataset_IBD_test.tsv")
write_tsv(To_test,"TABLES/Dataset_IBD_test.tsv")


