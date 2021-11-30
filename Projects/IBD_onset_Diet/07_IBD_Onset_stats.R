library(tidyverse)
wideScreen <- function(howWide=160) {
  options(width=as.integer(howWide))
}
wideScreen()

#1. Get samples with disase info
#ID      colitisulcerosa_followup_adu_q_1_2      crohns_followup_adu_q_1_2       colitisulcerosa_followup_adu_q_1_3      crohns_followup_adu_q_1_3
IBD_developers = read_tsv("IBD_developers.tsv") %>% mutate(Status_UC = ifelse(colitisulcerosa_followup_adu_q_1_2==1 | colitisulcerosa_followup_adu_q_1_3==1 | colitisulcerosa_followup_adu_q_1_1.b==1 | colitisulcerosa_followup_adu_q_1_1.c==1  , 1, 0), Status_CD=ifelse(crohns_followup_adu_q_1_2==1 | crohns_followup_adu_q_1_3==1 |crohns_followup_adu_q_1_1.b==1 | crohns_followup_adu_q_1_1.c==1 , 1 , 0) ) %>% select(-c(colitisulcerosa_followup_adu_q_1_2, crohns_followup_adu_q_1_2, colitisulcerosa_followup_adu_q_1_3, crohns_followup_adu_q_1_3,colitisulcerosa_followup_adu_q_1_1.b, colitisulcerosa_followup_adu_q_1_1.c, crohns_followup_adu_q_1_1.c, crohns_followup_adu_q_1_1.b)) %>% mutate(IBD = Status_UC+Status_CD)

IBD_developers %>% filter(Status_UC== 1) %>% dim()
IBD_developers %>% filter(Status_CD== 1) %>% dim()

#Status_UC, Status_CD

#Healthy Samples
Healthy = read_tsv("IBD_Healthy.tsv")



#2. Change names to the old IDs
Link = read_tsv("/groups/umcg-lifelines/prm03/releases/pheno_lifelines/v1/phenotype_linkage_file_project_pseudo_id.txt")
#PROJECT_PSEUDO_ID       PSEUDOIDEXT
Link %>% mutate(ID = PROJECT_PSEUDO_ID) %>% select(-PROJECT_PSEUDO_ID) -> Link
Link %>% filter(ID %in% c(IBD_developers$ID, Healthy$ID)) -> Link
left_join(Link, IBD_developers, by ="ID") %>% mutate(Participant = PSEUDOIDEXT) %>% select(-c(PSEUDOIDEXT, ID)) -> IBD_developers

#People who develop both are removed
IBD_developers[is.na(IBD_developers)] = 0
IBD_developers %>% filter(!IBD == 2) -> To_test



#3. Add Covariates
#Age_n Sex_n  BMI_n Smoking_n Smoking_ever       ID X7    Smoking
read_tsv("../TABLES/DIET_SCORES.tsv") -> Diet_score
left_join(To_test, Diet_score, by="Participant") -> To_test
To_test = To_test[!duplicated(To_test$Participant),]


#New covariates
read_tsv("../TABLES/Metadata2_LLD.tsv") -> Covariates
Covariates %>% mutate(Participant = as.numeric(ID)) %>% select(-ID) -> Covariates
left_join(To_test, Covariates, by="Participant") -> To_test
#Add PCs
read_tsv("../TABLES/PCs.tsv") -> PCs
colnames(PCs)[6] = "Participant"
left_join(To_test, PCs, by="Participant") -> To_test
#Add robust PCs
read_tsv("../TABLES/PCs_robust.tsv") -> PCs_r
colnames(PCs_r)[6] = "Participant"
left_join(To_test, PCs_r, by="Participant") -> To_test
#Add Food items
read_tsv("../TABLES/Table_participant_GDAG.tsv_PCA_input") -> Food_groups
col_remove = str_split(str="SumOfkcal       SumOfkJ SumOfeiwittot   SumOfeiwitplant SumOfeiwitdier  SumOfvettot     SumOfkhtot      SumOfmodis      SumOfpolys     SumOfalcohol     SumOffree_sug   SumOfadd_sug    SumOfglucose    SumOffructose   SumOflactose    SumOfmaltose    SumOfsucrose    SumOfGlyc_ind   SumOfGlyc_load  leeftijd  GESLACHT        GEWICHT LENGTE  BMI     BMR     EI_BMR",pattern=" ")
col_remove = unique(col_remove[[1]][col_remove[[1]] != ""])
Food_groups %>% select(-col_remove)   -> Food_groups
left_join(To_test, Food_groups) -> To_test

##########Filter BMR
read_tsv("../TABLES/Samples_BMR_filtered.tsv") -> To_keep


To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())

print("Filtering samples based on cal intake")
To_test %>% filter(Participant %in% To_keep$Participant) -> To_test


#Confint calculates profile confidence intervals
Logistic_model = function(To_test_m){
	To_test_m %>% drop_na() -> To_test
	#PCA model

	Model_PCs = glm(Status ~ RC1 + RC2 + RC3 + RC4 + RC5,To_test, family=binomial(link="logit"))
	print(summary(Model_PCs))
	confint(Model_PCs) -> CONF
	as.data.frame(CONF)[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model_PCs)$coefficients)[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> Results_PCA_noCov
		
	left_join(Results_PCA_noCov, CONFS) ->Results_PCA_noCov
	Results_PCA_noCov %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_PCA_noCov
	
	#Diet scores
	glm(Status ~ Protein ,To_test, family=binomial(link="logit")) -> Model1
	as.data.frame(confint(Model1))["Protein",] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model1)$coefficients)["Protein",] %>% rownames_to_column("Regressor") -> Results_potein_no_cov
	left_join(Results_potein_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) ->Results_potein_no_cov

	glm(Status ~ LL ,To_test, family=binomial(link="logit")) -> Model2
	as.data.frame(confint(Model2))["LL",] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model2)$coefficients)["LL",] %>% rownames_to_column("Regressor") -> Results_LL_no_cov
	left_join(Results_LL_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) ->Results_LL_no_cov
	
	
	glm(Status ~ HEI ,To_test, family=binomial(link="logit")) -> Model2
        as.data.frame(confint(Model2))["HEI",] %>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model2)$coefficients)["HEI",] %>% rownames_to_column("Regressor") -> Results_HEI_no_cov
        left_join(Results_HEI_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) ->Results_HEI_no_cov

        glm(Status ~ aMED ,To_test, family=binomial(link="logit")) -> Model2
        as.data.frame(confint(Model2))["aMED",] %>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model2)$coefficients)["aMED",] %>% rownames_to_column("Regressor") -> Results_aMED_no_cov
        left_join(Results_aMED_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) ->Results_aMED_no_cov

	
	glm(Status ~ mMED ,To_test, family=binomial(link="logit")) -> Model2
	as.data.frame(confint(Model2))["mMED",] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model2)$coefficients)["mMED",] %>% rownames_to_column("Regressor") -> Results_mMED_no_cov
	left_join(Results_mMED_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) ->Results_mMED_no_cov
	

	#Covariates
	glm(Status ~ Age_n + Sex_n + BMI_n + Smoking_n,To_test, family=binomial(link="logit")) -> Model3
	as.data.frame(confint(Model3))[c("Age_n","Sex_n","BMI_n","Smoking_n"),] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model3)$coefficients)[c("Age_n","Sex_n","BMI_n","Smoking_n"),] %>% rownames_to_column("Regressor") -> Results_Cov_no_cov
	left_join(Results_Cov_no_cov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_Cov_no_cov


	#Robust PCA
	glm(Status ~ V1 + V2 + V3 + V4 + V5,To_test, family=binomial(link="logit")) -> Model4
	as.data.frame(confint(Model4))[c("V1","V2","V3","V4","V5"),] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model4)$coefficients)[c("V1","V2","V3","V4","V5"),]  %>% rownames_to_column("Regressor") -> Results_PCA_robust_noCov
	left_join(Results_PCA_robust_noCov, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_PCA_robust_noCov


	rbind(rbind(rbind(rbind(rbind(rbind(rbind(Results_PCA_noCov,Results_potein_no_cov),Results_LL_no_cov),Results_HEI_no_cov),Results_aMED_no_cov),Results_mMED_no_cov),Results_Cov_no_cov),Results_PCA_robust_noCov) -> uncorrected_models

	print("including covariates")
	#Includeing covariates
	glm(Status ~ RC1 + RC2 + RC3 + RC4 + RC5 + Age_n + Sex_n + BMI_n + Smoking_n,To_test, family=binomial(link="logit")) -> Model5
	as.data.frame(confint(Model5))[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model5)$coefficients)[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> Results_PCA
	left_join(Results_PCA, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_PCA

	glm(Status ~ Protein + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model6
	as.data.frame(confint(Model6))["Protein",]%>% rownames_to_column("Regressor") -> CONFS
	as.data.frame(summary(Model6)$coefficients)["Protein",] %>% rownames_to_column("Regressor") -> Results_potein
	left_join(Results_potein, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_potein

	glm(Status ~ LL + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model6
        as.data.frame(confint(Model6))["LL",]%>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model6)$coefficients)["LL",] %>% rownames_to_column("Regressor") -> Results_LL
        left_join(Results_LL, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_LL

	
	glm(Status ~ HEI + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model6
        as.data.frame(confint(Model6))["HEI",]%>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model6)$coefficients)["HEI",] %>% rownames_to_column("Regressor") -> Results_HEI
        left_join(Results_HEI, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_HEI
	
	
        glm(Status ~ aMED + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model6
        as.data.frame(confint(Model6))["aMED",]%>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model6)$coefficients)["aMED",] %>% rownames_to_column("Regressor") -> Results_aMED
        left_join(Results_aMED, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_aMED
	
        glm(Status ~ mMED + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model6
        as.data.frame(confint(Model6))["mMED",]%>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model6)$coefficients)["mMED",] %>% rownames_to_column("Regressor") -> Results_mMED
        left_join(Results_mMED, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_mMED
	
        glm(Status ~ V1 + V2 + V3 + V4 + V5 + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")) -> Model7
        as.data.frame(confint(Model7))[c("V1","V2","V3","V4","V5"),] %>% rownames_to_column("Regressor") -> CONFS
        as.data.frame(summary(Model7)$coefficients)[c("V1","V2","V3","V4","V5"),]  %>% rownames_to_column("Regressor") -> Results_PCA_robust
        left_join(Results_PCA_robust, CONFS)  %>% mutate(Odds = exp(as.numeric(Estimate)), Odds_5= exp(as.numeric(`2.5 %`)), Odds_95= exp(as.numeric(`97.5 %`))) -> Results_PCA_robust	
		
	
	rbind(rbind(rbind(rbind(rbind(rbind(Results_PCA,Results_potein),Results_LL),Results_HEI),Results_aMED),Results_mMED),Results_PCA_robust) -> Corrected_models

	return(list(uncorrected_models,Corrected_models))	

}
#smoking_n2 includes both Current smokers and smokers ever
To_test %>% mutate(Smoking_n = Smoking) -> To_test
To_test %>% select(-c("Smoking_ever","X7")) -> To_test

print(To_test)
To_test %>% filter(! is.na(BMI_n)) -> test_2
test_2 %>% filter(! is.na(Smoking_n)) -> test_3
test_3 %>% filter(! is.na(LL)) -> test_4
test_4  %>% filter(!is.na(Age_n)) -> test_5
test_5 %>% filter(!is.na(Sex_n)) -> test_6
dim(To_test) ; dim(test_2) ; dim(test_3) ; dim(test_4) ; dim(test_5) ; dim(test_6)

Make_quantiles = function(x){
  quantile(x) -> Qs
  y = as.vector(x)
  y[x <= Qs[2]] = 1 
  y[x> Qs[2] & x <= Qs[3]] = 2 
  y[x> Qs[3] & x <= Qs[4]] = 3 
  y[x > Qs[4]] = 4 
  return(y)
}

To_test_q = To_test
for (i in c("RC1","RC2","RC3","RC4","RC5", "V1","V2","V3","V4","V5") ){
	i = which(colnames(To_test_q) == i)
	print(i)
        To_test_q[i] = Make_quantiles(as_vector(To_test_q[i]))
}




print("Testing CD and UC")
To_test %>% filter(Status_CD == 1) %>% mutate(Status=Status_CD) -> CD_to_Test
To_test %>% filter(Status_UC == 1) %>% mutate(Status=Status_UC)  -> UC_to_Test
To_test %>% filter(Status_UC == 0) %>% mutate(Status=0)  -> Healthy

To_test  %>% group_by(Status_CD) %>% summarise(n())

rbind(CD_to_Test, Healthy) -> CD_to_Test
rbind(UC_to_Test, Healthy) -> UC_to_Test

To_test_q %>% filter(Status_CD == 1) %>% mutate(Status=Status_CD) -> CD_to_Test_q ; To_test_q %>% filter(Status_UC == 1) %>% mutate(Status=Status_UC)  -> UC_to_Test_q ; To_test_q %>% filter(Status_UC == 0) %>% mutate(Status=0)  -> Healthy_q
rbind(CD_to_Test_q, Healthy_q) -> CD_to_Test_q ; rbind(UC_to_Test_q, Healthy_q) -> UC_to_Test_q


Late = F
###Late onset?
if (Late == T){
#Just for late onset
CD_to_Test %>% filter(Age_n>40) -> CD_to_Test
UC_to_Test %>% filter(Age_n>40) -> UC_to_Test

CD_to_Test %>% group_by(Status) %>% summarise(n())
UC_to_Test %>% group_by(Status) %>% summarise(n())
}


TEST = F
if (TEST == T){
CD_to_Test %>% group_by(Status) %>% summarise(n())
UC_to_Test %>% group_by(Status) %>% summarise(n())

write_tsv(CD_to_Test, "Regression_input_CD.tsv")
write_tsv(UC_to_Test, "Regression_input_UC.tsv")



Logistic_model(CD_to_Test) -> CD_results
print(CD_results)
Logistic_model(UC_to_Test) -> UC_results
print(UC_results)


write_tsv(x=CD_results[[1]], "RESULTS/CD_associations_uncorrected.tsv")
write_tsv(x=CD_results[[2]], "RESULTS/CD_associations_corrected.tsv")

write_tsv(x=UC_results[[1]],"RESULTS/UC_associations_uncorrected.tsv")
write_tsv(x=UC_results[[2]], "RESULTS/UC_associations_corrected.tsv")


###Run quantile analysis
Logistic_model(CD_to_Test_q) -> CD_results_q
print(CD_results_q)


Logistic_model(UC_to_Test_q) -> UC_results_q
print(UC_results_q)

write_tsv(x=CD_results_q[[1]], "RESULTS/CD_associations_QUANTILES_uncorrected.tsv")
write_tsv(x=CD_results_q[[2]], "RESULTS/CD_associations_QUANTILES_corrected.tsv")

write_tsv(x=UC_results_q[[1]],"RESULTS/UC_associations_QUANTILES_uncorrected.tsv")
write_tsv(x=UC_results_q[[2]], "RESULTS/UC_associations_QUANTILES_corrected.tsv")
}



options("width"=200)


Final_food_group = tibble()
Final_food_group_Corrected = tibble()
for (Food_Group in colnames(To_test)[23:length( colnames(To_test))]){
        Regressor  <- as.vector(as_vector(select(To_test , Food_Group)))
        To_test %>% mutate(Regressor = Regressor) -> Test_food_group
        Test_food_group %>% filter(! is.na(Status_CD)) %>% mutate(Status = Status_CD) %>% select(-c(Status_CD,Status_UC)) -> CD_to_Test
        Test_food_group %>% filter(! is.na(Status_UC)) %>% mutate(Status = Status_UC) %>% select(-c(Status_CD,Status_UC))  -> UC_to_Test        
        
        as.data.frame(summary(glm(Status ~ Regressor,CD_to_Test, family=binomial(link="logit")))$coefficients)["Regressor",] %>% mutate(Regressor = Food_Group) %>% mutate(Status = "CD_dev") -> R1
        as.data.frame(summary(glm(Status ~ Regressor,UC_to_Test, family=binomial(link="logit")))$coefficients)["Regressor",] %>% mutate(Regressor = Food_Group) %>% mutate(Status = "UC_dev") -> R2
        Final_food_group = rbind(Final_food_group, R1)
        Final_food_group = rbind(Final_food_group, R2)

        as.data.frame(summary(glm(Status ~ Regressor +Age_n + Sex_n + BMI_n + Smoking_n ,CD_to_Test, family=binomial(link="logit")))$coefficients)["Regressor",] %>% mutate(Regressor = Food_Group) %>% mutate(Status = "CD_dev") -> R1
         as.data.frame(summary(glm(Status ~ Regressor +Age_n + Sex_n + BMI_n + Smoking_n,UC_to_Test, family=binomial(link="logit")))$coefficients)["Regressor",] %>% mutate(Regressor = Food_Group) %>% mutate(Status = "UC_dev") -> R2
        Final_food_group_Corrected =  rbind(rbind(Final_food_group_Corrected, R1), R2)

}

Final_food_group %>% mutate(FDR = p.adjust(`Pr(>|z|)`, "fdr")) -> Food_group_results
write_tsv(x=Food_group_results, "RESULTS/ALL_Food_groups_associations.tsv")
Final_food_group_Corrected %>% mutate(FDR = p.adjust(`Pr(>|z|)`, "fdr")) -> Final_food_group_Corrected
write_tsv(x=Final_food_group_Corrected, "RESULTS/ALL_Food_groups_associations_corrected.tsv")







