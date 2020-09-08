library(tidyverse)

#####Read dataset
To_test = read_tsv("TABLES/Dataset_IBD_test.tsv")
##########Filter BMR
read_tsv("TABLES/Samples_BMR_filtered.tsv") -> To_keep
#########

To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())

print("Filtering samples based on cal intake")
To_test %>% filter(! ID %in% To_keep$Participant) %>% group_by(Status_UC, Status_CD) %>% summarise(n())

To_test %>% filter(ID %in% To_keep$Participant) -> To_test



Logistic_model = function(To_test_m){
	To_test_m %>% drop_na() -> To_test
	#PCA model
	as.data.frame(summary(glm(Status ~ RC1 + RC2 + RC3 + RC4 + RC5,To_test, family=binomial(link="logit")))$coefficients)[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> Results_PCA_noCov
	#Diet scores
	as.data.frame(summary(glm(Status ~ Protein ,To_test, family=binomial(link="logit")))$coefficients)["Protein",] %>% rownames_to_column("Regressor") -> Results_potein_no_cov
	as.data.frame(summary(glm(Status ~ LL ,To_test, family=binomial(link="logit")))$coefficients)["LL",] %>% rownames_to_column("Regressor") -> Results_LL_no_cov
	as.data.frame(summary(glm(Status ~ HEI ,To_test, family=binomial(link="logit")))$coefficients)["HEI",] %>% rownames_to_column("Regressor") -> Results_HEI_no_cov
	as.data.frame(summary(glm(Status ~ aMED ,To_test, family=binomial(link="logit")))$coefficients)["aMED",] %>% rownames_to_column("Regressor") -> Results_aMED_no_cov
	as.data.frame(summary(glm(Status ~ mMED ,To_test, family=binomial(link="logit")))$coefficients)["mMED",] %>% rownames_to_column("Regressor") -> Results_mMED_no_cov
	#Covariates
	as.data.frame(summary(glm(Status ~ Age_n + Sex_n + BMI_n + Smoking_n,To_test, family=binomial(link="logit")))$coefficients)[c("Age_n","Sex_n","BMI_n","Smoking_n"),] %>% rownames_to_column("Regressor") -> Results_Cov_no_cov


	#Robust PCA
	as.data.frame(summary(glm(Status ~ V1 + V2 + V3 + V4 + V5,To_test, family=binomial(link="logit")))$coefficients)[c("V1","V2","V3","V4","V5"),] %>% rownames_to_column("Regressor") -> Results_PCA_robust_noCov


	rbind(rbind(rbind(rbind(rbind(rbind(rbind(Results_PCA_noCov,Results_potein_no_cov),Results_LL_no_cov),Results_HEI_no_cov),Results_aMED_no_cov),Results_mMED_no_cov),Results_Cov_no_cov),Results_PCA_robust_noCov) -> uncorrected_models


	as.data.frame(summary(glm(Status ~ RC1 + RC2 + RC3 + RC4 + RC5 + Age_n + Sex_n + BMI_n + Smoking_n,To_test, family=binomial(link="logit")))$coefficients)[c("RC1","RC2","RC3","RC4","RC5"),] %>% rownames_to_column("Regressor") -> Results_PCA
	as.data.frame(summary(glm(Status ~ Protein + Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")))$coefficients)["Protein",] %>% rownames_to_column("Regressor") -> Results_potein
        as.data.frame(summary(glm(Status ~ LL+Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")))$coefficients)["LL",] %>% rownames_to_column("Regressor") -> Results_LL
        as.data.frame(summary(glm(Status ~ HEI+Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")))$coefficients)["HEI",] %>% rownames_to_column("Regressor") -> Results_HEI
        as.data.frame(summary(glm(Status ~ aMED+Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")))$coefficients)["aMED",] %>% rownames_to_column("Regressor") -> Results_aMED
	as.data.frame(summary(glm(Status ~ mMED+Age_n + Sex_n + BMI_n + Smoking_n ,To_test, family=binomial(link="logit")))$coefficients)["mMED",] %>% rownames_to_column("Regressor") -> Results_mMED
	as.data.frame(summary(glm(Status ~ V1 + V2 + V3 + V4 + V5,To_test, family=binomial(link="logit")))$coefficients)[c("V1","V2","V3","V4","V5"),] %>% rownames_to_column("Regressor") -> Results_PCA_robust	
	
	rbind(rbind(rbind(rbind(rbind(rbind(Results_PCA,Results_potein),Results_LL),Results_HEI),Results_aMED),Results_mMED),Results_PCA_robust) -> Corrected_models

	return(list(uncorrected_models,Corrected_models))	

}

To_test %>% mutate(Smoking_n = ifelse(Smoking_n == "No", 0, 1)) -> To_test


print("After filtering")
To_test %>% group_by(Status_UC, Status_CD) %>% summarise(n())

print("Testing CD and UC")
To_test %>% filter(! is.na(Status_CD)) %>% mutate(Status = Status_CD) %>% select(-c(Status_CD,Status_UC)) -> CD_to_Test
To_test %>% filter(! is.na(Status_UC)) %>% mutate(Status = Status_UC) %>% select(-c(Status_CD,Status_UC))  -> UC_to_Test


Logistic_model(CD_to_Test) -> CD_results
Logistic_model(UC_to_Test) -> UC_results

write_tsv(x=CD_results[[1]], "RESULTS/CD_associations_uncorrected.tsv")
write_tsv(x=CD_results[[2]], "RESULTS/CD_associations_corrected.tsv")

write_tsv(x=UC_results[[1]],"RESULTS/UC_associations_uncorrected.tsv")
write_tsv(x=UC_results[[2]], "RESULTS/UC_associations_corrected.tsv")


#print(CD_results)
#print(UC_results)
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
write_tsv(x=Food_group_results, "RESULTS/Food_groups_associations.tsv")
Final_food_group_Corrected %>% mutate(FDR = p.adjust(`Pr(>|z|)`, "fdr")) -> Final_food_group_Corrected
write_tsv(x=Final_food_group_Corrected, "RESULTS/Food_groups_associations_corrected.tsv")



