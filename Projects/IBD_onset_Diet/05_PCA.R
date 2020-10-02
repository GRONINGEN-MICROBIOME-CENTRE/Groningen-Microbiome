#Script for making a PCA

library(tidyverse)
library(psych)

set.seed(99)


#ID      Age_n   Sex_n   BMI_n   Smoking_n
Meta = read_tsv("TABLES/Metadata_LLD.tsv")
Meta %>% mutate(Participant = ID) -> Meta
DF = read_tsv("TABLES/Table_participant_GDAG.tsv_PCA_input")

#DF = read_tsv("TABLES/Table_participant_KCAL.tsv_PCA_input")
left_join(DF,select(Meta, -Sex_n), by="Participant") -> Ex_data
Ex_data %>% mutate(Sex_n = GESLACHT) -> Ex_data

Data_filter = function(Exploration_data = Ex_data, Filter = 3){
	#Filters
	Exploration_data %>% filter(is.na(Sex_n)) %>% summarise(n())
	Exploration_data %>% drop_na(Sex_n) %>% mutate(SumOfkcal =as.numeric(SumOfkcal)) -> Exploration_data
	quantile(filter(Exploration_data, Sex_n=="Male")$SumOfkcal,probs=c(0.025, 0.975)) -> Male_bounds
	quantile(filter(Exploration_data, Sex_n=="Female")$SumOfkcal,probs=c(0.025, 0.975)) -> Female_bounds
	#Filter 1: Fixed lower bound, high upper bound
	if (Filter == 1){
	Exploration_data %>% filter(! (Sex_n == "Female" &  SumOfkcal < 500) ) %>% filter(! (Sex_n == "Male" &  SumOfkcal < 800) ) -> Lower_bound_filtered
	Lower_bound_filtered %>% filter(! (Sex_n == "Female" &  SumOfkcal > Female_bounds[2]) ) %>%  filter(! (Sex_n == "Male" &  SumOfkcal > Male_bounds[2]) ) -> Filtered
	} else if (Filter == 2){
	#Filter 2: Filter on up and low quantiles
	Exploration_data %>% filter(! (Sex_n == "Female" &  SumOfkcal < Female_bounds[1])) %>% filter(! (Sex_n == "Female" &  SumOfkcal > Female_bounds[2]) ) %>% filter(! (Sex_n == "Male" &  SumOfkcal < Male_bounds[1])) %>% filter(! (Sex_n == "Male" &  SumOfkcal > Male_bounds[2]) ) -> Filtered
	} else if (Filter == 3){
	Exploration_data %>% drop_na(BMR) %>% mutate(BMR =as.numeric(BMR)) -> Exploration_data
	quantile(filter(Exploration_data, Sex_n=="Male")$BMR,probs=c(0.025, 0.975)) -> Male_bounds
	quantile(filter(Exploration_data, Sex_n=="Female")$BMR,probs=c(0.025, 0.975)) -> Female_bounds
	 Exploration_data %>% filter(! (Sex_n == "Female" & BMR < Female_bounds[1])) %>% filter(! (Sex_n == "Female" &  BMR  > Female_bounds[2]) ) %>% filter(! (Sex_n == "Male" & BMR  < Male_bounds[1])) %>% filter(! (Sex_n == "Male" &  BMR  > Male_bounds[2]) ) -> Filtered	
	}

	return(Filtered)

}


Summaries = F
if (Summaries == T){
	Exploration_data %>% ggplot(aes(x=Age_n, col=Sex_n)) + geom_density() + theme_bw() -> PLOT1
	print("Age sum")
	summary(Exploration_data$Age_n)
	print("Sex sum")
	summary(as.factor(Exploration_data$Sex_n))
	print("KCAL sum")
	summary(Exploration_data$SumOfkcal)
	Exploration_data %>% ggplot(aes(x=SumOfkcal, col=Sex_n)) + geom_density() + theme_bw() -> PLOT

	ggsave(plot=PLOT, filename= "Distribution_SumOfkcal.pdf")
	ggsave(plot=PLOT1, filename= "Distribution_age.pdf")

	Exploration_data %>% group_by(Sex_n) %>% summarise(mean(Age_n), mean(SumOfkcal), sd(Age_n), sd(SumOfkcal), max(Age_n), max(SumOfkcal), min(Age_n), min(SumOfkcal)) -> Ssummary
	write_tsv(Ssummary, "Metadata_summary.tsv")

	q()

}


col_remove = c("SUMOFALCOHOL","SUMOFEIWITDIER","SUMOFEIWITPLANT","SUMOFEIWITTOT","SumOfkcal","SUMOFKHTOT","SUMOFKJ","SUMOFVETTOT")
col_remove = str_split(str="SumOfkcal       SumOfkJ SumOfeiwittot   SumOfeiwitplant SumOfeiwitdier  SumOfvettot     SumOfkhtot      SumOfmodis      SumOfpolys     SumOfalcohol     SumOffree_sug   SumOfadd_sug    SumOfglucose    SumOffructose   SumOflactose    SumOfmaltose    SumOfsucrose    SumOfGlyc_ind   SumOfGlyc_load  leeftijd  GESLACHT        GEWICHT LENGTE  BMI     BMR     EI_BMR",pattern=" ")
col_remove = unique(col_remove[[1]][col_remove[[1]] != ""])

DF %>% select(-col_remove) -> DF


Compute_PCA = function(DF, ncomp, Do_others=F){

	apply(DF, 2, FUN= scale) %>% as_tibble() -> DF
	if (Do_others  == T){
		print(" provide KMO and Barlett’s test of sphericity significance from the models")
		print("correlation matrix")
		cor_matrix = cor(DF)
		print(round(cor_matrix, 3))

		t_bart_1 <- psych::cortest.bartlett(cor_matrix, n=nrow(DF))
		print("Barlett results")
		print(t_bart_1)

		print(" KMO (Kaiser-Meyer-Olkin)")
		t_kmo <- psych::KMO(cor_matrix)
		print(t_kmo)
		#summary(t_kmo$MSAi)
		print("MSAi below 0.5")
		print(t_kmo$MSAi[t_kmo$MSAi<0.5])
	}
	# package psych does PCA and varimax in one go.
	res_pca_rotated <- psych::principal(DF, rotate="varimax", nfactors=ncomp, scores=TRUE)
	print(res_pca_rotated)
	return(res_pca_rotated)
}



Filtered = Data_filter(Filter=1)
write_tsv(Filtered , "TABLES/Samples_BMR_filtered.tsv")
DF %>% filter(Participant %in% Filtered$Participant) %>% select(-Participant) -> DF_filtered1
print("PCA on Filter #1")
#Compute_PCA(DF_filtered1,2, T)
Compute_PCA(DF_filtered1,5,Do_others=T) -> PCs
PCs <-loadings(PCs)



apply(DF_filtered1, 2, FUN= scale) %*% PCs -> RMS

as_tibble(RMS) %>% mutate(ID = filter(DF, Participant %in% Filtered$Participant)$Participant) -> RMS

write.table(RMS,file='TABLES/PCs.tsv', quote=FALSE, sep='\t',row.names = F)



#Robust PCA
apply(DF_filtered1, 2, FUN= scale) %>% as_tibble() -> DF_f
library(rrcov)
.libPaths("/groups/umcg-lifelines/tmp01/users/umcg-sandreusanchez/Diet_IBD/R_library")
library(pracma)

res_pca_2 <- PcaHubert(DF_f, k=5, alpha =0.75, scale=FALSE, center = FALSE)
#Rotate loadings
sdev <- sqrt(res_pca_2$eigenvalues) #compute SD of PCs from eigenvalues
raw_loadings_2 <- res_pca_2$loadings[,1:5] %*% diag(sdev, 5, 5)
rotated_loadings <- varimax(raw_loadings_2, normalize = FALSE)$loadings
inv_loadings <- t(pracma::pinv(rotated_loadings))
res_pca_rot_scores <- scale(DF_f) %*% inv_loadings
as_tibble(res_pca_rot_scores) %>% mutate(ID = filter(DF, Participant %in% Filtered$Participant)$Participant) -> res_pca_rot_scores
write.table(res_pca_rot_scores,file='TABLES/PCs_robust.tsv', quote=FALSE, sep='\t',row.names = F)

print(inv_loadings)

# Variance (Proportion)
r3 <- diag(cor(DF_f, use = "pairwise"))
communalities3 <- rowSums(rotated_loadings^2)
uniquenesses3 <- r3 - communalities3
vx3 <- colSums(rotated_loadings^2)
vtotal3 <- sum(communalities3 + uniquenesses3)
# Proportion Var
vx3/vtotal3
#Cumulative Var
cumsum(vx3/vtotal3)
#Proportion Explained
vx3/sum(vx3)
#Cumulative Proportion
cumsum(vx3/sum(vx3))





#####Do PCA separately for CD/UC and healthy
Info_Disease = read_tsv("TABLES/Dataset_IBD_test.tsv") #After script 7
Info_Disease %>% filter(Status_CD == 1) -> CD_samples
Info_Disease %>% filter(Status_UC == 1) -> UC_samples
Info_Disease %>% filter(Status_UC == 0 & Status_CD == 0) -> Healthy_samples

DF %>% filter(Participant %in% CD_samples$ID) %>% filter(Participant %in% Filtered$Participant) %>% select(-Participant) -> DF_filtered_CD
DF %>% filter(Participant %in% UC_samples$ID) %>% filter(Participant %in% Filtered$Participant) %>% select(-Participant) -> DF_filtered_UC
DF %>% filter(Participant %in% Healthy_samples$ID) %>% filter(Participant %in% Filtered$Participant) %>% select(-Participant) -> DF_filtered_health


print("PCA CD dev")
Compute_PCA(DF_filtered_CD,5) -> PCA_CD
print("PCA UC dev")
Compute_PCA(DF_filtered_UC, 5) -> PCA_UC
print("PCA healthy")
Compute_PCA(DF_filtered_health,5) -> PCA_Healthy
#######




