library(tidyverse)
wideScreen <- function(howWide=160) {
  options(width=as.integer(howWide))
}
wideScreen()

Family = F
if (Family == T){
	D = read_tsv("../TABLES/Family_disease.tsv")
	D %>% filter(! Disease %in% c("crohns_siblings_fam_q_1", "ulcerativecolitis_siblings_fam_q_1", "crohns_children_fam_q_1", "ulcerativecolitis_children_fam_q_1") ) -> D #Remopve sibligns because they are weird
	spread(D, Disease, `Baseline assessment, questionnaire 1`) -> D2
	D2$ID -> IDs
	print(D2)
	D2 %>% select(-ID) %>% as.matrix()  -> D3
	storage.mode(D3) <- "numeric"
	D3[D3 == 3] = NA
	D3[D3 == 1] = 0
	apply(D3,1 , FUN= function(x){ sum(x, na.rm = T) } ) -> D3

	D2 %>% mutate( presence = D3) %>% mutate(Family_case = ifelse(presence > 0, 1, 0)) -> D2
	D2 -> Family

	write_tsv(D2, "TABLES/Family_disease_final.tsv")
}

##Change metadata file

D2 = read_tsv("../TABLES/Metadata_LLD.tsv",na = c("", "NA", "NaN")) 
#Smoking_n       Smoking_ever
print(D2)
D2 %>% group_by(Smoking_n) %>% summarise(n())
D2 %>% mutate(Smoking = ifelse(is.na(Smoking_n), NA , ifelse(Smoking_n == "Yes", 2, ifelse(is.na(Smoking_ever), 0, ifelse(Smoking_ever  == "Yes", 1, 0))))) -> D2
D2 %>% group_by(Smoking) %>% summarise(n())

print(D2)

write_tsv(D2, "../TABLES/Metadata2_LLD.tsv")



###############################################################
#####################From here, old version ###################
###############################################################
q()

x = as.numeric(D2$Smoking_n == "Yes")
y = as.numeric(D2$Smoking_ever == "Yes")



print(D2)
q()
Smoking = rowSums( cbind (x,y), na.rm=TRUE) >= 1

Smoking[is.na(x)] = NA
Smoking[Smoking] = "Yes"
Smoking[Smoking == F] = "No"

D2 %>% mutate(Smoking_n2 = Smoking) -> D2
#D2 %>% filter(Smoking_n != Smoking_n2) %>% print()

write_tsv(D2, "TABLES/Metadata2_LLD.tsv")


##Change test file
D2 = read_tsv("TABLES/Data_v2.tsv")
x = as.numeric(D2$Smoking_n == "Yes") ; y = as.numeric(D2$Smoking_ever == "Yes")
Smoking = rowSums( cbind (x,y), na.rm=TRUE) >= 1
Smoking[is.na(x)] = NA ; Smoking[Smoking] = "Yes" ; Smoking[Smoking == F] = "No"
D2 %>% mutate(Smoking_n2 = Smoking) -> D2 

filter(Family, ID %in% D2$ID) %>% arrange(ID) -> Family
D2 %>% arrange(ID) -> D2

left_join(D2, select(Family, c(ID, Family_case )), by="ID") -> D2

write_tsv(D2, "TABLES/Data_v3.tsv")
