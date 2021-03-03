# By: Weersma Group, UMCG (2020)
# =================================================

# =================================================
# Diet consistency analysis
# =================================================
# Script compares dietarry questionnaries between
# fecal sampling (T=0) and previously collected (T-5 years)
# diet questionnaries with higher number of questions 
# to establish diet consistency over 5-year time period
#
# Details:
# ============
# This script opens a tsv file where there are 4 columns (sample ID, Diet code, Baseline answer,  Follow-up answer)
# Note: Each row should be a different combination sample-diet_item
# It removes all rows with some NAs (missing Baseline or Followup answer) and count how many answers were identical between both timepoints
# and then calculate euclydian distance between time-points

library(tidyverse)

euc.dist <- function(x1, x2){ B = abs(as.numeric(x1) - as.numeric(x2)) ; B = B[!is.na(B)] ; return(mean(B)) }
euc.dist_wrong  <- function(x1, x2){ B = abs(as.numeric(x1) - as.numeric(x2)) ; B = B[!is.na(B)] ; B = B[!B==0] ; return(mean(B)) }

#Code for speeding up process
#Data = read_tsv("Matched_data.tsv")
#Data %>% filter(! is.na(Followup)) %>% filter(! is.na(Baseline)) -> Data_filtered
#print(dim(Data_filtered))
#write_tsv(Data_filtered,"Matched_data_filtered.tsv")

#Data_filtered = read_tsv("Matched_data_filtered.tsv")
Data_filtered = read_tsv("Matched_data.tsv")


Data_filtered[Data_filtered$Diet == "FOOD1",] -> Q1
Q1[Q1$Baseline != "No",]$ID -> IDs_BS
Q1[Q1$Followup != "No",]$ID -> IDs_FA



Require_4_scale = c("BB2Q1A", "BB1Q10A", "BB2Q14A")
Scale_4 = function(Vec){
	
	Vec = toupper(Vec)
	Vec[Vec=="ALWAYS"]  = 4 ; Vec[Vec=="NEVER"]  = 1 ; Vec[Vec=="OFTEN"]  = 3 ; Vec[Vec=="SOMETIMES"]  = 2 
	Vec[Vec=="NOT THIS MONTH"]  = 1 ; Vec[Vec=="1 DAY A MONTH"]  = 1 ; Vec[Vec=="2-3 DAYS A MONTH"]  = 2; Vec[Vec=="1 DAY A WEEK"]  = 2  ; Vec[Vec=="2-3 DAYS A WEEK"]  = 3 ; Vec[Vec=="4-5 DAYS A WEEK"] =3 ; Vec[Vec=="6-7 DAYS A WEEK"]  = 4
	
	return(Vec)
}

Phenos_dep = c("FOOD4E", "FOOD2A", "FOOD4A", "FOOD2G", "FOOD2D", "FOOD2C", "FOOD2B", "FOOD2F", "FOOD2E", "FOOD4B", "FOOD4C", "FOOD4D")
FD = function(x){
	if (x["ID"] %in% IDs_BS && is.na(x["Baseline"])){ x["Baseline"] = "No"}
	return(x)
}
FD2 = function(x){
	 if (x["ID"] %in% IDs_FA && is.na(x["Folowup"])){ x["Followup"] = "No"}
	return(x)
}


Make_diet_factors = function(Vector){
  
  if (!is.na(as.numeric(Vector[1]))){
    Vector = as.numeric(Vector)
    return(Vector)
  } 
 
  Vector[Vector=="ALWAYS"] =  4
  Vector[Vector=="OFTEN"] =  3
  Vector[Vector=="SOMETIMES"] =  2
  Vector[Vector=="NEVER"] =  1
  
  Vector[Vector=="NO"] = 1
  Vector[Vector=="YES"] =  2
  
  Vector[Vector=="YES, ALWAYS"] = 3
  Vector[Vector=="YES, SOMETIMES"] = 2
  
  
  Vector[Vector=="NOT THIS MONTH"]  = 1 ; Vector[Vector=="1 DAY A MONTH"]  = 1 ; Vector[Vector=="2-3 DAYS A MONTH"]  = 2; Vector[Vector=="1 DAY A WEEK"]  = 2  ; Vector[Vector=="2-3 DAYS A WEEK"]  = 3 ; Vector[Vector=="4-5 DAYS A WEEK"] =3 ; Vector[Vector=="6-7 DAYS A WEEK"]  = 4
  
  
  return(Vector)
  
}

To_do =  unique(Data_filtered$Diet)
To_do = c("BB2Q14A",   "BB2Q8B3" ,  "BB1Q11B1",  "BB2Q8B4" ,  "BB1Q11B7",  "BB1Q11B10", "FOOD2G"   , "BB1Q11B13",
"BB2Q4E1",   "BB1Q11B5",  "BB2Q4E2" ,  "BB2Q4E4" ,  "BB2Q8B7" ,  "FOOD32A1" , "FOOD32A2" , "BB2Q8B5"  ,
"FOOD36B",   "FOOD2A",    "BB1Q11B12", "BB1Q11B9",  "FOOD28E" ,  "BB1Q11B3" , "BB1Q11B11", "FOOD36A"  ,
"BB1Q11",  "FOOD32A3",  "BB1Q11B4",    "BB2Q4E5" ,  "BB2Q1A"  ,  "FOOD35"   , "FOOD33B"  , "FOOD2C"   ,
"FOOD4E",  "BB1Q10A" ,  "FOOD28D" ,    "FOOD5"   ,  "FOOD1"   ,  "FOOD5B"   , "FOOD4A"   , "FOOD33A"  ,
"BB2Q4E",  "FOOD2E"  ,  "FOOD2F"  ,    "FOOD3"   ,  "FOOD2D"  ,  "FOOD2B"   , "FOOD4B"   , "FOOD4C"   ,
"FOOD4D",  "BB1Q11A" ,  "BB1Q15A" ,    "BB1Q16A" ,  "BB1Q1A"  ,  "BB1Q4A"   , "BB2Q10A"  , "BB2Q12A"  ,
"BB2Q13",  "BB2Q15A" ,  "BB2Q18A" ,    "BB2Q3A"  ,  "BB2Q4A"  ,  "BB2Q5A"   , "BB2Q6A"   , "BB2Q7A"   ,
"BB2Q8A",  "BB3Q15A" ,  "BB3Q16A" ,    "BB3Q17A" ,  "BB3Q1A"  ,  "BB3Q4A"   , "BB3Q5A"   , "CHFOOD27" ,
"FOOD5A")

length(To_do)

Result = tibble()
for (i in To_do){
	#Go variable by variable
	print(i)
	Data_filtered %>% filter(Diet==i) -> Fil
	#Special phenotypes were a NA is a 'No' if the first question is NO 
	if (i %in% Phenos_dep){
		#apply(Fil, 1, FUN= FD) ->New_f
		#t(New_f) %>%  as_tibble() %>%  `colnames<-`(rownames(New_f)) -> Fil
		#apply(Fil, 1, FUN= FD2) -> New_f
		#t(New_f) %>%  as_tibble() %>%  `colnames<-`(rownames(New_f)) -> Fil
		Fil %>% mutate(Baseline = ifelse(ID %in% IDs_BS & is.na(Baseline), "No", Baseline )) -> Fil
		Fil %>% mutate(Followup = ifelse(ID %in% IDs_BS & is.na(Followup), "No", Followup )) -> Fil
	}
	#Remove all NA records
	Fil %>% drop_na() -> Fil
	if (dim(Fil)[1] == 0){ next}

	#MAtch grammar
	if ("1 day per month" %in% unique(Fil$Followup) | "1 day a month" %in% unique(Fil$Baseline) |  "1 day per week" %in% unique(Fil$Baseline)) {Fil$Baseline = str_replace(Fil$Followup,"per month", "a month") ; Fil$Followup = str_replace(Fil$Followup,"per month", "a month") }
	if (i %in% Require_4_scale){ Fil$Baseline = Scale_4(Fil$Baseline) ;  Fil$Followup = Scale_4(Fil$Followup)}
	
	
	Fil$Baseline = toupper(Fil$Baseline)
	Fil$Followup = toupper(Fil$Followup)
	Number_samples = dim(Fil)[1]
	
	#Make everything in 4-levels
	Fil %>% mutate(Followup_n= Make_diet_factors(Fil$Followup), Baseline_n= Make_diet_factors(Fil$Baseline)) -> Fil
	#print(unique(Fil$Followup_n))
	#print(unique(Fil$Baseline_n))
	Dist = euc.dist(Fil$Followup_n, Fil$Baseline_n)
	Dist_wrong = euc.dist_wrong(Fil$Followup_n, Fil$Baseline_n)	
	Constant = (sum(Fil$Baseline == Fil$Followup)/Number_samples) * 100
	#Distance between factors: Would need to order them somehow
	temporal = tibble("Diet" = i, "Number_samples"=Number_samples, "Constant(%)"=Constant,"EuDist"=Dist,"Distance_change"= Dist_wrong, "Answers_baseline"= paste(sort(unique(Fil$Baseline)), collapse=","), "Answers_followup"= paste(sort(unique(Fil$Followup)), collapse=","))
	
	Result = rbind(Result, temporal)
}
#write_tsv(Result,"Summary_diery.tsv")

Result %>% arrange(`Constant(%)`) -> Result

write_tsv(Result ,"Summary.tsv")


