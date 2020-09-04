# aMED Score LL CLuster, Dietary-scores_Mediterannean_Diet_Score
# aMED-Fung_et_al
library(tidyverse)
phenosdiet1 = read_tsv("TABLES/Table_participant_GDAG.tsv_script_eva_input")

#ffq data, phenosdiet2 #rows is samples, column is foodgroups   
MDS<-phenosdiet1

#Positive consumption above median gets 1 points, below median 0 points. 

#1 Vegetables,  "group_vegetables" includes; "vegetables_cooked_nobutter, vegetables_stirfried, vegetables_cooked_butter"
veg.total <- MDS$vegetables_cooked_nobutter + MDS$vegetables_stirfried + MDS$vegetables_cooked_butter




#veg.total<-MDS$group_vegetables
x<-median(veg.total)                    
MDS$veg.points <- ifelse(veg.total>= x, 1, 0)

#2 Legumes, only 1 column legumes in database; "legumes" 
leg.total<- MDS$legumes
x<- median(leg.total)
MDS$leg.points <- ifelse(leg.total>= x, 1, 0)

#3 fruits,  "fruitjuice","fruit" (group_fruits contains applesauce)
fruit.total<-MDS$fruit+MDS$fruitjuice  
x<-median(fruit.total)  
MDS$fruit.points<- ifelse(fruit.total>= x, 1, 0)

#4 nuts, "group_nuts" consists of peanutbutter, nut_d, snack_nut
#nut.total<-MDS$group_nuts
nut.total <- MDS$nut_d + MDS$snack_nut + MDS$peanutbutter
x<-median(nut.total)
MDS$nut.points<-ifelse(nut.total>= x, 1, 0)

#5 whole grains, "bread"
whgrains.total<-MDS$bread
x<-median(whgrains.total)
MDS$whgrains.points<-ifelse(whgrains.total>= x, 1, 0)

#6 fish, group_fish contains herring_salted, fish_white_fried, fish_lean, fish_fatty, fish_other, fish_prepared_fat
#fish.total<-MDS$group_fish
fish.total <- MDS$herring_salted + MDS$fish_white_fried + MDS$fish_lean + MDS$fish_fatty + MDS$fish_other + MDS$fish_prepared_fat
x<-median(fish.total)
MDS$fish.points<-ifelse(fish.total>= x, 1, 0)

#7 red and processed meats; containing: sausage_smoked, meats_other_cold_cuts, meat, beef_lean, beef_fat
#pork_lean, pork_fat, pork_processed, meats_fat, snack_meats
meat.total<- MDS$sausage_smoked_cold_cuts + MDS$meat + MDS$meats_fat+ MDS$meats_other_cold_cuts+ MDS$beef_fat+ MDS$beef_lean + 
  MDS$pork_fat + MDS$pork_lean + MDS$pork_processed + MDS$snack_meats
x<-median(meat.total)
MDS$meat.points<-ifelse(meat.total>= x, 0, 1)

#8 alcohol; group_alcohol consist of beer, wine_red, wine_white, wine_fort, spirits, other_alc_drinks
#  1 point if 10-20 g/d (1-2 glasses a day),  1 glas alcohol is 10 g / day 
#summary(MDS$group_alcohol)
MDS$group_alcohol <- MDS$beer + MDS$wine_red + MDS$wine_white + MDS$wine_fort + MDS$spirits + MDS$other_alc_drinks
MDS$alcohol.points<-ifelse((MDS$group_alcohol >=10 & MDS$group_alcohol <= 20), 1, 0 ) 

#9 M/S ratio no data 

#add the points together [max amount of 8 points, min amount of 0 points]
MDS$aMED <- MDS$veg.points + MDS$leg.points + MDS$fruit.points + MDS$nut.points + 
  MDS$whgrains.points + MDS$fish.points + MDS$meat.points + MDS$alcohol.points

MDS %>% select(Participant, aMED) -> To_write
write_tsv(To_write,"TABLES/aMED_score.tsv")

#ggplot(MDS,aes(x=aMED)) + geom_histogram() + theme_bw() -> PLOT
#ggsave(filename="Histogram_amed.pdf", plot=PLOT)

