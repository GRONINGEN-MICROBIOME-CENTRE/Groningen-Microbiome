#Creator: Eva Schuttert, mMED Score Khalili et al. 

#load data, rows is samples, column is foodgroups   
library(tidyverse)
MDS<- read_tsv("TABLES/Table_participant_GDAG.tsv_script_eva_input")

#Positive consumption above median gets 1 points, below median 0 points.Negative consumption the other way around.  
#1 Vegetables and  fruit; includes: "vegetables_cooked_nobutter, vegetables_stirfried, vegetables_cooked_butter
#                           fruits: "fruitjuice","fruit"
vegfruit.total <- MDS$vegetables_cooked_nobutter + MDS$vegetables_stirfried + MDS$vegetables_cooked_butter + MDS$fruit + MDS$fruitjuice 
x<-median(vegfruit.total)                    
MDS$vegfruit.points <- ifelse(vegfruit.total>= x, 1, 0)

#2 Legumes and nuts: legumes, peanutbutter, nut_d, snack_nut
legnut.total<- MDS$legumes + MDS$nut_d + MDS$snack_nut + MDS$peanutbutter
x<- median(legnut.total)
MDS$legnut.points <- ifelse(legnut.total>= x, 1, 0)

#3 non-refined or high fibre grains, "bread"
whgrains.total<-MDS$bread
x<-median(whgrains.total)
MDS$whgrains.points<-ifelse(whgrains.total>= x, 1, 0)

#4 fermented dairy products; 
ferdairy.total <- MDS$cheese_20 + MDS$cheese_40 + MDS$cheese_48 + MDS$cheese_other_cream_cheese + MDS$cheese_d + MDS$snack_cheese + 
                  MDS$buttermilk + MDS$yoghurt_drink_added_sugar + MDS$yoghurt_drink_other + MDS$custard_ff + MDS$yoghurt_ff + 
                  MDS$yoghurt_lf + MDS$yoghurt_lf_fruits + MDS$fromage_frais_fruits + MDS$other_custard_yoghurt_fromagefrais
x <- median(ferdairy.total)
MDS$ferdairy.points <- ifelse(ferdairy.total>= x, 1, 0)

#5 fish, group_fish contains herring_salted, fish_white_fried, fish_lean, fish_fatty, fish_other, fish_prepared_fat
#fish.total<-MDS$group_fish
fish.total <- MDS$herring_salted + MDS$fish_white_fried + MDS$fish_lean + MDS$fish_fatty + MDS$fish_other + MDS$fish_prepared_fat
x<-median(fish.total)
MDS$fish.points<-ifelse(fish.total>= x, 1, 0)

#6 red and processed meats; containing: sausage_smoked, meats_other_cold_cuts, meat, beef_lean, beef_fat
#pork_lean, pork_fat, pork_processed, meats_fat, snack_meats
meat.total<- MDS$sausage_smoked_cold_cuts + MDS$meat + MDS$meats_fat+ MDS$meats_other_cold_cuts+ MDS$beef_fat+ MDS$beef_lean + 
  MDS$pork_fat + MDS$pork_lean + MDS$pork_processed + MDS$snack_meats
x<-median(meat.total)
MDS$meat.points<-ifelse(meat.total>= x, 0, 1)

#7 olive oil / rapeseed oil -> no data yet, maybe use fatty acids? 
oils.total<-MDS$margarine_lfb + MDS$saladdressing_f + MDS$saladdressing_w + MDS$mayonaise + MDS$mayonaise_snack
x<-median(oils.total)
MDS$oil.points<-ifelse(oils.total>= x, 1, 0)

#8 alcohol; group_alcohol consist of beer, wine_red, wine_white, wine_fort, spirits, other_alc_drinks
#  1 point if 5-15 g/d (1-2 glasses a day),  1 glas alcohol is 10 g / day 
#summary(MDS$group_alcohol)
MDS$group_alcohol <- MDS$beer + MDS$wine_red + MDS$wine_white + MDS$wine_fort + MDS$spirits + MDS$other_alc_drinks
MDS$alcohol.points<-ifelse((MDS$group_alcohol >=5 & MDS$group_alcohol <= 15), 1, 0 ) 

#add the points together [max amount of 7 points, min amount of 0 points]
MDS$mMED <- MDS$vegfruit.points + MDS$legnut.points  +  MDS$whgrains.points + MDS$ferdairy.points + 
 + MDS$fish.points + MDS$meat.points + MDS$oil.points + MDS$alcohol.points

MDS %>% select(Participant, mMED) -> To_write
write_tsv(To_write,"TABLES/mMED_score.tsv")


