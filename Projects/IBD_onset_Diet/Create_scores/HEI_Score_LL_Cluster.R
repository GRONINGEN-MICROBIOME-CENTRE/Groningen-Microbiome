# Healthy Eating Index, HEI 2010, LL Cluster 
# adapted from Guenther et al. 2013

library(tidyverse)
phenosdiet1 = read_tsv("TABLES/Table_participant_GDAG.tsv_script_eva_input")
phenosdiet2 = read_tsv("TABLES/Table_participant_KCAL.tsv_script_eva_input")
#data FFQ gram/day 
HEI<-phenosdiet1
HEI2<-phenosdiet2

# more points if higher intake: Total Fruit, whole fruit, total vegetables, greens and beans,
#                               whole grains, dairy, Total Protein Foods, Seafood and Plant Proteins   
# more points if lower intake: Fatty acids, refined grains, sodium, empty calories 


#Components need to be divided by calories dived by 1000 to get per 1000 calorie amount
Kcal<- HEI$SumOfkcal/1000

# ++++ ADEQUACY COMPONENTS ++++ #

###Whole fruit###Only whole fruit and does not make use of eg fruitpies etc, in our data "fruit" 
HEI$wholefruit<-HEI$fruit/Kcal
HEI$wholefruit.points<- ntile(HEI$wholefruit, 6)-1  #divide into quintiles 

###Total Fruit ### "fruit" "fruitjuice"
HEI$totalfruit<- (HEI$fruit+HEI$fruitjuice)/Kcal
HEI$totalfruit.points<- ntile(HEI$totalfruit, 6)-1 

###Total Vegetables ### "group_vegetables" consisting of: "vegetables_cooked_nobutter" "vegetables_stirfried" "vegetables_cooked_butter"
#HEI$totalvegetables<- HEI$group_vegetables/Kcal
HEI$totalvegetables<- (HEI$vegetables_cooked_nobutter + HEI$vegetables_cooked_butter + HEI$vegetables_stirfried)/Kcal
HEI$totalvegetables.points<-ntile(HEI$totalvegetables,6)-1

###Greans and Beans ### "legumes"
HEI$greensandbeans<- HEI$legumes/Kcal
HEI$greensandbeans.points<- ntile(HEI$greensandbeans,6)-1

###wholegrains ### "bread"
HEI$wholegrains<-HEI$bread/Kcal
HEI$wholegrains.points<- ntile(HEI$wholegrains,11)-1

###dairy ### group cheese(cheese_20, cheese_40, cheese_48, cheese_other_cream_cheese, cheese_d, snack_cheese)
#            group dairy (breakfast_drink, milk_whole, milk_semiskimmed
#milk_skimmed, buttermilk, chocolatemilk, yoghurt_drink_added_sugar, yoghurt_drink_other, custard_ff
#yoghurt_ff, yoghurt_lf, yoghurt_lf_fruits, fromage_frais_fruits, porridge, other_custard_yoghurt_fromagefrais
#icecream_dairy, whipped_cream, coffeecreamer_ff, coffeecreamer_hf, coffeecreamer_p, milk_coffee) #21 foodsitems 
HEI$group_cheese<-(HEI$cheese_20 + HEI$cheese_40 + HEI$cheese_48 + HEI$cheese_other_cream_cheese + HEI$cheese_d + HEI$snack_cheese)
HEI$group_dairy <- (HEI$breakfast_drink + HEI$milk_whole + HEI$milk_semiskimmed + HEI$milk_skimmed + HEI$buttermilk +
                  HEI$chocolatemilk + HEI$yoghurt_drink_added_sugar + HEI$yoghurt_drink_other + HEI$custard_ff + 
                  HEI$yoghurt_ff + HEI$yoghurt_lf + HEI$yoghurt_lf_fruits + HEI$fromage_frais_fruits + HEI$porridge + 
                  HEI$other_custard_yoghurt_fromagefrais + HEI$icecream_dairy + HEI$whipped_cream + HEI$coffeecreamer_ff + 
                  HEI$coffeecreamer_hf + HEI$coffeecreamer_p + HEI$milk_coffee)
HEI$dairy<-(HEI$group_cheese + HEI$group_dairy)/Kcal
HEI$dairy.points<- ntile(HEI$dairy,11)-1

###Total Protein Foods ### group_meat, group_egg + #greensandbeans +# seafood and plant protein 
#group_meat; sausage_smoked_cold_cuts, meats_other_cold_cuts, meat, beef_lean, beef_fat, pork_lean, pork_fat
#           pork_processed, meats_fat, chicken, other_meat_poultry, gravy, snack_meats
HEI$group_meat <- (HEI$sausage_smoked_cold_cuts + HEI$meats_other_cold_cuts + HEI$meat + HEI$beef_lean + HEI$beef_fat + HEI$pork_lean + 
                  HEI$pork_fat + HEI$pork_processed + HEI$meats_fat + HEI$chicken + HEI$other_meat_poultry + HEI$gravy + HEI$snack_meats)
HEI$group_eggs <- HEI$egg_baked + HEI$egg_cooked
HEI$group_fish <- HEI$herring_salted + HEI$fish_white_fried + HEI$fish_lean + HEI$fish_fatty + HEI$fish_other + HEI$fish_prepared_fat 
HEI$group_nuts <- HEI$nut_d + HEI$snack_nut + HEI$peanutbutter


HEI$TPF<- (HEI$group_meat + HEI$group_eggs + HEI$legumes + HEI$group_fish + HEI$group_nuts)/Kcal 
HEI$TPF.points <-ntile(HEI$TPF,6)-1

###Seafood and plant proteins ### fish and nuts 
HEI$SFPP<- (HEI$group_fish+ HEI$group_nuts)/Kcal 
HEI$SFPP.points<-ntile(HEI$SFPP,6)-1

###Fatty acids ### no data 

###Refined Grains### "porridge" "pasta" ""rice" "cereals" "crackers" "rolls" 
HEI$refinedgrains<- (HEI$porridge + HEI$pasta + HEI$rice + HEI$cereals + HEI$crackers + HEI$rolls)/Kcal
HEI$refinedgrains.points<- ntile(desc(HEI$refinedgrains), 11) -1

###Sodium ### no data 

###Empty Calories### 
#This category is in kcal/day instead of gram/day 
#alcohol, coffee, tea, nonalc_drinks, pastry, prepared meal, sauces, savoury_snacks, sugar_sweets, sweet spreads, dairy 
HEI$group_alcohol <- HEI2$beer + HEI2$wine_red + HEI2$wine_white + HEI2$wine_fort + HEI2$spirits + HEI2$other_alc_drinks
HEI$group_nonalc_drinks <- rowSums(HEI2[,c("softdrink_sugar", "softdrink_no_sugar", "beer_af")]) 
HEI$group_pastry <- rowSums(HEI2[,c( "pasta", "biscuits_s", "cake", "pastry", "spiced_cake")])
HEI$group_prepared_meal <- rowSums(HEI2[,c("ready_meals_chinese_indian", "meals_fast_food", "ready_meals_other","pizza" )]) 
HEI$group_sauces <- (HEI2$sauce_hot + HEI2$mayonaise + HEI2$`Halvanaise_mayonnaise_non-red sauces` + HEI2$`Halvanaise_mayonnaise_non-red sauces_snack` + 
                    HEI2$saladdressing_f + HEI2$saladdressing_w + HEI2$mayonaise_snack)
HEI$group_savoury_snacks <- rowSums(HEI2[,c("sandwichspread_butter","snack_savoury_hot", "crisps", "salad_toast")]) 
HEI$group_sugar_sweets <- rowSums(HEI2[,c("chocolatespreads",  "spreads_sweet" ,"sugar_yoghurt", "sugar_coffee" ,
                                         "sugar_tea" , "candybars","chocolate", "sweets")])
HEI$group_dairy <- (HEI2$breakfast_drink + HEI2$milk_whole + HEI2$milk_semiskimmed + HEI2$milk_skimmed + HEI2$buttermilk +
                      HEI2$chocolatemilk + HEI2$yoghurt_drink_added_sugar + HEI2$yoghurt_drink_other + HEI2$custard_ff + 
                      HEI2$yoghurt_ff + HEI2$yoghurt_lf + HEI2$yoghurt_lf_fruits + HEI2$fromage_frais_fruits + HEI2$porridge + 
                      HEI2$other_custard_yoghurt_fromagefrais + HEI2$icecream_dairy + HEI2$whipped_cream + HEI2$coffeecreamer_ff + 
                      HEI2$coffeecreamer_hf + HEI2$coffeecreamer_p + HEI2$milk_coffee)  
HEI$emptycalories<- (HEI$group_alcohol + HEI2$coffee + HEI2$tea + HEI$group_nonalc_drinks
                     + HEI$group_pastry + HEI$group_prepared_meal + HEI$group_sauces
                     + HEI$group_savoury_snacks + HEI$group_sugar_sweets + HEI$group_dairy)


#Anything with over 50% of total energy from empty calories gets 0, 
#Anything with under 19% of total energy from empty calories gets 20 points. 
#everything inbetween should be divided in strata 

totalKCal<- HEI$SumOfkcal
HEI$emptycaloriesperc <- HEI$emptycalories / HEI$SumOfkcal
summary(HEI$emptycaloriesperc)
#from here I could not do with my database, so I am not sure if it works. 
x<- HEI$emptycaloriesperc
HEI$emptycalorie.points<- ifelse(x >= 0.5, 0, 
                              ifelse(x < 0.5 & x >= 0.4845, 1, 
                                     ifelse(x < 0.4845 & x>=0.469, 2,
                                            ifelse(x < 0.469 & x>= 0.4535 , 3,
                                                   ifelse(x < 0.4535 & x >=0.438, 4,
                                                          ifelse(x < 0.438 & x>= 0.4225, 5, 
                                                                 ifelse(x < 0.4225 & x>=0.407, 6,
                                                                        ifelse(x < 0.407 & x>= 0.3915, 7,
                                                                               ifelse(x < 0.3915 & x>= 0.376, 8,
                                                                                      ifelse(x < 0.376 & x>0.3605, 9,
                                                                                             ifelse(x < 0.3605 & x >= 0.345, 10, 
                                                                                                    ifelse(x < 0.345 & x>= 0.3295, 11,
                                                                                                           ifelse(x < 0.3295 & x>= 0.314 , 12,
                                                                                                                  ifelse(x < 0.314 & x >= 0.2985, 13,
                                                                                                                         ifelse(x < 0.2985 & x>= 0.283, 14, 
                                                                                                                                ifelse(x < 0.283 & x>=0.2675, 15,
                                                                                                                                       ifelse(x < 0.2675 & x>=0.252, 16,
                                                                                                                                              ifelse(x < 0.252 & x>= 0.2365, 17,
                                                                                                                                                     ifelse(x < 0.2365 & x >=0.221, 18,
                                                                                                                                                            ifelse(x < 0.221 & x>0.20559, 19, 20))))))))))))))))))))



#this is the way I calculated it, with dividing it into strata, like the other categories. 
#HEI$emptycalorie.points<- ntile(desc(HEI$emptycaloriesperc), 21)-1

###Total Score ### 
all.components <- cbind(HEI$totalfruit.points, HEI$wholefruit.points, HEI$totalvegetables.points, HEI$greensandbeans.points, 
                        HEI$wholegrains.points, HEI$dairy.points, HEI$TPF.points, HEI$SFPP.points, HEI$refinedgrains.points, HEI$emptycalorie.points)

HEI$hei<- (HEI$totalfruit.points + HEI$wholefruit.points + HEI$totalvegetables.points + HEI$greensandbeans.points +HEI$wholegrains.points+ HEI$dairy.points+ HEI$TPF.points+HEI$SFPP.points+ HEI$refinedgrains.points + HEI$emptycalorie.points)

HEI %>% select(Participant, hei) -> To_write
write_tsv(To_write,"TABLES/HEI_score.tsv")

#ggplot(HEI,aes(x=hei)) + geom_histogram() + theme_bw() -> PLOT
#ggsave(filename="Histogram_hei.pdf", plot=PLOT)

 #there is a maximum score of 80 and a minimum of 0, LLDEEP database had min 15 and max 70. 
