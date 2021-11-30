library(tidyverse)



options("width"=200)

CD = read_tsv("Regression_input_CD.tsv")
UC = read_tsv("Regression_input_UC.tsv")


To_test = rbind(CD, UC) %>% distinct() 

Food_items = read_tsv("../TABLES/Table_participant_GDAG.tsv_script_eva_input")
Food_items %>% mutate(fruits = fruitjuice + fruit, vegetables = vegetables_cooked_nobutter+ vegetables_stirfried + vegetables_cooked_butter, nuts=peanutbutter + nut_d + snack_nut, grains = bread, 
fish= herring_salted+fish_white_fried+fish_lean+fish_fatty+fish_other+fish_prepared_fat, meat=sausage_smoked_cold_cuts+meats_other_cold_cuts+meat+beef_lean+beef_fat+pork_processed+snack_meats+pork_fat+pork_lean+meats_fat, alcohol = wine_red+wine_white+wine_fort+spirits+other_alc_drinks) %>% select(Participant, fruits, nuts, grains, fish, meat, alcohol, vegetables, legumes ) -> Food_items


left_join(To_test, Food_items, by="Participant") -> To_test

Final_food_group = tibble()
Final_food_group_Corrected = tibble()
for (Food_Group in colnames(select(Food_items, -Participant)) ){
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


