#Protein Score LL Cluster 

#Protein score is sum of 1) percentage of energy from total protein 
#                        2) percentage of energy from plant to animal protein ratio  

#Protein score - quantiles for protein to all calories 
library(tidyverse)
proteinscore = read_tsv("TABLES/Table_participant_GDAG.tsv_script_eva_input")

#Make database proteinscore, rows is samples, columns is foodgroups 
#proteinscore<-phenosdiet[,c("SUMOFEIWITTOT", "SUMOFEIWITPLANT", "SUMOFEIWITDIER","SUMOFKCAL") ]

# 1. total protein score 
# Protein (kcal) to all Kcal = energy-%, calculate  by multiplying by 4
proteinscore$ProttoKcal <- ((proteinscore$SumOfeiwittot*4)/(proteinscore$SumOfkcal))*100 #= identifcal to Prot_en

#2. ntile gives scorepoint 
proteinscore$TotProteinscore <- ntile(proteinscore$ProttoKcal, 11) -1   #score from zero to 10, so 11 bins

#3. protein plant to protein animal 
proteinscore$PlanttoAnimal <- (proteinscore$SumOfeiwitplant*4)/(proteinscore$SumOfeiwitdier*4)  

#5. protein animal to protein pant
proteinscore$ScorePlantAnimalProt <- ntile(proteinscore$PlanttoAnimal, 11) - 1

#6. sum protein score and plant protein score
proteinscore$ScoreProtein <- proteinscore$TotProteinscore + proteinscore$ScorePlantAnimalProt
 
#make histogram of distribution of points [max 20 points(2*10), min 0 points] 
proteinscore %>% select(Participant, ScoreProtein) -> To_write
write_tsv(To_write,"TABLES/Protein_score.tsv")

#ggplot(proteinscore,aes(x=ScoreProtein)) + geom_histogram() + theme_bw() -> PLOT
#ggsave(filename="Histogram_protein.pdf", plot=PLOT)
