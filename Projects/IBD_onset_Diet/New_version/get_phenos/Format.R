library(tidyverse)
options("width"=200)

#ID      colitisulcerosa_followup_adu_q_1_1.b    crohns_followup_adu_q_1_1.b     crohns_presence_adu_q_1_1       colitisulcerosa_presence_adu_q_1_1      colitisulcerosa_followup_adu_q_1_2      crohns_followup_adu_q_1_2       colitisulcerosa_followup_adu_q_1_1.c    crohns_followup_adu_q_1_1.c     colitisulcerosa_followup_adu_q_1_3      crohns_followup_adu_q_1_3
Data = read_tsv("IBD_dev.tsv")

Data[is.na(Data)] = 0
Data[Data==2] = 0
Data %>% summarise_if(is.numeric, sum) %>% t()

print("Total Data size")
print(dim(Data))

Healthy = Data
print("Filtering IBD only")
Data %>% filter(crohns_followup_adu_q_1_2  == 1 | colitisulcerosa_followup_adu_q_1_2 == 1 | crohns_followup_adu_q_1_3  == 1 | colitisulcerosa_followup_adu_q_1_3 == 1 | colitisulcerosa_followup_adu_q_1_1.b == 1 |  crohns_followup_adu_q_1_1.b == 1 |  colitisulcerosa_followup_adu_q_1_1.c ==1 | crohns_followup_adu_q_1_1.c == 1 )  -> Data
Data[Data== 2] = 0
print(dim(Data))

print("Filtering IBD people in BL")
Data %>% filter( ! crohns_presence_adu_q_1_1 == 1) %>% filter(! colitisulcerosa_presence_adu_q_1_1 == 1) -> Data
Data %>% select( -c(crohns_presence_adu_q_1_1, colitisulcerosa_presence_adu_q_1_1)) -> Data
print(dim(Data))

Healthy %>%  filter( ! crohns_presence_adu_q_1_1 == 1) %>% filter(! colitisulcerosa_presence_adu_q_1_1 == 1) -> Healthy

Healthy %>% filter(! ID %in%  Data$ID ) -> Healthy
write_tsv(Healthy, "IBD_Healthy.tsv")

print("How many developers per timepoint")
Data %>% summarise_if(is.numeric, sum) %>% t() 
print("Saving")
write_tsv(Data, "IBD_developers.tsv")
q()


Data %>% mutate(crohns_presence_adu_q_1 = ifelse(crohns_presence_adu_q_1 == 1, 1, NA), colitisulcerosa_presence_adu_q_1 = ifelse(colitisulcerosa_presence_adu_q_1 == 1, 1, NA), colitisulcerosa_followup_adu_q_1 = ifelse(colitisulcerosa_followup_adu_q_1 == 2, 0, colitisulcerosa_followup_adu_q_1),  crohns_followup_adu_q_1 = ifelse(crohns_followup_adu_q_1 == 2, 0, crohns_followup_adu_q_1) ) -> Data

Data %>% filter(! (crohns_presence_adu_q_1 == 1  |   colitisulcerosa_presence_adu_q_1 == 1) ) %>%
 filter(crohns_followup_adu_q_1  == 1 | colitisulcerosa_followup_adu_q_1 == 1) -> Data
print(Data)

