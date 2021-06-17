library(tidyverse)
multi = read_tsv("../Resilio Sync/SNV_sync/data/Benchmark_multiref.tsv")%>% select(-X9) %>% mutate(Specificity = TN/(TN+FN), Sensitivity = TP/(TP+FN), Precision=TP/(TP+FP) ) 
uni = read_tsv("../Resilio Sync/SNV_sync/data/Benchmark_uniref.tsv")%>% select(-X9) %>% mutate(Specificity = TN/(TN+FN), Sensitivity = TP/(TP+FN), Precision=TP/(TP+FP) ) 
Stats = read_tsv("../Resilio Sync/SNV_sync/table/Supplementary_table_1.tsv") 

Stats %>%mutate(Taxa = Species ) %>% select(Taxa, `Number contigs`, N50 ) -> Stats2
left_join(Stats2, select(multi, c(Taxa,Sensitivity, Precision, Tool,Benchmark_Design)), by="Taxa") -> Stats_multi
left_join(Stats2, select(uni, c(Taxa,Sensitivity, Precision,Tool,Benchmark_Design)), by="Taxa") -> Stats_uni


Stats_multi%>% mutate(N50 = as.vector(scale(N50)), `Number contigs` = as.vector(scale(`Number contigs`))) -> Stats_multi
Stats_uni%>% mutate(N50 = as.vector(scale(N50)), `Number contigs` = as.vector(scale(`Number contigs`))) -> Stats_uni
rbind(Stats_multi, Stats_uni) -> Combined_stats

###Effect of Number contigs
#senstivity
lm(Sensitivity ~  `Number contigs`*Tool*Benchmark_Design  , Combined_stats) -> m2
lm(Sensitivity ~  `Number contigs`+Tool*Benchmark_Design , Combined_stats) -> m1
lm(Sensitivity ~  Tool * Benchmark_Design , Combined_stats) -> m0
anova(m0,m1) ; anova(m1,m2)
#Precision
lm(Precision ~  `Number contigs`*Tool*Benchmark_Design  , Combined_stats) -> m2_p
lm(Precision ~  `Number contigs`+Tool*Benchmark_Design , Combined_stats) -> m1_p
lm(Precision ~  Tool * Benchmark_Design , Combined_stats) -> m0_p
anova(m0_p,m1_p) ; anova(m1_p,m2_p)

###Effect N50
lm(Sensitivity ~  N50*Tool*Benchmark_Design  , Combined_stats) -> m2.2
lm(Sensitivity ~  N50+Tool*Benchmark_Design , Combined_stats) -> m1.2
anova(m1.2,m0)
#Precision
lm(Precision ~  N50*Tool*Benchmark_Design  , Combined_stats) -> m2.2_p
lm(Precision ~  N50+Tool*Benchmark_Design , Combined_stats) -> m1.2_p
anova(m0_p,m1.2_p)
