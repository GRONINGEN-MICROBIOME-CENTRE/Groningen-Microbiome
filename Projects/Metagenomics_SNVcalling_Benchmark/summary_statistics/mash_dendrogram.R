library(tidyverse)
library(reshape2)
library(ape)
# setwd to source file location
mash<-read.delim("../00.data/mash_s10000_k17.txt ", sep = "\t", quote = "", header = F)[1:2116,]
as_tibble(mash) %>% filter(! V3 == 0) %>% summarise(M = mean(V3), Max= max(V3), Min=min(V3) )
mash$V1 <- mash$V1 %>% 
  str_replace_all("/data/p290416/project/SNV_benchmark/test_genome_anno/SEQS/","") %>%
  str_replace_all("\\.fa", "")
mash$V2 <- mash$V2 %>% 
  str_replace_all("/data/p290416/project/SNV_benchmark/test_genome_anno/SEQS/","") %>%
  str_replace_all("\\.fa", "")

melt(mash[,1:3])

mash_mat <- mash[,1:3] %>% 
  reshape(idvar = "V1", timevar = "V2", direction = "wide") %>%
  data.frame(row.names = "V1")
colnames(mash_mat) <- str_replace_all(colnames(mash_mat), "V3\\.","")



mash_tree = nj(as.dist(mash_mat))

pdf("05.mash/mash_s10000_k17.pdf",width = 6, height = 10)
plot(mash_tree)
dev.off()
