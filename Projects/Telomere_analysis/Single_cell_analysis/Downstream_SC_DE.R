library(tidyverse)
library(clusterProfiler)
library('org.Hs.eg.db')
library(biomaRt)

setwd("~/Resilio Sync/LLD phenotypes for Sergio/scripts/Pathway_enrichment/")



#Second file

Subsampling_test = function(Test_Data, n=20){
  cHR_counts = table(Test_Data$DE, Test_Data$chromosome_name)
  Test_Data %>% filter(DE == 0) -> To_sample
  set.seed(100)
  Keep= c()
  for (Chromosome in  colnames(cHR_counts)){
    if (! Chromosome %in% filter(Test_Data, DE=T)$chromosome_name  ){ next }
    To_sample %>% filter(chromosome_name == Chromosome) -> Chr_samples
    Chr_samples[sample(nrow(Chr_samples), round(cHR_counts[2,][Chromosome])*n),] -> Sampled
    Keep = c(Keep, Sampled$entrezgene_id)
  }
  To_sample %>% filter(! entrezgene_id %in% Keep) -> To_remove
  Test_Data %>% filter(! entrezgene_id %in% To_remove$entrezgene_id) -> Chr_matched
  return(Chr_matched)
}

Just_distance = function(DA){
  DA %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG #filter(p_val_adj > 0.05) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG
  
  read_tsv("BG_missing.txt", col_names = F) -> Correct
  left_join(BG, Correct) -> BG
  BG %>% mutate(X1 = ifelse( is.na(X2), X1, X2 )) -> BG
  BG_g = mapIds(org.Hs.eg.db, BG$X1, 'ENTREZID', 'SYMBOL')
  
  #1. Mapping to genomic positions
  biolist <- as.data.frame(listMarts())
  ensembl=useMart("ensembl")
  esemblist <- as.data.frame(listDatasets(ensembl)) #Human genes (GRCh38.p13) GRCh38.p13
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  
  t2g<-getBM(attributes=c('entrezgene_id','chromosome_name','start_position','end_position'), mart = ensembl)
  my_ids_bg <- data.frame(entrezgene_id=BG_g)
  my_ids4 <- merge(my_ids_bg, t2g, by= 'entrezgene_id')
  my_ids4 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids4
  my_ids4 %>% filter(chromosome_name %in% c("X", as.character(seq(1:22)))) -> my_ids4
  my_ids4[!duplicated(my_ids4$entrezgene_id),] -> my_ids4
  
  
  #3. Get chromosome length as a proxy for telomere
  Chr_length = tibble(chromosome_name= c(as.character(seq(1:22)), "X"), L=c(248956422, 242193529, 198295559,190214555, 181538259,
                                                                            170805979, 159345973 ,145138636, 138394717, 133797422, 135086622,
                                                                            133275309, 114364328, 107043718, 101991189, 90338345,
                                                                            83257441, 80373285, 58617616, 64444167,46709983 ,50818468, 156040895) )
  left_join(my_ids4, Chr_length) -> Test_Data
  Test_Data %>% mutate(Distance_tel = pmin( start_position, (L-end_position) )) -> Test_Data
  
  #4. Add names
  BG %>% mutate(entrezgene_id= BG_g) %>% as_tibble() %>%  dplyr::select(-X2) -> All
  rbind(All, tibble(X1=c("MT-ND5","ATP5MJ"), entrezgene_id = c("4540", "9556"))) -> All
  left_join(Test_Data, All) -> Test_Data
  
  Test_Data %>% filter(chromosome_name %in% c(as.character(seq(1:22)), "X")) -> Test_Data
  
  Sig = DA %>% filter(p_val_adj < 0.05)
  Test_Data %>% mutate(status = ifelse(X1 %in% filter(Sig, logFC>0)$gene_symbol , "up", ifelse(X1 %in% filter(Sig, logFC<0)$gene_symbol , "down", "no" ))) -> Test_Data
  mutate(Test_Data, DE= ifelse(status=="no", 0, 1)) -> Test_Data
  Test_Data %>% mutate(Subtelomeric = ifelse(Distance_tel< 4*10^6, 1, 0 )) -> Test_Data
  Test_Data %>% mutate(Long_tel_effect = ifelse(Distance_tel < 10^7, 1, 0 )) -> Test_Data
  
  
  #5. Plots
  Test_Data$status = ordered(as.factor(Test_Data$status), c("no", "down", "up"))
  #Test_Data %>% mutate(DE = ifelse(status == "no", F, T ) ) -> Test_Data
  
  Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=status)) + geom_density(alpha=0.5) + theme_bw() +
    geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(3))  -> Distr_fig
  
  Test_Data %>% drop_na() %>% ggplot(aes(y=log10(Distance_tel), x=status, fill=status)) + 
    ggforce::geom_sina(alpha=0.3, col="grey")+ geom_boxplot(alpha=0.5,outlier.shape = NA) + theme_bw() + geom_hline(yintercept = log10(4*10^6),linetype="dotted")  + geom_hline(yintercept = log10(10^7)) +
    scale_fill_manual(values = viridis::viridis(3)) + coord_flip() + labs(fill= "DE") +ylab("log10(Distance to telomere (bp))") + xlab("Gene group") -> Box_fig
  
  Fig = Distr_fig + Box_fig
  return(list(Fig, Test_Data))
  
}

Analysis = function(DA, Subsampling=F, n = 20, Significance=0.05){
  DA %>% filter(p_val_adj < Significance) %>% filter(logFC > 0) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> Up
  DA %>% filter(p_val_adj < Significance) %>% filter(logFC < 0) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> Down
  DA %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG #filter(p_val_adj > 0.05) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG
  
  #N=3 missing LINC00152, AIM1 and TMEM66
  Up = Up %>% mutate(X1 = ifelse(X1=="C7orf55-LUC7L2" ,"FMC1-LUC7L2", X1),
                     X1 = ifelse(X1=="C17orf89" ,"NDUFAF8", X1),
                     X1 = ifelse(X1=="ATP5L" ,"ATP5MG", X1),
                     X1 = ifelse(X1== "LINC00152", "CYTOR", X1), X1 = ifelse(X1=="TMEM66" ,"SARAF", X1),
                     X1 = ifelse(X1== "AIM1", "CRYBG1", X1),
                     X1 = ifelse(X1 == "GNB2L1", "RACK1",X1))
  Down = Down %>% mutate(X1 = ifelse(X1=="H3F3B" ,"H3-3B", X1),
                         X1 = ifelse(X1== "LINC00152", "CYTOR", X1), X1 = ifelse(X1=="TMEM66" ,"SARAF", X1),
                         X1 = ifelse(X1=="C14orf2" ,"ATP5MJ", X1), X1 = ifelse(X1== "AIM1", "CRYBG1", X1),
                         X1 = ifelse(X1=="SEPT7" ,"SEPTIN7", X1))
  #BG = read_tsv("v2_bg.txt", col_names = F)
  
  read_tsv("BG_missing.txt", col_names = F) -> Correct
  left_join(BG, Correct) -> BG
  BG %>% mutate(X1 = ifelse( is.na(X2), X1, X2 )) -> BG
  
  Up_g  = mapIds(org.Hs.eg.db, Up$X1, 'ENTREZID', 'SYMBOL')
  Down_g  = mapIds(org.Hs.eg.db, Down$X1, 'ENTREZID', 'SYMBOL')
  BG_g = mapIds(org.Hs.eg.db, BG$X1, 'ENTREZID', 'SYMBOL')
  
  #1. Mapping to genomic positions
  biolist <- as.data.frame(listMarts())
  ensembl=useMart("ensembl")
  esemblist <- as.data.frame(listDatasets(ensembl)) #Human genes (GRCh38.p13) GRCh38.p13
  ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
  filters = listFilters(ensembl)
  attributes = listAttributes(ensembl)
  
  t2g<-getBM(attributes=c('entrezgene_id','chromosome_name','start_position','end_position'), mart = ensembl)
  
  
  my_ids <- data.frame(entrezgene_id=Up_g)
  my_ids2 <- merge(my_ids, t2g, by= 'entrezgene_id')
  my_ids2 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids2
  
  
  my_ids_down <- data.frame(entrezgene_id=Down_g)
  my_ids3 <- merge(my_ids_down, t2g, by= 'entrezgene_id')
  my_ids3 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids3
  
  my_ids_bg <- data.frame(entrezgene_id=BG_g)
  my_ids4 <- merge(my_ids_bg, t2g, by= 'entrezgene_id')
  my_ids4 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids4
  my_ids4 %>% filter(chromosome_name %in% c("X", as.character(seq(1:22)))) -> my_ids4
  my_ids4[!duplicated(my_ids4$entrezgene_id),] -> my_ids4
  
  
  #2. Add all of them
  rbind(mutate(my_ids2, status="up"), mutate(my_ids3, status="down")) -> Affected
  mutate(my_ids4, status="no") %>% filter(!entrezgene_id %in% Affected$entrezgene_id) -> nO_affected
  
  #3. Get chromosome length as a proxy for telomere
  Chr_length = tibble(chromosome_name= c(as.character(seq(1:22)), "X"), L=c(248956422, 242193529, 198295559,190214555, 181538259,
                                                                            170805979, 159345973 ,145138636, 138394717, 133797422, 135086622,
                                                                            133275309, 114364328, 107043718, 101991189, 90338345,
                                                                            83257441, 80373285, 58617616, 64444167,46709983 ,50818468, 156040895) )
  left_join(rbind(Affected, nO_affected), Chr_length) -> Test_Data
  Test_Data %>% mutate(Distance_tel = pmin( start_position, (L-end_position) )) -> Test_Data
  
  #4. Add names
  Up %>% mutate(entrezgene_id= Up_g) -> Upregulated
  Down %>% mutate(entrezgene_id= Down_g) -> Downregulated
  BG %>% mutate(entrezgene_id= BG_g) %>% as_tibble() %>%  dplyr::select(-X2) -> All
  rbind(All, tibble(X1=c("MT-ND5","ATP5MJ"), entrezgene_id = c("4540", "9556"))) -> All
  left_join(Test_Data, All) -> Test_Data
  
  Test_Data %>% filter(chromosome_name %in% c(as.character(seq(1:22)), "X")) -> Test_Data
  mutate(Test_Data, DE= ifelse(status=="no", 0, 1)) -> Test_Data
  Test_Data %>% mutate(Subtelomeric = ifelse(Distance_tel< 4*10^6, 1, 0 )) -> Test_Data
  Test_Data %>% mutate(Long_tel_effect = ifelse(Distance_tel < 10^7, 1, 0 )) -> Test_Data
  
  
  #5. Plots
  Test_Data$status = ordered(as.factor(Test_Data$status), c("no", "down", "up"))
  Test_Data %>% mutate(DE = ifelse(status == "no", F, T ) ) -> Test_Data
  
  Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=status)) + geom_density(alpha=0.5) + theme_bw() +
    geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(3))  #viridis::scale_fill_viridis()
  
  Test_Data %>% drop_na() %>% ggplot(aes(y=log10(Distance_tel), x=status, fill=status)) + 
    ggforce::geom_sina(alpha=0.3, col="grey")+ geom_boxplot(alpha=0.5,outlier.shape = NA) + theme_bw() + geom_hline(yintercept = log10(4*10^6),linetype="dotted")  + geom_hline(yintercept = log10(10^7)) +
    scale_fill_manual(values = viridis::viridis(3)) + coord_flip() + labs(fill= "DE") +ylab("log10(Distance to telomere (bp))") + xlab("Gene group") -> Box_fig
  
  Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=as.factor(DE))) + geom_density(alpha=0.5) + theme_bw() +
    geom_vline(xintercept = log10(4*10^6),linetype="dotted")  + geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(2)) + xlab("log10(Distance to telomere (bp))") +
    labs(fill = "DE") -> Distr_fig
  
  library(patchwork)
  Distr_fig + Box_fig + plot_annotation(tag_levels = 'A') -> Fig
  
  #Statistical test
  summary(glm(DE ~ log10(Distance_tel), Test_Data , family= binomial(link="logit"))) -> Sum1 # logistic
  kruskal.test(Test_Data$Distance_tel, Test_Data$status) -> Sum2
  dunn.test::dunn.test(Test_Data$Distance_tel, Test_Data$status) -> Sum3
  
  filter(Test_Data, ! status == "up" ) -> Test_down ; fisher.test(table(Test_down$DE, Test_down$Subtelomeric) ) -> Sum4
  filter(Test_Data, ! status == "down" ) -> Test_up ; fisher.test(table(Test_up$DE, Test_up$Subtelomeric) ) -> Sum5
  
  #fisher.test(table(Test_Data$DE, Test_Data$Subtelomeric) )
  
  #Subsampling
  if (Subsampling == T){
    Chr_matched = Subsampling_test(Test_Data, n = n)
    table(Chr_matched$DE, Chr_matched$chromosome_name)
    wilcox.test(Chr_matched$Distance_tel ~ Chr_matched$DE) -> Sum6 #Non parametric
    kruskal.test(Chr_matched$Distance_tel, Chr_matched$status) -> Sum7
    dunn.test::dunn.test(Chr_matched$Distance_tel, Chr_matched$status) -> Sum8
    
    Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=status)) + geom_density(alpha=0.5) + theme_bw() +
      geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(3))  -> subfig2
    
    Test_Data %>% drop_na() %>% ggplot(aes(y=log10(Distance_tel), x=status, fill=status)) + 
      ggforce::geom_sina(alpha=0.3, col="grey")+ geom_boxplot(alpha=0.5,outlier.shape = NA) + theme_bw() + geom_hline(yintercept = log10(4*10^6),linetype="dotted")  + geom_hline(yintercept = log10(10^7)) +
      scale_fill_manual(values = viridis::viridis(3)) + coord_flip() + labs(fill= "DE") +ylab("log10(Distance to telomere (bp))") + xlab("Gene group") -> subfig3
    
    Fig2 = subfig2 + subfig3
    
    filter(Test_Data, ! status == "up" ) -> Test_down ; fisher.test(table(Test_down$DE, Test_down$Subtelomeric) ) -> Sum9
    filter(Test_Data, ! status == "down" ) -> Test_up ; fisher.test(table(Test_up$DE, Test_up$Subtelomeric) ) -> Sum10
    
    Final = list(Test_Data, Sum1, Sum2, Sum3, Sum4, Sum5, Fig, Fig2,Sum6, Sum7, Sum8, Sum9, Sum10 )
  } else {
    
    Final = list(Test_Data, Sum1, Sum2, Sum3,Sum4, Sum5, Fig )
  }
  
    
  
  return(Final)
}


dea.list_formatted <- readRDS('dea.list_formatted.rds') #read DEA output list
dea.list <- dea.list_formatted$v2 #pick only v2 results
dea.list_notna <- lapply(dea.list, function(x) x[!is.na(x$logFC),]) #remove NA estimates

N= names(dea.list_notna)
for (cell_n in length(N) ){
  Name = N[cell_n]
  DF = as_tibble(as.data.frame(dea.list_notna[cell_n]))
  colnames(DF) = c("gene_symbol", "logFC", "p_val",  "p_val_adj")
  
  if (Name %in% c("CD8Tcells", "Tmemorycells", "Tnaivecells")){ 
    Result_cell = Just_distance(DF)
    next 
  } else if( cell_n == 1){
    Result_cell = Analysis(DF, Subsampling = T)
  } else {
    Result_cell = Analysis(DF, Subsampling = T, n = 11)
  }

  as.data.frame(Result_cell[[2]]$coefficients)["log10(Distance_tel)","Pr(>|z|)" ] #For DE
  Enrichment_P = Result_cell[[6]]$p.value #For up
  Enrichment_P_corrected = Result_cell[[13]]$p.value #For up
}

#N=3 missing LINC00152, AIM1 and TMEM66



########
read_tsv("~/../Downloads/dea.T_cells_naive.notna.txt") -> broad
broad %>% mutate(gene_symbol = ifelse(gene_symbol == "C1orf228", "ARMH1", gene_symbol)) -> broad
Analysis(broad) -> Test_Data2
left_join(Test_Data2, Correct) %>% mutate(`Original name` = ifelse( is.na(`Original name`), X1, `Original name` ) ) -> Test_Data2
left_join(Test_Data2, Correct) %>% mutate(`Original name` = ifelse( `Original name`=="ARMH1","C1orf228", `Original name` ) ) -> Test_Data2
Test_Data2 %>% filter(! duplicated(entrezgene_id)) -> Test_Data2
write_tsv(Test_Data2, "All_gene_info_broad.tsv")
#######



























DA = read_tsv("~/Downloads/dea_ct.txt")
DA %>% filter(p_val_adj < 0.05) %>% filter(logFC > 0) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> Up
DA %>% filter(p_val_adj < 0.05) %>% filter(logFC < 0) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> Down
DA %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG #filter(p_val_adj > 0.05) %>% mutate(X1 = gene_symbol ) %>% dplyr::select(X1) -> BG

Up = Up %>% mutate(X1 = ifelse(X1=="C7orf55-LUC7L2" ,"FMC1-LUC7L2", X1),
                   X1 = ifelse(X1=="C17orf89" ,"NDUFAF8", X1),
                   X1 = ifelse(X1=="ATP5L" ,"ATP5MG", X1),
                   X1 = ifelse(X1 == "GNB2L1", "RACK1",X1))
Down = Down %>% mutate(X1 = ifelse(X1=="H3F3B" ,"H3-3B", X1),
                       X1 = ifelse(X1=="C14orf2" ,"ATP5MJ", X1),
                       X1 = ifelse(X1=="SEPT7" ,"SEPTIN7", X1))
#BG = read_tsv("v2_bg.txt", col_names = F)

read_tsv("BG_missing.txt", col_names = F) -> Correct
left_join(BG, Correct) -> BG
BG %>% mutate(X1 = ifelse( is.na(X2), X1, X2 )) -> BG

Up_g  = mapIds(org.Hs.eg.db, Up$X1, 'ENTREZID', 'SYMBOL')
Down_g  = mapIds(org.Hs.eg.db, Down$X1, 'ENTREZID', 'SYMBOL') ; Down_g["MT-ND5"] = "4540"  ; Down_g["ATP5MJ"] = "9556"
BG_g = mapIds(org.Hs.eg.db, BG$X1, 'ENTREZID', 'SYMBOL')



#########################
###Enrichment analysis###
#########################

##Try with BG 

#############Kegg############
search_kegg_organism('hsa', by='kegg_code')
kk <- enrichKEGG(gene = Up_g, organism = 'hsa', universe = unique(BG_g), pAdjustMethod = "fdr",qvalueCutoff = 0.5,pvalueCutoff = 0.5)
kk2 <- enrichKEGG(gene = Down_g, organism = 'hsa', universe = unique(BG_g), pAdjustMethod = "fdr",qvalueCutoff = 0.5,pvalueCutoff = 0.5)
as_tibble(kk2) ; as_tibble(kk)
###########Go biological process ############
Go_1 = enrichGO(gene  = Up_g,
                universe =  unique(BG_g), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff  = 0.5, qvalueCutoff  = 0.5, readable      = TRUE)
GO_2  = enrichGO(gene  = Down_g,
                 universe =  unique(BG_g), OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff  = 0.5, qvalueCutoff  = 0.5, readable      = TRUE)
as_tibble(Go_1) ; as_tibble(GO_2)
######################

All_DE = c(Up_g, Down_g)
as_tibble(enrichKEGG(gene = All_DE, organism = 'hsa', universe = BG_g, pAdjustMethod = "fdr",qvalueCutoff = 0.5,pvalueCutoff = 0.5))
as_tibble(enrichGO(gene  = All_DE, universe =  BG_g, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff  = 0.5, qvalueCutoff  = 0.5, readable      = TRUE))


#without BG set
####KEGG#############
kk <- enrichKEGG(gene = Up_g, organism = 'hsa', pAdjustMethod = "fdr")
as_tibble(kk)
kk2 <- enrichKEGG(gene = Down_g, organism = 'hsa', pAdjustMethod = "fdr")
as_tibble(kk2)
#########GO############
Go_1 = enrichGO(gene  = Up_g,
                OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable      = TRUE)
GO_2  = enrichGO(gene  = Down_g,
                 OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff  = 0.01, qvalueCutoff  = 0.05, readable      = TRUE)
as_tibble(Go_1) ; as_tibble(GO_2)




####Are DE genes predominant in subtelomeric regions?

#1. Mapping to genomic positions
biolist <- as.data.frame(listMarts())
ensembl=useMart("ensembl")
esemblist <- as.data.frame(listDatasets(ensembl)) #Human genes (GRCh38.p13) GRCh38.p13
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
attributes = listAttributes(ensembl)

t2g<-getBM(attributes=c('entrezgene_id','chromosome_name','start_position','end_position'), mart = ensembl)


my_ids <- data.frame(entrezgene_id=Up_g)
my_ids2 <- merge(my_ids, t2g, by= 'entrezgene_id')
my_ids2 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids2
#add SNH69 which cannot be found
rbind(my_ids2, tibble(entrezgene_id="735301", chromosome_name="16", start_position=1964959, end_position=1965509)) -> my_ids2

my_ids_down <- data.frame(entrezgene_id=Down_g)
my_ids3 <- merge(my_ids_down, t2g, by= 'entrezgene_id')
my_ids3 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids3

my_ids_bg <- data.frame(entrezgene_id=BG_g)
my_ids4 <- merge(my_ids_bg, t2g, by= 'entrezgene_id')
my_ids4 %>% as_tibble()  %>% filter(!is.na(entrezgene_id)) -> my_ids4
my_ids4 %>% filter(chromosome_name %in% c("X", as.character(seq(1:22)))) -> my_ids4
my_ids4[!duplicated(my_ids4$entrezgene_id),] -> my_ids4


#2. Add all of them
rbind(mutate(my_ids2, status="up"), mutate(my_ids3, status="down")) -> Affected
mutate(my_ids4, status="no") %>% filter(!entrezgene_id %in% Affected$entrezgene_id) -> nO_affected

#3. Get chromosome length as a proxy for telomere
Chr_length = tibble(chromosome_name= c(as.character(seq(1:22)), "X"), L=c(248956422, 242193529, 198295559,190214555, 181538259,
                                                                          170805979, 159345973 ,145138636, 138394717, 133797422, 135086622,
                                                                          133275309, 114364328, 107043718, 101991189, 90338345,
                                                                          83257441, 80373285, 58617616, 64444167,46709983 ,50818468, 156040895) )
left_join(rbind(Affected, nO_affected), Chr_length) -> Test_Data
Test_Data %>% mutate(Distance_tel = pmin( start_position, (L-end_position) )) -> Test_Data

#4. Add names
Up %>% mutate(entrezgene_id= Up_g) -> Upregulated
Down %>% mutate(entrezgene_id= Down_g) -> Downregulated
BG %>% mutate(entrezgene_id= BG_g) %>% as_tibble() %>%  dplyr::select(-X2) -> All
rbind(All, tibble(X1=c("MT-ND5","ATP5MJ"), entrezgene_id = c("4540", "9556"))) -> All
left_join(Test_Data, All) -> Test_Data

Test_Data %>% filter(chromosome_name %in% c(as.character(seq(1:22)), "X")) -> Test_Data
mutate(Test_Data, DE= ifelse(status=="no", 0, 1)) -> Test_Data
Test_Data %>% mutate(Subtelomeric = ifelse(Distance_tel< 4*10^6, 1, 0 )) -> Test_Data
Test_Data %>% mutate(Long_tel_effect = ifelse(Distance_tel < 10^7, 1, 0 )) -> Test_Data
Test_Data %>% filter(! duplicated(entrezgene_id)) -> Test_Data

#5. Plots
Test_Data$status = ordered(as.factor(Test_Data$status), c("no", "down", "up"))

#####Distribution by chr
#Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=status)) + geom_density(alpha=0.5) + theme_bw() +
#  facet_wrap(~chromosome_name) + geom_vline(x=4*10^6)
#####Distribution all chr
Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=status)) + geom_density(alpha=0.5) + theme_bw() +
  geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(3))  #viridis::scale_fill_viridis()

Test_Data %>% drop_na() %>% ggplot(aes(y=log10(Distance_tel), x=status, fill=status)) + 
  ggforce::geom_sina(alpha=0.3, col="grey")+ geom_boxplot(alpha=0.5,outlier.shape = NA) + theme_bw() + geom_hline(yintercept = log10(4*10^6),linetype="dotted")  + geom_hline(yintercept = log10(10^7)) +
  scale_fill_manual(values = viridis::viridis(3)) + coord_flip() + labs(fill= "DE") +ylab("log10(Distance to telomere (bp))") + xlab("Gene group") -> Box_fig

Test_Data %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=as.factor(DE))) + geom_density(alpha=0.5) + theme_bw() +
  geom_vline(xintercept = log10(4*10^6),linetype="dotted")  + geom_vline(xintercept = log10(10^7)) +scale_fill_manual(values = viridis::viridis(2)) + xlab("log10(Distance to telomere (bp))") +
  labs(fill = "DE") -> Distr_fig

library(patchwork)
Distr_fig + Box_fig + plot_annotation(tag_levels = 'A')

#6. Statistical test
#Logistic - group specific
summary(glm(as.factor(status) ~ log10(Distance_tel), filter(Test_Data, ! status == "down" ), family= binomial(link="logit"))) 
summary(glm(as.factor(status) ~ log10(Distance_tel), filter(Test_Data, ! status == "up" ), family= binomial(link="logit")))
#Group test
kruskal.test(Test_Data$Distance_tel, Test_Data$status)
dunn.test::dunn.test(Test_Data$Distance_tel, Test_Data$status)

#All DE test
wilcox.test(Test_Data$Distance_tel ~ Test_Data$DE) #Non parametric
summary(glm(DE ~ log10(Distance_tel), Test_Data , family= binomial(link="logit"))) # logistic


#Enrichment test
filter(Test_Data, ! status == "up" ) -> Test_down ; fisher.test(table(Test_down$DE, Test_down$Subtelomeric) )
filter(Test_Data, ! status == "down" ) -> Test_up ; fisher.test(table(Test_up$DE, Test_up$Subtelomeric) )
fisher.test(table(Test_Data$DE, Test_Data$Subtelomeric))
chisq.test(Test_Data$DE, Test_Data$Subtelomeric)


#Matching BG so that the chromosome proportion is the same with DE
cHR_counts = table(Test_Data$DE, Test_Data$chromosome_name)
Test_Data %>% filter(DE == 0) -> To_sample
set.seed(100)
Keep= c()
for (Chromosome in  colnames(cHR_counts)){
  To_sample %>% filter(chromosome_name == Chromosome) -> Chr_samples
  Chr_samples[sample(nrow(Chr_samples), round(cHR_counts[2,][Chromosome])*11),] -> Sampled
  Keep = c(Keep, Sampled$entrezgene_id)
}
To_sample %>% filter(! entrezgene_id %in% Keep) -> To_remove
Test_Data %>% filter(! entrezgene_id %in% To_remove$entrezgene_id) -> Chr_matched
table(Chr_matched$DE, Chr_matched$chromosome_name)
wilcox.test(Chr_matched$Distance_tel ~ Chr_matched$DE) #Non parametric
kruskal.test(Chr_matched$Distance_tel, Chr_matched$status)
dunn.test::dunn.test(Chr_matched$Distance_tel, Chr_matched$status)


Chr_matched %>% drop_na() %>% ggplot(aes(x=log10(Distance_tel), fill=as.factor(DE))) + geom_density(alpha=0.5) + theme_bw() +
  geom_vline(xintercept = log10(4*10^6))



colnames(Correct) = c("Original name", "X1")
left_join(Test_Data, Correct) %>% mutate(`Original name` = ifelse( is.na(`Original name`), X1, `Original name` ) ) -> Test_Data
write_tsv(Test_Data,"All_gene_info.tsv")



Test_Data = read_tsv("All_gene_info.tsv")
Test_Data %>% filter(Distance_tel < 4*10^6)









