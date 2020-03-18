#### find taxa that did not run in the gwas analyisis


library(data.table)
library(tidyverse)

gwas.list<-fread("/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa/gwas.list",data.table = F,header=F)
mb.list<-fread("/groups/umcg-ugli/tmp01/umcg-elopera/testDAG//features_and_phenotypes/all.taxa.v2.list",data.table = F,header=F)
#mb.list<-fread("/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/2nd.list",data.table = F,header=F) 
head(gwas.list)
head(mb.list)
nrow(mb.list)
nrow(gwas.list)
gwas.list$V2<-"1"
#gwas.list<-gwas.list[c(227:nrow(gwas.list)),]


both.lists<-left_join(mb.list,gwas.list,by="V1")
nrow(both.lists)
head(both.lists)
which(is.na(both.lists$V2))
tableNA<-both.lists[is.na(both.lists$V2),"V1" ]
nrow(tableNA)
tail(tableNA)
left_out<-data.frame(cbind(which(is.na(both.lists$V2))+226,tableNA))
head(left_out)
class(left_out)
nrow(left_out)

both.lists[261,]
colnames(left_out)<-c("row.in.list","taxaname")
left_out$possible_reason<-"NA"
left_out$possible_reason<-cbind(c(rep("no_convergence",19),rep("cluster_error",nrow(left_out)-19)))
head(left_out,20)

write.table(left_out,"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/left_out_v2.list",quote = F,row.names = F,col.names = T)


#243-226
#252-226
#both.lists[c(17:26),]
#write.table(both.lists[c(17:26),"V1"],"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/3rd.list",quote = F,row.names = F,col.names = F)


#mb #96 = k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri

gwas.list[c(95),]
head(gwas.list)
gwas.list3<-fread("/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa/gwas.list3",data.table = F,header=F)
gwas.list3$V2<-"1"
nrow(gwas.list3)
mb.toremove<-gwas.list$V1[c(1:96)]
toclump<-gwas.list3[-which(gwas.list3$V1 %in% mb.toremove),]
write.table(gwas.list3[-which(gwas.list3$V1 %in% mb.toremove),],"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa/to.clump.list",quote = F,row.names = F,col.names = F)
length(toclump)

toclump[toclump %like% "k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri"]
