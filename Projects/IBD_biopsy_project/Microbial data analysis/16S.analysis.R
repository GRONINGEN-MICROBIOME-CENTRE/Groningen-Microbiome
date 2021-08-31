# ===============================================================================================
# this script is to analyze 16S from dada2 output
# ===============================================================================================
library(ggplot2)
source("Microbiome.function.R")
library(vegan)
library(ggsci)
library(stringr)
library(data.table)
library(crayon)
library(randomcoloR)
library(cosinor2)
library(reshape2)
library(umapr)
library(tidyverse)
library(uwot)
library(MMUPHin)
library(magrittr)
library(dplyr)
library("viridis") 
# ===============================================================================================
# dada2 statistic report
# ===============================================================================================

dada2_run1=read.table("DadaOutput/Run_1///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run1$SampleID=paste(rownames(dada2_run1))
dada2_run1$batch="run_1"

dada2_run2=read.table("DadaOutput/Run_2///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run2$SampleID=paste(rownames(dada2_run2))
dada2_run2$batch="run_2"

dada2_run3=read.table("DadaOutput/Run_3///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run3$SampleID=paste(rownames(dada2_run3))
dada2_run3$batch="run_3"

dada2_run4=read.table("DadaOutput/Run_4//dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run4$SampleID=paste(rownames(dada2_run4))
dada2_run4$batch="run_4"

dada2_run5=read.table("DadaOutput/Run_5///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run5$SampleID=paste(rownames(dada2_run5))
dada2_run5$batch="run_5"

dada2_run6=read.table("DadaOutput/Run_6///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run6$SampleID=paste(rownames(dada2_run6))
dada2_run6$batch="run_6"

dada2_run7=read.table("DadaOutput/Run_7//dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run7$SampleID=paste(rownames(dada2_run7))
dada2_run7$batch="run_7"

dada2_run8=read.table("DadaOutput/Run_8///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run8$SampleID=paste(rownames(dada2_run8))
dada2_run8$batch="run_8"

dada2_run9=read.table("DadaOutput/Run_9///dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)
dada2_run9$SampleID=paste(rownames(dada2_run9))
dada2_run9$batch="run_9"

dada2_runall=rbind(dada2_run1,dada2_run2,dada2_run3,dada2_run4,dada2_run5,dada2_run6,dada2_run7,dada2_run8,dada2_run9)
dada2_runall=dada2_runall[,c("SampleID","batch","nonchim","input")]

# ===============================================================================================
# aglined to kingdom reads count
# ===============================================================================================

dada2_align_run1=read.table("DadaOutput/Run_1///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run1=as.data.frame(t(dada2_align_run1))
dada2_align_run1=dada2_align_run1[,"Bacteria",drop=F]
dada2_align_run1$SampleID=paste(rownames(dada2_align_run1),"run1",sep = "_")
dada2_align_run1$batch="run_1"

dada2_align_run2=read.table("DadaOutput/Run_2///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run2=as.data.frame(t(dada2_align_run2))
dada2_align_run2=dada2_align_run2[,"Bacteria",drop=F]
dada2_align_run2$SampleID=paste(rownames(dada2_align_run2),"run2",sep = "_")
dada2_align_run2$batch="run_2"

dada2_align_run3=read.table("DadaOutput/Run_3///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run3=as.data.frame(t(dada2_align_run3))
dada2_align_run3=dada2_align_run3[,"Bacteria",drop=F]
dada2_align_run3$SampleID=paste(rownames(dada2_align_run3),"run3",sep = "_")
dada2_align_run3$batch="run_3"

dada2_align_run4=read.table("DadaOutput/Run_4///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run4=as.data.frame(t(dada2_align_run4))
dada2_align_run4=dada2_align_run4[,"Bacteria",drop=F]
dada2_align_run4$SampleID=paste(rownames(dada2_align_run4),"run4",sep = "_")
dada2_align_run4$batch="run_4"

dada2_align_run5=read.table("DadaOutput/Run_5///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run5=as.data.frame(t(dada2_align_run5))
dada2_align_run5=dada2_align_run5[,"Bacteria",drop=F]
dada2_align_run5$SampleID=paste(rownames(dada2_align_run5),"run5",sep = "_")
dada2_align_run5$batch="run_5"

dada2_align_run6=read.table("DadaOutput/Run_6///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run6=as.data.frame(t(dada2_align_run6))
dada2_align_run6=dada2_align_run6[,"Bacteria",drop=F]
dada2_align_run6$SampleID=paste(rownames(dada2_align_run6),"run6",sep = "_")
dada2_align_run6$batch="run_6"

dada2_align_run7=read.table("DadaOutput/Run_7///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run7=as.data.frame(t(dada2_align_run7))
dada2_align_run7=dada2_align_run7[,"Bacteria",drop=F]
dada2_align_run7$SampleID=paste(rownames(dada2_align_run7),"run7",sep = "_")
dada2_align_run7$batch="run_7"

dada2_align_run8=read.table("DadaOutput/Run_8///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run8=as.data.frame(t(dada2_align_run8))
dada2_align_run8=dada2_align_run8[,"Bacteria",drop=F]
dada2_align_run8$SampleID=paste(rownames(dada2_align_run8),"run8",sep = "_")
dada2_align_run8$batch="run_8"

dada2_align_run9=read.table("DadaOutput/Run_9///Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run9=as.data.frame(t(dada2_align_run9))
dada2_align_run9=dada2_align_run9[,"Bacteria",drop=F]
dada2_align_run9$SampleID=paste(rownames(dada2_align_run9),"run9",sep = "_")
dada2_align_run9$batch="run_9"

dada2_align_all=rbind(dada2_align_run1,dada2_align_run2,dada2_align_run3,dada2_align_run4,dada2_align_run5,dada2_align_run6,dada2_align_run7,dada2_align_run8,dada2_align_run9)
dada2_align_all$AlignedReads=dada2_align_all$Bacteria
dada2_align_all=dada2_align_all[order(dada2_align_all$SampleID),]
dada2_runall=dada2_runall[order(dada2_runall$SampleID),]
dada2_align_all$Unaligned=dada2_runall$nonchim-dada2_align_all$AlignedReads
dada2_align_all$ratio=dada2_align_all$AlignedReads/(dada2_align_all$AlignedReads+dada2_align_all$Unaligned)
dada2_align_all$ratio[is.na(dada2_align_all$ratio)]=0

# ===============================================================================================
# combine with sequencing metadata
# ===============================================================================================

DNAqc=read.table("DadaOutput/MiSeq datasheet MDL run  merged.txt",sep = "\t",header = T,stringsAsFactors = F)
dada2_align_all=merge(dada2_align_all,DNAqc,by="SampleID",all=F)

# ===============================================================================================
# reads count plot
# ===============================================================================================

dada2_align_all_a=dada2_align_all[,c("SampleID","Qbit","AlignedReads","batch","ratio")]
dada2_align_all_b=dada2_align_all[,c("SampleID","Qbit","Unaligned","batch","ratio")]
colnames(dada2_align_all_a)[3]="Reads"
colnames(dada2_align_all_b)[3]="Reads"
dada2_align_all_a$Group="Aligned Reads"
dada2_align_all_b$Group="Unaligned Reads"

dada2_align_all=rbind(dada2_align_all_a,dada2_align_all_b)
dada2_align_all=dada2_align_all[order(dada2_align_all$ratio),]
sample_order=as.factor(dada2_align_all$SampleID[dada2_align_all$Group=="Aligned Reads"])
dada2_align_all$SampleID=factor(dada2_align_all$SampleID,levels = sample_order)

ggplot(data=dada2_align_all) +
  geom_bar(mapping=aes(x=SampleID, y=Reads, fill=Group),stat="identity")+
  geom_line(mapping=aes(x=SampleID, y=ratio*100000),size=0.5,group = 1)+
  facet_wrap(~batch, ncol =3, scales = "free_x", drop = TRUE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 3))+
  scale_y_continuous(sec.axis = sec_axis(trans= ~./100000)) +
  scale_fill_npg()+guides(fill=F)
ggsave("Plot/dada2.reads.count.pdf",width = 12,height = 6)
write.table(dada2_align_all,file = "OutputTable/AlignedReads.QC.txt",sep = "\t",row.names = F,quote = F)

# ===============================================================================================
# before decontamination, genus
# ===============================================================================================

dada2_genus_run1=read.table("DadaOutput/Run_1/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run2=read.table("DadaOutput/Run_2/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run3=read.table("DadaOutput/Run_3/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run4=read.table("DadaOutput/Run_4/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run5=read.table("DadaOutput/Run_5/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run6=read.table("DadaOutput/Run_6/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run7=read.table("DadaOutput/Run_7/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run8=read.table("DadaOutput/Run_8/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)
dada2_genus_run9=read.table("DadaOutput/Run_9/Genus.level.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1)

colnames(dada2_genus_run1)=paste(colnames(dada2_genus_run1),"run1",sep = "_")
colnames(dada2_genus_run2)=paste(colnames(dada2_genus_run2),"run2",sep = "_")
colnames(dada2_genus_run3)=paste(colnames(dada2_genus_run3),"run3",sep = "_")
colnames(dada2_genus_run4)=paste(colnames(dada2_genus_run4),"run4",sep = "_")
colnames(dada2_genus_run5)=paste(colnames(dada2_genus_run5),"run5",sep = "_")
colnames(dada2_genus_run6)=paste(colnames(dada2_genus_run6),"run6",sep = "_")
colnames(dada2_genus_run7)=paste(colnames(dada2_genus_run7),"run7",sep = "_")
colnames(dada2_genus_run8)=paste(colnames(dada2_genus_run8),"run8",sep = "_")
colnames(dada2_genus_run9)=paste(colnames(dada2_genus_run9),"run9",sep = "_")
dada2_genus_run9=dada2_genus_run9[,colnames(dada2_genus_run9) %in% DNAqc$SampleID]

l=list(dada2_genus_run1,dada2_genus_run2,dada2_genus_run3,dada2_genus_run4,dada2_genus_run5,dada2_genus_run6,dada2_genus_run7,dada2_genus_run8,dada2_genus_run9)
dada2_genus_all=Reduce(merge, lapply(l, function(x) data.frame(x, rn = row.names(x))))
rownames(dada2_genus_all)=dada2_genus_all$rn
dada2_genus_all$rn=NULL
dada2_genus_all=as.data.frame(t(dada2_genus_all))

# ===============================================================================================
# before decontamination, the major genus composition, low-quality sample decision
# ===============================================================================================

dada2_genus_run1=as.data.frame(t(dada2_genus_run1))
dada2_genus_run2=as.data.frame(t(dada2_genus_run2))
dada2_genus_run3=as.data.frame(t(dada2_genus_run3))
dada2_genus_run4=as.data.frame(t(dada2_genus_run4))
dada2_genus_run5=as.data.frame(t(dada2_genus_run5))
dada2_genus_run6=as.data.frame(t(dada2_genus_run6))
dada2_genus_run7=as.data.frame(t(dada2_genus_run7))
dada2_genus_run8=as.data.frame(t(dada2_genus_run8))
dada2_genus_run9=as.data.frame(t(dada2_genus_run9))
dada2_genus_run9=dada2_genus_run9[rowSums(dada2_genus_run9)>0,]

dada2_genus_run1=apply(dada2_genus_run1,1,function(x){
  x=x/sum(x)
})
dada2_genus_run1=as.data.frame(t(dada2_genus_run1))
dada2_genus_run1=Filter(dada2_genus_run1,relative = 0.01,present = 0.1)
dada2_genus_run1=dada2_genus_run1[,colSums(dada2_genus_run1)>0]

dada2_genus_run2=apply(dada2_genus_run2,1,function(x){
  x=x/sum(x)
})
dada2_genus_run2=as.data.frame(t(dada2_genus_run2))
dada2_genus_run2=Filter(dada2_genus_run2,relative = 0.01,present = 0.1)
dada2_genus_run2=dada2_genus_run2[,colSums(dada2_genus_run2)>0]

dada2_genus_run3=apply(dada2_genus_run3,1,function(x){
  x=x/sum(x)
})
dada2_genus_run3=as.data.frame(t(dada2_genus_run3))
dada2_genus_run3=Filter(dada2_genus_run3,relative = 0.01,present = 0.1)
dada2_genus_run3=dada2_genus_run3[,colSums(dada2_genus_run3)>0]

dada2_genus_run4=apply(dada2_genus_run4,1,function(x){
  x=x/sum(x)
})
dada2_genus_run4=as.data.frame(t(dada2_genus_run4))
dada2_genus_run4=Filter(dada2_genus_run4,relative = 0.01,present = 0.1)
dada2_genus_run4=dada2_genus_run4[,colSums(dada2_genus_run4)>0]

dada2_genus_run5=apply(dada2_genus_run5,1,function(x){
  x=x/sum(x)
})
dada2_genus_run5=as.data.frame(t(dada2_genus_run5))
dada2_genus_run5=Filter(dada2_genus_run5,relative = 0.01,present = 0.1)
dada2_genus_run5=dada2_genus_run5[,colSums(dada2_genus_run5)>0]

dada2_genus_run6=apply(dada2_genus_run6,1,function(x){
  x=x/sum(x)
})
dada2_genus_run6=as.data.frame(t(dada2_genus_run6))
dada2_genus_run6=Filter(dada2_genus_run6,relative = 0.01,present = 0.1)
dada2_genus_run6=dada2_genus_run6[,colSums(dada2_genus_run6)>0]
dada2_genus_run6=dada2_genus_run6[rowSums(dada2_genus_run6)>0,]

dada2_genus_run7=apply(dada2_genus_run7,1,function(x){
  x=x/sum(x)
})
dada2_genus_run7=as.data.frame(t(dada2_genus_run7))
dada2_genus_run7=Filter(dada2_genus_run7,relative = 0.01,present = 0.1)
dada2_genus_run7=dada2_genus_run7[,colSums(dada2_genus_run7)>0]
dada2_genus_run7=dada2_genus_run7[rowSums(dada2_genus_run7)>0,]

dada2_genus_run8=apply(dada2_genus_run8,1,function(x){
  x=x/sum(x)
})
dada2_genus_run8=as.data.frame(t(dada2_genus_run8))
dada2_genus_run8=Filter(dada2_genus_run8,relative = 0.01,present = 0.1)
dada2_genus_run8=dada2_genus_run8[,colSums(dada2_genus_run8)>0]

dada2_genus_run9=apply(dada2_genus_run9,1,function(x){
  x=x/sum(x)
})
dada2_genus_run9=as.data.frame(t(dada2_genus_run9))
dada2_genus_run9=Filter(dada2_genus_run9,relative = 0.01,present = 0.1)
dada2_genus_run9=dada2_genus_run9[,colSums(dada2_genus_run9)>0]
dada2_genus_run9=dada2_genus_run9[rowSums(dada2_genus_run9)>0,]

run1_genus_count=data.frame(SampleID=rownames(dada2_genus_run1),core_genus=rowSums(dada2_genus_run1>0))
run2_genus_count=data.frame(SampleID=rownames(dada2_genus_run2),core_genus=rowSums(dada2_genus_run2>0))
run3_genus_count=data.frame(SampleID=rownames(dada2_genus_run3),core_genus=rowSums(dada2_genus_run3>0))
run4_genus_count=data.frame(SampleID=rownames(dada2_genus_run4),core_genus=rowSums(dada2_genus_run4>0))
run5_genus_count=data.frame(SampleID=rownames(dada2_genus_run5),core_genus=rowSums(dada2_genus_run5>0))
run6_genus_count=data.frame(SampleID=rownames(dada2_genus_run6),core_genus=rowSums(dada2_genus_run6>0))
run7_genus_count=data.frame(SampleID=rownames(dada2_genus_run7),core_genus=rowSums(dada2_genus_run7>0))
run8_genus_count=data.frame(SampleID=rownames(dada2_genus_run8),core_genus=rowSums(dada2_genus_run8>0))
run9_genus_count=data.frame(SampleID=rownames(dada2_genus_run9),core_genus=rowSums(dada2_genus_run9>0))

core_genus=rbind(run1_genus_count,run2_genus_count,run3_genus_count,run4_genus_count,run5_genus_count,run6_genus_count,run7_genus_count,run8_genus_count,run9_genus_count)
core_genus=merge(core_genus,dada2_align_all,by="SampleID",all=F)
ggplot(core_genus, aes(x=log2(Reads), y=core_genus)) +
  geom_point(size=2, shape=23)
core_genus$bins=findInterval(core_genus$Reads, seq(0, 100000, by=500))
core_genus=core_genus[core_genus$Group=="Aligned Reads",]
ggplot(core_genus, aes(x=bins, y=core_genus, group=bins,color=bins)) + 
  geom_violin(width=3,trim=FALSE)+
  geom_boxplot(width=0.5, fill="white")+scale_fill_npg()+theme_bw()+
  xlab("")+xlim(0,100)
ggsave("Plot/InputReads.vs.core.genus.pdf",width = 15,height = 3)

# ===============================================================================================
# before decontamination, remove samples, by aligned reads <2,000 (see plot above)
# ===============================================================================================

# by mapping rate
poor_samples=dada2_align_all[dada2_align_all$Group=="Aligned Reads",]
poor_samples$SampleID=as.character(poor_samples$SampleID)
poor_samples=unique(as.character(poor_samples$SampleID[poor_samples$Reads<2000]))

# by Qbit
poor_samples=append(poor_samples,as.character(dada2_align_all$SampleID[dada2_align_all$Qbit<1]))
poor_samples=unique(poor_samples)

# check if NC was removed by the cretiria above
DNAqc$PoorSample=NA
DNAqc$PoorSample[DNAqc$SampleID %in% poor_samples]="PoorSample"
DNAqc$PoorSample[!DNAqc$SampleID %in% poor_samples]="GoodSample"

# ===============================================================================================
# decontamination, combine all negative controls, remove ASVs <10% each batch
# ===============================================================================================

library(phyloseq)
library(decontam); packageVersion("decontam")
library(ggplot2)
library(DECIPHER)
library(dada2)
summarize_taxa <- function(counts, taxonomy) {
  require('plyr')
  if(is.matrix(taxonomy)) {
    alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    tapply(counts, taxonomy, sum)
  }
}

taxa1=read.table("DadaOutput/Run_1//SV.taxa.txt",header = T,sep = "\t")
taxa1=as.matrix(taxa1)
seqtab.nochim1=read.table("./DadaOutput/Run_1/SV.table.txt",sep = "\t",header = T)
seqtab.nochim1=seqtab.nochim1[,colSums(seqtab.nochim1>2)>nrow(seqtab.nochim1)*0.1]
seqtab.nochim1=as.matrix(seqtab.nochim1)
rownames(seqtab.nochim1)=paste(rownames(seqtab.nochim1),"run1",sep = "_")
samdf1=DNAqc[DNAqc$SampleID %like% "run1",]
rownames(samdf1)=samdf1$SampleID
ps1 = phyloseq(otu_table(seqtab.nochim1, taxa_are_rows=FALSE),sample_data(samdf1),tax_table(taxa1))

taxa2=read.table("DadaOutput/Run_2/SV.taxa.txt",header = T,sep = "\t")
taxa2=as.matrix(taxa2)
seqtab.nochim2=read.table("./DadaOutput/Run_2/SV.table.txt",sep = "\t",header = T)
seqtab.nochim2=seqtab.nochim2[,colSums(seqtab.nochim2>2)>nrow(seqtab.nochim2)*0.1]
seqtab.nochim2=as.matrix(seqtab.nochim2)
rownames(seqtab.nochim2)=paste(rownames(seqtab.nochim2),"run2",sep = "_")
samdf2=DNAqc[DNAqc$SampleID %like% "run2",]
rownames(samdf2)=samdf2$SampleID
ps2 = phyloseq(otu_table(seqtab.nochim2, taxa_are_rows=FALSE),sample_data(samdf2),tax_table(taxa2))

taxa3=read.table("DadaOutput/Run_3/SV.taxa.txt",header = T,sep = "\t")
taxa3=as.matrix(taxa3)
seqtab.nochim3=read.table("./DadaOutput/Run_3/SV.table.txt",sep = "\t",header = T)
seqtab.nochim3=seqtab.nochim3[,colSums(seqtab.nochim3>2)>nrow(seqtab.nochim3)*0.1]
seqtab.nochim3=as.matrix(seqtab.nochim3)
rownames(seqtab.nochim3)=paste(rownames(seqtab.nochim3),"run3",sep = "_")
samdf3=DNAqc[DNAqc$SampleID %like% "run3",]
rownames(samdf3)=samdf3$SampleID
ps3 = phyloseq(otu_table(seqtab.nochim3, taxa_are_rows=FALSE),sample_data(samdf3),tax_table(taxa3))

taxa4=read.table("DadaOutput/Run_4/SV.taxa.txt",header = T,sep = "\t")
taxa4=as.matrix(taxa4)
seqtab.nochim4=read.table("./DadaOutput/Run_4/SV.table.txt",sep = "\t",header = T)
seqtab.nochim4=seqtab.nochim4[,colSums(seqtab.nochim4>2)>nrow(seqtab.nochim4)*0.1]
seqtab.nochim4=as.matrix(seqtab.nochim4)
rownames(seqtab.nochim4)=paste(rownames(seqtab.nochim4),"run4",sep = "_")
samdf4=DNAqc[DNAqc$SampleID %like% "run4",]
rownames(samdf4)=samdf4$SampleID
ps4 = phyloseq(otu_table(seqtab.nochim4, taxa_are_rows=FALSE),sample_data(samdf4),tax_table(taxa4))

taxa5=read.table("DadaOutput/Run_5/SV.taxa.txt",header = T,sep = "\t")
taxa5=as.matrix(taxa5)
seqtab.nochim5=read.table("./DadaOutput/Run_5/SV.table.txt",sep = "\t",header = T)
seqtab.nochim5=seqtab.nochim5[,colSums(seqtab.nochim5>2)>nrow(seqtab.nochim5)*0.1]
seqtab.nochim5=as.matrix(seqtab.nochim5)
rownames(seqtab.nochim5)=paste(rownames(seqtab.nochim5),"run5",sep = "_")
samdf5=DNAqc[DNAqc$SampleID %like% "run5",]
rownames(samdf5)=samdf5$SampleID
ps5 = phyloseq(otu_table(seqtab.nochim5, taxa_are_rows=FALSE),sample_data(samdf5),tax_table(taxa5))

taxa6=read.table("DadaOutput/Run_6/SV.taxa.txt",header = T,sep = "\t")
taxa6=as.matrix(taxa6)
seqtab.nochim6=read.table("./DadaOutput/Run_6/SV.table.txt",sep = "\t",header = T)
seqtab.nochim6=seqtab.nochim6[,colSums(seqtab.nochim6>2)>nrow(seqtab.nochim6)*0.1]
seqtab.nochim6=as.matrix(seqtab.nochim6)
rownames(seqtab.nochim6)=paste(rownames(seqtab.nochim6),"run6",sep = "_")
samdf6=DNAqc[DNAqc$SampleID %like% "run6",]
rownames(samdf6)=samdf6$SampleID
ps6 = phyloseq(otu_table(seqtab.nochim6, taxa_are_rows=FALSE),sample_data(samdf6),tax_table(taxa6))

taxa7=read.table("DadaOutput/Run_7/SV.taxa.txt",header = T,sep = "\t")
taxa7=as.matrix(taxa7)
seqtab.nochim7=read.table("./DadaOutput/Run_7/SV.table.txt",sep = "\t",header = T)
seqtab.nochim7=seqtab.nochim7[,colSums(seqtab.nochim7>2)>nrow(seqtab.nochim7)*0.1]
seqtab.nochim7=as.matrix(seqtab.nochim7)
rownames(seqtab.nochim7)=paste(rownames(seqtab.nochim7),"run7",sep = "_")
samdf7=DNAqc[DNAqc$SampleID %like% "run7",]
rownames(samdf7)=samdf7$SampleID
ps7 = phyloseq(otu_table(seqtab.nochim7, taxa_are_rows=FALSE),sample_data(samdf7),tax_table(taxa7))

taxa8=read.table("DadaOutput/Run_8/SV.taxa.txt",header = T,sep = "\t")
taxa8=as.matrix(taxa8)
seqtab.nochim8=read.table("./DadaOutput/Run_8/SV.table.txt",sep = "\t",header = T)
seqtab.nochim8=seqtab.nochim8[,colSums(seqtab.nochim8>2)>nrow(seqtab.nochim8)*0.1]
seqtab.nochim8=as.matrix(seqtab.nochim8)
rownames(seqtab.nochim8)=paste(rownames(seqtab.nochim8),"run8",sep = "_")
samdf8=DNAqc[DNAqc$SampleID %like% "run8",]
rownames(samdf8)=samdf8$SampleID
ps8 = phyloseq(otu_table(seqtab.nochim8, taxa_are_rows=FALSE),sample_data(samdf8),tax_table(taxa8))

taxa9=read.table("DadaOutput/Run_9/SV.taxa.txt",header = T,sep = "\t")
taxa9=as.matrix(taxa9)
seqtab.nochim9=read.table("./DadaOutput/Run_9/SV.table.txt",sep = "\t",header = T)
seqtab.nochim9=seqtab.nochim9[,colSums(seqtab.nochim1>2)>nrow(seqtab.nochim9)*0.1]
seqtab.nochim9=as.matrix(seqtab.nochim9)
rownames(seqtab.nochim9)=paste(rownames(seqtab.nochim9),"run9",sep = "_")
samdf9=DNAqc[DNAqc$SampleID %like% "run9",]
rownames(samdf9)=samdf9$SampleID
ps9 = phyloseq(otu_table(seqtab.nochim9, taxa_are_rows=FALSE),sample_data(samdf9),tax_table(taxa9))

ps_merge=merge_phyloseq(ps1,ps2,ps3,ps4,ps5,ps6,ps7,ps8,ps9)
aa=as.data.frame(otu_table(ps_merge))
write.table(aa,"FilterASV///ASV.table.txt",quote = F,sep = "\t")
svtable=as.matrix(aa)
rownames(svtable)=rownames(aa)
svtable=t(svtable)
txtable=as.matrix(as.data.frame((tax_table(ps_merge))))
mm=summarize_taxa(svtable, txtable)
write.table(as.data.frame(mm$Kingdom),file = "FilterASV//Kingdom.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Phylum),file = "FilterASV//Phylum.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Class),file = "FilterASV//Class.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Order),file = "FilterASV//Order.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Family),file = "FilterASV//Family.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Genus),file = "FilterASV//Genus.txt",sep = "\t",quote = F)

# identify contamination by prevelence
ps_merge <- subset_samples(ps_merge,!sample_names(ps_merge) %in% c("S67_run3","S1_run4"))
ps_merge <- subset_samples(ps_merge,sample_names(ps_merge) %in% DNAqc$SampleID[DNAqc$PoorSample=="GoodSample"])
sample_data(ps_merge)$is.neg <- sample_data(ps_merge)$Sample_or_Control == "Control"

contamdf.prev05 <- isContaminant(ps_merge, method="prevalence", neg="is.neg", threshold=0.5)
table(contamdf.prev05$contaminant)
ps_merge.pa <- transform_sample_counts(ps_merge, function(abund) 1*(abund>0))
ps_merge.pa.neg <- prune_samples(sample_data(ps_merge.pa)$Sample_or_Control == "Control", ps_merge.pa)
ps_merge.pa.pos <- prune_samples(sample_data(ps_merge.pa)$Sample_or_Control == "TrueSample", ps_merge.pa)
df.pa <- data.frame(pa.pos=taxa_sums(ps_merge.pa.pos), pa.neg=taxa_sums(ps_merge.pa.neg),
                    contaminant=contamdf.prev05$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")
ggsave("Plot/decontamination.presence.pdf")
hist(contamdf.prev05$p,breaks = 200)

# compare contamination taxa
ps_merge.contam <- prune_taxa(contamdf.prev05$contaminant, ps_merge)
aa=as.data.frame(otu_table(ps_merge.contam))
svtable=as.matrix(aa)
rownames(svtable)=rownames(aa)
svtable=t(svtable)
txtable=as.matrix(as.data.frame((tax_table(ps_merge.contam))))
mm=summarize_taxa(svtable, txtable)
contamination=as.data.frame(mm$Genus)
rownames(contamination)
contamination=as.data.frame(mm$Family)

# remove contamination ASV
ps.noncontam <- prune_taxa(!contamdf.prev05$contaminant, ps_merge)
aa=as.data.frame(otu_table(ps.noncontam))
write.table(aa,"Decontamination_dada2//ASV.2000reads.table.txt",quote = F,sep = "\t")
svtable=as.matrix(aa)
rownames(svtable)=rownames(aa)
svtable=t(svtable)
txtable=as.matrix(as.data.frame((tax_table(ps.noncontam))))
mm=summarize_taxa(svtable, txtable)

# remove known contaminates
write.table(as.data.frame(mm$Kingdom),file = "Decontamination_dada2//Kingdom.2000reads.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Phylum),file = "Decontamination_dada2//Phylum.2000reads.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Class),file = "Decontamination_dada2//Class.2000reads.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Order),file = "Decontamination_dada2//Order.2000reads.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Family),file = "Decontamination_dada2//Family.2000reads.txt",sep = "\t",quote = F)
write.table(as.data.frame(mm$Genus),file = "Decontamination_dada2//Genus.2000reads.txt",sep = "\t",quote = F)

# ===============================================================================================
# post decontamination, import phenotype
# ===============================================================================================

phenotype=read.table("Metadaata.16sRNAseq.update.mostcomplete.txt",sep = "\t",fill = T,stringsAsFactors = F,header = T,check.names = F)
phenotype$Location=NA
phenotype$Location[phenotype$`Location (corrected)`=="ascendens"]="colon"
phenotype$Location[phenotype$`Location (corrected)`=="coecum"]="colon"
phenotype$Location[phenotype$`Location (corrected)`=="descendens"]="colon"
phenotype$Location[phenotype$`Location (corrected)`=="sigmoid"]="colon"
phenotype$Location[phenotype$`Location (corrected)`=="transversum"]="colon"
phenotype$Location[phenotype$`Location (corrected)`=="ileum"]="ileum"
phenotype$Location[phenotype$`Location (corrected)`=="rectum"]="colon"

phenotype=phenotype[,colnames(phenotype) %in% c("Research ID","16SSampleID","biopsy_number",
                                                "Cohort","Location (corrected)","Location",
                                                "Inflammation","age_at_biopsy","sex",
                                                "BMI","smoking_DB","resec_part_colon","resec_ileocec","resec_part_small","biopsy_date")]

# ===============================================================================================
# post decontamination, beta and alpha diversity, genus level
# ===============================================================================================

dada2_genus_decon=read.table("Decontamination_dada2//Genus.2000reads.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
dada2_genus_decon=as.data.frame(t(dada2_genus_decon))
dada2_genus_decon=dada2_genus_decon[rowSums(dada2_genus_decon)>0,]
# shannon
shannon=as.data.frame(diversity((dada2_genus_decon), index="shannon"))
colnames(shannon)="shannon"

shannon=merge(shannon,DNAqc,by.x="row.names",by.y="SampleID",all = F)
shannon=merge(shannon,phenotype,by.x="Row.names",by.y="16SSampleID",all=F)
shannon$Batch=NA
shannon$Batch[shannon$Row.names %like% "run1"]="Run1"
shannon$Batch[shannon$Row.names %like% "run2"]="Run2"
shannon$Batch[shannon$Row.names %like% "run3"]="Run3"
shannon$Batch[shannon$Row.names %like% "run4"]="Run4"
shannon$Batch[shannon$Row.names %like% "run5"]="Run5"
shannon$Batch[shannon$Row.names %like% "run6"]="Run6"
shannon$Batch[shannon$Row.names %like% "run7"]="Run7"
shannon$Batch[shannon$Row.names %like% "run8"]="Run8"
shannon$Batch[shannon$Row.names %like% "run9"]="Run9"

# batch and storage time effect
ggplot(shannon, aes(x=Batch, y=shannon, fill=Batch)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+scale_fill_npg()+theme_bw()+
  xlab("")

shannon=shannon[!is.na(shannon$biopsy_date),]
shannon$biopsy_date=str_split_fixed(shannon$biopsy_date, "/", 2)[,1]
shannon$biopsy_date=as.numeric(shannon$biopsy_date)
shannon_tmp=shannon[shannon$biopsy_date>2000 & shannon$biopsy_date<2020,]
shannon_tmp=shannon_tmp[!is.na(shannon_tmp$biopsy_date),]
ggplot(shannon_tmp, aes(x=biopsy_date, y=shannon)) + 
  geom_point()+theme_bw()+geom_smooth(method='lm')+
  xlab("")
cor.test(shannon_tmp$biopsy_date,shannon_tmp$shannon,method = "spearman")

# bray-curtis
beta_diversity=vegdist((dada2_genus_decon),method = "bray")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates
pcoas=merge(pcoas,DNAqc,by.x="Sample",by.y="SampleID",all = F)
pcoas=merge(pcoas,phenotype,by.x="Sample",by.y="16SSampleID",all=F)
pcoas=merge(pcoas,dada2_align_all_a,by.x="Sample","SampleID",all=F)

pcoas$Batch=NA
pcoas$Batch[pcoas$Sample %like% "run1"]="Run1"
pcoas$Batch[pcoas$Sample %like% "run2"]="Run2"
pcoas$Batch[pcoas$Sample %like% "run3"]="Run3"
pcoas$Batch[pcoas$Sample %like% "run4"]="Run4"
pcoas$Batch[pcoas$Sample %like% "run5"]="Run5"
pcoas$Batch[pcoas$Sample %like% "run6"]="Run6"
pcoas$Batch[pcoas$Sample %like% "run7"]="Run7"
pcoas$Batch[pcoas$Sample %like% "run8"]="Run8"
pcoas$Batch[pcoas$Sample %like% "run9"]="Run9"

# cohort effect
ggplot (pcoas, aes(-PCoA1,PCoA2,color=Cohort,fill=Cohort)) + 
  geom_point() + theme_bw() +stat_ellipse(geom = "polygon", alpha = 0.1)+
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + 
  ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) +
  guides(size=F)+
  scale_color_npg()

distmat <- as.matrix(beta_diversity)
cohort1.dis=distmat[rownames(distmat) %in% pcoas$Sample[pcoas$Cohort=="DDTx"],colnames(distmat) %in% pcoas$Sample[pcoas$Cohort=="DDTx"]]
cohort2.dis=distmat[rownames(distmat) %in% pcoas$Sample[pcoas$Cohort=="1000IBD selection"],colnames(distmat) %in% pcoas$Sample[pcoas$Cohort=="1000IBD selection"]]
cohort3.dis=distmat[rownames(distmat) %in% pcoas$Sample[pcoas$Cohort=="1000IBD selection surgery"],colnames(distmat) %in% pcoas$Sample[pcoas$Cohort=="1000IBD selection surgery"]]
cohort4.dis=distmat[rownames(distmat) %in% pcoas$Sample[pcoas$Cohort=="Control"],colnames(distmat) %in% pcoas$Sample[pcoas$Cohort=="Control"]]

distmat1=data.frame(as.table(cohort1.dis))[lower.tri(cohort1.dis, diag = F), ]
distmat2=data.frame(as.table(cohort2.dis))[lower.tri(cohort2.dis, diag = F), ]
distmat3=data.frame(as.table(cohort3.dis))[lower.tri(cohort3.dis, diag = F), ]
distmat4=data.frame(as.table(cohort4.dis))[lower.tri(cohort4.dis, diag = F), ]
distmat1$group="DDTx"
distmat2$group="1000IBD selection"
distmat3$group="1000IBD selection surgery"
distmat4$group="Control"
distmat=rbind(distmat1,distmat2,distmat3,distmat4)

ggplot(distmat, aes(x=group, y=Freq,fill=group)) + 
  geom_boxplot(color="black")+scale_fill_jama()+theme_bw()+guides(color=F)
ggsave("Plot/dissimilarity.genus.pdf",width = 8,height = 6)
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="1000IBD selection"])
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="1000IBD selection surgery"])
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="DDTx"])

# ===============================================================================================
# post filtering, add PCoA arrows, genus level
# ===============================================================================================
dada2_genus_decon=as.data.frame(t(dada2_genus_decon))
beta_diversity=vegdist((dada2_genus_decon),method = "bray")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates
pcoas=merge(pcoas,DNAqc,by.x="Sample",by.y="SampleID",all = F)
pcoas=merge(pcoas,phenotype,by.x="Sample",by.y="16SSampleID",all=F)
pcoas=merge(pcoas,dada2_align_all_a,by.x="Sample","SampleID",all=F)

## PC1 correlations
PC1_BC_cors = apply(dada2_genus_decon,2,function(x){
  cor(pcoas$PCoA1,x,method = "spearman")})

PC1_BC_pvalues = p.adjust(apply(dada2_genus_decon,2,function(x){
  cor.test(pcoas$PCoA1,x,method = "spearman")$p.value}),method = "fdr")

PC1_pvalues = PC1_BC_pvalues[names(PC1_BC_cors)]
PC1_cors = data.frame("PC1_cors" = PC1_BC_cors,"PC1_pvalues" = PC1_pvalues)


## PC2 correlations
PC2_BC_cors = apply(dada2_genus_decon,2,function(x){
  cor(pcoas$PCoA2,x,method = "spearman")})

PC2_BC_pvalues = p.adjust(apply(dada2_genus_decon,2,function(x){
  cor.test(pcoas$PCoA2,x,method = "spearman")$p.value}),method = "fdr")

PC2_pvalues = PC2_BC_pvalues[names(PC2_BC_cors)]
PC2_cors = data.frame("PC2_cors" = PC2_BC_cors,"PC2_pvalues" = PC2_pvalues)

PC_cors = data.frame(cbind(PC1_cors,PC2_cors))
PC_cors = PC_cors[(PC_cors$PC1_pvalues <0.05 & PC_cors$PC2_pvalues <0.05),]

ggplot(pcoas,aes(x = PCoA1/max(abs(PCoA1)), y = PCoA2/max(abs(PCoA2))))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.08539776, yend = -0.08543441),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.21812126, yend = -0.31309078),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.29063994, yend = -0.85313237),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.98583329, yend = 0.29192570),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.11134703, yend = 0.08391757),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_text(x=0.08539776,y=-0.08543441,label="Akkermansia",color="black",fontface="italic")+
  geom_text(x=0.21812126,y=-0.31309078,label="Bifidobacterium",color="black",fontface="italic")+
  geom_text(x=0.29063994,y=-0.85313237,label="Faecalibacterium",color="black",fontface="italic")+
  geom_text(x=-0.98583329,y=0.29192570,label="Bacteroides",color="black",fontface="italic")+
  geom_text(x=-0.11134703,y=0.08391757,label="Escherichia/Shigella",color="black",fontface="italic")+
  theme_classic()
ggsave("Plot/Arrows.pdf")

# ===============================================================================================
# post filtering, genus level, batch composition
# ===============================================================================================

dada2_genus_decon=read.table("Decontamination_dada2//Genus.2000reads.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
dada2_genus_decon=as.data.frame(t(dada2_genus_decon))
dada2_genus_decon=dada2_genus_decon[rownames(dada2_genus_decon) %in% DNAqc$SampleID[DNAqc$Sample_or_Control=="TrueSample"],]
dada2_genus_decon=dada2_genus_decon[,colSums(dada2_genus_decon)>0]
dada2_genus_decon=dada2_genus_decon[rowSums(dada2_genus_decon)>0,]

dada2_genus_decon=apply(dada2_genus_decon,1,function(x){
  x=x/sum(x)
})
dada2_genus_decon=as.data.frame(t(dada2_genus_decon))
dada2_genus_decon[is.na(dada2_genus_decon)]=0
dada2_genus_decon=Filter(dada2_genus_decon,relative = 0.001,present = 0.1)
dada2_genus_decon=dada2_genus_decon[,colSums(dada2_genus_decon)>0]
dada2_genus_decon=dada2_genus_decon[rowSums(dada2_genus_decon)>0,]

dada2_genus_plot=CompositionTable(dada2_genus_decon,10)
dada2_genus_plot=merge(dada2_genus_plot,DNAqc,by.x="ID",by.y = "SampleID")
dada2_genus_plot$Batch=NA
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run1"]="run1"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run2"]="run2"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run3"]="run3"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run4"]="run4"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run5"]="run5"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run6"]="run6"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run7"]="run7"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run8"]="run8"
dada2_genus_plot$Batch[dada2_genus_plot$ID %like% "run9"]="run9"

dada2_run1_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run1")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run1_decon_aver$AverageAbundance[dada2_run1_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run1"])
}
dada2_run2_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run2")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run2_decon_aver$AverageAbundance[dada2_run2_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run2"])
}
dada2_run3_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run3")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run3_decon_aver$AverageAbundance[dada2_run3_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run3"])
}
dada2_run4_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run4")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run4_decon_aver$AverageAbundance[dada2_run4_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run4"])
}
dada2_run5_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run5")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run5_decon_aver$AverageAbundance[dada2_run5_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run5"])
}
dada2_run6_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run6")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run6_decon_aver$AverageAbundance[dada2_run6_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run6"])
}
dada2_run7_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run7")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run7_decon_aver$AverageAbundance[dada2_run7_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run7"])
}
dada2_run8_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run8")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run8_decon_aver$AverageAbundance[dada2_run8_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run8"])
}
dada2_run9_decon_aver=data.frame(Taxa=unique(as.character(dada2_genus_plot$Level)),AverageAbundance=NA,Batch="run9")
for(i in unique(as.character(dada2_genus_plot$Level))){
  dada2_run9_decon_aver$AverageAbundance[dada2_run9_decon_aver$Taxa==i]=mean(dada2_genus_plot$Relative[dada2_genus_plot$Level==i & dada2_genus_plot$Batch=="run9"])
}
dada2_genus_plot_average=rbind(dada2_run1_decon_aver,
                                dada2_run2_decon_aver,
                                dada2_run3_decon_aver,
                                dada2_run4_decon_aver,
                                dada2_run5_decon_aver,
                                dada2_run6_decon_aver,
                                dada2_run7_decon_aver,
                                dada2_run8_decon_aver,
                                dada2_run9_decon_aver)
ggplot(dada2_genus_plot_average, aes(x=Batch, y=AverageAbundance,fill = Taxa))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_npg()+
  ylab("Relative abundance")+
  theme_bw()+xlab("")+
  theme(plot.subtitle = element_text(vjust = 1), 
        plot.caption = element_text(vjust = 1), 
        axis.line = element_line(colour = 'white'), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line.x = element_line(colour = 'white'),
        axis.title.x=element_blank(),
        axis.text.x = element_text(hjust = 1,size = 15),
        axis.ticks.x=element_blank())+theme(legend.position="bottom")
ggsave("Plot/Batch.composition.genus.decon.pdf",width = 6,height = 4)
