# =========================================================================================
# By: Weersma Group, UMCG (2020)
# PLOTTER FOR DMP microbiome cluster anlaysis
# =========================================================================================

library(ggplot2)
library(ggridges)
library(ggsci)
library(RColorBrewer)
library(randomcoloR)
library(ggpubr)
library(MASS)
library(tidyverse)
library(caret)

CompositionTable <- function(x,n){ # change the metaphlan result to a composition table, select top n most abundant features
  require(foreach)
#  x[is.na(x)]=0
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  
  return(composition_table)
}

# genus level plot =========
# ==================================
# binominal distribution
# ==================================
dag3=read.table("DAG3_metaphlan_bacteria_archaea_filtered.txt",sep = "\t",
                header = T,check.names = F,stringsAsFactors = F)
rownames(dag3)=dag3$ID
dag3$ID=NULL

dag3_genus=dag3[,grep("g__",colnames(dag3))]
dag3_genus=dag3_genus[,grep("s__",colnames(dag3_genus),invert = T)]
colnames(dag3_genus)=lapply(colnames(dag3_genus),function(x){
  strsplit(x,"g__")[[1]][2]
})
dag3_genus_plot=CompositionTable(dag3_genus,10)
dag3_genus_plot$Level=factor(dag3_genus_plot$Level,levels = c("Bifidobacterium","Parabacteroides","Eubacterium",
                                                              "Bacteroidales_noname","Sutterella","Subdoligranulum",
                                                              "Prevotella","Faecalibacterium",
                                                              "Bacteroides","Alistipes"))
dag3_genus_plot[dag3_genus_plot==0]=NA
dag3_genus_plot=na.omit(dag3_genus_plot)
dag3_genus_plot$Relative=-log2(dag3_genus_plot$Relative)
set.seed(10)
n <- 10
palette <- distinctColorPalette(n)
ggplot(dag3_genus_plot, aes(x = Relative, y = Level,fill=Level,color=Level)) + theme_classic()+
  geom_density_ridges(alpha=0.8,scale = 2)+scale_fill_manual(values = palette)+scale_color_manual(values = palette)+ylab("")
ggsave("./Plot/Genus.distribution.pdf")

# ==================================
# technical factors
# ==================================
dag3_meta=read.table("DAG3.metadata.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
dag3_meta=dag3_meta[,c("DAG3_sampleID","postclean.reads","postclean.gc","conc.ng.ul","vol.ul","age",
                       "gender","BMI")]
rownames(dag3_meta)=dag3_meta$DAG3_sampleID
batch=read.table("dag3_batches_samples.txt",sep = "\t",header = F,stringsAsFactors = F)
batch=batch[!duplicated(batch$V2),]
rownames(batch)=batch$V2
colnames(batch)=c("Batch","Name")
dag3_meta=merge(dag3_meta,batch,by="row.names",all=F)
rownames(dag3_meta)=dag3_meta$Row.names
dag3_meta$Row.names=NULL

techniq=merge(dag3_meta,dag3,by="row.names",all=F)
techniq[techniq=="F"]=0
techniq[techniq=="M"]=1
techniq$gender=as.numeric(techniq$gender)

# Prevotella plot
myplots <- list()
count = 1
for(i in 3:9){
  tech=colnames(techniq)[i]
  df=data.frame(x=(techniq[,"Prevotella"]),y=techniq[,tech])
  mm=(cor.test(df$x,df$y,method = "spearman"))
  pvalue=round(mm$p.value,5)
  cor=round(mm$estimate,5)
  myplots[[count]]=ggplot(df, aes(x=x, y=y)) + 
    geom_point(color="#E18727FF",alpha=0.5)+theme_bw()+
    geom_smooth(method=lm,color="#0072B5FF")+xlab("Prevotella")+ylab(tech)+
    labs(title = paste("Pvalue is ",pvalue),subtitle = paste("CorrelationCoefficient is ",cor))+theme(plot.title = element_text(size=10),plot.subtitle = element_text(size = 10))
  count = count + 1
}
ggarrange(plot(myplots[[1]]),plot(myplots[[2]]),plot(myplots[[3]]),plot(myplots[[4]]),ncol = 2,nrow = 2)
ggsave("Plot/Prevotella.tech.pdf")

# Bacteroides plot
myplots <- list()
count = 1
for(i in 3:9){
  tech=colnames(techniq)[i]
  df=data.frame(x=(techniq[,"Bacteroides"]),y=techniq[,tech])
  mm=(cor.test(df$x,df$y,method = "spearman"))
  pvalue=round(mm$p.value,5)
  cor=round(mm$estimate,5)
  myplots[[count]]=ggplot(df, aes(x=x, y=y)) + 
    geom_point(color="#20854EFF",alpha=0.5)+theme_bw()+
    geom_smooth(method=lm,color="#0072B5FF")+xlab("Bacteroides")+ylab(tech)+
    labs(title = paste("Pvalue is ",pvalue),subtitle = paste("CorrelationCoefficient is ",cor))+theme(plot.title = element_text(size=10),plot.subtitle = element_text(size = 10))
  count = count + 1
}
ggarrange(plot(myplots[[1]]),plot(myplots[[2]]),plot(myplots[[3]]),plot(myplots[[4]]),ncol = 2,nrow = 2)
ggsave("Plot/Bacteoides.tech.pdf")

# ==================================
# two enterotypes
# ==================================
pcoa=read.table("PAM.PCoA.genus.txt",stringsAsFactors = F)
enterotype=read.table("PAM.enterotype.genus.txt",stringsAsFactors = F)
genus_enterotype=cbind(pcoa,enterotype)
rownames(genus_enterotype)=rownames(enterotype)
rownames(techniq)=techniq$Row.names
techniq$Row.names=NULL
genus_enterotype=merge(genus_enterotype,techniq,by="row.names",all=F)
genus_enterotype$Enterotype=as.factor(genus_enterotype$Enterotype)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(200))
p0=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =Enterotype),size=1) + theme_bw()+scale_color_npg()+ggtitle("Two Enterotypes")
p01=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Prevotella),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")
p02=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Bacteroides),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")
p03=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Alistipes),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")

p1=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =`postclean.reads`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p2=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =`postclean.gc`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p3=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =`conc.ng.ul`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p4=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =`vol.ul`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p5=ggplot (genus_enterotype, aes(V1,V2)) + geom_point(aes(colour =`Batch`),size=1) + theme_bw()+scale_color_npg()+ggtitle("Total batch")
p6=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch1",], aes(V1,V2)) + geom_point(colour ="#BC3C29FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch1")
p7=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch2",], aes(V1,V2)) + geom_point(colour ="#0072B5FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch2")
p8=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch3",], aes(V1,V2)) + geom_point(colour ="#E18727FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch3")
p9=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch4",], aes(V1,V2)) + geom_point(colour ="#20854EFF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch4")
p10=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch5",], aes(V1,V2)) + geom_point(colour ="#7876B1FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch5")
p11=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch6",], aes(V1,V2)) + geom_point(colour ="#6F99ADFF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch6")
p12=ggplot (genus_enterotype[genus_enterotype$Batch=="dag3_batch7",], aes(V1,V2)) + geom_point(colour ="#FFDC91FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch7")
# plot two enterotypes
ggarrange(p0,p01,p02,p03)
ggsave("Plot/Enterotype.genus.pdf",width = 14,height = 10)
# plot reads, gc content concerntration and volunm
ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)
ggsave("Plot/Enterotype.genus.tech.pdf",width = 10,height = 12)
# plot batch effects
ggarrange(p5,p6,p7,p8,p9,p10,p11,p12,ncol = 3,nrow = 3)
ggsave("Plot/Enterotype.genus.batch.pdf",width = 10,height = 12)


# ==================================
# compare two enterotypes
# ==================================
enterotype$ID=rownames(enterotype)
dag3_genus_plot=CompositionTable(dag3_genus,10)
dag3_genus_plot$Level=factor(dag3_genus_plot$Level,levels = c("Bifidobacterium","Parabacteroides","Eubacterium",
                                                              "Bacteroidales_noname","Sutterella","Subdoligranulum",
                                                              "Prevotella","Faecalibacterium",
                                                              "Bacteroides","Alistipes"))
dag3_genus_plot=merge(dag3_genus_plot,enterotype,by="ID",all=T)
dag3_genus_plot$Enterotype=as.factor(dag3_genus_plot$Enterotype)
ggplot(dag3_genus_plot, aes(Level,Relative,fill=Enterotype)) + 
  geom_boxplot(position=position_dodge(1),outlier.shape = NA)+
  theme_bw()+
  scale_color_npg()+
  scale_fill_npg()+
  guides(color=FALSE)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size = 6,hjust = 1))+
  xlab("Genera")+ylab("Relative Abundance")+coord_flip()
ggsave("Plot/Enterotype.genus.compare.pdf",width = 5,height = 8)
compare=compare_means(Relative~Enterotype,data=dag3_genus_plot,method = "wilcox.test",group.by = "Level")
write.table(compare,"Plot/Enterotype.genus.compare.txt",sep = "\t",quote = F,row.names = F)

dag3_genus_lda=merge(enterotype,dag3_genus,by="row.names",all=F)
rownames(dag3_genus_lda)=dag3_genus_lda$Row.names
dag3_genus_lda$Row.names=NULL
dag3_genus_lda$ID=NULL
dag3_genus_lda=dag3_genus_lda[,c("Enterotype","Bifidobacterium","Parabacteroides","Eubacterium",
                                 "Bacteroidales_noname","Sutterella","Subdoligranulum",
                                 "Prevotella","Faecalibacterium",
                                 "Bacteroides","Alistipes")]
compare_spearman = foreach(i=2:ncol(dag3_genus_lda),.combine = rbind) %do%  {
  mm=cor.test(dag3_genus_lda$Enterotype,dag3_genus_lda[,i],method = "spearman")
  return.string = data.frame(Level=colnames(dag3_genus_lda)[i],Coefficient=mm$estimate,Pvalue=mm$p.value)
}


# species level =============================================================================
# ==================================
# binominal distribution
# ==================================
dag3=read.table("DAG3_metaphlan_bacteria_archaea_filtered.txt",sep = "\t",
                header = T,check.names = F,stringsAsFactors = F)
rownames(dag3)=dag3$ID
dag3$ID=NULL

dag3_species=dag3[,grep("s__",colnames(dag3))]
dag3_species=dag3_species[,grep("t__",colnames(dag3_species),invert = T)]
colnames(dag3_species)=lapply(colnames(dag3_species),function(x){
  strsplit(x,"s__")[[1]][2]
})
dag3_species_plot=CompositionTable(dag3_species,10)
dag3_species_plot$Level=factor(dag3_species_plot$Level,levels = c("Bacteroides_vulgatus","Bacteroidales_bacterium_ph8","Alistipes_shahii",
                                                              "Sutterella_wadsworthensis","Alistipes_onderdonkii","Bacteroides_uniformis",
                                                              "Subdoligranulum_unclassified","Prevotella_copri",
                                                              "Faecalibacterium_prausnitzii","Alistipes_putredinis"))
dag3_species_plot[dag3_species_plot==0]=NA
dag3_species_plot=na.omit(dag3_species_plot)
dag3_species_plot$Relative=-log2(dag3_species_plot$Relative)
set.seed(10)
n <- 10
palette <- distinctColorPalette(n)
ggplot(dag3_species_plot, aes(x = Relative, y = Level,fill=Level,color=Level)) + theme_classic()+
  geom_density_ridges(alpha=0.8,scale = 2)+scale_fill_manual(values = palette)+scale_color_manual(values = palette)+ylab("")
ggsave("./Plot/species.distribution.pdf")

# ==================================
# techniqual factors
# ==================================
dag3_meta=read.table("DAG3.metadata.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
dag3_meta=dag3_meta[,c("DAG3_sampleID","postclean.reads","postclean.gc","conc.ng.ul","vol.ul","age",
                       "gender","BMI")]
rownames(dag3_meta)=dag3_meta$DAG3_sampleID
batch=read.table("dag3_batches_samples.txt",sep = "\t",header = F,stringsAsFactors = F)
batch=batch[!duplicated(batch$V2),]
rownames(batch)=batch$V2
colnames(batch)=c("Batch","Name")
dag3_meta=merge(dag3_meta,batch,by="row.names",all=F)
rownames(dag3_meta)=dag3_meta$Row.names
dag3_meta$Row.names=NULL

techniq=merge(dag3_meta,dag3_species,by="row.names",all=F)
techniq[techniq=="F"]=0
techniq[techniq=="M"]=1
techniq$gender=as.numeric(techniq$gender)

# Prevotella_copri
myplots <- list()
count = 1
for(i in 3:9){
  tech=colnames(techniq)[i]
  df=data.frame(x=(techniq[,"Prevotella_copri"]),y=techniq[,tech])
  mm=(cor.test(df$x,df$y,method = "spearman"))
  pvalue=round(mm$p.value,5)
  cor=round(mm$estimate,5)
  myplots[[count]]=ggplot(df, aes(x=x, y=y)) + 
    geom_point(color="#E18727FF",alpha=0.5)+theme_bw()+
    geom_smooth(method=lm,color="#0072B5FF")+xlab("Prevotella_copri")+ylab(tech)+
    labs(title = paste("Pvalue is ",pvalue),subtitle = paste("CorrelationCoefficient is ",cor))+theme(plot.title = element_text(size=10),plot.subtitle = element_text(size = 10))
  count = count + 1
}
ggarrange(plot(myplots[[1]]),plot(myplots[[2]]),plot(myplots[[3]]),plot(myplots[[4]]),ncol = 2,nrow = 2)
ggsave("Plot/Prevotella_copri.tech.pdf")

# Bacteroides_uniformis
myplots <- list()
count = 1
for(i in 3:9){
  tech=colnames(techniq)[i]
  df=data.frame(x=(techniq[,"Bacteroides_uniformis"]),y=techniq[,tech])
  mm=(cor.test(df$x,df$y,method = "spearman"))
  pvalue=round(mm$p.value,5)
  cor=round(mm$estimate,5)
  myplots[[count]]=ggplot(df, aes(x=x, y=y)) + 
    geom_point(color="#20854EFF",alpha=0.5)+theme_bw()+
    geom_smooth(method=lm,color="#0072B5FF")+xlab("Bacteroides_uniformis")+ylab(tech)+
    labs(title = paste("Pvalue is ",pvalue),subtitle = paste("CorrelationCoefficient is ",cor))+theme(plot.title = element_text(size=10),plot.subtitle = element_text(size = 10))
  count = count + 1
}
ggarrange(plot(myplots[[1]]),plot(myplots[[2]]),plot(myplots[[3]]),plot(myplots[[4]]),ncol = 2,nrow = 2)
ggsave("Plot/Bacteroides_uniformis.tech.pdf")

# ==================================
# two enterotypes
# ==================================
pcoa=read.table("PAM.PCoA.species.txt",stringsAsFactors = F)
enterotype=read.table("PAM.enterotype.species.txt",stringsAsFactors = F)
species_enterotype=cbind(pcoa,enterotype)
rownames(species_enterotype)=rownames(enterotype)
rownames(techniq)=techniq$Row.names
techniq$Row.names=NULL
species_enterotype=merge(species_enterotype,techniq,by="row.names",all=F)
species_enterotype$Enterotype=as.factor(species_enterotype$Enterotype)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))
sc <- scale_colour_gradientn(colours = myPalette(200))
p0=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Enterotype),size=1) + theme_bw()+scale_color_npg()+ggtitle("Two Enterotypes")
p01=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Prevotella_copri),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")
p02=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Bacteroides_uniformis),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")
p03=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =Alistipes_putredinis),size=1) + theme_bw()+sc+ggtitle("Two Enterotypes")

p1=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =`postclean.reads`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p2=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =`postclean.gc`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p3=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =`conc.ng.ul`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p4=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =`vol.ul`),size=1) + theme_bw()+sc+
  theme(legend.position = 'top')
p5=ggplot (species_enterotype, aes(V1,V2)) + geom_point(aes(colour =`Batch`),size=1) + theme_bw()+scale_color_npg()+ggtitle("Total batch")+guides(color=F)
p6=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch1",], aes(V1,V2)) + geom_point(colour ="#BC3C29FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch1")
p7=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch2",], aes(V1,V2)) + geom_point(colour ="#0072B5FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch2")
p8=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch3",], aes(V1,V2)) + geom_point(colour ="#E18727FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch3")
p9=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch4",], aes(V1,V2)) + geom_point(colour ="#20854EFF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch4")
p10=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch5",], aes(V1,V2)) + geom_point(colour ="#7876B1FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch5")
p11=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch6",], aes(V1,V2)) + geom_point(colour ="#6F99ADFF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch6")
p12=ggplot (species_enterotype[species_enterotype$Batch=="dag3_batch7",], aes(V1,V2)) + geom_point(colour ="#FFDC91FF",size=1) + theme_bw()+scale_color_npg()+ggtitle("dag3_batch7")
# plot two enterotypes
ggarrange(p0,p01,p02,p03)
ggsave("Plot/Enterotype.species.pdf",width = 14,height = 10)
# plot reads, gc content concerntration and volunm
ggarrange(p1,p2,p3,p4,ncol = 2,nrow = 2)
ggsave("Plot/Enterotype.species.tech.pdf",width = 10,height = 12)
# plot batch effects
ggarrange(p5,p6,p7,p8,p9,p10,p11,p12,ncol = 3,nrow = 3)
ggsave("Plot/Enterotype.species.batch.pdf",width = 10,height = 12)


# ==================================
# compare two enterotypes
# ==================================
enterotype$ID=rownames(enterotype)
dag3_species_plot=CompositionTable(dag3_species,10)
dag3_species_plot$Level=factor(dag3_species_plot$Level,levels = c("Bacteroides_vulgatus","Bacteroidales_bacterium_ph8","Alistipes_shahii",
                                                                  "Sutterella_wadsworthensis","Alistipes_onderdonkii","Bacteroides_uniformis",
                                                                  "Subdoligranulum_unclassified","Prevotella_copri",
                                                                  "Faecalibacterium_prausnitzii","Alistipes_putredinis"))
dag3_species_plot=merge(dag3_species_plot,enterotype,by="ID",all=T)
dag3_species_plot$Enterotype=as.factor(dag3_species_plot$Enterotype)
ggplot(dag3_species_plot, aes(Level,Relative,fill=Enterotype)) + 
  geom_boxplot(position=position_dodge(1),outlier.shape = NA)+
  theme_bw()+
  scale_color_npg()+
  scale_fill_npg()+
  guides(color=FALSE)+guides(fill=FALSE)+
  theme(axis.text.x = element_text(size = 6,hjust = 1))+
  xlab("Genera")+ylab("Relative Abundance")+coord_flip()
ggsave("Plot/Enterotype.species.compare.pdf",width = 5,height = 8)
compare=compare_means(Relative~Enterotype,data=dag3_species_plot,method = "wilcox.test",group.by = "Level")
write.table(compare,"Plot/Enterotype.species.compare.txt",sep = "\t",quote = F,row.names = F)

dag3_species_lda=merge(enterotype,dag3_species,by="row.names",all=F)
rownames(dag3_species_lda)=dag3_species_lda$Row.names
dag3_species_lda$Row.names=NULL
dag3_species_lda$ID=NULL
dag3_species_lda=dag3_species_lda[,c("Enterotype","Bacteroides_vulgatus","Bacteroidales_bacterium_ph8","Alistipes_shahii",
                                     "Sutterella_wadsworthensis","Alistipes_onderdonkii","Bacteroides_uniformis",
                                     "Subdoligranulum_unclassified","Prevotella_copri",
                                     "Faecalibacterium_prausnitzii","Alistipes_putredinis")]
compare_spearman = foreach(i=2:ncol(dag3_species_lda),.combine = rbind) %do%  {
  mm=cor.test(dag3_species_lda$Enterotype,dag3_species_lda[,i],method = "spearman")
  return.string = data.frame(Level=colnames(dag3_species_lda)[i],Coefficient=mm$estimate,Pvalue=mm$p.value)
}
















