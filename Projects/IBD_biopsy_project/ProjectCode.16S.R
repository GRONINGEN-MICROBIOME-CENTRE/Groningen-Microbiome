# this script is for analysis
library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(vegan)
library(RColorBrewer)
library(ggalluvial)
require(cowplot)
library(ggpubr)
library("dplyr")
library(edgeR)
library(lme4)
library(limma)
library(mixOmics)
library(robCompositions)
library(microbiome)
source("Microbiome.function.R")
library(lmerTest)

# *******===========================================================================================********
#                                   import data, analysis starts
#
#        genus.2000reads.txt: decontaminated, filtered, not harmonized with phenotype  
#        input.genus.txt: ready to analysis, harmonized with phenotype, before clr
#        CLR.genus.txt: read to analysis, after clr
#
# *******===========================================================================================********

covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
family=read.table("Decontamination_dada2/Input.family.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=read.table("Decontamination_dada2/Input.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
phylum=read.table("Decontamination_dada2/Input.phylum.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

# ===========================================================================================
# composition plot, core bacteria
# ===========================================================================================
phylum=phylum[,colSums(phylum)>0]
phylum=phylum[rowSums(phylum)>0,]
phylum=as.data.frame(apply(phylum,1,function(x){
  x=x/sum(x)
}))
phylum=as.data.frame(t(phylum))
phylum=phylum[,colSums(phylum>0)>0.05*nrow(phylum)]

set.seed(122)
Tax_s_order<-names(sort(apply(phylum, 2, mean),decreasing = T))[1:8]
palette=TaxaColor(Tax_s_order,level = "phyla")

dada2_phylum_plot=CompositionTable(phylum,8)
dada2_phylum_plot$Level<-factor(dada2_phylum_plot$Level,levels = rev(Tax_s_order))
dada2_phylum_plot_CD=dada2_phylum_plot[dada2_phylum_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"],]
dada2_phylum_plot_UC=dada2_phylum_plot[dada2_phylum_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"],]
dada2_phylum_plot_Control=dada2_phylum_plot[dada2_phylum_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="Control"],]

dada2_phylum_plot_CD$Location[dada2_phylum_plot_CD$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
dada2_phylum_plot_CD$Location[dada2_phylum_plot_CD$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

dada2_phylum_plot_UC$Location[dada2_phylum_plot_UC$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
dada2_phylum_plot_UC$Location[dada2_phylum_plot_UC$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

dada2_phylum_plot_Control$Location[dada2_phylum_plot_Control$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
dada2_phylum_plot_Control$Location[dada2_phylum_plot_Control$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

ggplot(dada2_phylum_plot_CD, aes(x=ID, y=Relative,fill = Level))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
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
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),
        axis.ticks.x=element_blank())+theme(legend.position="bottom")+
  facet_wrap(~Location,  scales = "free_x", drop = TRUE)+xlab("")
ggsave("OutputPlot//Composition.CD.phylum.pdf",width = 8,height = 4)

ggplot(dada2_phylum_plot_UC, aes(x=ID, y=Relative,fill = Level))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
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
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),
        axis.ticks.x=element_blank())+theme(legend.position="bottom")+
  facet_wrap(~Location,  scales = "free_x", drop = TRUE)+xlab("")
ggsave("OutputPlot//Composition.UC.phylum.pdf",width = 8,height = 4)

ggplot(dada2_phylum_plot_Control, aes(x=ID, y=Relative,fill = Level))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
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
        axis.text.x = element_text(angle = 90, hjust = 1,size = 5),
        axis.ticks.x=element_blank())+theme(legend.position="bottom")+
  facet_wrap(~Location,  scales = "free_x", drop = TRUE)+xlab("")
ggsave("OutputPlot//Composition.Control.phylum.pdf",width = 8,height = 4)

genus=genus[,colSums(genus)>0]
genus=genus[rowSums(genus)>0,]
genus=as.data.frame(apply(genus,1,function(x){
  x=x/sum(x)
}))
genus=as.data.frame(t(genus))
genus=genus[,colSums(genus>0)>0.1*nrow(genus)]

genus_plot=CompositionTable(genus,10)
Tax_s_order<-names(sort(apply(genus, 2, mean),decreasing = T))[1:10]
Tax_s_order=append(Tax_s_order,"Others")
palette=TaxaColor(Tax_s_order,level = "genus")
genus_plot_CD=genus_plot[genus_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"],]
genus_plot_UC=genus_plot[genus_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"],]
genus_plot_Control=genus_plot[genus_plot$ID %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="Control"],]

genus_plot_CD$Location[genus_plot_CD$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
genus_plot_CD$Location[genus_plot_CD$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

genus_plot_UC$Location[genus_plot_UC$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
genus_plot_UC$Location[genus_plot_UC$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

genus_plot_Control$Location[genus_plot_Control$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum"]]="ileum"
genus_plot_Control$Location[genus_plot_Control$ID %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon"]]="colon"

genus_plot_CD_ileum=AverageTable(genus_plot_CD[genus_plot_CD$Location=="ileum",])
colnames(genus_plot_CD_ileum)[1]="AverageAbundance"
genus_plot_CD_ileum$Location="ileum"
genus_plot_CD_colon=AverageTable(genus_plot_CD[genus_plot_CD$Location=="colon",])
colnames(genus_plot_CD_colon)[1]="AverageAbundance"
genus_plot_CD_colon$Location="colon"
genus_plot_CD=rbind(genus_plot_CD_ileum,genus_plot_CD_colon)

ggplot(genus_plot_CD, aes(x=Location, y=AverageAbundance,fill = Taxa))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
  ylab("Relative abundance")+
  theme_bw()
ggsave("OutputPlot//Composition.CD.genus.pdf",width = 4,height = 2)

genus_plot_UC_ileum=AverageTable(genus_plot_UC[genus_plot_UC$Location=="ileum",])
colnames(genus_plot_UC_ileum)[1]="AverageAbundance"
genus_plot_UC_ileum$Location="ileum"
genus_plot_UC_colon=AverageTable(genus_plot_UC[genus_plot_UC$Location=="colon",])
colnames(genus_plot_UC_colon)[1]="AverageAbundance"
genus_plot_UC_colon$Location="colon"
genus_plot_UC=rbind(genus_plot_UC_ileum,genus_plot_UC_colon)

ggplot(genus_plot_UC, aes(x=Location, y=AverageAbundance,fill = Taxa))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
  ylab("Relative abundance")+
  theme_bw()
ggsave("OutputPlot//Composition.UC.genus.pdf",width = 4,height = 2)

genus_plot_Control_ileum=AverageTable(genus_plot_Control[genus_plot_Control$Location=="ileum",])
colnames(genus_plot_Control_ileum)[1]="AverageAbundance"
genus_plot_Control_ileum$Location="ileum"
genus_plot_Control_colon=AverageTable(genus_plot_Control[genus_plot_Control$Location=="colon",])
colnames(genus_plot_Control_colon)[1]="AverageAbundance"
genus_plot_Control_colon$Location="colon"
genus_plot_Control=rbind(genus_plot_Control_ileum,genus_plot_Control_colon)

ggplot(genus_plot_Control, aes(x=Location, y=AverageAbundance,fill = Taxa))+
  geom_bar(position = "stack",stat="identity")+
  scale_fill_manual(breaks = Tax_s_order,
                    values = palette)+
  ylab("Relative abundance")+
  theme_bw()
ggsave("OutputPlot//Composition.Control.genus.pdf",width = 4,height = 2)

genus_CD_colon=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon" & covariate_bac$Diagnosis=="CD"],]
genus_CD_ileum=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum" & covariate_bac$Diagnosis=="CD"],]
genus_UC_colon=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon" & covariate_bac$Diagnosis=="UC"],]
genus_UC_ileum=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum" & covariate_bac$Diagnosis=="UC"],]
genus_control_colon=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="colon" & covariate_bac$Diagnosis=="Control"],]
genus_control_ileum=genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Location_rough=="ileum" & covariate_bac$Diagnosis=="Control"],]

summary_statistics_metadata(genus_control_ileum)

# ===========================================================================================
# bacteria PCA, shannon/PCA with surgery
# ===========================================================================================

# shannon diversity 
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F)
genus=genus[rowSums(genus)>0,colSums(genus>0)>(0.1*nrow(genus))]
shannon <- as.data.frame(diversity(t(genus), index="shannon"))
shannon=merge(shannon,covariate_bac,by="row.names",all=F)
shannon=shannon[(shannon$Diagnosis!=""),]
shannon$Diagnosis=factor(shannon$Diagnosis,levels = c("CD","UC","Control"))
ggplot(shannon, aes(x=Diagnosis, y=shannon, fill=Diagnosis)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_manual(values=c("#EE8866","#AA4499","#44BB99"))+
  xlab("")+
  theme_classic() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') 
ggsave("OutputPlot/Shannono.Diagnosis.pdf",width = 3,height = 2.5)

shannon$Split=paste(shannon$Diagnosis, shannon$Location_rough)
ggplot(shannon, aes(x=Split, y=shannon, fill=Split)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_npg()+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top')+
  xlab("")
ggsave("OutputPlot/test.pdf",width = 6,height = 3)
kruskal.test(shannon ~ Diagnosis, data = shannon)
wilcox.test(shannon$shannon[shannon$Split=="CD ileum"],shannon$shannon[shannon$Split=="Control colon"])
wilcox.test(shannon$shannon[shannon$Diagnosis=="UC"],shannon$shannon[shannon$Diagnosis=="Control"])
wilcox.test(shannon$shannon[shannon$Diagnosis=="CD"],shannon$shannon[shannon$Diagnosis=="UC"])
write.table(shannon,"OutputTable/Shannon.genus.txt",sep = "\t",quote = F,row.names = F)

# PCA
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
beta_diversity=vegdist((genus),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,covariate_bac,all = F,by="row.names")
pca_analysis[pca_analysis==""]=NA
pca_analysis$resec_ileocec=as.factor(pca_analysis$resec_ileocec)
pca_analysis=pca_analysis[!is.na(pca_analysis$Diagnosis),]
ggplot (pca_analysis, aes(V1,V2,color=Diagnosis)) + 
  geom_point(size=2) + theme_classic() +
  guides(size=F)+
  scale_color_npg()+
  stat_ellipse(type = "norm")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') +xlab("PCA1")+ylab("PCA2")
ggscatter(pca_analysis, x = "V1", y = "V2",
          color = "Diagnosis", palette = c("#EE8866","#44BB99","#AA4499"),
          ellipse = TRUE, 
          mean.point = TRUE,
          star.plot = TRUE,alpha=0.2,xlab = "PC1",ylab = "PC2")
ggsave("OutputPlot//PCA.bac.diagnosis.pdf",width = 3.6,height = 2.8)

distmat <- as.matrix(beta_diversity)
cohort1.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="CD"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="CD"]]
cohort2.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="UC"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="UC"]]
cohort3.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="Control"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$Diagnosis=="Control"]]

distmat1=data.frame(as.table(cohort1.dis))[lower.tri(cohort1.dis, diag = F), ]
distmat2=data.frame(as.table(cohort2.dis))[lower.tri(cohort2.dis, diag = F), ]
distmat3=data.frame(as.table(cohort3.dis))[lower.tri(cohort3.dis, diag = F), ]

distmat1$group="CD"
distmat2$group="UC"
distmat3$group="Control"
distmat=rbind(distmat1,distmat2,distmat3)

distmat$group=factor(distmat$group,levels = c("CD","UC","Control"))
ggplot(distmat, aes(x=group, y=Freq,fill=group)) + 
  geom_boxplot(color="black",size=0.3)+scale_fill_manual(values=c("#EE8866","#AA4499","#44BB99"))+theme_classic()+guides(fill=F)
ggsave("OutputPlot/dissimilarity.genus.pdf",width = 3.2,height = 2.2)
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="CD"])
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="UC"])
wilcox.test(distmat$Freq[distmat$group=="CD"],distmat$Freq[distmat$group=="UC"])

## PC1 correlations
PC1_BC_cors = apply(genus,2,function(x){
  cor(pca_analysis$V1,x,method = "spearman")})

PC1_BC_pvalues = p.adjust(apply(genus,2,function(x){
  cor.test(pca_analysis$V1,x,method = "spearman")$p.value}),method = "fdr")

PC1_pvalues = PC1_BC_pvalues[names(PC1_BC_cors)]
PC1_cors = data.frame("PC1_cors" = PC1_BC_cors,"PC1_pvalues" = PC1_pvalues)


## PC2 correlations
PC2_BC_cors = apply(genus,2,function(x){
  cor(pca_analysis$V2,x,method = "spearman")})

PC2_BC_pvalues = p.adjust(apply(genus,2,function(x){
  cor.test(pca_analysis$V2,x,method = "spearman")$p.value}),method = "fdr")

PC2_pvalues = PC2_BC_pvalues[names(PC2_BC_cors)]
PC2_cors = data.frame("PC2_cors" = PC2_BC_cors,"PC2_pvalues" = PC2_pvalues)

PC_cors = data.frame(cbind(PC1_cors,PC2_cors))
PC_cors = PC_cors[(PC_cors$PC1_pvalues <0.05 & PC_cors$PC2_pvalues <0.05),]
ggplot(pca_analysis,aes(x = V1, y = V2))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = -0.97439692, yend = 0.4439545),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.50872905, yend = -0.9060646),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.36393732, yend = -0.4866545),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.35510964, yend = -0.3966603),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.21104951, yend = -0.2429104),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.20018007, yend = -0.3081907),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.11328872, yend = -0.2513477),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.08572697, yend = 0.1271985),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.07994902, yend = -0.4135981),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
  geom_segment(color="black",linetype=1,aes(x = 0, y = 0, xend = 0.07526764, yend = 0.1891876),size=.7, arrow = arrow(length = unit(0.5, "cm")))+
#  geom_text(x=-0.98109247,y=0.39593729,label="Bacteroides",color="black",fontface="italic")+
#  geom_text(x=0.41847382,y=-0.85353803,label="Faecalibacterium",color="black",fontface="italic")+
#  geom_text(x=0.30945125,y=-0.44160344,label="Agathobacter",color="black",fontface="italic")+
#  geom_text(x=0.28998603,y=-0.50774242,label="Blautia",color="black",fontface="italic")+
#  geom_text(x=0.21115063,y=-0.30497296,label="Bifidobacterium",color="black",fontface="italic")+
#  geom_text(x=0.08043040,y=0.11420504,label="Escherichia/Shigella",color="black",fontface="italic")+
  theme_classic()
ggsave("OutputPlot/Arrows.pdf",width = 8,height = 8)

# ===========================================================================================
# bacteria surgery within disease
# ===========================================================================================

covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]
genus=genus[!rownames(genus) %in% remove,]
genus=genus[!rownames(genus) %like% "1056",]
genus=genus[!rownames(genus) %like% "680",]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(genus),]
genus=genus[rownames(genus) %in% rownames(covariate_bac),]

covariate_bac$Surgery=NA
covariate_bac$Surgery[covariate_bac$resec_part_colon==1]="Part.colon"
covariate_bac$Surgery[covariate_bac$resec_ileocec==1]="Ileocec"
covariate_bac$Surgery[covariate_bac$resec_part_small==1]="Part.small"
covariate_bac$Surgery[is.na(covariate_bac$Surgery)]="No.surgery"
covariate_bac$Surgery[covariate_bac$resec_ileocec==1 & covariate_bac$resec_part_colon==1]="Multi.surgery"
covariate_bac$Surgery[covariate_bac$resec_part_small==1 & covariate_bac$resec_part_colon==1]="Multi.surgery"
covariate_bac$Surgery[covariate_bac$resec_part_small==1 & covariate_bac$resec_ileocec==1]="Multi.surgery"

bacteria=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
shannon=as.data.frame(diversity(t(bacteria), index="shannon"))
colnames(shannon)="shannon"
shannon=merge(shannon,covariate_bac,by="row.names",all = F)
shannon$Surgery=factor(shannon$Surgery,levels = c("Multi.surgery","Ileocec","Part.colon","No.surgery","Part.small"))

# CD, part colon, ileocec and part small intestine
beta_diversity=vegdist((genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"],]),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=4))
pca_analysis=merge(pca_analysis,covariate_bac,all = F,by="row.names")
pca_analysis$resec_ileocec=as.factor(pca_analysis$resec_ileocec)
ggplot (pca_analysis, aes(V1,V2,color=Surgery)) + 
  geom_point() + theme_bw() +
  guides(size=F)+
  scale_color_npg()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")

shannon_CD=shannon[shannon$Diagnosis=="CD",]
ggplot(shannon_CD, aes(x=Surgery, y=shannon, fill=Surgery)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_jama()+
  xlab("")+theme_bw()+guides(fill=F)
ggsave("OutputPlot/CD.surgery.pdf",width = 8,height = 3)
wilcox.test(shannon_CD$shannon[shannon_CD$Surgery=="Part.small"],shannon_CD$shannon[shannon_CD$Surgery=="No.surgery"])

# UC, part colon
beta_diversity=vegdist((genus[rownames(genus) %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"],]),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=4))
pca_analysis=merge(pca_analysis,covariate_bac,all = F,by="row.names")
pca_analysis$resec_part_colon=as.factor(pca_analysis$resec_part_colon)

ggplot (pca_analysis, aes(V1,V2,color=resec_part_colon)) + 
  geom_point() + theme_bw() +
  guides(size=F)+
  scale_color_npg()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")

shannon_UC=shannon[shannon$Diagnosis=="UC",]
shannon_UC=shannon_UC[shannon_UC$Surgery!="Ileocec" & shannon_UC$Surgery!="Multi.surgery",]
ggplot(shannon_UC, aes(x=Surgery, y=shannon, fill=Surgery)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_jama()+
  xlab("")+theme_bw()+guides(fill=F)
ggsave("OutputPlot/UC.surgery.pdf",width = 4,height = 3)
wilcox.test(shannon_UC$shannon[shannon_UC$Surgery=="Part.colon"],shannon_UC$shannon[shannon_UC$Surgery=="No.surgery"])

# compare part colon effect at CD and UC
bacteria=read.table("OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria_CD=bacteria[rownames(bacteria) %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="CD"],]
bacteria_CD=bacteria_CD[rownames(bacteria_CD) %in% rownames(covariate_bac)[covariate_bac$Surgery=="Part.colon" | covariate_bac$Surgery=="No.surgery"],]
bacteria_UC=bacteria[rownames(bacteria) %in% rownames(covariate_bac)[covariate_bac$Diagnosis=="UC"],]
bacteria_UC=bacteria_UC[rownames(bacteria_UC) %in% rownames(covariate_bac)[covariate_bac$Surgery=="Part.colon" | covariate_bac$Surgery=="No.surgery"],]

result_CD = foreach(i=1:ncol(bacteria_CD),.combine = rbind) %do%  {
  tmp.bacteria=colnames(bacteria_CD)[i]
  x=bacteria_CD[,i,drop=F]
  tmp.data=merge(x,covariate_bac,by="row.names",all=F)
  tmp.data$Row.names=NULL
  tmp.data$Surgery[tmp.data$Surgery=="Part.colon"]=1
  tmp.data$Surgery[tmp.data$Surgery=="No.surgery"]=0
  tmp.data$Surgery=as.numeric(tmp.data$Surgery)
  
  fit<- lm(tmp.data[,tmp.bacteria]~Surgery+age_at_biopsy+sex+BMI,data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit)))
    cat(tmp.bacteria,"is done","\n")
    return.string=data.frame(bacteria=tmp.bacteria,surgery.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Surgery"],
                             surgery.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Surgery"],
                             surgery.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Surgery"])

}
result_UC = foreach(i=1:ncol(bacteria_UC),.combine = rbind) %do%  {
  tmp.bacteria=colnames(bacteria_UC)[i]
  x=bacteria_UC[,i,drop=F]
  tmp.data=merge(x,covariate_bac,by="row.names",all=F)
  tmp.data$Row.names=NULL
  tmp.data$Surgery[tmp.data$Surgery=="Part.colon"]=1
  tmp.data$Surgery[tmp.data$Surgery=="No.surgery"]=0
  tmp.data$Surgery=as.numeric(tmp.data$Surgery)
  
  fit<- lm(tmp.data[,tmp.bacteria]~Surgery+age_at_biopsy+sex+BMI,data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit)))
  cat(tmp.bacteria,"is done","\n")
  
  return.string=data.frame(bacteria=tmp.bacteria,surgery.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Surgery"],
                           surgery.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Surgery"],
                           surgery.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Surgery"])
  
}
result_CD=result_CD[order(result_CD$surgery.P),]
result_UC=result_UC[order(result_UC$surgery.P),]
result_CD$FDR=p.adjust(result_CD$surgery.P,method = "BH")
result_UC$FDR=p.adjust(result_UC$surgery.P,method = "BH")

# ===========================================================================================
# bacteria abundance transformation, clr
# ===========================================================================================

# fucntions are from Johannes, rows are samples and colums are taxa, and the input is abundance table
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
family=read.table("Decontamination_dada2/Input.family.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
phylum=read.table("Decontamination_dada2/Input.phylum.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
class=read.table("Decontamination_dada2/Input.class.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
ASV=read.table("Decontamination_dada2/Input.ASV.txt",sep = " ",row.names = 1,header = T,check.names = F,stringsAsFactors = F,quote = "")

genus=genus[,colSums(genus)>0]
genus=genus[rowSums(genus)>0,]
genus=apply(genus,1,function(x){
  x=x/sum(x)
})
genus=as.data.frame(t(genus))

family=family[,colSums(family)>0]
family=family[rowSums(family)>0,]
family=apply(family,1,function(x){
  x=x/sum(x)
})
family=as.data.frame(t(family))

phylum=phylum[,colSums(phylum)>0]
phylum=phylum[rowSums(phylum)>0,]
phylum=apply(phylum,1,function(x){
  x=x/sum(x)
})
phylum=as.data.frame(t(phylum))

class=class[,colSums(class)>0]
class=class[rowSums(class)>0,]
class=apply(class,1,function(x){
  x=x/sum(x)
})
class=as.data.frame(t(class))

orders=orders[,colSums(orders)>0]
orders=orders[rowSums(orders)>0,]
orders=apply(orders,1,function(x){
  x=x/sum(x)
})
orders=as.data.frame(t(orders))

genus_clr <- zCompositions::cmultRepl(genus, method="CZM", label=0)
genus_clr = genus_clr[,colSums(genus>0)>nrow(genus) * 0.1]

family_clr <- zCompositions::cmultRepl(family, method="CZM", label=0)
family_clr = family_clr[,colSums(family>0)>nrow(family) * 0.1]

order_clr <- zCompositions::cmultRepl(orders, method="CZM", label=0)
order_clr = order_clr[,colSums(orders>0)>nrow(orders) * 0.1]

phylum_clr <- zCompositions::cmultRepl(phylum, method="CZM", label=0)
phylum_clr = phylum_clr[,colSums(phylum>0)>nrow(phylum) * 0.1]

class_clr <- zCompositions::cmultRepl(class, method="CZM", label=0)
class_clr = class_clr[,colSums(class>0)>nrow(class) * 0.1]

write.table(genus_clr,"OutputTable/CLR.genus.txt",row.names = T,quote = F,sep = "\t")
write.table(family_clr,"OutputTable/CLR.family.txt",row.names = T,quote = F,sep = "\t")
write.table(order_clr,"OutputTable/CLR.order.txt",row.names = T,quote = F,sep = "\t")
write.table(phylum_clr,"OutputTable/CLR.phylum.txt",row.names = T,quote = F,sep = "\t")
write.table(class_clr,"OutputTable/CLR.class.txt",row.names = T,quote = F,sep = "\t")
write.table(cbind(genus_clr,family_clr,class_clr,order_clr,phylum_clr),"OutputTable/CLR.bacteria.txt",row.names = T,quote = F,sep = "\t")

covariate_bac=covariate_bac[!is.na(covariate_bac$Diagnosis),]
covariate_bac=covariate_bac[!is.na(covariate_bac$Inflammation),]
covariate_bac$BMI[is.na(covariate_bac$BMI)]=median(covariate_bac$BMI[!is.na(covariate_bac$BMI)])
covariate_bac$smoking_DB[is.na(covariate_bac$smoking_DB)]=median(covariate_bac$smoking_DB[!is.na(covariate_bac$smoking_DB)])

write.table(covariate_bac,"OutputTable/Covariate_bac.organized.txt",sep = "\t",row.names = T,quote = F)

# ===========================================================================================
# bacteria abundance characterization, inflammation/disease
# ===========================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(bacteria),]

# CD vs. HC at ileum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
covariate1=covariate1[,c("Diagnosis","Inflammation","ResearchID","age_at_biopsy","sex","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result1 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>100){
    fit<-lm(tmp.data[,tmp.genus]~Diagnosis,data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit)))
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,diagnosis.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.beta=coefs.fit$Estimate[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.se=coefs.fit$Std..Error[rownames(coefs.fit)=="DiagnosisControl"])
  }else{return.string=data.frame(bacteria=tmp.genus,diagnosis.P=NA,
                                 diagnosis.beta=NA,
                                 diagnosis.se=NA)}

}

# CD vs. HC at colon
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
covariate1=covariate1[,c("Diagnosis","Inflammation","ResearchID","age_at_biopsy","sex","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result2 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>100){
    fit<-lm(tmp.data[,tmp.genus]~Diagnosis,data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit)))
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,diagnosis.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.beta=coefs.fit$Estimate[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.se=coefs.fit$Std..Error[rownames(coefs.fit)=="DiagnosisControl"])
  }else{return.string=data.frame(bacteria=tmp.genus,diagnosis.P=NA,
                                 diagnosis.beta=NA,
                                 diagnosis.se=NA)}
  
}


# UC vs. HC at colon
covariate1=covariate_bac[covariate_bac$Diagnosis=="UC" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
covariate1=covariate1[,c("Diagnosis","Inflammation","ResearchID","age_at_biopsy","sex","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
covariate1$Diagnosis=factor(covariate1$Diagnosis,levels = c("UC","Control"))
result3 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>100){
    fit<-lm(tmp.data[,tmp.genus]~Diagnosis,data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit)))
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,diagnosis.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.beta=coefs.fit$Estimate[rownames(coefs.fit)=="DiagnosisControl"],
                             diagnosis.se=coefs.fit$Std..Error[rownames(coefs.fit)=="DiagnosisControl"])
  }else{return.string=data.frame(bacteria=tmp.genus,diagnosis.P=NA,
                                 diagnosis.beta=NA,
                                 diagnosis.se=NA)}
  
}

result1$group=1
result2$group=2
result3$group=3

result1$FDR=p.adjust(result1$diagnosis.P)
result2$FDR=p.adjust(result2$diagnosis.P)
result3$FDR=p.adjust(result3$diagnosis.P)
result=rbind(result1,result2,result3)

tmp.data$Diagnosis=as.factor(tmp.data$Inflammation)
ggplot(tmp.data, aes(x=Inflammation, y=tmp.data[,1],fill=Inflammation)) + 
  geom_boxplot(color="black")+scale_fill_nejm()+theme_bw()+guides(color=F)
ggsave("OutputPlot/Faecalibacterium.UC.colon.pdf")

write.table(result,"OutputTable/Comparison.composition.txt",sep = "\t",row.names = F)

#===========================================================================================
# bacteria presence/absennce characterization, inflammation
# ===========================================================================================
bacteria=read.table("Decontamination_dada2/Input.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria[bacteria>0]=1
bacteria=bacteria[,colSums(bacteria!=0)>0.1*nrow(bacteria)]

# CD vs. HC at ileum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]
covariate1=covariate1[,c("Inflammation","ResearchID","age_at_biopsy","sex","BMI","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result1 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tryCatch(  {
  tmp.genus=colnames(bacteria1)[i]
  tmp.data=merge(bacteria1[,i,drop=F],covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  #tmp.data=tmp.data[tmp.data$Inflammation!=1,]
  
  tmp.data$sex=scale(tmp.data$sex)
  tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
  tmp.data$Inflammation=scale(tmp.data$Inflammation)
  #tmp.data$Inflammation=as.factor(tmp.data$Inflammation)
  
  fit0 <- glmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+(1|Batch),data=tmp.data,family = binomial(link='logit'))
  #resid=residuals(fit0)
  #tmp.data$feature=resid
  #fit1<-lm(feature~Inflammation,data=tmp.data)
  #coefs.fit <- data.frame(coef(summary(fit1)))
  cat(tmp.genus,"is done","\n")
  
  #fit0 <- glm(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI,data=tmp.data,family = binomial(link='logit'))
  #fit2 <- lmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+(1|ResearchID),data=tmp.data)
  coefs.fit <- data.frame(coef(summary(fit0)))
  
  count.data=table(tmp.data$Inflammation,tmp.data[,1])
  ptab <- as.data.frame(prop.table(count.data, margin=1))
  ptab=ptab[ptab$Var2==1,]

  return.string=data.frame(bacteria=tmp.genus,inflammation.P=coefs.fit$Pr...z..[rownames(coefs.fit)=="Inflammation"],
                           inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                           inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"],
                           #inflammation1=ptab$Freq[ptab$Var1==1],
                           inflammation2=ptab$Freq[1],
                           inflammation3=ptab$Freq[2]
                           )}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}

# CD vs. HC at colon
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
covariate1=covariate1[,c("Inflammation","ResearchID","age_at_biopsy","sex","BMI","smoking_DB","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result2 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tryCatch(  {
    tmp.genus=colnames(bacteria1)[i]
    tmp.data=merge(bacteria1[,i,drop=F],covariate1,by="row.names",all=F)
    tmp.data$Row.names=NULL
    
    #tmp.data=tmp.data[tmp.data$Inflammation!=1,]
    
    tmp.data$sex=scale(tmp.data$sex)
    tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
    tmp.data$Inflammation=scale(tmp.data$Inflammation)
    
    cat(tmp.genus,"is done","\n")
    
    fit2 <- glmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+(1|Batch),data=tmp.data,family = binomial(link='logit'))
    coefs.fit <- data.frame(coef(summary(fit2)))
    
    count.data=table(tmp.data$Inflammation,tmp.data[,1])
    ptab <- as.data.frame(prop.table(count.data, margin=1))
    ptab=ptab[ptab$Var2==1,]
    
    return.string=data.frame(bacteria=tmp.genus,inflammation.P=coefs.fit$Pr...z..[rownames(coefs.fit)=="Inflammation"],
                             inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                             inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"],
                             #inflammation1=ptab$Freq[ptab$Var1==1],
                             inflammation2=ptab$Freq[1],
                             inflammation3=ptab$Freq[2]
    )}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}


# UC vs. HC at colon
covariate1=covariate_bac[covariate_bac$Diagnosis=="UC" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="colon",]
covariate1=covariate1[,c("Inflammation","ResearchID","age_at_biopsy","sex","BMI","smoking_DB","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result3 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tryCatch(  {
    tmp.genus=colnames(bacteria1)[i]
    tmp.data=merge(bacteria1[,i,drop=F],covariate1,by="row.names",all=F)
    tmp.data$Row.names=NULL
    
    #tmp.data=tmp.data[tmp.data$Inflammation!=1,]
    
    tmp.data$sex=scale(tmp.data$sex)
    tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
    tmp.data$Inflammation=scale(tmp.data$Inflammation)
    
    cat(tmp.genus,"is done","\n")
    
    fit2 <- glmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+(1|Batch),data=tmp.data,family = binomial(link='logit'))
    coefs.fit <- data.frame(coef(summary(fit2)))
    
    count.data=table(tmp.data$Inflammation,tmp.data[,1])
    ptab <- as.data.frame(prop.table(count.data, margin=1))
    ptab=ptab[ptab$Var2==1,]
    
    return.string=data.frame(bacteria=tmp.genus,inflammation.P=coefs.fit$Pr...z..[rownames(coefs.fit)=="Inflammation"],
                             inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                             inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"],
                             #inflammation1=ptab$Freq[ptab$Var1==1],
                             inflammation2=ptab$Freq[1],
                             inflammation3=ptab$Freq[2]
    )}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}


result1$group=1
result2$group=2
result3$group=3
result=rbind(result1,result2,result3)
result$FDR=p.adjust(result$inflammation.P)

tmp.data$Inflammation=as.numeric(as.character(tmp.data$Inflammation))
tmp.data$Acidaminococcaceae=as.factor(tmp.data$`Acidaminococcaceae`)
ggplot(tmp.data,aes(x = (Inflammation), group = (`Acidaminococcaceae`), fill = `Acidaminococcaceae`)) + 
  geom_bar(position = "fill")+theme_bw()+scale_fill_aaas()
ggsave("OutputPlot/Ruminococcaceae_UCG-002.absence.pdf",width = 5,height = 4)

# combined inflamed vs. non-inflamed
covariate1=covariate_bac[covariate_bac$Diagnosis=="UC" | covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[,c("Inflammation","ResearchID","age_at_biopsy","sex","Batch","Location_rough","BMI")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result0 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tryCatch(  {
    tmp.genus=colnames(bacteria1)[i]
    tmp.data=merge(bacteria1[,i,drop=F],covariate1,by="row.names",all=F)
    tmp.data$Row.names=NULL
    
    tmp.data=tmp.data[tmp.data$Inflammation!=1,]
    
    tmp.data$sex=scale(tmp.data$sex)
    tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
    tmp.data$Inflammation=scale(tmp.data$Inflammation)
    
    cat(tmp.genus,"is done","\n")
    
    fit2 <- lmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+Location_rough+(1|Batch),data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit2)))
    
    count.data=table(tmp.data$Inflammation,tmp.data[,1])
    ptab <- as.data.frame(prop.table(count.data, margin=1))
    ptab=ptab[ptab$Var2==1,]
    
    return.string=data.frame(bacteria=tmp.genus,inflammation.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Inflammation"],
                             inflammation.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Inflammation"],
                             inflammation.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Inflammation"],
                             #inflammation1=ptab$Freq[ptab$Var1==1],
                             inflammation2=ptab$Freq[1],
                             inflammation3=ptab$Freq[2]
    )}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}
result0$FDR=p.adjust(result0$inflammation.P)

tmp.data$Inflammation=as.numeric(as.character(tmp.data$Inflammation))
tmp.data$Blautia=as.factor(tmp.data$Blautia)
ggplot(tmp.data,aes(x = (Inflammation), group = (Blautia), fill = Blautia)) + 
  geom_bar(position = "fill")+theme_bw()+scale_fill_aaas()
ggsave("OutputPlot/Blautia.absence.combined.pdf",width = 5,height = 4)


# ===========================================================================================
# bacteria abundance characterization, location
# ===========================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)

# CD, colon vs ileum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD",]
covariate1=covariate1[,c("Location_rough","ResearchID","age_at_biopsy","sex","BMI","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result4 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  tmp.data$sex=scale(tmp.data$sex)
  tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
  
  if(nrow(tmp.data)>100){
    fit<-lmer(tmp.data[,tmp.genus]~Location_rough+age_at_biopsy+sex+BMI+(1 | Batch),data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit)))
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,location.P=coefs.fit$Pr...t..[rownames(coefs.fit)=="Location_roughileum"],
                             location.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Location_roughileum"],
                             location.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Location_roughileum"])
  }else{return.string=data.frame(bacteria=tmp.genus,location.P=NA,
                                 location.beta=NA,
                                 location.se=NA)}
  
}
result4$FDR=p.adjust(result4$location.P,method = "BH")

tmp.data$Location_rough=as.factor(tmp.data$Location_rough)
ggplot(tmp.data, aes(x=Location_rough, y=tmp.data[,1],fill=Location_rough)) + 
  geom_boxplot(color="black")+scale_fill_nejm()+theme_bw()+guides(color=F)



#===========================================================================================
# bacteria presence/absennce characterization, location
# ===========================================================================================
bacteria=read.table("Decontamination_dada2/Input.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria[bacteria>0]=1
bacteria=bacteria[,colSums(bacteria!=0)>0.1*nrow(bacteria)]

# CD, colon vs lieum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[,c("Location_rough","ResearchID","age_at_biopsy","sex","BMI","Batch")]
bacteria1=bacteria[rownames(bacteria) %in% rownames(covariate1),]
covariate1=covariate1[rownames(covariate1) %in% rownames(bacteria1),]

bacteria1=bacteria1[order(rownames(bacteria1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result5 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tryCatch(  {
    tmp.genus=colnames(bacteria1)[i]
    tmp.data=merge(bacteria1[,i,drop=F],covariate1,by="row.names",all=F)
    tmp.data$Row.names=NULL
    
    tmp.data$sex=scale(tmp.data$sex)
    tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
    
    #tmp.data$Inflammation=as.factor(tmp.data$Inflammation)
    
    fit0 <- glmer(tmp.data[,tmp.genus]~Location_rough+age_at_biopsy+sex+BMI+(1|Batch),data=tmp.data,family = binomial(link='logit'))
    #resid=residuals(fit0)
    #tmp.data$feature=resid
    #fit1<-lm(feature~Inflammation,data=tmp.data)
    coefs.fit <- data.frame(coef(summary(fit0)))
    cat(tmp.genus,"is done","\n")
    
    #fit2 <- lmer(tmp.data[,tmp.genus]~Location_rough+age_at_biopsy+sex+BMI+(1|ResearchID)+(1|Batch),data=tmp.data)
    #coefs.fit <- data.frame(coef(summary(fit2)))

    return.string=data.frame(bacteria=tmp.genus,location.P=coefs.fit$Pr...z..[rownames(coefs.fit)=="Location_roughileum"],
                             location.beta=coefs.fit$Estimate[rownames(coefs.fit)=="Location_roughileum"],
                             location.se=coefs.fit$Std..Error[rownames(coefs.fit)=="Location_roughileum"]
    )}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}

result5$FDR=p.adjust(result5$location.P,method = "BH")




# ===========================================================================================
# bacteria abundance characterization, diagnosis
# ===========================================================================================

# CD and UC vs. HC
covariate1=covariate_bac[covariate_bac$Diagnosis!="1000IBD_surgery",]
covariate1$Inflammation[covariate1$Inflammation=="Yes"]=3
covariate1$Inflammation[covariate1$Inflammation=="No" & covariate1$Diagnosis=="UC"]=2
covariate1$Inflammation[covariate1$Inflammation=="No" & covariate1$Diagnosis=="Control"]=1
covariate1=covariate1[,c("Diagnosis","Inflammation","ResearchID","age_at_biopsy","sex","BMI","smoking_DB","Batch","Location_rough","Cohort")]
covariate1$Batch=as.factor(covariate1$Batch)
covariate1$Inflammation=as.numeric(covariate1$Inflammation)
covariate1$ResearchID=as.factor(covariate1$ResearchID)
covariate1$age_at_biopsy[is.na(covariate1$age_at_biopsy)]=median(covariate1$age_at_biopsy[!is.na(covariate1$age_at_biopsy)])
covariate1$BMI[is.na(covariate1$BMI)]=median(covariate1$BMI[!is.na(covariate1$BMI)])
covariate1$smoking_DB[is.na(covariate1$smoking_DB)]=median(covariate1$smoking_DB[!is.na(covariate1$smoking_DB)])
covariate1$Location_rough=as.factor(covariate1$Location_rough)
covariate1$Diagnosis=as.factor(covariate1$Diagnosis)
covariate1$Cohort[covariate1$Cohort=="1000IBD"]=1
covariate1$Cohort[covariate1$Cohort=="Control"]=0
covariate1$Cohort=as.numeric(covariate1$Cohort)
which(is.na(covariate1))
genus1=genus.clr[rownames(genus.clr) %in% rownames(covariate1),]

genus1=genus1[order(rownames(genus1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result0 = foreach(i=1:ncol(genus1),.combine = rbind) %do%  {
  tmp.genus=colnames(genus1)[i]
  x=genus1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>30){
    fit<-lmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+smoking_DB+Batch+Location_rough+(1|ResearchID),data=tmp.data)
    residuals=as.data.frame(resid(fit))
    mm1=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="CD"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
    mm2=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="UC"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
    nn=wilcox.test(residuals$`resid(fit)`[tmp.data$Cohort==1],residuals$`resid(fit)`[tmp.data$Cohort==0])
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,diagnosis.CD=mm1$p.value,diagnosis.UC=mm2$p.value,cohort.p=nn$p.value)
  }else{return.string=data.frame(bacteria=tmp.genus,
                                 diagnosis.CD=NA,diagnosis.UC=NA,cohort.p=NA)}
}

# ===========================================================================================
# bacteria absence/presence characterization, diagnosis
# ===========================================================================================

# CD and UC vs. HC
covariate1=covariate_bac[covariate_bac$Diagnosis!="1000IBD_surgery",]
covariate1$Inflammation[covariate1$Inflammation=="Yes"]=3
covariate1$Inflammation[covariate1$Inflammation=="No" & covariate1$Diagnosis=="UC"]=2
covariate1$Inflammation[covariate1$Inflammation=="No" & covariate1$Diagnosis=="Control"]=1
covariate1=covariate1[,c("Diagnosis","Inflammation","ResearchID","age_at_biopsy","sex","BMI","smoking_DB","Batch","Location_rough","Cohort")]
covariate1$Batch=as.factor(covariate1$Batch)
covariate1$Inflammation=as.numeric(covariate1$Inflammation)
covariate1$ResearchID=as.factor(covariate1$ResearchID)
covariate1$age_at_biopsy[is.na(covariate1$age_at_biopsy)]=median(covariate1$age_at_biopsy[!is.na(covariate1$age_at_biopsy)])
covariate1$BMI[is.na(covariate1$BMI)]=median(covariate1$BMI[!is.na(covariate1$BMI)])
covariate1$smoking_DB[is.na(covariate1$smoking_DB)]=median(covariate1$smoking_DB[!is.na(covariate1$smoking_DB)])
covariate1$Location_rough=as.factor(covariate1$Location_rough)
covariate1$Cohort[covariate1$Cohort=="1000IBD"]=1
covariate1$Cohort[covariate1$Cohort=="Control"]=0
covariate1$Cohort=as.numeric(covariate1$Cohort)
which(is.na(covariate1))

genus1=genus
genus1[genus1>0]=1
genus1=genus1[rownames(genus1) %in% rownames(covariate1),]

genus1=genus1[order(rownames(genus1)),]
covariate1=covariate1[order(rownames(covariate1)),]
result0 = foreach(i=1:ncol(genus1),.combine = rbind) %do%  {
  tryCatch(  {
  tmp.genus=colnames(genus1)[i]
  tmp.data=merge(genus1[,i,drop=F],covariate1,by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  tmp.data$Row.names=NULL
  tmp.data$age_at_biopsy=scale(tmp.data$age_at_biopsy)
  tmp.data$BMI=scale(tmp.data$BMI)
  fit<-glmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+smoking_DB+Batch+(1|ResearchID),data=tmp.data,family = binomial)
  residuals=as.data.frame(resid(fit))
  mm1=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="CD"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
  mm2=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="UC"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
  nn=wilcox.test(residuals$`resid(fit)`[tmp.data$Cohort==1],residuals$`resid(fit)`[tmp.data$Cohort==0])
  cat(tmp.genus,"is done","\n")
  
  return.string=data.frame(bacteria=tmp.genus,diagnosis.CD=mm1$p.value,diagnosis.UC=mm2$p.value,cohort.p=nn$p.value)}, error=function(e){cat("Error",conditionMessage(e), "\n")})
}
result0$diagnosis.CD.FDR=p.adjust(result0$diagnosis.CD)
result0$diagnosis.UC.FDR=p.adjust(result0$diagnosis.UC)
result0$cohort.p.FDR=p.adjust(result0$cohort.p)

tmp.CD=genus1[rownames(genus1) %in% rownames(covariate1)[covariate1$Diagnosis=="CD"],]
tmp.CD.pre=data.frame(Pre=colSums(tmp.CD))
tmp.CD.pre$Pre=tmp.CD.pre$Pre/nrow(tmp.CD)
  
tmp.UC=genus1[rownames(genus1) %in% rownames(covariate1)[covariate1$Diagnosis=="UC"],]
tmp.UC.pre=data.frame(Pre=colSums(tmp.UC))
tmp.UC.pre$Pre=tmp.UC.pre$Pre/nrow(tmp.UC)

tmp.IBD=genus1[rownames(genus1) %in% rownames(covariate1)[covariate1$Diagnosis=="CD" | covariate1$Diagnosis=="UC"],]
tmp.IBD.pre=data.frame(Pre=colSums(tmp.IBD))
tmp.IBD.pre$Pre=tmp.IBD.pre$Pre/nrow(tmp.IBD)

tmp.HC=genus1[rownames(genus1) %in% rownames(covariate1)[covariate1$Diagnosis=="Control"],]
tmp.HC.pre=data.frame(Pre=colSums(tmp.HC))
tmp.HC.pre$Pre=tmp.HC.pre$Pre/nrow(tmp.HC)

tmp.CD.pre$Group="CD"
tmp.UC.pre$Group="UC"
tmp.IBD.pre$Group="IBD"
tmp.HC.pre$Group="HC"
tmp.CD.pre$bacteria=rownames(tmp.CD.pre)
tmp.UC.pre$bacteria=rownames(tmp.UC.pre)
tmp.IBD.pre$bacteria=rownames(tmp.IBD.pre)
tmp.HC.pre$bacteria=rownames(tmp.HC.pre)
prevalence=rbind(tmp.CD.pre,tmp.UC.pre,tmp.IBD.pre,tmp.HC.pre)
prevalence$Group=factor(prevalence$Group,levels = c("HC","IBD","CD","UC"))

rownames(tmp.CD.pre)=tmp.CD.pre$bacteria
tmp.CD.pre$bacteria=NULL
tmp.CD.pre$Group=NULL
colnames(tmp.CD.pre)[1]="CD"

rownames(tmp.UC.pre)=tmp.UC.pre$bacteria
tmp.UC.pre$bacteria=NULL
tmp.UC.pre$Group=NULL
colnames(tmp.UC.pre)[1]="UC"

rownames(tmp.IBD.pre)=tmp.IBD.pre$bacteria
tmp.IBD.pre$bacteria=NULL
tmp.IBD.pre$Group=NULL
colnames(tmp.IBD.pre)[1]="IBD"

rownames(tmp.HC.pre)=tmp.HC.pre$bacteria
tmp.HC.pre$bacteria=NULL
tmp.HC.pre$Group=NULL
colnames(tmp.HC.pre)[1]="HC"

prevalence.tmp=cbind(tmp.CD.pre,tmp.UC.pre,tmp.IBD.pre,tmp.HC.pre)
prevalence.dist <- dist(prevalence.tmp)
prevalence.hclust <- hclust(prevalence.dist, method = "complete")
pdf("clust.pdf")
plot(prevalence.hclust, labels = FALSE)
dev.off()

sample_order=rownames(prevalence.tmp)[prevalence.hclust$order]
prevalence$bacteria=factor(prevalence$bacteria,levels= as.factor(sample_order))
ggplot(prevalence, aes(x=Group, y=bacteria, fill=Pre)) + 
  geom_tile(colour="white") + 
  scale_fill_gradient2(low="#456BB3", high = "#F26A55", mid = "grey", midpoint = 0.5) + 
  theme_classic() + theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = "", y= "") + 
  coord_fixed(ratio=1)
ggsave("OutputPlot/bacteria.absence.diagnosis.all.pdf",width = 9,height = 6)

tmp.data1=tmp.data
tmp.data1=tmp.data1[tmp.data1$Diagnosis!="Control",]
tmp.data1$Diagnosis="IBD"
tmp.data=rbind(tmp.data,tmp.data1)
tmp.data$Diagnosis=(factor(tmp.data$Diagnosis,levels = c("Control","IBD","CD","UC")))
tmp.data$Streptococcus=as.factor(tmp.data$Streptococcus)
ggplot(tmp.data,aes(x = (Diagnosis), group = (Streptococcus), fill = Streptococcus)) + 
  geom_bar(position = "fill")+theme_bw()+scale_fill_aaas()
ggsave("OutputPlot/Streptococcus..diagnosis.combined.pdf",width = 6,height = 4)


# ===========================================================================================
# bacteria data correction
# ===========================================================================================

basic_factors=covariate_bac[,c("age_at_biopsy","sex","BMI","smoking_DB","ResearchID","Batch","Inflammation","Location_rough")]
basic_factors$BMI[is.na(basic_factors$BMI)]=median(basic_factors$BMI[!is.na(basic_factors$BMI)])
basic_factors$Inflammation[is.na(basic_factors$Inflammation)]=median(basic_factors$Inflammation[!is.na(basic_factors$Inflammation)])
basic_factors$smoking_DB[is.na(basic_factors$smoking_DB)]=median(basic_factors$smoking_DB[!is.na(basic_factors$smoking_DB)])
bacteria.clr=bacteria.clr[rownames(bacteria.clr) %in% rownames(basic_factors),]
bacteria.clr=bacteria.clr[order(rownames(bacteria.clr)),]
basic_factors=basic_factors[order(rownames(basic_factors)),]
basic_factors$ResearchID=as.factor(basic_factors$ResearchID)
basic_factors$Batch=as.factor(basic_factors$Batch)

bacteria.clr[bacteria.clr==0]=NA
bacteria.clr.correct = apply(bacteria.clr,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="bacteria"
  
  fit=lmer(bacteria~age_at_biopsy+sex+BMI+smoking_DB+Batch+Inflammation+Location_rough+(1|ResearchID),data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

bacteria.binary=bacteria
bacteria.binary[bacteria.binary>0]=1
bacteria.binary=bacteria.binary[rownames(bacteria.binary) %in% rownames(basic_factors),]
bacteria.binary=bacteria.binary[order(rownames(bacteria.binary)),]
basic_factors=basic_factors[order(rownames(basic_factors)),]

bacteria.binary.correct = apply(bacteria.binary,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="bacteria"
  subset.data$age_at_biopsy=scale(subset.data$age_at_biopsy)
  subset.data$BMI=scale(subset.data$BMI)
  
  fit1=glmer(bacteria~age_at_biopsy+sex+BMI+smoking_DB+(1|ResearchID),data=subset.data,family = binomial)
  x.resid = resid(fit1)
  subset.data$bacteria=x.resid
  
  fit<-lm(bacteria~Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

plot(bacteria.clr[,53],bacteria.clr.correct[,53])
plot(bacteria.binary[,23],bacteria.binary.correct[,23])
cor.test(bacteria.clr[,53],bacteria.clr.correct[,53])

output.bacteria.clr=bacteria.clr.correct[rownames(bacteria.clr.correct) %in% rownames(covariate_bac)[covariate_bac$Cohort!="1000IBD_surgery"],]
output.bacteria.binary=bacteria.binary.correct[rownames(bacteria.binary.correct) %in% rownames(covariate_bac)[covariate_bac$Cohort!="1000IBD_surgery"],]

write.table(output.bacteria.clr,"OutputTable/bacteria.clr.correction.txt",sep = "\t",row.names = T,quote = F)
write.table(output.bacteria.binary,"OutputTable/bacteria.binary.correction.txt",sep = "\t",row.names = T,quote = F)

output.bacteria.clr=read.table("OutputTable/bacteria.clr.correction.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = F)
output.bacteria.binary=read.table("OutputTable/bacteria.binary.correction.txt",row.names = 1,sep = "\t",check.names = F,stringsAsFactors = )

output.bacteria.clr=as.data.frame(t(output.bacteria.clr))
output.bacteria.binary=as.data.frame(t(output.bacteria.binary))

corrected_data.clr = cbind(rownames(output.bacteria.clr),output.bacteria.clr)
colnames(corrected_data.clr)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data.clr),
                   Gene = rownames(corrected_data.clr),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data.clr, file = "OutputTable/16S_numeric.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "OutputTable/16S.annot",sep="\t",row.names=F,quote = F)

corrected_data.binary = cbind(rownames(output.bacteria.binary),output.bacteria.binary)
colnames(corrected_data.binary)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(corrected_data.binary),
                   Gene = rownames(corrected_data.binary),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"

write.table(corrected_data.binary, file = "OutputTable/16S_binary.txt",sep="\t",row.names = F,quote = F)

bacteria=read.table("Decontamination_dada2/Input.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
bacteria[bacteria>0]=1
bacteria=bacteria[,colSums(bacteria)>0.05*nrow(bacteria)]
write.table(bacteria, file = "OutputTable/16S_binary.txt",sep="\t",row.names = T,quote = F)

bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
write.table(bacteria, file = "OutputTable/16S_numeric.txt",sep="\t",row.names = T,quote = F)

# ===========================================================================================
# gene expression correction
# ===========================================================================================

covariate_rna=read.table("Covariate.rna.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genes=read.table("OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)

basic_factors=covariate_rna[,c("Diagnosis","age_at_biopsy","sex","BMI","ResearchID","Batch","Inflammation","Location_rough")]
basic_factors$Inflammation[basic_factors$Inflammation=="Yes"]=3
basic_factors$Inflammation[basic_factors$Inflammation=="No" & basic_factors$Diagnosis=="UC"]=2
basic_factors$Inflammation[basic_factors$Inflammation=="No" & basic_factors$Diagnosis=="CD"]=2
basic_factors$Inflammation[basic_factors$Inflammation=="No" & basic_factors$Diagnosis=="Control"]=1
basic_factors$Inflammation=as.numeric(basic_factors$Inflammation)
basic_factors$BMI[is.na(basic_factors$BMI)]=median(basic_factors$BMI[!is.na(basic_factors$BMI)])
basic_factors$Inflammation[is.na(basic_factors$Inflammation)]=median(basic_factors$Inflammation[!is.na(basic_factors$Inflammation)])
basic_factors$ResearchID=NULL
basic_factors$Location_rough=as.factor(basic_factors$Location_rough)
basic_factors$Batch[basic_factors$Batch=="batch1"]=1
basic_factors$Batch[basic_factors$Batch=="batch2"]=2
basic_factors$Batch=as.numeric(basic_factors$Batch)

basic_factors_CD=basic_factors[basic_factors$Diagnosis=="CD",]
basic_factors_UC=basic_factors[basic_factors$Diagnosis=="UC",]
basic_factors_Control=basic_factors[basic_factors$Diagnosis=="Control",]

# overall correction
genes_all=genes[rownames(genes) %in% rownames(basic_factors),]
basic_factors_all=basic_factors[rownames(basic_factors) %in% rownames(genes_all),]
genes_all=genes_all[order(rownames(genes_all)),]
basic_factors_all=basic_factors_all[order(rownames(basic_factors_all)),]

genes_all[genes_all==0]=NA
genes_all.correct = apply(genes_all,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors_all[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="gene"
  
  fit=lm(gene~age_at_biopsy+sex+BMI+Inflammation+Location_rough+Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

write.table(genes_all.correct,file = "RNAseq.CLR.all.correction.txt",sep = "\t",quote = F,row.names = T)

# CD correction
genes_CD=genes[rownames(genes) %in% rownames(basic_factors_CD),]
basic_factors_CD=basic_factors_CD[rownames(basic_factors_CD) %in% rownames(genes_CD),]
genes_CD=genes_CD[order(rownames(genes_CD)),]
basic_factors_CD=basic_factors_CD[order(rownames(basic_factors_CD)),]

genes_CD[genes_CD==0]=NA
genes_CD.correct = apply(genes_CD,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors_CD[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="gene"
  
  fit=lm(gene~age_at_biopsy+sex+BMI+Inflammation+Location_rough+Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

write.table(genes_CD.correct,file = "RNAseq.CLR.CD.correction.txt",sep = "\t",quote = F,row.names = T)

# UC correction

genes_UC=genes[rownames(genes) %in% rownames(basic_factors_UC),]
basic_factors_UC=basic_factors_UC[rownames(basic_factors_UC) %in% rownames(genes_UC),]
genes_UC=genes_UC[order(rownames(genes_UC)),]
basic_factors_UC=basic_factors_UC[order(rownames(basic_factors_UC)),]

genes_UC[genes_UC==0]=NA
genes_UC.correct = apply(genes_UC,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors_UC[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="gene"
  
  fit=lm(gene~age_at_biopsy+sex+BMI+Inflammation+Location_rough+Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

write.table(genes_UC.correct,file = "RNAseq.CLR.UC.correction.txt",sep = "\t",quote = F,row.names = T)

# Control correction

genes_Control=genes[rownames(genes) %in% rownames(basic_factors_Control),]
basic_factors_Control=basic_factors_Control[rownames(basic_factors_Control) %in% rownames(genes_Control),]
genes_Control=genes_Control[order(rownames(genes_Control)),]
basic_factors_Control=basic_factors_Control[order(rownames(basic_factors_Control)),]

genes_Control[genes_Control==0]=NA
genes_Control.correct = apply(genes_Control,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = basic_factors_Control[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="gene"
  
  fit=lm(gene~age_at_biopsy+sex+BMI+Location_rough+Batch,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})

write.table(genes_Control.correct,file = "RNAseq.CLR.Control.correction.txt",sep = "\t",quote = F,row.names = T)



# ===========================================================================================
# host-gene interaction
# ===========================================================================================

genes=read.table("RNAseq.correction.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
bacteria1=read.table("bacteria.binary.correction.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)
bacteria2=read.table("bacteria.clr.correction.txt",stringsAsFactors = F,header = T,sep = "\t",row.names = 1)

result1=matrix(nrow = 1,ncol = 4)
result1=as.data.frame(result1)
colnames(result1)=c("Bacteria","Genes","Pvalue","Estimate")
for (i in 1:ncol(bacteria1)){

    tmp.genus=colnames(bacteria1)[i]
    tmp.data=(bacteria1[,i,drop=F])
    tmp.data=na.omit(tmp.data)
    tmp.gene=genes[rownames(genes) %in% rownames(tmp.data),]
    tmp.data=tmp.data[rownames(tmp.data) %in% rownames(tmp.gene),,drop=F]
    tmp.gene=tmp.gene[order(rownames(tmp.gene)),]
    tmp.data=tmp.data[order(rownames(tmp.data)),,drop=F]
    
    tmp.result=matrix(nrow = ncol(tmp.gene),ncol = 4)
    tmp.result=as.data.frame(tmp.result)
    colnames(tmp.result)=c("Bacteria","Genes","Pvalue","Estimate")
    tmp.result$Bacteria=tmp.genus
    for(n in 1:ncol(tmp.gene)){
      tmp.result$Genes[n]=colnames(tmp.gene)[n]
      mm=cor.test(tmp.data[,1],tmp.gene[,n])
      tmp.result$Pvalue[n]=mm$p.value
      tmp.result$Estimate[n]=mm$estimate
      
      cat(tmp.genus,"-------",colnames(tmp.gene)[n],"is done","\n")
    }
    
    result1=rbind(result1,tmp.result)
    cat(tmp.genus,"DONE=======","\n")
}
result1=result1[-1,]
result1$FDR=p.adjust(result1$Pvalue)

write.table(result1,file = "Binary.correlation.txt",sep = "\t",row.names = F,quote = F)

result2=matrix(nrow = 1,ncol = 4)
result2=as.data.frame(result2)
colnames(result2)=c("Bacteria","Genes","Pvalue","Estimate")
for (i in 1:ncol(bacteria2)){
  
  tmp.genus=colnames(bacteria2)[i]
  tmp.data=(bacteria2[,i,drop=F])
  tmp.data=na.omit(tmp.data)
  tmp.gene=genes[rownames(genes) %in% rownames(tmp.data),]
  tmp.data=tmp.data[rownames(tmp.data) %in% rownames(tmp.gene),,drop=F]
  tmp.gene=tmp.gene[order(rownames(tmp.gene)),]
  tmp.data=tmp.data[order(rownames(tmp.data)),,drop=F]
  
  tmp.result=matrix(nrow = ncol(tmp.gene),ncol = 4)
  tmp.result=as.data.frame(tmp.result)
  colnames(tmp.result)=c("Bacteria","Genes","Pvalue","Estimate")
  tmp.result$Bacteria=tmp.genus
  for(n in 1:ncol(tmp.gene)){
    tmp.result$Genes[n]=colnames(tmp.gene)[n]
    mm=cor.test(tmp.data[,1],tmp.gene[,n])
    tmp.result$Pvalue[n]=mm$p.value
    tmp.result$Estimate[n]=mm$estimate
    
    cat(tmp.genus,"-------",colnames(tmp.gene)[n],"is done","\n")
  }
  
  result2=rbind(result2,tmp.result)
  cat(tmp.genus,"DONE=======","\n")
}
result2=result2[-1,]
result2$FDR=p.adjust(result2$Pvalue)

write.table(result2,file = "Quantitative.correlation.txt",sep = "\t",row.names = F,quote = F)

# ===========================================================================================
# host-gene interaction, target approach (only inflammed-related genes)
# ===========================================================================================
library(reshape2)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation_genes=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
inflammation=inflammation[inflammation$Gene %in% inflammation_genes,]
inflammation$UpDown=NA
inflammation$UpDown[inflammation$inflammation.beta>0]="Up"
inflammation$UpDown[inflammation$inflammation.beta<0]="Down"

inflammation=inflammation[!duplicated(inflammation$Gene),]
annotation=inflammation[,c("Gene","UpDown")]
rownames(annotation)=annotation$Gene
annotation$Gene=NULL

associations_CD=read.table("OutputTable/Quantitative.correlation.CD.txt",sep = "\t",header = T,stringsAsFactors = F)
associations_CD=associations_CD[associations_CD$Genes %in% inflammation_genes,]
associations_CD=associations_CD[associations_CD$Bacteria %in% colnames(genus),]
associations_CD$FDR=p.adjust(associations_CD$Pvalue,method = "BH")

associations_UC=read.table("OutputTable/Quantitative.correlation.UC.txt",sep = "\t",header = T,stringsAsFactors = F)
associations_UC=associations_UC[associations_UC$Genes %in% inflammation_genes,]
associations_UC=associations_UC[associations_UC$Bacteria %in% colnames(genus),]
associations_UC$FDR=p.adjust(associations_UC$Pvalue,method = "BH")

associations_Control=read.table("OutputTable/Quantitative.correlation.Control.txt",sep = "\t",header = T,stringsAsFactors = F)
associations_Control=associations_Control[associations_Control$Genes %in% inflammation_genes,]
associations_Control=associations_Control[associations_Control$Bacteria %in% colnames(genus),]
associations_Control$FDR=p.adjust(associations_Control$Pvalue,method = "BH")

overlap_gene=intersect(associations_CD$Genes[associations_CD$FDR<0.1],intersect(associations_UC$Genes[associations_UC$FDR<0.1],
                                                                                associations_Control$Genes[associations_Control$FDR<0.1]))
overlap_bacteria=intersect(associations_CD$Bacteria[associations_CD$FDR<0.1],intersect(associations_UC$Bacteria[associations_UC$FDR<0.1],
                                                                                associations_Control$Bacteria[associations_Control$FDR<0.1]))

associations_Control$Disease="Control"
associations_CD$Disease="CD"
associations_UC$Disease="UC"
edge=rbind(associations_CD[associations_CD$Bacteria==overlap_bacteria[1] & associations_CD$FDR<0.1,],
           associations_UC[associations_UC$Bacteria==overlap_bacteria[1] & associations_UC$FDR<0.1,],
           associations_Control[associations_Control$Bacteria==overlap_bacteria[1] & associations_Control$FDR<0.1,])
edge_list=edge[,c("Bacteria","Genes","Disease")]
colnames(edge_list)[1:2]=c("from","to")
node_list=data.frame(id=unique(c(edge_list$from,edge_list$to)),type=NA)
node_list$Type[node_list$id==overlap_bacteria[1]]="Bacteria"
node_list$Type[node_list$id!=overlap_bacteria[1]]="Genes"

library(igraph)
net = graph_from_data_frame(edge_list,vertices = node_list,directed = T)
plot(net)

colrs.v = c(Bacteria = "lightblue",Genes = "gold") #node colours
V(net)$color = colrs.v[V(net)$Type]

colrs.e = c(CD = "red", UC = "blue",Control="green") #edge colours
E(net)$color = colrs.e[E(net)$Disease] 

# meta analysis, CD, UC, control
library(metap)
genus=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
genus_CD=genus[rownames(genus) %in% rownames(basic_factors_CD),]
genus_CD=genus_CD[,colSums(genus_CD!=0)>(0.1*nrow(genus_CD))]
genus_UC=genus[rownames(genus) %in% rownames(basic_factors_UC),]
genus_UC=genus_UC[,colSums(genus_UC!=0)>(0.1*nrow(genus_UC))]
genus_Control=genus[rownames(genus) %in% rownames(basic_factors_Control),]
genus_Control=genus_Control[,colSums(genus_Control!=0)>(0.1*nrow(genus_Control))]
select=c(intersect(colnames(genus_Control),intersect(colnames(genus_CD),colnames(genus_UC))))
associations_CD=associations_CD[associations_CD$Bacteria %in% select,]
associations_UC=associations_UC[associations_UC$Bacteria %in% select,]
associations_Control=associations_Control[associations_Control$Bacteria %in% select,]

nrow(associations_Control[associations_Control$FDR<0.1,])
associations_CD_sig=associations_CD[associations_CD$Pvalue<0.05,]
associations_UC_sig=associations_UC[associations_UC$Pvalue<0.05,]
associations_Control_sig=associations_Control[associations_Control$Pvalue<0.05,]
associations_CD_sig$pair=paste(associations_CD_sig$Bacteria,associations_CD_sig$Genes,sep="_")
associations_UC_sig$pair=paste(associations_UC_sig$Bacteria,associations_UC_sig$Genes,sep="_")
associations_Control_sig$pair=paste(associations_Control_sig$Bacteria,associations_Control_sig$Genes,sep="_")

associations_CD$pair=paste(associations_CD$Bacteria,associations_CD$Genes,sep="_")
associations_UC$pair=paste(associations_UC$Bacteria,associations_UC$Genes,sep="_")
associations_Control$pair=paste(associations_Control$Bacteria,associations_Control$Genes,sep="_")
result_meta_association = foreach(i=c(1:nrow(associations_Control)),.combine = rbind) %do% {
  tmp.gene=associations_Control$pair[i]
  cat(i,"++++++",sep = "\n")
  
  p_vector=c(associations_CD$Pvalue[associations_CD$pair==tmp.gene],
             associations_UC$Pvalue[associations_UC$pair==tmp.gene],
             associations_Control$Pvalue[associations_Control$pair==tmp.gene])
  w_vector=c(nrow(genes_CD.correct),nrow(genes_UC.correct),nrow(genes_Control.correct))
  
  metaP=sumz(p_vector,weights = w_vector)
  metaP=metaP$p
  
  return.string=data.frame(Pair=tmp.gene,
                           CD.beta=associations_CD$Estimate[associations_CD$pair==tmp.gene],
                           UC.beta=associations_UC$Estimate[associations_UC$pair==tmp.gene],
                           Control.beta=associations_Control$Estimate[associations_Control$pair==tmp.gene],
                           CD.pvalue=associations_CD$Pvalue[associations_CD$pair==tmp.gene],
                           UC.pvalue=associations_UC$Pvalue[associations_UC$pair==tmp.gene],
                           Control.pvalue=associations_Control$Pvalue[associations_Control$pair==tmp.gene],
                           metaP=metaP)
  
}
aa=data.frame(str_split_fixed(result_meta_association$Pair, "_", 2))
result_meta_association$Bacteria=aa$X1
result_meta_association$Genes=aa$X2
result_meta_association$FDR=p.adjust(result_meta_association$metaP,method = "BH")
result_meta_association1=result_meta_association[result_meta_association$CD.beta<0 & result_meta_association$UC.beta<0 & result_meta_association$Control.beta<0,]
result_meta_association2=result_meta_association[result_meta_association$CD.beta>0 & result_meta_association$UC.beta>0 & result_meta_association$Control.beta>0,]
result_meta_association_clean=rbind(result_meta_association1,result_meta_association2)
result_meta_association_clean=result_meta_association_clean[order(result_meta_association_clean$FDR),]

associations_cor=as.data.frame(t(acast(associations, Bacteria~Genes, value.var="Estimate")))
associations_sig=as.data.frame(t(acast(associations, Bacteria~Genes, value.var="FDR")))
associations_cor=associations_cor[,colnames(associations_cor) %in% colnames(associations_sig)]
associations_cor=associations_cor[,colnames(associations_cor) %in% colnames(genus)]
associations_sig=associations_sig[,colnames(associations_sig) %in% colnames(genus)]

associations_cor=result_meta_association_clean[1:20,c("Control.beta","CD.beta","UC.beta","Pair")]
rownames(associations_cor)=associations_cor$Pair
associations_cor$Pair=NULL
associations_sig=result_meta_association_clean[1:20,c("Control.pvalue","CD.pvalue","UC.pvalue","Pair")]
rownames(associations_sig)=associations_sig$Pair
associations_sig$Pair=NULL

library(pheatmap)
heatmap.data=(associations_cor)
heatmap.data.lable=(associations_sig)
heatmap.data.lable[heatmap.data.lable<0.05]="*"
heatmap.data.lable[heatmap.data.lable>0.05]=""
my_colour = list(
  UpDown = c(Up = "#5AAE61", Down = "#F67E4B")
)

pdf("OutputPlot//host-microbe.association.pdf",width = 5,height = 8)
pheatmap(heatmap.data,cluster_cols = T, cluster_rows = T,
         show_rownames=T, show_colnames=T,
         cellheight = 20,cellwidth = 20,fontsize_number=10,fontsize_row=10,fontsize_col = 10,
         #annotation_col  =annotation,
         #annotation_colors = my_colour,
         display_numbers = heatmap.data.lable,color = colorRampPalette(c("#BB5566",mid="white", "#0077BB"))(100))

dev.off()

# lineaer plot
RNAseq1=read.table("RNAseq.CLR.Control.correction.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
RNAseq2=read.table("RNAseq.CLR.CD.correction.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
RNAseq3=read.table("RNAseq.CLR.UC.correction.txt",row.names = 1,header = T,stringsAsFactors = F,sep = "\t")
genus=read.table("OutputTable/CLR.genus.txt",sep="\t",row.names = 1,header = T)

# Alistipes_CHST13
bac="Alistipes"
bac=genus[,bac,drop=F]
geneset="CHST13"
geneset=RNAseq1[,colnames(RNAseq1) %in% geneset,drop=F]

tmp.data=merge(geneset,bac,by="row.names",all=F)
ggplot(tmp.data, aes(Alistipes, CHST13,fill="#BB5566")) +
  geom_point(shape = 21,color = "black", size = 3)+
  geom_smooth(method = lm,color="black")+
  theme_bw()+ theme(legend.position="bottom")+guides(fill=F)
ggsave("OutputPlot/Alistipes.Control.pdf",width = 3,height = 3)

# ===========================================================================================
# alpha diversity and genes
# ===========================================================================================
library(reshape2)
bacteria=read.table("Decontamination_dada2/Input.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F)
bacteria=bacteria[,colSums(bacteria>0)>0.05*nrow(bacteria)]
shannon=as.data.frame(diversity(t(bacteria), index="shannon"))
colnames(shannon)="shannon"

genes_CD.correct=read.table("RNAseq.CLR.CD.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genes_UC.correct=read.table("RNAseq.CLR.UC.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genes_Control.correct=read.table("RNAseq.CLR.Control.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genes_all.correct=read.table("RNAseq.CLR.all.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
genes_CD.correct=genes_CD.correct[,colnames(genes_CD.correct) %in% inflammation]
genes_UC.correct=genes_UC.correct[,colnames(genes_UC.correct) %in% inflammation]
genes_Control.correct=genes_Control.correct[,colnames(genes_Control.correct) %in% inflammation]
genes_all.correct=genes_all.correct[,colnames(genes_all.correct) %in% inflammation]

# all
result0=matrix(ncol = 3,nrow = ncol(genes_all.correct))
result0=as.data.frame(result0)
colnames(result0)=c("Gene","Estimate","Pvalue")
for(i in 1:ncol(genes_all.correct)){
  tmp.gene=colnames(genes_all.correct)[i]
  tmp.data=merge(shannon,genes_all.correct[,tmp.gene,drop=F],all=F,by="row.names")
  mm=cor.test(tmp.data[,2],tmp.data[,3])
  result0$Gene[i]=tmp.gene
  result0$Estimate[i]=mm$estimate
  result0$Pvalue[i]=mm$p.value
  
  cat(i,"+++++","\n")
}
result0$FDR=p.adjust(result0$Pvalue,method = "BH")

# CD
result1=matrix(ncol = 3,nrow = ncol(genes_CD.correct))
result1=as.data.frame(result1)
colnames(result1)=c("Gene","Estimate","Pvalue")
for(i in 1:ncol(genes_CD.correct)){
  tmp.gene=colnames(genes_CD.correct)[i]
  tmp.data=merge(shannon,genes_CD.correct[,tmp.gene,drop=F],all=F,by="row.names")
  mm=lm(tmp.data[,2]~tmp.data[,3],data = tmp.data)
  
  tmp.coef=as.data.frame(summary(mm)$coef)
  tmp.coef$CIup=tmp.coef$Estimate[2]+1.96*tmp.coef$`Std. Error`[2]
  tmp.coef$CIdown=tmp.coef$Estimate[2]-1.96*tmp.coef$`Std. Error`[2]
  
  result1$Gene[i]=tmp.gene
  result1$Estimate[i]=tmp.coef$Estimate[2]
  result1$Pvalue[i]=tmp.coef$`Pr(>|t|)`[2]
  result1$CIup[i]=tmp.coef$CIup[2]
  result1$CIdown[i]=tmp.coef$CIdown[2]
  
  cat(i,"+++++","\n")
}
result1$FDR=p.adjust(result1$Pvalue,method = "BH")

# UC
result2=matrix(ncol = 3,nrow = ncol(genes_UC.correct))
result2=as.data.frame(result2)
colnames(result2)=c("Gene","Estimate","Pvalue")
for(i in 1:ncol(genes_UC.correct)){
  tmp.gene=colnames(genes_UC.correct)[i]
  tmp.data=merge(shannon,genes_UC.correct[,tmp.gene,drop=F],all=F,by="row.names")
  mm=lm(tmp.data[,2]~tmp.data[,3],data = tmp.data)
  
  tmp.coef=as.data.frame(summary(mm)$coef)
  tmp.coef$CIup=tmp.coef$Estimate[2]+1.96*tmp.coef$`Std. Error`[2]
  tmp.coef$CIdown=tmp.coef$Estimate[2]-1.96*tmp.coef$`Std. Error`[2]
  
  result2$Gene[i]=tmp.gene
  result2$Estimate[i]=tmp.coef$Estimate[2]
  result2$Pvalue[i]=tmp.coef$`Pr(>|t|)`[2]
  result2$CIup[i]=tmp.coef$CIup[2]
  result2$CIdown[i]=tmp.coef$CIdown[2]
  
  cat(i,"+++++","\n")
}
result2$FDR=p.adjust(result2$Pvalue,method = "BH")

# Control
result3=matrix(ncol = 3,nrow = ncol(genes_Control.correct))
result3=as.data.frame(result3)
colnames(result3)=c("Gene","Estimate","Pvalue")
for(i in 1:ncol(genes_Control.correct)){
  tmp.gene=colnames(genes_Control.correct)[i]
  tmp.data=merge(shannon,genes_Control.correct[,tmp.gene,drop=F],all=F,by="row.names")
  mm=lm(tmp.data[,2]~tmp.data[,3],data = tmp.data)
  
  tmp.coef=as.data.frame(summary(mm)$coef)
  tmp.coef$CIup=tmp.coef$Estimate[2]+1.96*tmp.coef$`Std. Error`[2]
  tmp.coef$CIdown=tmp.coef$Estimate[2]-1.96*tmp.coef$`Std. Error`[2]
  
  result3$Gene[i]=tmp.gene
  result3$Estimate[i]=tmp.coef$Estimate[2]
  result3$Pvalue[i]=tmp.coef$`Pr(>|t|)`[2]
  result3$CIup[i]=tmp.coef$CIup[2]
  result3$CIdown[i]=tmp.coef$CIdown[2]
  
  cat(i,"+++++","\n")
}
result3$FDR=p.adjust(result3$Pvalue,method = "BH")

# plot
geneset=c("REG4","DEFB4A","MT1F","MT1M","MT1G","MT1H","ACOT9","MUC2","LPCAT2","MT1E")
geneset=genes_Control.correct[,colnames(genes_Control.correct) %in% geneset,drop=F]
tmp.data=merge(shannon,geneset,by="row.names",all=F)
tmp.data=melt(tmp.data, id.vars=c("Row.names", "shannon"))
tmp.data$variable=factor(tmp.data$variable,levels = c("REG4","DEFB4A","MT1F","MT1M","MT1G","MT1H","ACOT9","MUC2","LPCAT2","MT1E"))
ggplot(tmp.data, aes(shannon, value,fill=variable,color="grey")) +
  geom_point(shape = 21, 
             color = "black", size = 2)+
  geom_smooth(method = lm)+
  facet_wrap( .~ variable, scales="free_y",nrow = 2)+scale_color_jama()+
  theme_bw()+ theme(legend.position="bottom")+guides(fill=F)
ggsave("OutputPlot/alaphadiversity.gene.pdf",width = 8,height = 5)

# meta CD, UC and control
library(metap)
result_meta = foreach(i=c(1:nrow(result0)),.combine = rbind) %do% {
  tmp.gene=result0$Gene[i]
  cat(i,"++++++",sep = "\n")
  
  p_vector=c(result1$Pvalue[i],result2$Pvalue[i],result3$Pvalue[i])
  w_vector=c(nrow(genes_CD.correct),nrow(genes_UC.correct),nrow(genes_Control.correct))
  
  metaP=sumz(p_vector,weights = w_vector)
  metaP=metaP$p
  
  return.string=data.frame(Gene=tmp.gene,CD.beta=result1$Estimate[i],CD.pvalue=result1$Pvalue[i],
                           UC.beta=result2$Estimate[i],UC.pvalue=result2$Pvalue[i],
                           Control.beta=result3$Estimate[i],Control.pvalue=result3$Pvalue[i],
                           metaP=metaP)
  
}
result_meta=result_meta[result_meta$CD.pvalue<0.05 & result_meta$UC.pvalue<0.05 & result_meta$Control.pvalue<0.05,]
result_meta$FDR=p.adjust(result_meta$metaP,method = 'BH')

result1$Disease="CD"
result2$Disease="UC"
result3$Disease="Control"

combine=rbind(result1[result1$Gene %in% c("ITGA6","MET"),],result2[result2$Gene %in% c("ITGA6","MET"),],result3[result3$Gene %in% c("ITGA6","MET"),])
combine$Disease=factor(combine$Disease,levels = c("Control","CD","UC"))
ggplot(combine, aes(x=Disease, y=Estimate, ymin=CIdown, ymax=CIup,col=Disease,fill=Disease)) + 
  geom_hline(yintercept=0, lty=2) +
  geom_point(size=3, shape=21, colour="white", stroke = 0.5,position=position_dodge(width = 0.5)) +
  geom_pointrange(aes(col=Disease),position=position_dodge(width = 0.5))+
  coord_flip() +scale_color_jama()+scale_fill_jama()+theme_bw()+facet_wrap(~Gene,  scales = "free_x", drop = TRUE)+guides(color=F)+guides(fill=F)
ggsave("OutputPlot/meta.diversity.pdf",width = 3,height = 3)




















