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
kruskal.test(shannon ~ Diagnosis, data = shannon)
wilcox.test(shannon$shannon[shannon$Split=="CD ileum"],shannon$shannon[shannon$Split=="Control colon"])
wilcox.test(shannon$shannon[shannon$Diagnosis=="UC"],shannon$shannon[shannon$Diagnosis=="Control"])
wilcox.test(shannon$shannon[shannon$Diagnosis=="CD"],shannon$shannon[shannon$Diagnosis=="UC"])

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
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="CD"])
wilcox.test(distmat$Freq[distmat$group=="Control"],distmat$Freq[distmat$group=="UC"])
wilcox.test(distmat$Freq[distmat$group=="CD"],distmat$Freq[distmat$group=="UC"])

# ===========================================================================================
# bacteria surgery within disease
# ===========================================================================================

covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]

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

# ===========================================================================================
# bacteria abundance characterization, inflammation/disease
# ===========================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)

# CD vs. HC at ileum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD" | covariate_bac$Diagnosis=="Control",]
covariate1=covariate1[covariate1$Location_rough=="ileum",]

result1 = foreach(i=1:ncol(bacteria1),.combine = rbind) %do%  {
  tmp.genus=colnames(bacteria1)[i]
  x=bacteria1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>100){
    fit<-glmer(tmp.data[,tmp.genus]~Diagnosis + (1|ResearchID),data=tmp.data)
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
    fit<-glmer(tmp.data[,tmp.genus]~Diagnosis +(1|ResearchID),data=tmp.data)
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
    fit<-glmer(tmp.data[,tmp.genus]~Diagnosis +(1|ResearchID),data=tmp.data)
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

# ===========================================================================================
# bacteria abundance characterization, location
# ===========================================================================================
bacteria=read.table("OutputTable/CLR.bacteria.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)

# CD, colon vs ileum
covariate1=covariate_bac[covariate_bac$Diagnosis=="CD",]

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
    fit<-glmer(tmp.data[,tmp.genus]~Location_rough+age_at_biopsy+sex+BMI+(1|ResearchID),data=tmp.data)
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

# ===========================================================================================
# bacteria abundance characterization, diagnosis
# ===========================================================================================

# CD and UC vs. HC
result0 = foreach(i=1:ncol(genus1),.combine = rbind) %do%  {
  tmp.genus=colnames(genus1)[i]
  x=genus1[,i,drop=F]
  x[x==0]=NA
  x=na.omit(x)
  tmp.data=merge(x,covariate1,by="row.names",all=F)
  tmp.data=na.omit(tmp.data)
  tmp.data$Row.names=NULL
  
  if(nrow(tmp.data)>30){
    fit<-glmer(tmp.data[,tmp.genus]~Inflammation+age_at_biopsy+sex+BMI+smoking_DB+Batch+Location_rough+(1|ResearchID),data=tmp.data)
    residuals=as.data.frame(resid(fit))
    mm1=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="CD"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
    mm2=wilcox.test(residuals$`resid(fit)`[tmp.data$Diagnosis=="UC"],residuals$`resid(fit)`[tmp.data$Diagnosis=="Control"])
    nn=wilcox.test(residuals$`resid(fit)`[tmp.data$Cohort==1],residuals$`resid(fit)`[tmp.data$Cohort==0])
    cat(tmp.genus,"is done","\n")
    return.string=data.frame(bacteria=tmp.genus,diagnosis.CD=mm1$p.value,diagnosis.UC=mm2$p.value,cohort.p=nn$p.value)
  }else{return.string=data.frame(bacteria=tmp.genus,
                                 diagnosis.CD=NA,diagnosis.UC=NA,cohort.p=NA)}
}


