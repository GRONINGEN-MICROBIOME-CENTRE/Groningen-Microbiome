# ===============================================================================================
# this script is to analyze HMP2 16S from dada2 output
# author: Shixian Hu
# this script inludes count reads, low-quality sample removal, rare ASV removeal, decontamination
# low-quality sample : reads <2000
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
library("ape")
library(dendextend)

# ===============================================================================================
# dada2 statistic report
# ===============================================================================================

dada2_run=read.table("HMP2/dada2.read.statistics.txt",sep = "\t",row.names = 1,check.names = F)

# ===============================================================================================
# aglined to kingdom reads count
# ===============================================================================================

dada2_align_run=read.table("HMP2/Kingdom.level.txt",sep = "\t",header = T,check.names = F)
dada2_align_run=as.data.frame(t(dada2_align_run))
dada2_align_run=dada2_align_run[,"Bacteria",drop=F]
dada2_align_run=cbind(dada2_align_run,dada2_run)
colnames(dada2_align_run)[1]="AlignedReads"
dada2_align_run$Unaligned=dada2_align_run$nonchim-dada2_align_run$AlignedReads
dada2_align_run$ratio=dada2_align_run$AlignedReads/(dada2_align_run$AlignedReads+dada2_align_run$Unaligned)
dada2_align_run$SampleID=rownames(dada2_align_run)
dada2_align_run=dada2_align_run[order(dada2_align_run$AlignedReads,decreasing = F),]
dada2_align_run$SampleID=factor(dada2_align_run$SampleID,levels = dada2_align_run$SampleID)

ggplot(data=dada2_align_run) +
  geom_bar(mapping=aes(x=SampleID, y=AlignedReads),stat="identity")+
  geom_line(mapping=aes(x=SampleID, y=ratio*100000),size=0.5,group = 1)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1,size = 3))+
  scale_y_continuous(sec.axis = sec_axis(trans= ~./100000)) 

# ===============================================================================================
# before decontamination, remove samples, by aligned reads <2,000
# ===============================================================================================

# by mapping rate
dada2_align_run$SampleID=as.character(dada2_align_run$SampleID)
poor_samples=unique(as.character(dada2_align_run$SampleID[dada2_align_run$AlignedReads<2000]))

# ===============================================================================================
# import phenotype
# ===============================================================================================

phenotype=read.table("HMP2/Coupling.16S.txt",sep = "\t",fill = T,stringsAsFactors = F,header = T,check.names = F)
phenotype$Location=NA
phenotype$Location[phenotype$biopsy_location!="Ileum"]="Colon"
phenotype$Location[phenotype$biopsy_location=="Ileum"]="Ileum"

# ===============================================================================================
# import taxa and CLR
# ===============================================================================================
genus=read.table("HMP2/Genus.level.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
genus=as.data.frame(t(genus))
genus=genus[,colSums(genus)>0]
genus=genus[rowSums(genus)>0,]
genus=apply(genus,1,function(x){
  x=x/sum(x)
})
genus=as.data.frame(t(genus))

genus_clr <- zCompositions::cmultRepl(genus, method="CZM", label=0)
genus_clr = compositions::clr(genus_clr)
genus_clr=as.data.frame(genus_clr)

phylum=read.table("HMP2/phylum.level.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
phylum=as.data.frame(t(phylum))
phylum=phylum[,colSums(phylum)>0]
phylum=phylum[rowSums(phylum)>0,]
phylum=apply(phylum,1,function(x){
  x=x/sum(x)
})
phylum=as.data.frame(t(phylum))

phylum_clr <- zCompositions::cmultRepl(phylum, method="CZM", label=0)
phylum_clr = compositions::clr(phylum_clr)
phylum_clr=as.data.frame(phylum_clr)

class=read.table("HMP2/class.level.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
class=as.data.frame(t(class))
class=class[,colSums(class)>0]
class=class[rowSums(class)>0,]
class=apply(class,1,function(x){
  x=x/sum(x)
})
class=as.data.frame(t(class))

class_clr <- zCompositions::cmultRepl(class, method="CZM", label=0)
class_clr = compositions::clr(class_clr)
class_clr=as.data.frame(class_clr)

order=read.table("HMP2/order.level.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
order=as.data.frame(t(order))
order=order[,colSums(order)>0]
order=order[rowSums(order)>0,]
order=apply(order,1,function(x){
  x=x/sum(x)
})
order=as.data.frame(t(order))

order_clr <- zCompositions::cmultRepl(order, method="CZM", label=0)
order_clr = compositions::clr(order_clr)
order_clr=as.data.frame(order_clr)

family=read.table("HMP2/family.level.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F)
family=as.data.frame(t(family))
family=family[,colSums(family)>0]
family=family[rowSums(family)>0,]
family=apply(family,1,function(x){
  x=x/sum(x)
})
family=as.data.frame(t(family))

family_clr <- zCompositions::cmultRepl(family, method="CZM", label=0)
family_clr = compositions::clr(family_clr)
family_clr=as.data.frame(family_clr)

# ===============================================================================================
# harmpnize all data, filter 10%
# ===============================================================================================
#genus
genus=genus[,colnames(genus) %in% colnames(genus_clr)]
dada2_align_run=dada2_align_run[!dada2_align_run$SampleID %in% poor_samples,]
dada2_align_run$SampleID=gsub("SS","S",dada2_align_run$SampleID)
phenotype=phenotype[phenotype$Run %in% dada2_align_run$SampleID,]

rownames(genus)=gsub("SS","S",rownames(genus))
rownames(genus_clr)=gsub("SS","S",rownames(genus_clr))

genus=genus[rownames(genus) %in% phenotype$Run,]
genus_clr=genus_clr[rownames(genus_clr) %in% phenotype$Run,]

genus = genus[,colSums(genus>0)>nrow(genus) * 0.1]
genus_clr=genus_clr[,colnames(genus_clr) %in% colnames(genus)]

write.table(genus,file = "HMP2/Input.genus.txt",row.names = T,quote = F,sep = "\t")
write.table(genus_clr,file = "HMP2/Input.genus_CLR.txt",row.names = T,quote = F,sep = "\t")
write.table(phenotype,file = "HMP2/Input.phenotype.txt",row.names = T,quote = F,sep = "\t")

# phylum
phylum=phylum[,colnames(phylum) %in% colnames(phylum_clr)]
dada2_align_run=dada2_align_run[!dada2_align_run$SampleID %in% poor_samples,]
dada2_align_run$SampleID=gsub("SS","S",dada2_align_run$SampleID)
phenotype=phenotype[phenotype$Run %in% dada2_align_run$SampleID,]

rownames(phylum)=gsub("SS","S",rownames(phylum))
rownames(phylum_clr)=gsub("SS","S",rownames(phylum_clr))

phylum=phylum[rownames(phylum) %in% phenotype$Run,]
phylum_clr=phylum_clr[rownames(phylum_clr) %in% phenotype$Run,]

phylum = phylum[,colSums(phylum>0)>nrow(phylum) * 0.1]
phylum_clr=phylum_clr[,colnames(phylum_clr) %in% colnames(phylum)]

write.table(phylum,file = "HMP2/Input.phylum.txt",row.names = T,quote = F,sep = "\t")
write.table(phylum_clr,file = "HMP2/Input.phylum_CLR.txt",row.names = T,quote = F,sep = "\t")
write.table(phenotype,file = "HMP2/Input.phenotype.txt",row.names = T,quote = F,sep = "\t")

# family
family=family[,colnames(family) %in% colnames(family_clr)]
dada2_align_run=dada2_align_run[!dada2_align_run$SampleID %in% poor_samples,]
dada2_align_run$SampleID=gsub("SS","S",dada2_align_run$SampleID)
phenotype=phenotype[phenotype$Run %in% dada2_align_run$SampleID,]

rownames(family)=gsub("SS","S",rownames(family))
rownames(family_clr)=gsub("SS","S",rownames(family_clr))

family=family[rownames(family) %in% phenotype$Run,]
family_clr=family_clr[rownames(family_clr) %in% phenotype$Run,]

family = family[,colSums(family>0)>nrow(family) * 0.1]
family_clr=family_clr[,colnames(family_clr) %in% colnames(family)]

write.table(family,file = "HMP2/Input.family.txt",row.names = T,quote = F,sep = "\t")
write.table(family_clr,file = "HMP2/Input.family_CLR.txt",row.names = T,quote = F,sep = "\t")

# order
order=order[,colnames(order) %in% colnames(order_clr)]
dada2_align_run=dada2_align_run[!dada2_align_run$SampleID %in% poor_samples,]
dada2_align_run$SampleID=gsub("SS","S",dada2_align_run$SampleID)
phenotype=phenotype[phenotype$Run %in% dada2_align_run$SampleID,]

rownames(order)=gsub("SS","S",rownames(order))
rownames(order_clr)=gsub("SS","S",rownames(order_clr))

order=order[rownames(order) %in% phenotype$Run,]
order_clr=order_clr[rownames(order_clr) %in% phenotype$Run,]

order = order[,colSums(order>0)>nrow(order) * 0.1]
order_clr=order_clr[,colnames(order_clr) %in% colnames(order)]

write.table(order,file = "HMP2/Input.order.txt",row.names = T,quote = F,sep = "\t")
write.table(order_clr,file = "HMP2/Input.order_CLR.txt",row.names = T,quote = F,sep = "\t")

# class
class=class[,colnames(class) %in% colnames(class_clr)]
dada2_align_run=dada2_align_run[!dada2_align_run$SampleID %in% poor_samples,]
dada2_align_run$SampleID=gsub("SS","S",dada2_align_run$SampleID)
phenotype=phenotype[phenotype$Run %in% dada2_align_run$SampleID,]

rownames(class)=gsub("SS","S",rownames(class))
rownames(class_clr)=gsub("SS","S",rownames(class_clr))

class=class[rownames(class) %in% phenotype$Run,]
class_clr=class_clr[rownames(class_clr) %in% phenotype$Run,]

class = class[,colSums(class>0)>nrow(class) * 0.1]
class_clr=class_clr[,colnames(class_clr) %in% colnames(class)]

write.table(class,file = "HMP2/Input.class.txt",row.names = T,quote = F,sep = "\t")
write.table(class_clr,file = "HMP2/Input.class_CLR.txt",row.names = T,quote = F,sep = "\t")

# merge all
identical(rownames(genus_clr), rownames(family_clr))
bacteria=cbind(phylum_clr,order_clr,class_clr,family_clr,genus_clr)
write.table(bacteria,file = "HMP2/Input.bacteria_CLR.txt",row.names = T,quote = F,sep = "\t")

# ===============================================================================================
# starting analysis
# ===============================================================================================

phenotype=read.table("HMP2/Coupling.16S.txt",sep = "\t",fill = T,stringsAsFactors = F,header = T,check.names = F)
phenotype$Location=NA
phenotype$Location[phenotype$biopsy_location!="Ileum"]="Colon"
phenotype$Location[phenotype$biopsy_location=="Ileum"]="Ileum"

# ===============================================================================================
# beta and alpha diversity, genus level
# ===============================================================================================

genus=read.table("HMP2/Input.genus.txt")
genus_clr=read.table("HMP2/Input.genus_CLR.txt")

# shannon
shannon=as.data.frame(diversity(t(genus), index="shannon"))
colnames(shannon)="shannon"
shannon=merge(shannon,phenotype,by.x="row.names",by.y="Run",all=F)
ggplot(shannon, aes(x=diagnosis, y=shannon, fill=diagnosis)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+
  scale_fill_npg()+
  xlab("")+
  theme_bw() +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') 
wilcox.test(shannon$shannon[shannon$diagnosis=="CD"],shannon$shannon[shannon$diagnosis=="nonIBD"])
wilcox.test(shannon$shannon[shannon$diagnosis=="UC"],shannon$shannon[shannon$diagnosis=="nonIBD"])
wilcox.test(shannon$shannon[shannon$diagnosis=="CD"],shannon$shannon[shannon$diagnosis=="UC"])
ggsave("OutputPlot/HMP2.alpha.diversity.pdf",width = 5,height = 3)

# bray-curtis
beta_diversity=vegdist((genus_clr),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,phenotype,all = F,by.x="row.names",by.y="Run")
ggplot (pca_analysis, aes(V1,V2,color=diagnosis)) + 
  geom_point() + theme_bw() +
  guides(size=F)+
  scale_color_npg()+
  stat_ellipse(type = "norm")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),legend.position = 'top') +xlab("PCA1")+ylab("PCA2")
ggsave("OutputPlot//HMP2.PCA.bac.diagnosis.pdf",width = 5,height = 3)

distmat <- as.matrix(beta_diversity)
cohort1.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="CD"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="CD"]]
cohort2.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="UC"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="UC"]]
cohort3.dis=distmat[rownames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="nonIBD"],colnames(distmat) %in% pca_analysis$Row.names[pca_analysis$diagnosis=="nonIBD"]]

distmat1=data.frame(as.table(cohort1.dis))[lower.tri(cohort1.dis, diag = F), ]
distmat2=data.frame(as.table(cohort2.dis))[lower.tri(cohort2.dis, diag = F), ]
distmat3=data.frame(as.table(cohort3.dis))[lower.tri(cohort3.dis, diag = F), ]
distmat1$group="CD"
distmat2$group="UC"
distmat3$group="nonIBD"
distmat=rbind(distmat1,distmat2,distmat3)

ggplot(distmat, aes(x=group, y=Freq,fill=group)) + 
  geom_boxplot(color="black")+scale_fill_npg()+theme_bw()+guides(fill=F)
ggsave("OutputPlot//HMP2.dissimilarity.genus.pdf",width = 5,height = 3)
wilcox.test(distmat$Freq[distmat$group=="nonIBD"],distmat$Freq[distmat$group=="CD"])
wilcox.test(distmat$Freq[distmat$group=="nonIBD"],distmat$Freq[distmat$group=="UC"])
wilcox.test(distmat$Freq[distmat$group=="CD"],distmat$Freq[distmat$group=="UC"])

# ===============================================================================================
# paired analysis, highli personale
# ===============================================================================================

beta_diversity=vegdist((genus_clr),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,phenotype,all = F,by.x="row.names",by.y="Run")
rownames(pca_analysis)=pca_analysis$Row.names
dd <- dist(scale(genus_clr), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
dend <- as.dendrogram(hc)
pca_analysis$`Participant ID`=as.factor(pca_analysis$`Participant ID`)
groupCodes=pca_analysis$`Participant ID`
colorCodes <- distinctColorPalette(368)
labels_colors(dend) <- colorCodes[groupCodes][order.dendrogram(dend)]

pdf("OutputPlot/HMP2.Tree.plot.pdf",width = 12,height = 6)
#circlize_dendrogram(dend)
plot(dend)
dev.off()

# ===============================================================================================
# HALLA
# ===============================================================================================

halla_covariate=read.table("HMP2/Input.phenotype.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
halla_phylum=read.table("HMP2/Input.genus_CLR.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
rownames(halla_covariate)=halla_covariate$Run
halla_covariate$Run=NULL
halla_covariate$Project=NULL
halla_covariate$`External ID`=NULL
halla_covariate$`Participant ID`=NULL
halla_covariate$data_type=NULL
halla_covariate$`Research Project`=NULL
halla_covariate$biopsy_location=NULL

halla_covariate$CD.vs.control=NA
halla_covariate$CD.vs.control[halla_covariate$diagnosis=="CD"]=1
halla_covariate$CD.vs.control[halla_covariate$diagnosis=="nonIBD"]=0

halla_covariate$UC.vs.control=NA
halla_covariate$UC.vs.control[halla_covariate$diagnosis=="UC"]=1
halla_covariate$UC.vs.control[halla_covariate$diagnosis=="nonIBD"]=0

halla_covariate$CD.vs.UC=NA
halla_covariate$CD.vs.UC[halla_covariate$diagnosis=="UC"]=0
halla_covariate$CD.vs.UC[halla_covariate$diagnosis=="CD"]=1

halla_covariate$IBD.vs.control=NA
halla_covariate$IBD.vs.control[halla_covariate$diagnosis=="nonIBD"]=0
halla_covariate$IBD.vs.control[halla_covariate$diagnosis!="nonIBD"]=1

halla_covariate$diagnosis=NULL
halla_covariate$Antibiotics=NULL
halla_covariate$Chemotherapy=NULL
halla_covariate$`Immunosuppressants (e.g. oral corticosteroids)`=NULL
halla_covariate$Arthralgia=NULL
halla_covariate$hbi=NULL
halla_covariate$`CRP (mg/L)`=NULL
halla_covariate$baseline_montreal_location=NULL
halla_covariate$BMI=NULL

halla_covariate$Location[halla_covariate$Location=="Colon"]=1
halla_covariate$Location[halla_covariate$Location=="Ileum"]=0
halla_covariate$sex[halla_covariate$sex=="Female"]=1
halla_covariate$sex[halla_covariate$sex=="Male"]=0
halla_covariate$is_inflamed[halla_covariate$is_inflamed=="Yes"]=1
halla_covariate$is_inflamed[halla_covariate$is_inflamed=="No"]=0
halla_covariate$is_inflamed=as.numeric(halla_covariate$is_inflamed)
halla_covariate$Location=as.numeric(halla_covariate$Location)
halla_covariate$sex=as.numeric(halla_covariate$sex)

halla_covariate=as.data.frame(t(halla_covariate))
halla_phylum=as.data.frame(t(halla_phylum))

halla_covariate=halla_covariate[,order(colnames(halla_covariate))]
halla_phylum=halla_phylum[,order(colnames(halla_phylum))]

mystudy_phylum=read.table("OutputTable/CLR.genus.txt",sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
halla_phylum=halla_phylum[rownames(halla_phylum) %in% colnames(mystudy_phylum),]

write.table(halla_phylum,"HMP2/Halla.genus.txt",sep = "\t",row.names = T,quote = F)
write.table(halla_covariate,"HMP2/Halla.covariate.txt",sep = "\t",row.names = T,quote = F)

# ===============================================================================================
# dysbiosis
# ===============================================================================================
covariate=read.table("HMP2/Input.phenotype.txt",sep = "\t",check.names = F,stringsAsFactors = F,header = T)
genus=read.table("HMP2/Input.genus_CLR.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
mystudy_genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)

# calculate median distance score as dysbiosis
beta_diversity=vegdist(genus,method = "euclidean")
distmat <- (as.matrix(beta_diversity))
distmat=data.frame(as.table(distmat))
distmat_IBD=distmat[distmat$Var1 %in% covariate$Run[covariate$diagnosis=="CD"] | distmat$Var1 %in% covariate$Run[covariate$diagnosis=="UC"],]
distmat_IBD=distmat_IBD[distmat_IBD$Var2 %in% (covariate$Run)[covariate$diagnosis=="nonIBD"],]
distmat_control=distmat[distmat$Var1 %in% (covariate$Run)[covariate$diagnosis=="nonIBD"],]
distmat_control=distmat_control[distmat_control$Var2 %in% (covariate$Run)[covariate$diagnosis=="nonIBD"],]
distmat_control=distmat_control[distmat_control$Freq!=0,]

sample_IBD=as.character(unique(distmat_IBD$Var1))
dysbiosis_IBD = foreach(i=1:length(sample_IBD),.combine = rbind) %do%  {
  tmp.sample=sample_IBD[i]
  tmp.cal=distmat_IBD[distmat_IBD$Var1==tmp.sample,]
  tmp.median=median(tmp.cal$Freq)
  
  return.string=data.frame(ID=tmp.sample,Dysbiosis=tmp.median,Group="IBD")
}
sample_control=as.character(unique(distmat_control$Var1))
dysbiosis_control = foreach(i=1:length(sample_control),.combine = rbind) %do%  {
  tmp.sample=sample_control[i]
  tmp.cal=distmat_control[distmat_control$Var1==tmp.sample,]
  tmp.median=median(tmp.cal$Freq)
  
  return.string=data.frame(ID=tmp.sample,Dysbiosis=tmp.median,Group="Control")
}
dysbiosis_all=rbind(dysbiosis_control,dysbiosis_IBD)

pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,dysbiosis_all,all = F,by.x="row.names",by.y="ID")

ggplot (pca_analysis, aes(V1,-V2,color=Dysbiosis)) + 
  geom_point() + theme_bw() +
  guides(size=F)+scale_color_jcolors_contin("pal3", reverse = TRUE, bias = 2.25)+
  theme(legend.position = 'top')+xlab("PCA1")+ylab("PCA2")
ggsave("OutputPlot///HMP2.Dysbiosis.PCA.pdf",width = 5,height = 5)
write.table(dysbiosis_all,file = "OutputTable/HMP2.Dysbiosis.score.all.txt",sep = "\t",row.names = F,quote = F)

# dysbiosis vs. disease group
dysbiosis_all_group=dysbiosis_all
dysbiosis_all_group$Group=as.character(dysbiosis_all_group$Group)
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% (covariate$Run)[covariate$diagnosis=="CD"]]="CD"
dysbiosis_all_group$Group[dysbiosis_all_group$ID %in% (covariate$Run)[covariate$diagnosis=="UC"]]="UC"
dysbiosis_all_group=dysbiosis_all_group[dysbiosis_all_group$Group!="IBD",]
ggplot(dysbiosis_all_group, aes(x=Dysbiosis)) + geom_density(aes(group=Group, fill=Group), alpha=0.5) +
  scale_color_npg()+scale_fill_nejm()+theme_bw()
ggsave("OutputPlot/HMP2.Dysbiosis.density.pdf",width = 7,height = 4)

# dysbiosis vs. bacteria
dys_genus = foreach(i=1:ncol(genus),.combine = rbind) %do%  {
  tmp.genus=genus[,i,drop=F]
  tmp.data=merge(dysbiosis_all,tmp.genus,by.x="ID",by.y="row.names",all=F)
  tmp.test=cor.test(tmp.data$Dysbiosis,tmp.data[,4],method = "spearman")
  
  return.string=data.frame(Bacteria=colnames(genus)[i],Cor=tmp.test$estimate,Pvalue=tmp.test$p.value)
}
dys_genus=dys_genus[dys_genus$Bacteria %in% colnames(mystudy_genus),]
dys_genus$FDR=p.adjust(dys_genus$Pvalue,method = "BH")

bacteria="Escherichia/Shigella"
tmp.data=merge(dysbiosis_all,genus[,bacteria,drop=F],by.x="ID",by.y="row.names",all=F)
ggplot(tmp.data, aes(`Escherichia/Shigella`,Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2,fill="#DDAA33")+
  geom_smooth(method = lm,color="black")+theme_bw()

# heatmap
library(pheatmap)
mystudy_dys=read.table("OutputTable/Dys.genus.txt",sep = "\t",header = T,stringsAsFactors = F)
dys_genus=dys_genus[dys_genus$Bacteria %in% mystudy_dys$Bacteria,]

associations_cor=merge(mystudy_dys[,1:2],dys_genus[,1:2],by="Bacteria",all=F)
colnames(associations_cor)=c("Bacteria","Discovery","Replication")
associations_sig=merge(mystudy_dys[,c(1,4)],dys_genus[,c(1,4)],by="Bacteria",all=F)
colnames(associations_sig)=c("Bacteria","Discovery","Replication")
rownames(associations_cor)=associations_cor$Bacteria
associations_cor$Bacteria=NULL
rownames(associations_sig)=associations_sig$Bacteria
associations_sig$Bacteria=NULL
associations_sig$Replication[associations_sig$Discovery>0.05]=9
associations_sig$Replication[associations_cor$Discovery*associations_cor$Replication<0]=9

associations_cor=as.data.frame(t(associations_cor))
associations_sig=as.data.frame(t(associations_sig))

heatmap.data=(associations_cor)
heatmap.data.lable=(associations_sig)
heatmap.data.lable[heatmap.data.lable<0.05 & heatmap.data.lable>=0]="*"
heatmap.data.lable[heatmap.data.lable>0.05]=""

pdf("OutputPlot//host-microbe.association.pdf",width = 5,height = 10)
pheatmap(heatmap.data,cluster_cols = T, cluster_rows = F,
         show_rownames=T, show_colnames=T,
         cellheight = 8,cellwidth = 10,fontsize_number=12,fontsize_row=6,fontsize_col = 8,
         display_numbers = heatmap.data.lable,color = colorRampPalette(c("#4393C3",mid="white", "#FEDA8B"))(100))
dev.off()

library(circlize) 
pdf("OutputPlot/test.font.pdf")
col_fun1 = colorRamp2(c(-1, 0, 1), c("#0077BB",mid="white", "#997700"))
circos.heatmap(as.matrix(t(heatmap.data)),col = col_fun1,dend.side = "inside",track.height = 0.1,rownames.side = "outside")
circos.clear()
dev.off()

bacteria="Veillonella"
tmp.data=merge(dysbiosis_all,genus[,bacteria,drop=F],by.x="ID",by.y="row.names",all=F)
ggplot(tmp.data, aes(`Veillonella`,Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2,fill="#CC3311")+
  geom_smooth(method = lm,color="black")+theme_bw()
ggsave("OutputPlot/HMP2.Dysbiosis.bacteria.pdf",width = 3,height = 3)
consistent_genus=as.data.frame(t(heatmap.data.lable))
consistent_genus=rownames(consistent_genus)[consistent_genus$Replication=="*"]
write.table(consistent_genus,file = "OutputTable/1000IBD,HMP2.consistent.genus.txt",row.names = F,quote = F)

# ===============================================================================================
# RNAseq data CLR
# ===============================================================================================
rnaseq_phenotype=read.table("HMP2/Coupling.RNAseq.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t")
rnaseq=read.table("HMP2/GSE111889_host_tx_counts.tsv",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
mystudy_rnaseq=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
rnaseq=rnaseq[rownames(rnaseq) %in% colnames(mystudy_rnaseq),]
rnaseq=as.data.frame(t(rnaseq))
rnaseq_clr <- zCompositions::cmultRepl(rnaseq, method="CZM", label=0)
rnaseq_clr=compositions::clr(rnaseq_clr)
rnaseq_clr=as.data.frame(rnaseq_clr)

# ===============================================================================================
# RNAseq data correction
# ===============================================================================================
# RNA data correction (age, gender, BMI, smoking, bacth, inflammation, location (no diagnosis, location == diagnosis))
covariate_rna=rnaseq_phenotype[,c("biopsy_location","consent_age","is_inflamed","sex")]
rownames(covariate_rna)=rnaseq_phenotype$`External ID`
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Terminal ileum"]="Ileum"
covariate_rna=covariate_rna[covariate_rna$biopsy_location!="Non-inflamed",]
covariate_rna$biopsy_location[covariate_rna$biopsy_location!="Ileum"]="Colon"
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Ileum"]=0
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Colon"]=1
covariate_rna$is_inflamed[covariate_rna$is_inflamed=="Yes"]=1
covariate_rna$is_inflamed[covariate_rna$is_inflamed=="No"]=0
covariate_rna$sex[covariate_rna$sex=="Female"]=0
covariate_rna$sex[covariate_rna$sex=="Male"]=1
covariate_rna$biopsy_location=as.numeric(covariate_rna$biopsy_location)
covariate_rna$is_inflamed=as.numeric(covariate_rna$is_inflamed)
covariate_rna$sex=as.numeric(covariate_rna$sex)
covariate_rna=na.omit(covariate_rna)

rnaseq_clr=rnaseq_clr[rownames(rnaseq_clr) %in% rownames(covariate_rna),]
covariate_rna=covariate_rna[rownames(covariate_rna) %in% rownames(rnaseq_clr),]
rnaseq_clr=rnaseq_clr[order(rownames(rnaseq_clr)),]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]

pca=prcomp(rnaseq_clr,scale = TRUE)
fviz_eig(pca)

eigenvalue=get_eig(pca)
ind <- get_pca_ind(pca)
pca_matrix=as.data.frame(ind$coord)
pca_matrix=pca_matrix[,1:6]

pca_matrix=merge(pca_matrix,covariate_rna,by="row.names",all=F)
pca_matrix$is_inflamed=as.factor(pca_matrix$is_inflamed)
pca_matrix$biopsy_location=as.factor(pca_matrix$biopsy_location)

p1=ggscatter(pca_matrix, x = "Dim.1", y = "Dim.2",
             color = "biopsy_location", palette = "jco",
             shape = "biopsy_location",
             ellipse = TRUE, 
             mean.point = TRUE,
             star.plot = TRUE,alpha=0.5)
p2=ggscatter(pca_matrix, x = "Dim.1", y = "Dim.2",
             color = "is_inflamed", palette = c("#117733", "#88CCEE", "#BB5566"),
             shape = "is_inflamed",
             ellipse = TRUE, 
             mean.point = TRUE,
             star.plot = TRUE,alpha=0.5)
ggarrange(p1,p2)
ggsave("HMP2/PCA.CLR.RNAseq.pdf",width = 8,height = 4)

rnaseq_correct = apply(rnaseq_clr,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_rna[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="tmp.rnaseq"
  
  fit=lm(tmp.rnaseq~biopsy_location+sex+consent_age+is_inflamed,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
rnaseq_correct=as.data.frame(rnaseq_correct)
write.table(rnaseq_correct,"HMP2//CLR.HMP2.rnqseq.txt",sep = "\t",quote = F,row.names = T)

# ===============================================================================================
# RNAseq data + 16S data harmonize, keep only baceria the same in my dataset
# ===============================================================================================
rnaseq_correct=read.table("HMP2/CLR.HMP2.rnqseq.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
bacteria_clr=read.table("HMP2/Input.bacteria_CLR.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)

mystudy_bacteria=read.table(file ="OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T)
bacteria_clr=bacteria_clr[,colnames(bacteria_clr) %in% colnames(mystudy_bacteria)]

coupling=read.table("HMP2/Coupling.RNAseq+16S.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t")
coupling=coupling[coupling$RNAseqID %in% rownames(rnaseq_correct) & coupling$`16SID` %in% rownames(bacteria_clr),]
rnaseq_correct=rnaseq_correct[rownames(rnaseq_correct) %in% coupling$RNAseqID,]
bacteria_clr=bacteria_clr[rownames(bacteria_clr) %in% coupling$`16SID`,]

rnaseq_correct=rnaseq_correct[order(rownames(rnaseq_correct)),]
coupling=coupling[order((coupling$RNAseqID)),]
rownames(rnaseq_correct)=coupling$`16SID`

identical(rownames(rnaseq_correct),rownames(bacteria_clr))

# ===============================================================================================
# sparseCCA, mystudy inflammation-related genes + mystudy bacteria replication
# ===============================================================================================
mystudy_rnaseq=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation_gene=intersect(mystudy_rnaseq$Gene[mystudy_rnaseq$group==1],intersect(mystudy_rnaseq$Gene[mystudy_rnaseq$group==2],mystudy_rnaseq$Gene[mystudy_rnaseq$group==3]))
HMP2_inflammation_gene_correct=rnaseq_correct[,colnames(rnaseq_correct) %in% inflammation_gene]
HMP2_inflammation_gene_correct=HMP2_inflammation_gene_correct[order(rownames(HMP2_inflammation_gene_correct)),]
bacteria_clr=bacteria_clr[order(rownames(bacteria_clr)),]

identical(rownames(bacteria_clr),rownames(HMP2_inflammation_gene_correct))
write.table(HMP2_inflammation_gene_correct,file = "HMP2/HMP2_inflammation_gene_correct_sparseCCA.txt",row.names = T,sep = "\t",quote = F)
write.table(bacteria_clr,file = "HMP2/HMP2_bacteria_clr_sparseCCA.txt",row.names = T,sep = "\t",quote = F)

# ===============================================================================================
# individual gene-bacteria replication
# ===============================================================================================
library(foreach)
rnaseq_correct=read.table("HMP2/CLR.HMP2.rnqseq.txt",sep = "\t",header = T,check.names = F,stringsAsFactors = F,row.names = 1)
bacteria_clr=read.table("HMP2/Input.bacteria_CLR.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)

mystudy_bacteria=read.table(file ="OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T)
bacteria_clr=bacteria_clr[,colnames(bacteria_clr) %in% colnames(mystudy_bacteria)]

coupling=read.table("HMP2/Coupling.RNAseq+16S.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t")
coupling=coupling[coupling$RNAseqID %in% rownames(rnaseq_correct) & coupling$`16SID` %in% rownames(bacteria_clr),]
rnaseq_correct=rnaseq_correct[rownames(rnaseq_correct) %in% coupling$RNAseqID,]
bacteria_clr=bacteria_clr[rownames(bacteria_clr) %in% coupling$`16SID`,]

rnaseq_correct=rnaseq_correct[order(rownames(rnaseq_correct)),]
coupling=coupling[order((coupling$RNAseqID)),]
rownames(rnaseq_correct)=coupling$`16SID`

rnaseq_correct=rnaseq_correct[order(rownames(rnaseq_correct)),]
bacteria_clr=bacteria_clr[order(rownames(bacteria_clr)),]
identical(rownames(rnaseq_correct),rownames(bacteria_clr))

mystudy_pair=read.table("OutputTable/gene.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F)
mystudy_pair=mystudy_pair[mystudy_pair$FDR<0.05,]

rnaseq_correct=rnaseq_correct[,colnames(rnaseq_correct) %in% mystudy_pair$Gene]
bacteria_clr=bacteria_clr[,colnames(bacteria_clr) %in% mystudy_pair$Taxa]

HMP2_pair = foreach(i=1:nrow(mystudy_pair),.combine = rbind) %do%  {
  tmp.gene=mystudy_pair$Gene[i]
  tmp.bacteria=mystudy_pair$Taxa[i]
  
  if(tmp.gene %in% colnames(rnaseq_correct) & tmp.bacteria %in% colnames(bacteria_clr)){
    tmp.data=merge(rnaseq_correct[,tmp.gene,drop=F],bacteria_clr[,tmp.bacteria,drop=F],by="row.names",all=F)
    colnames(tmp.data)[2:3]=c("gene","bacteria")
    tmp.mod=lm(gene ~ bacteria, data = tmp.data)
    tmp.mod=as.data.frame(summary(tmp.mod)$coef)
    
    return.string=data.frame(gene=tmp.gene,bacteria=tmp.bacteria,Beta=tmp.mod$Estimate[rownames(tmp.mod)=="bacteria"],
                             SE=tmp.mod$`Std. Error`[rownames(tmp.mod)=="bacteria"],
                             Zscore=tmp.mod$`t value`[rownames(tmp.mod)=="bacteria"],
                             P=tmp.mod$`Pr(>|t|)`[rownames(tmp.mod)=="bacteria"])
  }else{
    return.string=data.frame(gene=tmp.gene,bacteria=tmp.bacteria,Beta=NA,SE=NA,Zscore=NA,P=NA)
  }
}

HMP2_pair=HMP2_pair[order(HMP2_pair$Zscore,decreasing = T),]
mystudy_pair=mystudy_pair[order(mystudy_pair$Zscore,decreasing = T),]
HMP2_pair$rank=1:nrow(HMP2_pair)
mystudy_pair$rank=1:nrow(mystudy_pair)

HMP2_pair$pair=paste(HMP2_pair$gene,HMP2_pair$bacteria)
mystudy_pair$pair=paste(mystudy_pair$Gene,mystudy_pair$Taxa)
HMP2_pair=HMP2_pair[order(HMP2_pair$pair),]
mystudy_pair=mystudy_pair[order(mystudy_pair$pair),]

cor.test(mystudy_pair$rank,HMP2_pair$rank,method = "spearman")
cor.test(mystudy_pair$Beta,HMP2_pair$Beta)

write.table(HMP2_pair,"HMP2/replication.bacteria.gene.txt",sep="\t",row.names = F,quote = F)

# ===============================================================================================
# classify dysbiosis
# ===============================================================================================
quantile(dysbiosis_all$Dysbiosis)
dysbiosis_quantitle=dysbiosis_all
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis<=31.52727]=0
dysbiosis_quantitle$Dysbiosis[dysbiosis_quantitle$Dysbiosis>31.52727]=1
write.table(dysbiosis_all,file = "OutputTable/HMP2.dysbiosis.group.txt",sep = "\t",row.names = F,quote = F)

# ===============================================================================================
# model: gene ~ shannon 
# ===============================================================================================
shannon <- as.data.frame(diversity(t(genus), index="shannon"))
modelfit = foreach(i=1:ncol(rnaseq_correct),.combine = rbind) %do%  {
  tmp.gene=colnames(rnaseq_correct)[i]
  tmp.data=merge(rnaseq_correct[,tmp.gene,drop=F],shannon,all=F,by="row.names")
  colnames(tmp.data)=c("sampleID","gene","shannon")
  
  tmp.model=lm(gene ~ shannon,data = tmp.data)
  tmp.model=as.data.frame(summary(tmp.model)$coefficients)
  cat(green(i,"+++",tmp.gene,"\n"))
  
  return.string=data.frame(Gene=tmp.gene,Bac="Shannon",Beta=tmp.model$Estimate[2],Pvalue=tmp.model$`Pr(>|t|)`[2])
}

# ===============================================================================================
# model: gene ~ bac + dysbiosis 
# ===============================================================================================
mystudy_nointeraction=read.table("OutputTable/Dysbiosis.Nointeraction.result.txt",header = T)
mystudy_nointeraction=mystudy_nointeraction[mystudy_nointeraction$FDR_taxa<0.05,]
rnaseq_correct_tmp=rnaseq_correct[,colnames(rnaseq_correct) %in% mystudy_nointeraction$Gene]

inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation_genes=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
rnaseq_correct_tmp=rnaseq_correct[,colnames(rnaseq_correct) %in% inflammation_genes]

modelfit0 = foreach(i=1:ncol(rnaseq_correct_tmp),.combine = rbind) %do%  {
  
  tmp.gene=colnames(rnaseq_correct_tmp)[i]
  tmp.data=merge(rnaseq_correct_tmp[,tmp.gene,drop=F],genus_clr,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(genus_clr),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_dys","P_taxa","P_dys")
  
  for(n in 1:ncol(genus_clr)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=lm(Gene ~ Bac + Dys,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$Beta_dys[n]=tmp.model$Estimate[rownames(tmp.model)=="Dys"]
    tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    tmp.result$P_dys[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Dys"]
    
  }
  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=as.data.frame(tmp.result)
}
match=as.data.frame(which( outer(mystudy_nointeraction$Gene, modelfit0$Gene, "==") & 
         outer(mystudy_nointeraction$Taxa, modelfit0$Taxa, "=="), 
       arr.ind=TRUE))
modelfit0=modelfit0[match$col,]
modelfit0$FDR=p.adjust(modelfit0$P_taxa,method = "BH")
write.table(modelfit0,file = "HMP2//HMP2.noInteraction.txt",row.names = F,quote = F)

# compare with two results, 1000IBD and HMP2
modelfit0$pair=paste(modelfit0$Gene,modelfit0$Taxa,sep = "_")
mystudy_nointeraction$pair=paste(mystudy_nointeraction$Gene,mystudy_nointeraction$Taxa,sep = "_")
compare=merge(mystudy_nointeraction,modelfit0,by="pair",all=F)
plot(compare$Beta_taxa.x,compare$Beta_taxa.y)
cor.test(compare$Beta_taxa.x,compare$Beta_taxa.y)

mystudy_nointeraction=mystudy_nointeraction[mystudy_nointeraction$pair %in% modelfit0$pair,]
modelfit0=modelfit0[modelfit0$pair %in% mystudy_nointeraction$pair,]

mystudy_nointeraction=mystudy_nointeraction[order(abs(mystudy_nointeraction$Beta_taxa),decreasing = T),]
modelfit0=modelfit0[order(abs(modelfit0$Beta_taxa),decreasing = T),]
mystudy_nointeraction$rank=1:nrow(mystudy_nointeraction)
modelfit0$rank=1:nrow(modelfit0)

mystudy_nointeraction=mystudy_nointeraction[1:500,]
compare=merge(mystudy_nointeraction,modelfit0,by="pair",all=F)
plot(compare$Beta_taxa.x,compare$Beta_taxa.y)
cor.test(compare$Beta_taxa.x,compare$Beta_taxa.y)

library("ggExtra")
tmp.genus="Intestinibacter"
tmp.gene="SLC25A29"
tmp.data=merge(dysbiosis_quantitle,genus_clr[,tmp.genus,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data=merge(tmp.data,rnaseq_correct[,tmp.gene,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data$Dysbiosis[tmp.data$Dysbiosis==0]="Dysbiosis score [0~75%]"
tmp.data$Dysbiosis[tmp.data$Dysbiosis==1]="Dysbiosis score [75~100%]"

p=ggplot(tmp.data) +
  geom_point(aes(Intestinibacter,SLC25A29,fill=Dysbiosis),shape = 21, color = "black", size = 2)+
  geom_smooth(aes(Intestinibacter,SLC25A29),method = lm,color="black",linetype="dashed") +
  theme(legend.position="bottom")+scale_fill_manual(values=c("#BB5566", "#004488"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
y_box = axis_canvas(p, axis = "y", coord_flip = TRUE) +
  geom_density(data = tmp.data, aes(x=SLC25A29,fill=Dysbiosis),alpha=0.5) +
  scale_fill_manual(values = c("#BB5566", "#004488"))+
  coord_flip()
combined_plot <- insert_yaxis_grob(p, y_box, position = "right")
x_box = axis_canvas(p, axis = "x") +
  geom_density(data = tmp.data, aes(x=Intestinibacter,fill=Dysbiosis),alpha=0.5) +
  scale_fill_manual(values = c("#BB5566", "#004488"))
combined_plot <- insert_xaxis_grob(combined_plot, x_box, position = "top")
ggdraw(combined_plot)
ggsave("HMP2//Intestinibacter.SLC25A29.pdf",width = 5,height = 4)

# ===============================================================================================
# model: gene ~ bac + dysbiosis + bac * dysbiosis
# ===============================================================================================
mystudy_interaction=read.table("OutputTable/Dysbiosis.interaction.result.txt",header = T)
mystudy_interaction=mystudy_interaction[mystudy_interaction$FDR_interact<0.05,]
rnaseq_correct_tmp=rnaseq_correct[,colnames(rnaseq_correct) %in% mystudy_interaction$Gene]
modelfit1 = foreach(i=1:ncol(rnaseq_correct_tmp),.combine = rbind) %do%  {
  
  tmp.gene=colnames(rnaseq_correct_tmp)[i]
  tmp.data=merge(rnaseq_correct_tmp[,tmp.gene,drop=F],genus_clr,by="row.names",all=F)
  tmp.data=merge(dysbiosis_all,tmp.data,by.y="Row.names",by.x="ID",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(genus_clr),ncol = 6))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_dys","P_taxa","P_dys")
  
  for(n in 1:ncol(genus_clr)){
    tmp.test=tmp.data[,c(1:4,n+4)]
    colnames(tmp.test)=c("ID","Dys","Group","Gene","Bac")
    tmp.model=lm(Gene ~ Bac + Dys + Bac*Dys,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=colnames(tmp.data)[n+4]
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$Beta_dys[n]=tmp.model$Estimate[rownames(tmp.model)=="Dys"]
    tmp.result$Beta_interact[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac:Dys"]
    tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    tmp.result$P_dys[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Dys"]
    tmp.result$P_interact[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac:Dys"]
    
  }
  cat(green(i,"+++",tmp.gene,"\n"))
  return.string=as.data.frame(tmp.result)
}
match=as.data.frame(which( outer(mystudy_interaction$Gene, modelfit1$Gene, "==") & 
                             outer(mystudy_interaction$Taxa, modelfit1$Taxa, "=="), 
                           arr.ind=TRUE))
modelfit1=modelfit1[match$col,]
modelfit1$FDR=p.adjust(modelfit1$P_interact,method = "BH")

tmp.genus="Agathobacter"
tmp.gene="GUCA2A"
tmp.data=merge(dysbiosis_quantitle,genus_clr[,tmp.genus,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data=merge(tmp.data,rnaseq_correct[,tmp.gene,drop=F],by.x="ID",by.y="row.names",all=F)
tmp.data$Dysbiosis[tmp.data$Dysbiosis==0]="Dysbiosis score [0~75%]"
tmp.data$Dysbiosis[tmp.data$Dysbiosis==1]="Dysbiosis score [75~100%]"
ggplot(tmp.data, aes(Agathobacter,GUCA2A,fill=Dysbiosis)) +
  geom_point(shape = 21, color = "black", size = 2)+
  geom_smooth(method = lm,color="black",linetype="dashed")+
  facet_wrap(~ Dysbiosis)+
  theme(legend.position="bottom")+scale_fill_manual(values=c("#BB5566", "#004488"))+
  theme(panel.background = element_rect(fill = "white", colour = "black"))
ggsave("HMP2//Agathobacter.GUCA2A.pdf",width = 4,height = 4)

# ===============================================================================================
# cell type
# ===============================================================================================
rnaseq=read.table("HMP2/GSE111889_host_tx_counts.tsv",row.names = 1,header = T,stringsAsFactors = F,check.names = F)
deconvolution=xCellAnalysis(rnaseq)
deconvolution=as.data.frame(deconvolution)
deconvolution=as.data.frame(t(deconvolution))
write.table(deconvolution,file = "HMP2//Deconvolution.txt",sep = "\t",row.names = T,quote = F)

# RNA data correction (age, gender, BMI, smoking, bacth, inflammation, location (no diagnosis, location == diagnosis))
deconvolution=read.table("HMP2/Deconvolution.txt",sep = "\t",row.names = 1,header = T)
rnaseq_phenotype=read.table("HMP2/Coupling.RNAseq.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t")
covariate_rna=rnaseq_phenotype[,c("biopsy_location","consent_age","is_inflamed","sex","diagnosis")]
rownames(covariate_rna)=rnaseq_phenotype$`External ID`
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Terminal ileum"]="Ileum"
covariate_rna=covariate_rna[covariate_rna$biopsy_location!="Non-inflamed",]
covariate_rna$biopsy_location[covariate_rna$biopsy_location!="Ileum"]="Colon"
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Ileum"]=0
covariate_rna$biopsy_location[covariate_rna$biopsy_location=="Colon"]=1
covariate_rna$is_inflamed[covariate_rna$is_inflamed=="Yes"]=1
covariate_rna$is_inflamed[covariate_rna$is_inflamed=="No"]=0
covariate_rna$sex[covariate_rna$sex=="Female"]=0
covariate_rna$sex[covariate_rna$sex=="Male"]=1
covariate_rna$biopsy_location=as.numeric(covariate_rna$biopsy_location)
covariate_rna$is_inflamed=as.numeric(covariate_rna$is_inflamed)
covariate_rna$sex=as.numeric(covariate_rna$sex)
covariate_rna=na.omit(covariate_rna)

deconvolution=deconvolution[rownames(deconvolution) %in% rownames(covariate_rna),]
covariate_rna=covariate_rna[rownames(covariate_rna) %in% rownames(deconvolution),]
deconvolution=deconvolution[order(rownames(deconvolution)),]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]

deconvolution_sub=deconvolution[,c("cDC",	"Macrophages.M1",	"NK.cells",	"pDC",	"Macrophages.M2",
                                   "CD4..naive.T.cells",	"CD4..Tcm",	"CD8..naive.T.cells",	"CD8..Tcm",	
                                   "Tgd.cells",	"Th2.cells",	"Tregs",	"Th1.cells",
                                   "NKT",	"CD8..Tem",	"CD4..Tem",	"Class.switched.memory.B.cells",
                                   "Plasma.cells",	"naive.B.cells",	"Memory.B.cells",	"Basophils",
                                   "Mast.cells",	"Neutrophils",	"Eosinophils",	"Endothelial.cells",
                                   "Epithelial.cells","Fibroblasts"),drop=F]

# PCA
beta_diversity=vegdist((deconvolution_sub),method = "euclidean")
pca_analysis=as.data.frame(cmdscale(beta_diversity,k=8))
pca_analysis=merge(pca_analysis,covariate_rna,all = F,by="row.names")

pca_analysis[pca_analysis==""]=NA
ggplot (pca_analysis, aes(V1,V2,fill=diagnosis,alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  stat_ellipse(type = "norm")+
  scale_fill_jama()+
  theme(legend.position = 'top')
ggsave("HMP2//Inflammation.diagnosis.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=as.factor(is_inflamed),alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  scale_fill_aaas()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')
ggsave("HMP2//Inflammation.cell.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=as.factor(biopsy_location),alpha=0.5)) + 
  geom_point(shape = 21, size = 2)+ theme_bw() +
  guides(size=F)+
  scale_fill_simpsons()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')
ggsave("HMP2/Location.cell.PCA.pdf",width = 3,height = 3)
ggplot (pca_analysis, aes(V1,V2,fill=Batch,alpha=0.5)) + 
  geom_point(shape = 21, size = 2) + theme_bw() +
  guides(size=F)+
  scale_fill_nejm()+
  stat_ellipse(type = "norm")+
  theme(legend.position = 'top')
ggsave("OutputPlot/Location.batch.PCA.pdf",width = 3,height = 3)

# compare
compare=matrix(nrow = ncol(deconvolution_sub),ncol = 2)
compare=as.data.frame(compare)
colnames(compare)=c("CellType","Pvalue")
for(i in 1:ncol(deconvolution_sub)){
  cell=colnames(deconvolution_sub)[i]
  inflamed=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$is_inflamed==1],]
  noninf=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$is_inflamed==0],]
  mm=wilcox.test(inflamed[,cell],noninf[,cell])
  Pvalue=mm$p.value
  compare$CellType[i]=cell
  compare$Pvalue[i]=Pvalue
}
compare$FDR=p.adjust(compare$Pvalue)
compare=compare[order(compare$FDR),]

inflamed=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$is_inflamed==1],]
noninf=deconvolution_sub[rownames(deconvolution_sub) %in% rownames(covariate_rna)[covariate_rna$is_inflamed==0],]
compare_all=data.frame(Score=NA,Group=NA,Cell=NA,Sample=NA)
for(i in 1:nrow(compare)){
  cell=compare$CellType[i]
  sub.in=inflamed[,cell,drop=F]
  sub.non=noninf[,cell,drop=F]
  sub.in$Group="Inflammation"
  sub.non$Group="Non_inflammation"
  sub.in$Cell=cell
  sub.non$Cell=cell
  colnames(sub.in)[1]="Score"
  colnames(sub.non)[1]="Score"
  sub.data=rbind(sub.in,sub.non)
  sub.data$Sample=rownames(sub.data)
  rownames(sub.data)=NULL
  compare_all=rbind(compare_all,sub.data)
}
compare_all=na.omit(compare_all)
compare_all$Cell=factor(compare_all$Cell,levels = rev(unique(compare$CellType)))

ggplot (compare_all,aes(x= Cell, y=Score, fill=Group))+
  geom_boxplot(alpha=0.8,outlier.shape = NA)+
  geom_point(aes(color = Group),alpha=0.1,position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.5),size=1)+
  theme_bw()+
  scale_color_manual(values=c("#BB5566","#004488"))+
  scale_fill_manual(values=c("#BB5566","#004488"))+
  theme(axis.text.x = element_text(angle = 90,hjust = 1))+
  xlab("CellTypes")+ylab("EnrichmentScore")+
  coord_flip()+guides(fill=F)+guides(color=F)
ggsave("HMP2//Cell.compare.pdf",width = 8,height = 8)

# ====================================== ########## ============================================
# ====================================== predictors ============================================
# ====================================== ########## ============================================
library (hdi)
library(tidyverse)
library (dplyr)
library(glmnet)
library(glmnetUtils)
library(crayon)
library(caret)
library(pROC)
library(plyr)
library(readr)
lasso_fit=function(feature,outcome){
  if(length(colnames(feature))==0){
    result=data.frame(Protein=colnames(outcome),Feature="No",Lasso.beta=0)
  }else if(length(colnames(feature))==1){
    model=lm(outcome[,1]~ feature[,1])
    beta = summary(model)$coefficients[2]
    result=data.frame(Protein=colnames(outcome),Feature=colnames(feature),Lasso.beta=beta)
  }else{
    cv=cv.glmnet(as.matrix(feature),as.matrix(outcome), alpha = 1, nfolds = 10, type.measure="mse",standardize=T)
    beta <- (coef(cv, s = "lambda.min"))
    beta=as.data.frame(beta[-1,1])
    beta$variable=rownames(beta)
    colnames(beta)[1]="beta"
    result=data.frame(Protein=colnames(outcome),Feature=beta$variable,Lasso.beta=beta$beta)
  }
}

covariate_rna=covariate_rna[,c("biopsy_location","consent_age","is_inflamed","sex")]
deconvolution_sub=deconvolution_sub[,colSums(deconvolution_sub!=0)>(0.8*nrow(deconvolution_sub))]
deconvolution_sub[deconvolution_sub==0]=NA
deconvolution_correct = apply(deconvolution_sub,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_rna[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="tmp.rnaseq"
  
  fit=lm(tmp.rnaseq~biopsy_location+sex+consent_age+is_inflamed,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
deconvolution_correct=as.data.frame(deconvolution_correct)

mystudy_genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",check.names = F,stringsAsFactors = F,row.names = 1)
genus_clr=read.table("HMP2/Input.genus_CLR.txt",sep = "\t",check.names = F,row.names = 1,header = T,stringsAsFactors = F)
genus_clr=genus_clr[,colnames(genus_clr) %in% colnames(mystudy_genus)]

coupling=read.table("HMP2/Coupling.RNAseq+16S.txt",header = T,check.names = F,stringsAsFactors = F,sep = "\t")
coupling=coupling[coupling$RNAseqID %in% rownames(deconvolution_correct) & coupling$`16SID` %in% rownames(genus_clr),]
deconvolution_correct=deconvolution_correct[rownames(deconvolution_correct) %in% coupling$RNAseqID,]
genus_clr=genus_clr[rownames(genus_clr) %in% coupling$`16SID`,]

deconvolution_correct=deconvolution_correct[order(rownames(deconvolution_correct)),]
coupling=coupling[order((coupling$RNAseqID)),]
rownames(deconvolution_correct)=coupling$`16SID`

# CD
deconvolution_sub.correct.CD=deconvolution_correct

set.seed(100)
factors=genus_clr
deconvolution_sub_train=deconvolution_sub.correct.CD

trainset= foreach(i=1:ncol(deconvolution_sub_train),.combine = rbind) %do%  {
  
  tmp.protein=colnames(deconvolution_sub_train)[i]
  tmp.feature=as.character(colnames(factors))
  
  cat(yellow(tmp.protein,"lasso running","\n"))
  
  feature=na.omit(factors[,tmp.feature,drop=F])
  outcome=na.omit(deconvolution_sub_train[,tmp.protein,drop=F])
  
  feature=feature[rownames(feature) %in% rownames(outcome),,drop=F]
  outcome=outcome[rownames(outcome) %in% rownames(feature),,drop=F]
  
  feature=feature[order(rownames(feature)),,drop=F]
  outcome=outcome[order(rownames(outcome)),,drop=F]
  
  result=lasso_fit(feature,outcome)
  return.string = result
}
trainset=trainset[trainset$Lasso.beta!=0,]

trainset.variation=foreach(i=1:length(unique(trainset$Protein)),.combine = rbind) %do%  {
  pro=as.character(unique(trainset$Protein)[i])
  cov=as.character(trainset$Feature[trainset$Protein==pro])
  
  tmp.pro=deconvolution_sub_train[,pro,drop=F]
  tmp.cov=factors[,colnames(factors) %in% cov,drop=F]
  tmp.pro=na.omit(tmp.pro)
  tmp.cov=na.omit(tmp.cov)
  tmp.pro=(tmp.pro[rownames(tmp.pro) %in% rownames(tmp.cov),,drop=F])
  tmp.cov=(tmp.cov[rownames(tmp.cov) %in% rownames(tmp.pro),,drop=F])
  tmp.pro=tmp.pro[order(rownames(tmp.pro)),,drop=F]
  tmp.cov=tmp.cov[order(rownames(tmp.cov)),,drop=F]
  
  tmp.glm=glm(tmp.pro[,1]~.,data=tmp.cov,family = gaussian)
  tmp.coef=as.data.frame(summary(tmp.glm)$coef)
  tmp.coef$FDR=p.adjust(tmp.coef$`Pr(>|t|)`)
  tmp.coef$CIup=tmp.coef$Estimate+1.96*tmp.coef$`Std. Error`
  tmp.coef$CIdown=tmp.coef$Estimate-1.96*tmp.coef$`Std. Error`
  tmp.av=anova(tmp.glm)
  tmp.av$Explain=NA
  for(j in 2:nrow(tmp.av)){
    tmp.av$Explain[j]=(tmp.av$`Resid. Dev`[j-1]-tmp.av$`Resid. Dev`[j])/tmp.av$`Resid. Dev`[1]
  }
  tmp.av=merge(tmp.av,tmp.coef,by="row.names",all=F)
  
  if(nrow(tmp.av)==0){
    return.string=data.frame(MultiVariate=NA,MultiVariate.FDR=NA,MultiVariate.explain=NA,Protein=pro,CIup=NA,CIdonw=NA,Estimate=NA)
  }else{
    return.string=data.frame(MultiVariate=tmp.av$Row.names,MultiVariate.FDR=tmp.av$FDR,
                             MultiVariate.explain=tmp.av$Explain,Protein=pro,CIup=tmp.av$CIup,CIdonw=tmp.av$CIdown,Estimate=tmp.av$Estimate)
  }
}

variation_CD= foreach(i=1:length(unique(trainset.variation$Protein)),.combine = rbind) %do%  {
  
  tmp.pro=as.character(unique(trainset.variation$Protein)[i])
  tmp.var=trainset.variation[trainset.variation$Protein==tmp.pro,]
  tmp.summary=sum(tmp.var$MultiVariate.explain)
  
  
  cat(red(tmp.pro,"done","\n"))
  
  return.string = data.frame(Cell=tmp.pro,Variation=tmp.summary)
}

# compare
mystudy_cellvariation=read.table("OutputTable/Cell.variation.txt",sep = "\t",header = T)
ggplot(mystudy_cellvariation[mystudy_cellvariation$Protein=="Endothelial cells",], aes(x=MultiVariate, y=MultiVariate.explain,fill=MultiVariate,na.rm = F)) +
  geom_bar(stat="identity")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90))+
  guides(fill=F)+coord_flip()+ylab("Variation of Endothelial cells")
ggsave("OutputPlot/test.pdf",width = 6,height = 8)

tmp.pro="CD8+ Tcm"
tmp.bac="Lachnospiraceae_ND3007_group"
tmp.data=merge(bacteria1[,tmp.bac,drop=F],deconvolution_sub[,tmp.pro,drop=F],all=F,by="row.names")
tmp.data$Disease=NA
tmp.data$Disease[tmp.data$Row.names %in% rownames(covariate_rna)[covariate_rna$Diagnosis=="CD"]]="CD"
tmp.data$Disease[tmp.data$Row.names %in% rownames(covariate_rna)[covariate_rna$Diagnosis=="UC"]]="UC"
tmp.data=na.omit(tmp.data)
tmp.data$Lachnospiraceae_ND3007_group=as.factor(tmp.data$Lachnospiraceae_ND3007_group)
ggplot(tmp.data[tmp.data$Disease=="CD",], aes(Lachnospiraceae_ND3007_group, `CD8+ Tcm`,fill=Lachnospiraceae_ND3007_group)) +
  geom_boxplot() + 
  theme_bw()+ theme(legend.position="bottom")+guides(fill=F)+
  #facet_wrap(~ Disease)+
  scale_fill_npg()
ggsave("OutputPlot/CD8+ Tcm.pdf",height = 3,width = 3)
















