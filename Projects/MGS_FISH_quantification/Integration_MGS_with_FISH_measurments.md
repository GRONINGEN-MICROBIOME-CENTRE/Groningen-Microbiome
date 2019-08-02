---
title: "Exploring microbial cell quantification"
author: "Arnau Vich & Alex Kur"
date: "11/13/2018"
output: html_document
---
  
  # Microbial quantification: total cell density effect 
  
  Load Packages 
--
  ```{r}
library(vegan)
library (plyr)
library (reshape2)
library (ggplot2)
library (ppcor)
library (ggforce)
library (dplyr)
library (foreach)
library (Hmisc)
library (ggord)
library (ggpubr)
library (psych)
library(patchwork)
library(phyloseq)
```

## Import datasets

```{r, eval=F}
setwd("~/Desktop/QMP/Running_all")
phenos_LLD=read.table("LLD_phenos.txt", header=T, row.names = 1, sep = "\t")
phenos_1000IBD=read.table("IBD_phenos.txt", header=T, row.names = 1, sep = "\t")
tax_LLD=read.table("LLD_taxonomy_unstrat_clean.txt", header=T, row.names = 1, sep = "\t", check.names = F)
tax_IBD=read.table("IBD_taxonomy_unstrat_clean.txt", header=T, row.names = 1, sep = "\t")
meta_sp=read.table("metaphlan_sp.txt", sep="\t", header=T, row.names = 1, check.names = F)
phenos=read.table("db3.txt", sep="\t", header =T, row.names = 1)
resection=read.table("resection.txt", sep="\t", header=T, row.names = 1)
```

##Clean datasets

```{r, eval=F}

tax_IBD=as.data.frame(tax_IBD)
tax_LLD=as.data.frame(tax_LLD)
tax_LLD_match=tax_LLD[rownames(LLD_phenos)]
tax_IBD_match=tax_IBD[rownames(IBD_phenos)]

LLD_sp=tax_LLD_match[grep("s__", rownames(tax_LLD_match)), ]
IBD_sp=tax_IBD_match[grep("s__", rownames(tax_IBD_match)), ]

LLD_sp=LLD_sp[!grepl("t__", rownames(LLD_sp)), ]
IBD_sp=IBD_sp[!grepl("t__", rownames(IBD_sp)), ]

row.names(LLD_sp)=gsub(".*s__","",row.names(LLD_sp))
row.names(IBD_sp)=gsub(".*s__","",row.names(IBD_sp))
row.names(meta_sp)=gsub(".*s__","",row.names(meta_sp))

```


# Get summary statistics functions from: https://github.com/WeersmaLabIBD/General_Tools/blob/master/Metadata_summary.R

```{r, eval=F}

phenos$calpro_cat="Low"
phenos[phenos$Lab1Calprotectin>200,]$calpro_cat="High"
phenos$colectomy_or_pouch="No"
phenos[phenos$Stoma=="Yes" | phenos$Colectomy=="Yes" ,]$colectomy_or_pouch="Yes"

wt=wilcox.test(log10(phenos$FISH_total)~phenos$colectomy_or_pouch, conf.int = T) #=>  p-value = 0.000485, CI= 0.2981498-1.0158741
abs(qnorm(zz$p.value)) / sqrt (127)

ggplot(phenos, aes(colectomy_or_pouch,log10(phenos$FISH_total) )) + geom_boxplot(outlier.size = NA) +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_bw() + xlab ("Colectomy or stoma") + ylab("log10(FISH counts)")

summary_statistics_metadata(phenos)
#remove stomas and colectomies (represent small intestine)
phenos2=subset(phenos,phenos$Stoma=="No")
phenos2=subset(phenos2,phenos2$Individual_ID!="p037")
phenos2=subset(phenos2,phenos2$Colectomy=="No")
phenos2$Stoma=NULL
phenos2$Colectomy=NULL
summary_statistics_metadata(phenos2)
phenos2$FISH_log=log(phenos2$FISH_total)
phenos2$RD_log=log(phenos2$PF_Reads)
meta_sp=meta_sp[rownames(phenos2)]
```

#Calculate Shannon index

```
phenos2$Shannon = diversity(t(meta_sp),index = "shannon")
phenos_1000IBD$Shannon = diversity(t(IBD_sp),index = "shannon")
phenos_LLD$Shannon = diversity(t(LLD_sp),index = "shannon")
```

# Test loads versus overall composition

```{r, eval=F}

phenos_v2=read.table("db3.2.txt", sep="\t", header =T, row.names = 1)
phenos_v3=subset(phenos_v2,phenos_v2$Stoma=="No")
phenos_v3=subset(phenos_v3,phenos_v3$Colectomy=="No")
phenos_v3$Stoma=NULL
phenos_v3$Colectomy=NULL
phenos_v3$FISH_log=log10(phenos_v3$FISH_total)
phenos_v3$RD_log=log10(phenos_v3$PF_Reads)
phenos_v3$PF_Reads=NULL
phenos_v3$FISH_total=NULL
phenos_v3$Colon=NULL
taxa_data=as.data.frame(t(meta_sp))
#Check that phenotypes and microbiome tables are concordant in the number of samples
taxa_data=taxa_data[rownames(phenos_v3),]
dist.matrix <- vegdist(taxa_data,method = "bray")
#Univariate variance analysis 
#Start an empty matrix
adonis_meta <- matrix(ncol = 7, nrow=ncol(phenos_v3))
#Calculate variance explained per phenotype after 10.000 permutations
for (i in 1:ncol(phenos_v3)){
  ad<-adonis(dist.matrix ~ phenos_v3[,i],permutations=10000)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
rownames(adonis_meta) = colnames(phenos_v3)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")

ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + ylab ("Explained variance") + xlab ("Factor")  + theme(text = element_text(size=8))

ad<-adonis(dist.matrix ~ DNA_conc+ilealvalveinsitu+PPI+VitaminB12+Min2Resection+IlealResection+Thiopurines+OnlyColonicDisease+Any_Resection+Surgery+ResectionColon+IlialResection+IndicationAbcess+FrequencyResection+FrequencySurgery+FrequencyResectionColon+FrequencyIlialresection+FrequencyAbcess+TimetoLastResection+TimetoLastColonicresection+TimetoLastIlealResection+TimetoLastAbcessIndication+MontrealL+MontrealB.1+MontrealL.1+BMI+Sex+BSEgroep+ChromeCreat1+ChromeUrine1+CRPDifferencePrePostRiboflavin+FFQAlcohol+FFQEiwitplantaardig+FFQkcal+FFQkJ+FFQKoolhydratentotaal+HBI+IBDQ1+IBDQ1Bowel+IBDQ1Emotional+IBDQ1Social+IBDQ1Systemic+Lab1ALAT+Lab1BSE+Lab1Creat+Lab1Hb+liqstool1+SmokingPackyears+RD_log+FISH_log, data=phenos_v3,permutations=10000)

aov_table <- ad$aov.tab
adonis_meta_all= as.data.frame(aov_table)
adonis_meta_all$FDR=p.adjust(adonis_meta_all$`Pr(>F)`, method = "BH")
adonis_meta_all$Signi="No"
adonis_meta_all$Signi[adonis_meta_all$FDR<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta_all)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
adonis_meta_all=head(adonis_meta_all,-2)
ggplot(adonis_meta_all, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + geom_bar(stat = "identity") + coord_flip() + theme_bw() + ylab ("Explained variance") + xlab ("Factor")  + theme(text = element_text(size=8))
```



Check and remove FISH outliers (No outliers present after removing stomas/colectomies)
----
  
```{r, eval=FALSE}
phenos_v4=phenos_v3[ phenos_v3$FISH_log > max(boxplot(phenos_v3$FISH_log)$out),]
phenos_vs_fish=phenos_v3
phenos_vs_fish=select(phenos_vs_fish,FISH_log,everything())
phenos_vs_fish=select(phenos_vs_fish,RD_log,everything())
phenos_vs_fish$calpro_cat="a_Low"
phenos_vs_fish[phenos_vs_fish$Lab1Calprotectin>200,]$calpro_cat="b_High"
my_results=matrix(nrow=100, ncol=4)
a=1
for (i in 3:ncol(phenos_vs_fish)){
  if (is.numeric(phenos_vs_fish[,i])){
    my_test=cor.test(phenos_vs_fish[,2], phenos_vs_fish[,i], method="spearman")
    my_results[a,4]=my_test$p.value
    my_results[a,2]=my_test$estimate
    my_results[a,1]=colnames(phenos_vs_fish)[i]
    a=a+1
  } else {
    my_test=pairwise.wilcox.test(phenos_vs_fish[,2], phenos_vs_fish[,i] ,p.adjust.method = "none")
    my_test2=melt(my_test$p.value)
    for (x in 1:nrow(my_test2)){
      my_results[a,4]=as.character(my_test2[x,3])
      my_results[a,3]=as.character(my_test2[x,2])
      my_results[a,2]=as.character(my_test2[x,1])
      my_results[a,1]=colnames(phenos_vs_fish)[i]
      a=a+1  
    }
  }
}
my_results2=as.data.frame(my_results)
my_results2=my_results2[c(1:77),]
my_results2=my_results2[-c(35),]
my_results2=my_results2[-c(38),]
my_results2$V4=as.numeric(as.character(my_results2$V4))
my_results2$FDR=p.adjust(my_results2$V4, method="BH")

ggplot(phenos_vs_fish, aes(calpro_cat,phenos_vs_fish$FISH_log )) + geom_boxplot(outlier.size = NA) +  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) + theme_bw() + xlab ("Fecal calprotectin levels") + ylab("log10(FISH counts)")
ggscatter(phenos_v3, x="DNA_conc", y="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 1, label.y = 10.6) + xlab("DNA concentration (log10)") + ylab ("Bacterial loads (log10(FISH))")
ggscatter(phenos_v3, x="log_cal", y="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 1.5, label.y = 10.6) + xlab("Fecal calprotectin levels") + ylab ("Bacterial loads (log10(FISH))")
ggscatter(phenos_v3, x="RD_log", y="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 6.6, label.y = 10.6) + xlab("Sequencing read depth (log10)") + ylab ("Bacterial loads (log10(FISH))")
ggscatter(phenos_v3, x="liqstool1", y="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 0, label.y = 10.6) + xlab("Number of liquid stool before sampling") + ylab ("Bacterial loads (log10(FISH))")

```


Repeat only for patients with a resection
----

```{r,eval = FALSE}
phenos_resec=subset(phenos_vs_fish, phenos_vs_fish$Any_Resection=="Yes")
phenos_resec$ilealvalve=0
phenos_resec[phenos_resec$ilealvalveinsitu=="yes",]$ilealvalve=1
phenos_resec2=phenos_resec[,c("RD_log","FISH_log" ,"FrequencyResection","FrequencySurgery", "FrequencyResectionColon", "FrequencyIlialresection","FrequencyAbcess","FrequencyFistula","TimetoLastResection","TimetoLastColonicresection","TimetoLastIlealResection","TimetoLastAbcessIndication","TimetoLastFistulaIndication","ilealvalve")]

my_results_resec=matrix(nrow=12, ncol=4)

a=1

for (i in 3:ncol(phenos_resec2)){
    my_test=cor.test(phenos_resec2[,2], phenos_resec2[,i], method="spearman")
    my_results_resec[a,4]=my_test$p.value
    my_results_resec[a,2]=my_test$estimate
    my_results_resec[a,1]=colnames(phenos_resec2)[i]
    a=a+1
}
my_results_resec=as.data.frame(my_results_resec)
my_results_resec$FDR=p.adjust(as.numeric(as.character(my_results_resec$V4)), method="BH")
my_results_resec$V3=NULL
ggscatter(phenos_resec, x="FrequencyIlialresection", y="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 0, label.y = 10.6) + xlab("Number of resections in the ileum") + ylab ("Bacterial loads (log10(FISH))")
```


Filter taxa present in at least 10 % of each of the 3 cohorts, and do inverse rank transformation
-----


```{r, eval=F}
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
tax_fish = t(meta_sp)
tax_LLD = t(LLD_sp)
tax_1000IBD = t(IBD_sp)

bugs2keep = intersect(colnames(tax_LLD)[colSums(tax_LLD>0)>0.1 * nrow(tax_LLD)],
                      intersect(colnames(tax_1000IBD)[colSums(tax_1000IBD>0)>0.1 * nrow(tax_1000IBD)],
                                colnames(tax_fish)[colSums(tax_fish>0)>0.1 * nrow(tax_fish)]))
bugs_fish = apply(tax_fish[,bugs2keep],2,invrank)[rownames(phenos_v3),]
bugs_1000IBD = apply(tax_1000IBD[,bugs2keep],2,invrank)[rownames(phenos_1000IBD),]
bugs_LLD = apply(tax_LLD[,bugs2keep],2,invrank)[rownames(phenos_LLD),]

```



#Plot correlation FISH vs RD per bacteria

```{r, eval=F}

phenos_subset=phenos2[,c("FISH_log","RD_log")]
# View(phenos_subset)
phenos_subset=merge(phenos_subset,bugs_fish, by="row.names")
row.names(phenos_subset)=phenos_subset$Row.names
phenos_subset$Row.names=NULL

species_tito=matrix( ncol=4,nrow=131)
c=1
 for (b in 3:ncol(phenos_subset)){
     
     #Correlation RD
     x1=cor.test(phenos_subset[,b],phenos_subset[,2], method="spearman")
     species_tito[c,1]=x1$p.value
     species_tito[c,2]=x1$estimate
     
     #Correlation FISH
     x1=cor.test(phenos_subset[,b],phenos_subset[,1], method="spearman")
     species_tito[c,3]=x1$p.value
     species_tito[c,4]=x1$estimate
     c=c+1
 }
species_tito=as.data.frame(species_tito) 
rownames(species_tito)=colnames(phenos_subset)[3:ncol(phenos_subset)]
colnames(species_tito)=c("RD_pvalue", "RD_beta", "FISH_pval", "FISH_beta")
mean_abundance=as.data.frame(colMeans(as.data.frame(t(meta_sp))))
colnames(mean_abundance)=("rel_abund")
species_tito=merge(species_tito,mean_abundance, by="row.names")
species_tito$fdr=p.adjust(species_tito$FISH_pval, method = "BH")


ggplot (species_tito, aes(FISH_beta, RD_beta))  + geom_point(aes(size=rel_abund))  + theme_bw() + geom_abline(intercept =0 , slope = 1, color="red") + xlim(c(-0.6,0.6)) + ylim (c(-0.6,0.6)) + xlab("Spearman coefficient bacteria ~ cell density (FISH)") + ylab ("Spearman coefficient bacteria ~ Sequencing depth")+ geom_text(aes(label=ifelse(abs(FISH_beta)>0.3,as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) 


#Plot also the relative abundance of FISH counts
ggscatter(phenos2, y="FISHFprau", x="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman", label.x = 1) + xlab("F.Prausnitzii (FISH rel.abundance)") + ylab ("Bacterial loads (log10(FISH))")
ggscatter(phenos2, y="FISHEcoli", x="FISH_log", add="reg.line", add.params = list(color="red", fill="lightgray"), conf.int = T) + stat_cor(method = "spearman") + ylab("E.coli (FISH rel.abundance)") + xlab ("Bacterial loads (log10(FISH))")

```


Prepare phenotypes
------

```
IBD_phenos=merge(IBD_phenos, resection, by="row.names")
rownames(IBD_phenos)=IBD_phenos$Row.names
IBD_phenos$Row.names=NULL

IBD_phenos$PPI2=2
IBD_phenos[ IBD_phenos$PPI=="Non_user", ]$PPI2=1
IBD_phenos$Sex2=2
IBD_phenos[ IBD_phenos$Sex=="male", ]$Sex2=1
IBD_phenos$ResectionIlealAny2=2
IBD_phenos$ResectionIlealAny[is.na(IBD_phenos$ResectionIlealAny)] <- "no"
IBD_phenos[ IBD_phenos$ResectionIlealAny=="no", ]$ResectionIlealAny2=1


LLD_phenos$PPI2=2
LLD_phenos[ LLD_phenos$PPI=="Non_user", ]$PPI2=1
LLD_phenos$Sex2=2
LLD_phenos[ LLD_phenos$Sex=="male", ]$Sex2=1

phenos2$PPI2=2
phenos2[ phenos2$PPI=="No", ]$PPI2=1
phenos2$Sex2=2
phenos2[ phenos2$Sex==0, ]$Sex2=1
phenos2$IlealResection2=2
phenos2[ phenos2$IlealResection=="No", ]$IlealResection2=1

phenos_1000IBD=IBD_phenos[,c( "Concentration", "Shannon", "Age", "Sex2", "PPI2","ResectionIlealAny2", "NumberOfResectionsAny")]
phenos_fish=phenos2[,c("FISH_log", "DNA_conc", "Shannon", "Age", "Sex2", "PPI2","IlealResection2", "FrequencyResection")]
phenos_LLD=LLD_phenos[,c("Concentration","Shannon", "Age", "Sex2", "PPI2")]

colnames(phenos_1000IBD)=c("DNA_con", "Shannon", "Age", "Sex","PPI", "ResectionIleal", "NumberResections")
colnames(phenos_fish)=c("FISH_log","DNA_con", "Shannon", "Age", "Sex","PPI", "ResectionIleal", "NumberResections")
colnames(phenos_LLD)=c("DNA_con", "Shannon", "Age", "Sex","PPI")
```




Predict loads
-------

```
phenos_LLD$DNA_con = scale(phenos_LLD$DNA_con)
phenos_LLD$Shannon = scale(phenos_LLD$Shannon)

phenos_1000IBD$DNA_con = scale(phenos_1000IBD$DNA_con)
phenos_1000IBD$Shannon = scale(phenos_1000IBD$Shannon)

phenos_fish$DNA_con = scale(phenos_fish$DNA_con)
phenos_fish$Shannon = scale(phenos_fish$Shannon)

fish_model = lm(FISH_log ~ Shannon + DNA_con, data = phenos_fish)
summary(fish_model)

fish_pred_train = predict(fish_model, newdata = phenos_fish)
fish_pred_LLD = predict(fish_model, newdata = phenos_LLD)
fish_pred_1000IBD = predict(fish_model, newdata = phenos_1000IBD)
fish_real_train = setNames(phenos_fish$FISH_log, row.names(phenos_fish))

```

## FISH-phenotype association in different cohorts

```{r,eval = F}
pheno_names = colnames(phenos_fish)
pheno_names = setdiff(pheno_names, c("DNA_conc", "FISH_log","Shannon"))
fish.pheno = foreach( i = pheno_names, .combine = rbind)%do%{
  cor.train.real = cor.test(phenos_fish$FISH_log, phenos_fish[,i],method = "spearman")
  cor.train.pred = cor.test(fish_pred_train, phenos_fish[,i],method = "spearman")
  cor.1000IBD.pred = cor.test(fish_pred_1000IBD, phenos_1000IBD[,i],method = "spearman")
  if(i %in% colnames(phenos_LLD)){
    cor.LLD.pred = cor.test(fish_pred_LLD, phenos_LLD[,i],method = "spearman")
  } else {
    cor.LLD.pred = list()
    cor.LLD.pred$estimate = NA
    cor.LLD.pred$p.value = NA
  }
  data.frame(phenotype = i,
             train.real.cor = cor.train.real$estimate,
             train.real.p = cor.train.real$p.value,
             train.pred.cor = cor.train.pred$estimate,
             train.pred.p = cor.train.pred$p.value,
             pred.1000IBD.pred.cor = cor.1000IBD.pred$estimate,
             pred.1000IBD.pred.p = cor.1000IBD.pred$p.value,
             LLD.pred.cor = cor.LLD.pred$estimate,
             LLD.pred.p = cor.LLD.pred$p.value
  )
  
}
fish.pheno

bugs_1000IBD=bugs_1000IBD[match(rownames(phenos_1000IBD), rownames(bugs_1000IBD)), ]
bugs_fish=bugs_fish[match(rownames(phenos_fish), rownames(bugs_fish)), ]
bugs_LLD=bugs_LLD[match(rownames(phenos_LLD), rownames(bugs_LLD)), ]


out.1000IBD = rcorr(as.matrix(cbind(phenos_1000IBD,bugs_1000IBD)),type = "spearman")
out.1000IBD$r = out.1000IBD$r[1:7,8:138]
out.1000IBD$P = out.1000IBD$P[1:7,8:138]

out.QMP = rcorr(as.matrix(cbind(phenos_fish[,-1],bugs_fish)),type = "spearman")
out.QMP$r = out.QMP$r[1:7,8:138]
out.QMP$P = out.QMP$P[1:7,8:138]

out.LLD = rcorr(as.matrix(cbind(phenos_LLD,bugs_LLD)),type = "spearman")
out.LLD$r = out.LLD$r[1:5,6:136]
out.LLD$P = out.LLD$P[1:5,6:136]

#constrained version 
out.constrained.1000IBD = list()
out.constrained.1000IBD$r = foreach(i = 1:ncol(bugs_1000IBD),.combine = cbind)%do%{
  foreach(j = 1:7,.combine = c) %do% {
  p1 = pcor(cbind(bugs_1000IBD[,i],phenos_1000IBD[,j],fish_pred_1000IBD),method ="spearman")
  p1$estimate[1,2]
}};colnames(out.constrained.1000IBD$r) = colnames(bugs_1000IBD);rownames(out.constrained.1000IBD$r) = colnames(phenos_1000IBD)

out.constrained.1000IBD$P = foreach(i = 1:ncol(bugs_1000IBD),.combine = cbind)%do%{
  foreach(j = 1:7,.combine = c) %do% {
  p1 = pcor(cbind(bugs_1000IBD[,i],phenos_1000IBD[,j],fish_pred_1000IBD),method ="spearman")
  p1$p.value[1,2]
}};colnames(out.constrained.1000IBD$P) = colnames(bugs_1000IBD);rownames(out.constrained.1000IBD$P) = colnames(phenos_1000IBD)


#LLD
out.constrained.LLD = list()
out.constrained.LLD$r = foreach(i = 1:ncol(bugs_LLD),.combine = cbind)%do%{
  foreach(j = 1:5,.combine = c) %do% {
  p1 = pcor(cbind(bugs_LLD[,i],phenos_LLD[,j],fish_pred_LLD),method ="spearman")
  p1$estimate[1,2]
}};colnames(out.constrained.LLD$r) = colnames(bugs_LLD);rownames(out.constrained.LLD$r) = colnames(phenos_LLD)

out.constrained.LLD$P = foreach(i = 1:ncol(bugs_LLD),.combine = cbind)%do%{
  foreach(j = 1:5,.combine = c) %do% {
  p1 = pcor(cbind(bugs_LLD[,i],phenos_LLD[,j],fish_pred_LLD),method ="spearman")
  p1$p.value[1,2]
}};colnames(out.constrained.LLD$P) = colnames(bugs_LLD);rownames(out.constrained.LLD$P) = colnames(phenos_LLD)


#train
out.constrained.QMP = list()
out.constrained.QMP$r = foreach(i = 1:ncol(bugs_fish),.combine = cbind)%do%{
  foreach(j = 2:8,.combine = c) %do% {
  p1 = pcor(cbind(bugs_fish[,i],phenos_fish[,j],fish_pred_train),method ="spearman")
  p1$estimate[1,2]
}};colnames(out.constrained.QMP$r) = colnames(bugs_LLD);rownames(out.constrained.QMP$r) = colnames(phenos_fish)[2:8]



out.constrained.QMP$P = foreach(i = 1:ncol(bugs_fish),.combine = cbind)%do%{
  foreach(j = 2:8,.combine = c) %do% {
  p1 = pcor(cbind(bugs_fish[,i],phenos_fish[,j],fish_pred_train),method ="spearman")
  p1$p.value[1,2]
}};colnames(out.constrained.QMP$P) = colnames(bugs_LLD);rownames(out.constrained.QMP$P) = colnames(phenos_fish)[2:8]


#Real

out.real.QMP = list()
out.real.QMP$r = foreach(i = 1:ncol(bugs_fish),.combine = cbind)%do%{
  foreach(j = 2:8,.combine = c) %do% {
  p1 = pcor(cbind(bugs_fish[,i],phenos_fish[,j],fish_real_train),method ="spearman")
  p1$estimate[1,2]
}};colnames(out.real.QMP$r) = colnames(bugs_LLD);rownames(out.real.QMP$r) = colnames(phenos_fish)[2:8]


out.real.QMP$P = foreach(i = 1:ncol(bugs_fish),.combine = cbind)%do%{
  foreach(j = 2:8,.combine = c) %do% {
  p1 = pcor(cbind(bugs_fish[,i],phenos_fish[,j],fish_real_train),method ="spearman")
  p1$p.value[1,2]
}};colnames(out.real.QMP$P) = colnames(bugs_LLD);rownames(out.real.QMP$P) = colnames(phenos_fish)[2:8]

test1=as.data.frame(t(out.real.QMP$P))
test1=test1[3:7]
colnames(test1)=c("Age_RU_real", "Sex_RU_real", "PPI_RU_real", "ResectionIleal_RU_real", "NumberResections_RU_real")
pvalues_df=test1

test1=as.data.frame(t(out.QMP$P))
test1=test1[3:7]
colnames(test1)=c("Age_RU", "Sex_RU", "PPI_RU", "ResectionIleal_RU", "NumberResections_RU")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL

test1=as.data.frame(t(out.constrained.QMP$P))
test1=test1[3:7]
colnames(test1)=c("Age_RU_pred", "Sex_RU_pred", "PPI_RU_pred", "ResectionIleal_RU_pred", "NumberResections_RU_pred")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL



test1=as.data.frame(t(out.1000IBD$P))
test1=test1[3:7]
colnames(test1)=c("Age_1000IBD", "Sex_1000IBD", "PPI_1000IBD", "ResectionIleal_1000IBD", "NumberResections_1000IBD")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL


test1=as.data.frame(t(out.constrained.1000IBD$P))
test1=test1[3:7]
colnames(test1)=c("Age_1000IBD_pred", "Sex_1000IBD_pred", "PPI_1000IBD_pred", "ResectionIleal_1000IBD_pred", "NumberResections_1000IBD_pred")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL

test1=as.data.frame(t(out.LLD$P))
test1=test1[3:5]
colnames(test1)=c("Age_LLD", "Sex_LLD", "PPI_LLD")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL

test1=as.data.frame(t(out.constrained.LLD$P))
test1=test1[3:5]
colnames(test1)=c("Age_LLD_pred", "Sex_LLD_pred", "PPI_LLD_pred")
pvalues_df=merge(pvalues_df, test1, by="row.names")
row.names(pvalues_df)=pvalues_df$Row.names
pvalues_df$Row.names=NULL

```

Create figures
====

```{r, eval=F}
pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$PPI_RU_real) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$PPI_RU) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$PPI_RU_real) > -log10( 0.05/131) & -log10(pvalues_df$PPI_RU) > -log10( 0.05/131),]$color="red"

p1= ggplot(pvalues_df, aes(-log10(PPI_RU), -log10(PPI_RU_real), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(PPI_RU)>-log10( 0.05/131) | -log10(PPI_RU_real)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")

pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$PPI_RU_pred) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$PPI_RU) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$PPI_RU_pred) > -log10( 0.05/131) & -log10(pvalues_df$PPI_RU) > -log10( 0.05/131),]$color="red"
 
p1_2= ggplot(pvalues_df, aes(-log10(PPI_RU), -log10(PPI_RU_pred), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(PPI_RU)>-log10( 0.05/131) | -log10(PPI_RU_pred)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed") 

pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$PPI_1000IBD_pred) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$PPI_1000IBD) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$PPI_1000IBD_pred) > -log10( 0.05/131) & -log10(pvalues_df$PPI_1000IBD) > -log10( 0.05/131),]$color="red"

p2= ggplot(pvalues_df, aes(-log10(PPI_1000IBD), -log10(PPI_1000IBD_pred), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(PPI_1000IBD)>-log10( 0.05/131) | -log10(PPI_1000IBD_pred)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")



pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$PPI_LLD_pred) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$PPI_LLD) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$PPI_LLD_pred) > -log10( 0.05/131) & -log10(pvalues_df$PPI_LLD) > -log10( 0.05/131),]$color="red"

p3= ggplot(pvalues_df, aes(-log10(PPI_LLD), -log10(PPI_LLD_pred), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(PPI_LLD)>-log10( 0.05/131) | -log10(PPI_LLD_pred)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")



pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU_real) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU_real) > -log10( 0.05/131) & -log10(pvalues_df$ResectionIleal_RU) > -log10( 0.05/131),]$color="red"

p4=ggplot(pvalues_df, aes(-log10(ResectionIleal_RU), -log10(ResectionIleal_RU_real), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(ResectionIleal_RU)>-log10( 0.05/131) | -log10(ResectionIleal_RU_real)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")


pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU_pred) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$ResectionIleal_RU_pred) > -log10( 0.05/131) & -log10(pvalues_df$ResectionIleal_RU) > -log10( 0.05/131),]$color="red"

p4_2=ggplot(pvalues_df, aes(-log10(ResectionIleal_RU), -log10(ResectionIleal_RU_pred), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(ResectionIleal_RU)>-log10( 0.05/131) | -log10(ResectionIleal_RU_pred)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")


pvalues_df$color="black"
pvalues_df[-log10(pvalues_df$ResectionIleal_1000IBD_pred) > -log10( 0.05/131) ,]$color="blue"
pvalues_df[-log10(pvalues_df$ResectionIleal_1000IBD) > -log10( 0.05/131),]$color="orange"
pvalues_df[-log10(pvalues_df$ResectionIleal_1000IBD_pred) > -log10( 0.05/131) & -log10(pvalues_df$ResectionIleal_1000IBD) > -log10( 0.05/131),]$color="red"

p5=ggplot(pvalues_df, aes(-log10(ResectionIleal_1000IBD), -log10(ResectionIleal_1000IBD_pred), label=rownames(pvalues_df))) + geom_point(colour=pvalues_df$color) + theme_bw() + geom_hline(yintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_vline(xintercept=-log10( 0.05/131), linetype="dashed", color = "red") + geom_text(aes(label=ifelse(-log10(ResectionIleal_1000IBD)>-log10( 0.05/131) | -log10(ResectionIleal_1000IBD_pred)>-log10( 0.05/131),as.character(rownames(pvalues_df)),'')),hjust=0,vjust=0) + xlim(0,20) + ylim(0,20) + geom_abline(intercept = 0, slope = 1, linetype="dashed")

(p1/p1_2/p2/p3)|(p4/p4_2/p5/plot_spacer())
```

Test p-value consistency before and after correcting for microbial loads per sample
====
```
diff_table=pvalues_df
diff_table[diff_table<0.05]=1
diff_table[diff_table!=1]=0
mcnemar.test(diff_table$Age_RU_real,diff_table$Age_RU)
mcnemar.test(diff_table$Age_LLD_pred,diff_table$Age_LLD)
mcnemar.test(diff_table$Age_IBD_pred,diff_table$Age_IBD)
mcnemar.test(diff_table$Age_1000IBD_pred,diff_table$Age_1000IBD)
mcnemar.test(diff_table$Sex_RU_real,diff_table$Sex_RU)
mcnemar.test(diff_table$ResectionIleal_RU_real,diff_table$ResectionIleal_RU)
mcnemar.test(diff_table$NumberResections_1000IBD,diff_table$NumberResections_1000IBD_pred)
mcnemar.test(diff_table$NumberResections_RU_real,diff_table$NumberResections_RU)
```
