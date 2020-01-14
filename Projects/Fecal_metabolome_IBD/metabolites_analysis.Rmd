---
title: "Metabolomics"
author: "Arnau Vich"
date: "01/2020"
output: html_document
---

Introduction
===

In this project we aim to investigate the relation between fecal metabolites and inflammatory bowel diseases phenotypes. The cohort consist of 750 samples of different patients, devided in 255 controls from the population cohort LifeLines (NL), 265 patients with Crohns disease and 198 patients with ulcerative colitis from the 1000IBD cohort (NL). 

Metabolites were measured using METABOLON platform combining untargeted metabolite detection with short chain fatty acids measurments. Untargeted metabolite detection consisted on an ultrahigh performance liquid chromatography tandem mass spectrocopy (UPLC-MS/MS): 2 acidic positive ion condition (C18 columns) + 1 Basic negative ion optimazed (C18 column) + 1 Basic negative ionization (Hilic column).

Metabolites will be combined with the following data layers: matching shot-gun sequencing metagenomics (from the same sample), host phenotypes at time of sampling (including disease phenotypes), diet (food frequency questionnaires) & host genetics (GSA and WES data). 

Load functions
===

```{r, warning=F,message=FALSE, echo=TRUE}
library(corrplot)
library(plotly)
library(pheatmap)
library(reshape2)
library(ggrepel)
library(dplyr)
library(eulerr)
library(vegan)
library(ape)
library(ggridges)
library(ggplot2)
library(psych)
library(ggord)
library(caret)
library(Rtsne)
library(heatmaply)
library(ggpubr)
library(ggbiplot)
library(coin)
library(compositions)
library(RColorBrewer)

#Invert rank transformation function. 
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
```

Import data
---

a) Metabolites: Here we use both raw values provided by Metabolon (metabolites AUC) and transformed values. Transformation consist of two main procedures: adjusting medians of each metabolite to 1 (making distributions more symetric) and imputing missing values (missing values are replaced by the minimum recorded value of each metabolite). Raw values will be used to calculate prevalence and filter metabolites with high percentange of NA's

b) Short Chain Fatty Acids (SFCA): Together with the untargetted metabolite assessment, SCFA concentrations were measured in the same fecal samples (μg/g)

c) Host phenotypes: Information about the host (age,sex,BMI,etc.). 

d) Experiment variables: Different technical variables from the untargetted metabolic experiments. We will use this information as a potential counfounding factors. 

e) Faecal characteristics: Two phenotypes that describe the intestinal transit time and the sample moisture: bowel movement a day and stool consistency (soft or hard). Bristol stool scale was not available for the patients with IBD, however, this charactersitics were estimated from the harvey bradshaw score (HBI) and other clinical scores.

f) Food: Diet information captured from food frequency questionnaires (FFQs). Each element has been adjusted for nutrient density.

g) Metagenomics taxa: MetaPhlan's profiles (v2.7)

```{r, warning, message=F, echo=T}
#transformed metabolites (raw data transformed to adjust distrubution, medians = 1)
all_metabolites <- read.delim("~/Documents/Project_Metabolome/Takeda/all_metabolites2.txt", row.names=1)

#raw values (AUC of peaks)
all_metabolites_raw <- read.delim("~/Documents/Project_Metabolome/Takeda/all_metabolites_raw2.txt", row.names=1)

#Metabolites annotation
annot <- read.delim("~/Documents/Project_Metabolome/Takeda/info_metabolites_v2.txt", row.names=1)
annot$SUPER.PATHWAY.1=NULL

#phenotypes

#SCFA concentrations
phenos_scfa <- read.delim("~/Documents/Project_Metabolome/input_scfa_short.txt")

#Selection of host phenotypes
phenos=read.delim("~/Documents/Project_Metabolome/phenos_sel.txt", row.names = 1)

#info on techinical parameters of the metabolite measuments (columns used, day, etc.)
tech=read.delim("~/Documents/Project_Metabolome/tech.txt")

#Stool frequency & binary BSS (soft/hard)
stool=read.delim("~/Documents/Project_Metabolome/intestinal_info.txt")

#Diet (Food frequency questionnaires adjusted for nutrient density)

food <- read.delim("~/Documents/Project_Metabolome/Foodgroups.txt")

# Medication

# Storage

#translate id's 
IDs_LLD <- read.delim("~/Documents/Project_Metabolome/Takeda/IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("~/Documents/Project_Metabolome/Takeda/IDs_IBD.txt", header=FALSE)

#Metagenomics taxa

IBD_taxa <- read.delim("~/Documents/Project_Metabolome/taxa/IBD_taxonomy_unstrat_clean.txt", row.names=1)
LLD_taxa <- read.delim("~/Documents/Project_Metabolome/taxa/LLD_taxonomy_unstrat_clean.txt", row.names=1)

```

Adjust all ids
---

Translate all id's to the same id (easier to deal with multiple tables)

```{r, warning=F,message=FALSE, echo=TRUE}

rownames(tech)=tech$CLIENT_IDENTIFIER
tech$CLIENT_IDENTIFIER=NULL

#Change Metabolon ids to UMCG ids to connect later to phenotypes 
all_raw=as.data.frame(t(all_metabolites_raw))
all=as.data.frame(t(all_metabolites))

IDs_IBD=IDs_IBD[,c(1,5)]
IDs_LLD$V2=NULL

# [IMPORTANT] Table with matching ids between UMCG and Metabolon
IDs=data.frame(rbind(as.matrix(IDs_IBD), as.matrix(IDs_LLD)))
rownames(IDs)=IDs$V1
IDs$V1=NULL

# Merge to replace the ids Metabolon => UMCG
all_new_ID=merge(IDs,all, by="row.names")
all_new_ID_raw=merge(IDs,all_raw, by="row.names")
tech_ID=merge(IDs,tech, by="row.names")
rownames(tech_ID)=tech_ID$V5
tech_ID$V5=NULL
rownames(all_new_ID_raw)=all_new_ID_raw$V5
all_new_ID_raw$Row.names=NULL
all_new_ID_raw$V5=NULL
row.names(all_new_ID)=all_new_ID$V5
all_new_ID$V5=NULL
all_new_ID$Row.names=NULL

#Keep only the samples with metabolite measurments
stool=stool[complete.cases(stool[ ,1]),]
rownames(stool)=stool$UMCGIBDResearchIDorLLDeepID
stool$UMCGIBDResearchIDorLLDeepID=NULL
IDs2=IDs
rownames(IDs2)=IDs$V5
stool_ID=merge(IDs2,stool, by="row.names")
stool_ID$V5=NULL
rownames(stool_ID)=stool_ID$Row.names
stool_ID$Row.names=NULL

#Change also filter the diet file



#Merge host phenotypes with technical information
phenos2=merge(tech_ID,phenos, by="row.names")
rownames(phenos2)=phenos2$Row.names
phenos2$Row.names=NULL

#Merge also the stool phenotypes

phenos3=merge(phenos2,stool_ID, by="row.names")
rownames(phenos3)=phenos3$Row.names
phenos3$Row.names=NULL
```


Data preparation and descriptive statsitics
===


Normalize data for testing and split per phenotype (Control / CD / UC)
---

Use inverse rank transformation

```{r, warning=F,message=FALSE, echo=TRUE}
#Transform to log 10 ( metabolites values)
#mclean=log10(all_new_ID)
#Tranform to inv-rank values

CD=rownames(phenos)[phenos$DiagnosisCurrent=="CD"]
UC=rownames(phenos)[phenos$DiagnosisCurrent=="UC"]
CT=rownames(phenos)[phenos$DiagnosisCurrent=="generalpopulation"]

all_new_ID_raw_cd=all_new_ID_raw[rownames(all_new_ID_raw)%in%CD,]
all_new_ID_raw_uc=all_new_ID_raw[rownames(all_new_ID_raw)%in%UC,]
all_new_ID_raw_ct=all_new_ID_raw[rownames(all_new_ID_raw)%in%CT,]

all_new_ID_cd=all_new_ID[rownames(all_new_ID)%in%CD,]
all_new_ID_uc=all_new_ID[rownames(all_new_ID)%in%UC,]
all_new_ID_ct=all_new_ID[rownames(all_new_ID)%in%CT,]

norm_metab=as.data.frame(apply(all_new_ID,2,invrank))
norm_metab_CD=as.data.frame(apply(all_new_ID_cd,2,invrank))
norm_metab_UC=as.data.frame(apply(all_new_ID_uc,2,invrank))
norm_metab_CT=as.data.frame(apply(all_new_ID_ct,2,invrank))

```


Summary statistics for each metabolite
---

Important. Summary statistics are not calculated if the variance of a metabolite in the transformed table (median=1 & imputed missing values) is equal to zero, meaning that all the values are equal. This is because either all values were missing or only few samples (1 or 2) had a non-zero values that was used for imputation of the rest of the samples. 

- Calculate the number of zeros and non-zeros values on metabolites raw values

- Calculate the min, max, mean and median values on metabolites raw values

- Calculate potential outliers values (>3*IQR) on metabolites raw values 

- Calculate min and max values on transformed values (median adjusted)

- Normality test (Shapiro-test) on inverse-rank transformed data. 

```{r, warning=F,message=FALSE, echo=TRUE}

summary_met=matrix(nrow=ncol(all_new_ID_raw), ncol=49)
for (i in 1:ncol(all_new_ID_raw)){
  summary_met[i,2]=sum(is.na(all_new_ID_raw[,i]))
  summary_met[i,3]=sum(!is.na(all_new_ID_raw[,i]))
  if (var(all_new_ID[,i]) == 0){
  summary_met[i,4]="NA"
  summary_met[i,5]="NA"
  summary_met[i,6]="NA"
  summary_met[i,7]="NA"
  summary_met[i,8]="NA"
  summary_met[i,9]="NA"
  summary_met[i,10]="NA"
  summary_met[i,11]="NA"
  summary_met[i,12]="NA"
  summary_met[i,13]="NA"
  }
  else{
  summary_met[i,4]=min(all_new_ID_raw[,i], na.rm = T)
  summary_met[i,5]=max(all_new_ID_raw[,i], na.rm = T)
  summary_met[i,6]=mean(all_new_ID_raw[,i], na.rm = T)
  summary_met[i,7]=median(all_new_ID_raw[,i], na.rm = T)
  summary_met[i,8]=sd(all_new_ID_raw[,i], na.rm = T)
  summary_met[i,9]=length(boxplot(all_new_ID_raw[,i], plot=FALSE, range = 3)$out)
  summary_met[i,10]=min(all_new_ID[,i], na.rm = T)
  summary_met[i,11]=max(all_new_ID[,i], na.rm = T)
  summary_met[i,12]=sd(all_new_ID[,i], na.rm = T)
  nor=shapiro.test(norm_metab[,i])
  summary_met[i,13]=nor$p.value>=0.05
  }
#Controls
  summary_met[i,14]=sum(is.na(all_new_ID_raw_ct[,i]))
  summary_met[i,15]=sum(!is.na(all_new_ID_raw_ct[,i]))
  if (var(all_new_ID_ct[,i]) == 0){
  summary_met[i,16]="NA"
  summary_met[i,17]="NA"
  summary_met[i,18]="NA"
  summary_met[i,19]="NA"
  summary_met[i,20]="NA"
  summary_met[i,21]="NA"
  summary_met[i,22]="NA"
  summary_met[i,23]="NA"
  summary_met[i,24]="NA"
  summary_met[i,25]="NA"
  }
  else{
  summary_met[i,16]=min(all_new_ID_raw_ct[,i], na.rm = T)
  summary_met[i,17]=max(all_new_ID_raw_ct[,i], na.rm = T)
  summary_met[i,18]=mean(all_new_ID_raw_ct[,i], na.rm = T)
  summary_met[i,19]=median(all_new_ID_raw_ct[,i], na.rm = T)
  summary_met[i,20]=sd(all_new_ID_raw_ct[,i], na.rm = T)
  summary_met[i,21]=length(boxplot(all_new_ID_raw_ct[,i], plot=FALSE, range = 3)$out)
  summary_met[i,22]=min(all_new_ID_ct[,i], na.rm = T)
  summary_met[i,23]=max(all_new_ID_ct[,i], na.rm = T)
  summary_met[i,24]=sd(all_new_ID_ct[,i], na.rm = T)
  nor=shapiro.test(norm_metab_CT[,i])
  summary_met[i,25]=nor$p.value>=0.05
  }
#CD
  summary_met[i,26]=sum(is.na(all_new_ID_raw_cd[,i]))
  summary_met[i,27]=sum(!is.na(all_new_ID_raw_cd[,i]))
  if (var(all_new_ID_cd[,i]) == 0){
  summary_met[i,28]="NA"
  summary_met[i,29]="NA"
  summary_met[i,30]="NA"
  summary_met[i,31]="NA"
  summary_met[i,32]="NA"
  summary_met[i,33]="NA"
  summary_met[i,34]="NA"
  summary_met[i,35]="NA"
  summary_met[i,36]="NA"
  summary_met[i,37]="NA"
  }
  else{
  summary_met[i,28]=min(all_new_ID_raw_cd[,i], na.rm = T)
  summary_met[i,29]=max(all_new_ID_raw_cd[,i], na.rm = T)
  summary_met[i,30]=mean(all_new_ID_raw_cd[,i], na.rm = T)
  summary_met[i,31]=median(all_new_ID_raw_cd[,i], na.rm = T)
  summary_met[i,32]=sd(all_new_ID_raw_cd[,i], na.rm = T)
  summary_met[i,33]=length(boxplot(all_new_ID_raw_cd[,i], plot=FALSE, range = 3)$out)
  summary_met[i,34]=min(all_new_ID_cd[,i], na.rm = T)
  summary_met[i,35]=max(all_new_ID_cd[,i], na.rm = T)
  summary_met[i,36]=sd(all_new_ID_cd[,i], na.rm = T)
  nor=shapiro.test(norm_metab_CD[,i])
  summary_met[i,37]=nor$p.value>=0.05
  }
#UC
  summary_met[i,38]=sum(is.na(all_new_ID_raw_uc[,i]))
  summary_met[i,39]=sum(!is.na(all_new_ID_raw_uc[,i]))
  if (var(all_new_ID_uc[,i]) == 0){
  summary_met[i,40]="NA"
  summary_met[i,41]="NA"
  summary_met[i,42]="NA"
  summary_met[i,43]="NA"
  summary_met[i,44]="NA"
  summary_met[i,45]="NA"
  summary_met[i,46]="NA"
  summary_met[i,47]="NA"
  summary_met[i,48]="NA"
  summary_met[i,49]="NA"
  }
  else{
  summary_met[i,40]=min(all_new_ID_raw_uc[,i], na.rm = T)
  summary_met[i,41]=max(all_new_ID_raw_uc[,i], na.rm = T)
  summary_met[i,42]=mean(all_new_ID_raw_uc[,i], na.rm = T)
  summary_met[i,43]=median(all_new_ID_raw_uc[,i], na.rm = T)
  summary_met[i,44]=sd(all_new_ID_raw_uc[,i], na.rm = T)
  summary_met[i,45]=length(boxplot(all_new_ID_raw_uc[,i], plot=FALSE, range = 3)$out)
  summary_met[i,46]=min(all_new_ID_uc[,i], na.rm = T)
  summary_met[i,47]=max(all_new_ID_uc[,i], na.rm = T)
  summary_met[i,48]=sd(all_new_ID_uc[,i], na.rm = T)
  nor=shapiro.test(norm_metab_UC[,i])
  summary_met[i,49]=nor$p.value>=0.05
  }
}
summary_met[,1]=colnames(all_new_ID_raw)
summary_met=as.data.frame(summary_met)
colnames(summary_met)=c("Metabolite", "NAs","Non_NAs", "Min", "Max", "Mean", "Median","SD", "Outliers_x3IQR", "Min_trans_data", "Max_trans_data", "SD_trans_data","Normal_distrib","NAs_CT","Non_NAs_CT", "Min_CT", "Max_CT", "Mean_CT", "Median_CT","SD_CT", "Outliers_x3IQR_CT", "Min_trans_data_CT", "Max_trans_data_CT","SD_CT", "Normal_distrib_CT", "NAs_CD","Non_NAs_CD", "Min_CD", "Max_CD", "Mean_CD", "Median_CD", "SD_CD","Outliers_x3IQR_CD", "Min_trans_data_CD", "Max_trans_data_CD", "SD_CD", "Normal_distrib_CD","NAs_UC","Non_NAs_UC", "Min_UC", "Max_UC", "Mean_UC", "Median_UC", "SD_UC","Outliers_x3IQR_UC", "Min_trans_data_UC", "Max_trans_data_UC", "SD_UC", "Normal_distrib_UC")

#summary_met$perc_missing=((as.numeric(as.character(summary_met$NAs))/nrow(summary_met))*100)
summary_short=summary_met[,c("Metabolite","NAs","NAs_CT","NAs_CD", "NAs_UC")]
summary_short$missing=((as.numeric(as.character(summary_short$NAs))/nrow(all_new_ID))*100)
summary_short$missing_CT=((as.numeric(as.character(summary_short$NAs_CT))/nrow(all_new_ID_ct))*100)
summary_short$missing_CD=((as.numeric(as.character(summary_short$NAs_CD))/nrow(all_new_ID_cd))*100)
summary_short$missing_UC=((as.numeric(as.character(summary_short$NAs_UC))/nrow(all_new_ID_raw_uc))*100)
summary_short$NAs=NULL
summary_short$NAs_CT=NULL
summary_short$NAs_CD=NULL
summary_short$NAs_UC=NULL
summary_short2=melt(summary_short)
```

Distribution of the missigness

```{r, warning=F,message=FALSE, echo=F}
aa=ggplot(summary_short, aes(missing))+ geom_histogram(color="darkblue", fill="lightblue",bins=10) + theme_classic() + xlab("% missing samples") + ylab("Number of metabolites")

ggplotly(aa)
```


Distribution of the missigness per cohort (CT=controls, CD=Crohns, UC=Ulcerative colitis)

```{r, warning=F,message=FALSE, echo=F}
ggplot(summary_short2, aes(x=value,y=variable, fill=variable))+ geom_density_ridges_gradient(stat = "binline", bins = 20, scale=3) + theme_classic() + xlab("% missing samples") + ylab("Number of metabolites")
```

Same but different visualization

```{r, warning=F,message=FALSE, echo=F}
bb=ggplot(summary_short, aes(reorder(Metabolite, -missing),missing))+ geom_bar(color="darkblue", fill="lightblue", stat = "identity") + theme_classic() + ylab("% missingness") + xlab("Metabolites") + theme(axis.text.x=element_blank())

ggplotly(bb)
```


Summary statistics per sample
---


```{r, warning=F,message=FALSE, echo=TRUE}

all_new_ID_raw2=as.data.frame(t(all_new_ID_raw))
summary_sample=matrix(nrow=ncol(all_new_ID_raw2), ncol=3)
for (i in 1:ncol(all_new_ID_raw2)){
  summary_sample[i,2]=sum(is.na(all_new_ID_raw2[,i]))
  summary_sample[i,3]=sum(!is.na(all_new_ID_raw2[,i]))
}
summary_sample[,1]=colnames(all_new_ID_raw2)
summary_sample=as.data.frame(summary_sample)
colnames(summary_sample)=c("Metabolite", "NAs","Non_NAs")
summary_sample$perc_missing=((as.numeric(as.character(summary_sample$NAs))/nrow(all_new_ID_raw2))*100)
rownames(summary_sample)=summary_sample$Metabolite
summary_sample$Metabolite=NULL
summary_sample2=merge(summary_sample,phenos,by="row.names")
```

Plot representing the distribution of metabolites based on its missigness 

```{r, warning=F,message=FALSE, echo=F}

cc=ggplot(summary_sample2, aes(x=perc_missing, fill=DiagnosisCurrent))+ geom_histogram(color="black" ) + theme_classic() + xlab("Percentage of samples missing") + ylab("Number of metabolites")

ggplotly(cc)

```


Remove metabolites absent in more than 85% of the data per each cohort
---

```{r, warning=F,message=FALSE, echo=TRUE}
#to_remove=subset(summary_met, summary_met$perc_missing>85)
to_remove_90=subset(summary_short,missing_CT>90 | missing_CD>90 | missing_UC>90)
removeCT=subset(summary_short,missing_CT>90)
removeCT=factor(removeCT$Metabolite)
removeCD=subset(summary_short,missing_CD>90)
removeCD=factor(removeCD$Metabolite)
removeUC=subset(summary_short,missing_UC>90)
removeUC=factor(removeUC$Metabolite)
```

Venn diagram metabolites missing >90% per cohort

```{r, warning=F,message=FALSE, echo=F}
plot(euler(list(A = removeCT, B = removeCD, C = removeUC)),quantities=T, label=c("Controls","CD", "UC"))
```

```{r, warning=F,message=FALSE, echo=TRUE}
#to_remove=subset(summary_met, summary_met$perc_missing>85)
to_remove_85=subset(summary_short,missing_CT>85 | missing_CD>85 | missing_UC>85)
removeCT=subset(summary_short,missing_CT>85)
removeCT=factor(removeCT$Metabolite)
removeCD=subset(summary_short,missing_CD>85)
removeCD=factor(removeCD$Metabolite)
removeUC=subset(summary_short,missing_UC>85)
removeUC=factor(removeUC$Metabolite)

```

Venn diagram metabolites missing >85% per cohort

```{r, warning=F,message=FALSE, echo=F}
plot(euler(list(A = removeCT, B = removeCD, C = removeUC)),quantities=T, label=c("Controls","CD", "UC"))
```


Removing metabolites in the normalized metabolite tables and merge phenotypes

```{r, warning=F,message=FALSE, echo=T}

# Select metabolites to remove
to_remove=subset(summary_short,missing_CT>85 & missing_CD>85 & missing_UC>85)
to_remove2=factor(to_remove$Metabolite)
`%ni%` <- Negate(`%in%`)

#Filter/remove selected metabolites
norm_CD_filt=subset(norm_metab_CD,select = names(norm_metab_CD) %ni% to_remove2)
norm_CT_filt=subset(norm_metab_CT,select = names(norm_metab_CT) %ni% to_remove2)
norm_UC_filt=subset(norm_metab_UC,select = names(norm_metab_UC) %ni% to_remove2)
norm_filt=subset(norm_metab,select = names(norm_metab) %ni% to_remove2)
annot_v2=subset(annot,rownames(annot) %in% names(norm_filt))

#Adjust diagnosis in the phenotype file

phenos3$Diag="Controls"
phenos3[phenos3$DiagnosisCurrent=="CD",]$Diag="CD"
phenos3[phenos3$DiagnosisCurrent=="UC",]$Diag="UC"
#phenos3[phenos3$DiagnosisCurrent=="ReconsideringDiagnosis",]$Diag="IBDU"
#phenos3[phenos3$DiagnosisCurrent=="IBDI",]$Diag="IBDU"

rownames(phenos_scfa)=phenos_scfa$Row.names
phenos_scfa$Diag="Controls"
phenos_scfa[phenos_scfa$DiagnosisCurrent=="CD",]$Diag="CD"
phenos_scfa[phenos_scfa$DiagnosisCurrent=="UC",]$Diag="UC"
phenos_scfa[phenos_scfa$DiagnosisCurrent=="ReconsideringDiagnosis",]$Diag="IBDU"
phenos_scfa[phenos_scfa$DiagnosisCurrent=="IBDI",]$Diag="IBDU"

#Merge phenotyes
CD_f_t=merge(phenos3,norm_CD_filt,by="row.names")
CT_f_t=merge(phenos3,norm_CT_filt,by="row.names")
UC_f_t=merge(phenos3,norm_UC_filt,by="row.names")
all_f_t=merge(phenos3,norm_filt,by="row.names")

```


Explore the correlations and clustering between metabolites
----

Correlations using Spearman correlation. 

```{r, warning=F,message=FALSE, echo=TRUE}

corCD=cor(norm_CD_filt,method = "spearman")
corCT=cor(norm_CT_filt,method = "spearman")
corUC=cor(norm_UC_filt,method = "spearman")

annot_v3=annot_v2
annot_v3$SUB.PATHWAY=NULL
```

Correlation plots, metabolites are ordered based on the hierarchical clustering


Metabolites in population controls

```{r, warning=F,message=FALSE, echo=F, fig.width=13, fig.height=10}
heatmaply_cor(cor(norm_CT_filt, method = "spearman"), showticklabels = c(FALSE, FALSE), main ="Metabolites in population controls", plot_method = 'plotly', scale="none", row_side_colors =annot_v3)
```


Metabolites in patients with CD


```{r, warning=F,message=FALSE, echo=F, fig.width=13, fig.height=10}
heatmaply_cor(cor(norm_CD_filt, method = "spearman"), showticklabels = c(FALSE, FALSE), main ="Metabolites in patients with CD", plot_method = 'plotly', scale="none", row_side_colors =annot_v3)
```


Metabolites in patients with UC


```{r, warning=F,message=FALSE, echo=F, fig.width=13, fig.height=10}
heatmaply_cor(cor(norm_UC_filt, method = "spearman"), showticklabels = c(FALSE, FALSE), main ="Metabolites in patients with UC", plot_method = 'plotly', scale="none", row_side_colors =annot_v3, margins=NA)
```


Heatmap including metabolites and few clinical phenotypes
---

[Subjective interpretation]

Looking at the top (horitzontal) dendogram we can see 2~4 clusters (left to right):

  1) A group of mainly patients with IBD with an enrichment of patients with CD and without ileocecal valve (due to resection). We see an enrichent of metabolites that cluster in the middle of the heatmap.
  
  2) Small cluster of patients with (mainly) UC and enrichment for stoma characterized by a lower abundance of metabolites clustering on the top part of the heatmap and higher abundance of those metabolites clustering in the middle part. 
  
  3) Clusters 3 & 4 are a mixture of patients with IBD and controls with the 4th cluster with a slightly enrichment of patients with UC and a lower abundance of the metabolites distributed in the lower section of the heatmap. 
  
  
  Diagnose: UC [n=198] / CD [n=265] / Controls [n=287]
  
  Stoma: No [n=684] / Yes [n=66]
  
  Ileocecal valve in-situ: Yes [n=331] / no [n=159] (Only cases!)
  
  Active diseases: Active [n=122] / NotActive [n=366] (Only cases! Described as combination of markers: Calprotectin levels + HBI/SSCAI)
  
  Binary bristol stool index: soft (n=207) vs hard (n=420)
  
  Bowel Movement a day: min=0 (n=57) - max=15 (n=3) / median=2 / mean=2.2 / missing=9
  
```{r, warning=F,message=FALSE, echo=TRUE}

#annotation_sel=mclean_p[,c("Diag","Box", "LC.COLUMN", "SPL.AMNT.GRAM", "RUN.DAY")]
#annotation_sel3=phenos_scfa[,c("Diag","X2.Methylbutyric.acid","Acetic.acid","Butyric.acid","Hexanoic.acid","Isobutyric.acid","Isovaleric.acid","Propionic.acid","Valeric.acid")]
#annotation_sel3[,2:9]=log10(annotation_sel3[,2:9])
annot_v4=all_f_t[,c("Row.names","Stoma","IleocecalValveInSitu","ActiveDisease","Diag","BinaryBSS","BowelMovementADayDef")]
rownames(annot_v4)=annot_v4$Row.names
annot_v4$Row.names=NULL
```


```{r, warning=F,message=FALSE, echo=F, fig.height=15, fig.align='center',fig.width=10}
pheatmap(t(norm_filt), show_colnames = F,show_rownames = F, scale = "none", annotation_row = annot_v3, annotation_col = annot_v4, fontsize = 7, annotation_colors = list ( Diag=c("CD"="firebrick2", "UC"="blue2", "Controls"="white"),ActiveDisease=c("Active"="firebrick2","NotActive"="white"),IleocecalValveInSitu=c("no"="firebrick2","yes"="white"), Stoma=c("Yes"="red","No"="white"),BinaryBSS=c("Hard"="gray79","Soft"="brown"), SuperPathway=c("Amino Acid"="lightblue", "Carbohydrate"="red", "Cofactors and Vitamins" = "gold2", "Energy"="gray79", "Lipid"="darkviolet", "Nucleotide" = "orange", "Partially Characterized Molecules" ="turquoise", "Peptide" = "limegreen", "Unknown"="black", "Xenobiotics"="tan1" )))
```


Explorative correlation (co-abundance / co-occurance) of taxa and metabolites
---

Prepare taxonomic data: keep only species level, remove taxa present <20% of the samples, transform data using arsin square-root transformation. In total we keep 749 samples and 139 bacterial species in this step.

One IBD sample (G105251=>IBDFEC0435) does not have MGS data
 
```{r, warning=F,message=FALSE, echo=T}
#Clean names
LLD_sp=LLD_taxa[grep("s__", rownames(LLD_taxa)), ]
IBD_sp=IBD_taxa[grep("s__", rownames(IBD_taxa)), ]
LLD_sp=LLD_sp[!grepl("t__", rownames(LLD_sp)), ]
IBD_sp=IBD_sp[!grepl("t__", rownames(IBD_sp)), ]
row.names(LLD_sp)=gsub(".*s__","",row.names(LLD_sp))
row.names(IBD_sp)=gsub(".*s__","",row.names(IBD_sp))

#Merge
taxa=as.data.frame(t(merge(IBD_sp,LLD_sp, by="row.names")))
taxa=merge(IBD_sp,LLD_sp, by="row.names")
row.names(taxa)=taxa$Row.names
taxa$Row.names=NULL
taxa=as.data.frame(t(taxa))

#Match samples
taxa2=subset(taxa,rownames(taxa) %in% rownames(norm_filt))

#Filter present in at least 20% of the samples
taxa2=taxa2/100
taxa2[taxa2=="NA"]=0
taxa2[taxa2=="NaN"]=0
taxa_filt=taxa2[,((colSums(taxa2 !=0) / nrow(taxa2)) *100 )>20]
taxa_filt=asin(sqrt(taxa_filt))

taxa_and_mb=merge(taxa_filt,norm_filt,by="row.names")
rownames(taxa_and_mb)=taxa_and_mb$Row.names
taxa_and_mb$Row.names=NULL

#Split per phenotype
rownames(CD_f_t)=CD_f_t$Row.names
CD_f_t$Row.names=NULL

rownames(CT_f_t)=CT_f_t$Row.names
CT_f_t$Row.names=NULL

rownames(UC_f_t)=UC_f_t$Row.names
UC_f_t$Row.names=NULL

tmb_CD=subset(taxa_and_mb,rownames(taxa_and_mb) %in% rownames(CD_f_t))
tmb_UC=subset(taxa_and_mb,rownames(taxa_and_mb) %in% rownames(UC_f_t))
tmb_CT=subset(taxa_and_mb,rownames(taxa_and_mb) %in% rownames(CT_f_t))
```



Spearman correlation between bacteria abundances and metabolite values. Later we can consider more elegant methods.

Correlations are performed per cohort. Combining disease with control cohorts could influence the correlations (disease co-occurring events)


```{r, warning=F,message=FALSE, echo=T}

#Simple loop to correlate bateria with metabolites

pvs_cd=matrix(nrow=116, ncol=1438)
rhos_cd=matrix(nrow=116, ncol=1438)
x=0
for (i in 1:116){
  x=x+1
  y=0
  for (b in 117:ncol(tmb_CD)){
    y=y+1
    a=cor.test(tmb_CD[,i],tmb_CD[,b],method = "spearman")
    pvs_cd[x,y]=a$p.value
    rhos_cd[x,y]=a$estimate
  }
}
rownames(pvs_cd)=colnames(tmb_CD)[1:116]
rownames(rhos_cd)=colnames(tmb_CD)[1:116]
colnames(pvs_cd)=colnames(tmb_CD)[117:ncol(tmb_CD)]
colnames(rhos_cd)=colnames(tmb_CD)[117:ncol(tmb_CD)]

mpvs_cd=melt(pvs_cd)
mrhos_cd=melt(rhos_cd)
mpvs_cd$rho=mrhos_cd$value
mpvs_cd$fdr=p.adjust(mpvs_cd$value,method = "bonferroni")
#mlite=subset(mpvs,mpvs$fdr<0.05)


##UC


pvs_uc=matrix(nrow=116, ncol=1438)
rhos_uc=matrix(nrow=116, ncol=1438)
x=0
for (i in 1:116){
  x=x+1
  y=0
  for (b in 117:ncol(tmb_UC)){
    y=y+1
    a=cor.test(tmb_UC[,i],tmb_UC[,b],method = "spearman")
    pvs_uc[x,y]=a$p.value
    rhos_uc[x,y]=a$estimate
  }
}
rownames(pvs_uc)=colnames(tmb_UC)[1:116]
rownames(rhos_uc)=colnames(tmb_UC)[1:116]
colnames(pvs_uc)=colnames(tmb_UC)[117:ncol(tmb_UC)]
colnames(rhos_uc)=colnames(tmb_UC)[117:ncol(tmb_UC)]

mpvs_uc=melt(pvs_uc)
mrhos_uc=melt(rhos_uc)
mpvs_uc$rho=mrhos_uc$value
mpvs_uc$fdr=p.adjust(mpvs_uc$value,method = "bonferroni")
#mlite=subset(mpvs,mpvs$fdr<0.05)

#Controls


pvs_ct=matrix(nrow=116, ncol=1438)
rhos_ct=matrix(nrow=116, ncol=1438)
x=0
for (i in 1:116){
  x=x+1
  y=0
  for (b in 117:ncol(tmb_CT)){
    y=y+1
    a=cor.test(tmb_CT[,i],tmb_CT[,b],method = "spearman")
    pvs_ct[x,y]=a$p.value
    rhos_ct[x,y]=a$estimate
  }
}
rownames(pvs_ct)=colnames(tmb_CT)[1:116]
rownames(rhos_ct)=colnames(tmb_CT)[1:116]
colnames(pvs_ct)=colnames(tmb_CT)[117:ncol(tmb_CT)]
colnames(rhos_ct)=colnames(tmb_CT)[117:ncol(tmb_CT)]

mpvs_ct=melt(pvs_ct)
mrhos_ct=melt(rhos_ct)
mpvs_ct$rho=mrhos_ct$value
mpvs_ct$fdr=p.adjust(mpvs_ct$value,method = "bonferroni")
#mlite=subset(mpvs,mpvs$fdr<0.05)


```


Correlation metabolites-microbes in population controls

```{r, warning=F,message=FALSE, echo=F, fig.height=15, fig.align='center',fig.width=10}

ab=ggplot(data = mpvs_ct, aes(Var2, Var1, fill = rho))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Spearman\nCorrelation\nControls") +theme_minimal()+ theme(axis.text.x=element_blank(),axis.text.y = element_text(size=5)) + ylab("Taxa") +xlab("Metabolite")

ggplotly(ab)
```


Correlation metabolites-microbes in patients with Crohn's diseases

```{r, warning=F,message=FALSE, echo=F, fig.height=15, fig.align='center',fig.width=10}
ac=ggplot(data = mpvs_cd, aes(Var2, Var1, fill = rho))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Spearman\nCorrelation\nCD cohort") +theme_minimal()+ theme(axis.text.x=element_blank(),axis.text.y = element_text(size=5)) + ylab("Taxa") +xlab("Metabolite")

ggplotly(ac)
```


Correlation metabolites-microbes in patients with Ulcerative colitis. 


```{r, warning=F,message=FALSE, echo=F, fig.height=15, fig.align='center',fig.width=10}
ad=ggplot(data = mpvs_uc, aes(Var2, Var1, fill = rho))+ geom_tile(color = "white")+ scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-0.7,0.7), space = "Lab", name="Spearman\nCorrelation\nUC cohort") +theme_minimal()+ theme(axis.text.x=element_blank(),axis.text.y = element_text(size=5)) + ylab("Taxa") +xlab("Metabolite")
ggplotly(ad)
```



Preliminar differential abundance analyses: diseases and host phenotypes
---
  
First we conduct analysis to identify disease trend.

Important! Notice that here we don't take into account any counfunding or correction factors.

However, we remove patients with pouch or stomas (n=66) due to their dysbiotic profile seen at taxa and metabolic levels.

Here we simply look at the differential abundance of metabolic features comparing them between groups. We perform a Wilcoxon-test (non-parametric t-test).

Fold changes are estimated from the raw values provided by Metabolon (AUC) [log2(mean cases / mean controls)].
  
```{r, warning=F,message=FALSE, echo=T}
#mclean_raw=mclean_pre[,!colnames(mclean_pre) %in% to_remove]
m_raw_filt=subset(all_new_ID_raw,select = names(all_new_ID_raw) %ni% to_remove2)
m_raw_filt_p=merge(phenos3,m_raw_filt,by="row.names")

CD_controls=subset(all_f_t,all_f_t$Diag!="UC")

#Remove pouches or stomas
CD_controls=subset(CD_controls,CD_controls$Stoma=="No")
CD_controls_raw=subset(m_raw_filt_p,m_raw_filt_p$Diag!="UC")

UC_controls=subset(all_f_t,all_f_t$Diag!="CD")
UC_controls=subset(UC_controls,UC_controls$Stoma=="No")
UC_controls_raw=subset(m_raw_filt_p,all_f_t$Diag!="CD")


# ncol-45 (45 = number of phenotypes in the input dataframe)
cd_cnt=matrix(nrow = ncol(CD_controls)-45, ncol = 13)
x=0

for (i in 46:ncol(CD_controls)){
  x=x+1
  #cd_cnt[x,2]=median(CD_controls[,i][CD_controls$Diag=="CD"])
  cd_cnt[x,2]=median(CD_controls_raw[,i][CD_controls_raw$Diag=="CD"], na.rm = T) #2
  cd_cnt[x,3]=median(UC_controls_raw[,i][UC_controls_raw$Diag=="UC"], na.rm = T) #3
  #cd_cnt[x,3]=median(CD_controls[,i][CD_controls$Diag!="CD"])
  cd_cnt[x,4]=median(CD_controls_raw[,i][CD_controls_raw$Diag!="CD"], na.rm = T) #4
  #cd_cnt[x,4]=mean(CD_controls[,i][CD_controls$Diag=="CD"])
  cd_cnt[x,5]=mean(CD_controls_raw[,i][CD_controls_raw$Diag=="CD"], na.rm = T) #5
  cd_cnt[x,6]=mean(UC_controls_raw[,i][UC_controls_raw$Diag=="UC"], na.rm = T) #6
  #cd_cnt[x,5]=mean(CD_controls[,i][CD_controls$Diag!="CD"])
  cd_cnt[x,7]=mean(CD_controls_raw[,i][CD_controls_raw$Diag!="CD"], na.rm = T) #7
  #Wilcoxon effect size
  
  #Fold change
  #CD
  cd_cnt[x,8]=as.numeric(cd_cnt[x,5])/as.numeric(cd_cnt[x,7])
  cd_cnt[x,9]=log2(as.numeric(as.character(cd_cnt[x,8])))
  #UC
  cd_cnt[x,10]=as.numeric(cd_cnt[x,6])/as.numeric(cd_cnt[x,7])
  cd_cnt[x,11]=log2(as.numeric(as.character(cd_cnt[x,10])))
  
  #Wilcoxon test
  #CD
  asoc=wilcox.test(CD_controls[,i][CD_controls$Diag=="CD"],CD_controls[,i][CD_controls$Diag!="CD"])
  cd_cnt[x,1]=colnames(CD_controls)[i]
  cd_cnt[x,12]=asoc$p.value
  #UC
  asoc=wilcox.test(UC_controls[,i][UC_controls$Diag=="UC"],UC_controls[,i][UC_controls$Diag!="CD"])
  #cd_cnt[x,1]=colnames(CD_controls)[i]
  cd_cnt[x,13]=asoc$p.value
}

cd_cnt=as.data.frame(cd_cnt)

colnames(cd_cnt)=c("Metabolite", "Median_CD", "Median_UC","Median_Controls", "Mean_CD","Mean_UC", "Mean_Controls", "Fold_change_CD",  "Log2_Fold_change_CD","Fold_change_UC","Log2_Fold_change_UC", "p_value_CD", "p_value_UC")

cd_cnt$FDR_CD=p.adjust(as.numeric(as.character(cd_cnt$p_value_CD)), method = "bonferroni")
cd_cnt$FDR_UC=p.adjust(as.numeric(as.character(cd_cnt$p_value_UC)), method = "bonferroni")
cd_cnt$log10_pvalue_CD=-log10(cd_cnt$FDR_CD)
cd_cnt$log10_pvalue_UC=-log10(cd_cnt$FDR_UC)

cd_cnt2=cd_cnt
#ibd_cnt2=ibd_cnt[complete.cases(ibd_cnt$FDR),]

cd_cnt2$sig_cd="ns"
cd_cnt2$sig_uc="ns"
cd_cnt2[cd_cnt2$FDR_CD<0.05,]$sig_cd="sug"
cd_cnt2[cd_cnt2$FDR_CD<0.05 & as.numeric(as.character(cd_cnt2$Log2_Fold_change_CD))>1,]$sig_cd="up"
cd_cnt2[cd_cnt2$FDR_CD<0.05 & as.numeric(as.character(cd_cnt2$Log2_Fold_change_CD))< -1,]$sig_cd="down"
cd_cnt2[cd_cnt2$FDR_UC<0.05,]$sig_uc="sug"
cd_cnt2[cd_cnt2$FDR_UC<0.05 & as.numeric(as.character(cd_cnt2$Log2_Fold_change_UC))>1,]$sig_uc="up"
cd_cnt2[cd_cnt2$FDR_UC<0.05 & as.numeric(as.character(cd_cnt2$Log2_Fold_change_UC))< -1,]$sig_uc="down"
```

```{r, warning=F,message=FALSE, echo=T}
table(cd_cnt2$sig_cd)
table(cd_cnt2$sig_uc)
```


[WARNING, this is just an univariate analysis]

Using the above thresholds we found that almost half of the metabolites are differentially abundant on patients with CD (787/1438) and ~10% in patients with UC (158/1438). 


Volcano plot of differentially abundant metabolites when comparing controls to patients with CD


```{r, warning=F,message=FALSE, echo=F, fig.align='center'}
ba=ggplot(cd_cnt2,aes(x=as.numeric(as.character(Log2_Fold_change_CD)), y=log10_pvalue_CD, colour=sig_cd,text=Metabolite)) + geom_point() +
   ggtitle("Metabolites abundance CD vs Controls") +
   xlab("log2 fold change") + 
   ylab("-log10 adjusted p-value") +
   theme_bw() +
   scale_color_manual(values = c( "blue","black","gray38","coral1")) + xlim (-6,6) + geom_hline(yintercept = 1.30, linetype = "dashed", color="purple") +
   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color="purple") 
ggplotly(ba)
```



Volcano plot of differentially abundant metabolites when comparing controls to patients with UC


```{r, warning=F,message=FALSE, echo=F, fig.align='center'}
bb=ggplot(cd_cnt2,aes(x=as.numeric(as.character(Log2_Fold_change_UC)), y=log10_pvalue_UC, colour=sig_uc,text=Metabolite)) + geom_point() +
   ggtitle("Metabolites abundance UC vs Controls") +
   xlab("log2 fold change") + 
   ylab("-log10 adjusted p-value") +
   theme_bw() +
   scale_color_manual(values = c( "blue","black","gray38","coral1")) + xlim (-6,6) + geom_hline(yintercept = 1.30, linetype = "dashed", color="purple") +
   geom_vline(xintercept = c(-1, 1), linetype = "dashed", color="purple")
ggplotly(bb)
```





