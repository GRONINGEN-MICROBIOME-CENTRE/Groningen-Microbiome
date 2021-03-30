#Calculated using BiG-MAP: https://github.com/medema-group/BiG-MAP & https://github.com/victoriapascal/gutsmash/tree/gutsmash

#database: gutSMASH on a combined collection of 4.083 genomes from the Culturable Genome Reference (CGR), Human Microbiome Project (HMP) and complete genomes belonging to Clostridiales species. This resulted in 11.509 predicted metabolic gene clusters which we used as input for BiG-MAP.family


1.Calculate metabolomic gene clusters
----
  
 Input= Post QC paired fastq files

```{python}

#!/bin/bash


#SBATCH --job-name=B2
#SBATCH --error=B2.err
#SBATCH --output=B2.out
#SBATCH --mem=40gb
#SBATCH --time=23:50:00
#SBATCH --cpus-per-task=10



cd /groups/umcg-gastrocol/tmp04/for_Arnau_gutSMASH

export PATH="/groups/umcg-gastrocol/tmp04/for_Arnau_gutSMASH/my_miniconda/bin:$PATH"
conda activate /groups/umcg-gastrocol/tmp04/for_Arnau_gutSMASH/BiG-MAP/BM_process

ml Biopython/1.65-foss-2015b-Python-3.4.1
ml Bowtie2/2.3.4.1-foss-2015b
ml Python/3.6.3-foss-2015b

python3 ./BiG-MAP/src/BiG-MAP.map.py -th 10 -I1 ./b2/*1.fastq -I2 ./b2/*2.fastq -O ./results_b2 -P BiG-MAP.pickle


```


2. Prepare data 
----


```{r}
#Folder with different outputs form BiG-MAP (+GutSMASH) from different computational batches. 
#Folder contains BiG-MAP.map.results.ALL.csv file

#Input and merge different batches
setwd("~/Desktop/Metabolomics_v2/2.Input/BGC/")

file_list <- list.files()
flag=1
for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (flag==1){
    bgc <- read.table(file, header=TRUE,  check.names = F, sep=",", row.names = 1, stringsAsFactors = F)
    flag=5
  }
  
  # if the merged dataset does exist, append to it
  else {
    temp_bgc <-read.table(file, header=TRUE,  check.names = F, sep=",", row.names = 1, stringsAsFactors = F)
    bgc<-merge(bgc, temp_bgc, by="row.names", all.x=T)
    row.names(bgc)=bgc$Row.names
    bgc$Row.names=NULL
    rm(temp_bgc)
  }
  
}

cc_rd=read.table("~/Documents/Project_Metabolome/1.Input_files/IDs_read_depth.txt", header = T, row.names = "IDs_MGS")
cc_rd$row.names=NULL


#Get coverage of the core regions
bgc_coverage=select(bgc, contains(".corecov"))
colnames(bgc_coverage)=gsub(".corecov","",colnames(bgc_coverage))
bgc_coverage2=merge(cc_rd,t(bgc_coverage), by="row.names")
rownames(bgc_coverage2)=bgc_coverage2$PID
bgc_coverage2$PID=NULL
bgc_coverage2$Row.names=NULL
bgc_coverage2$IDs_MGS=NULL
bgc_coverage2$PF_RD=NULL
bgc_coverage3=subset(bgc_coverage2, row.names(bgc_coverage2) %in% row.names(my_phenos))


#Set a filtering

bgc_coverage4=bgc_coverage3
bgc_mean_coverage=data.frame(BCG=colnames(bgc_coverage4), Mean=colMeans(bgc_coverage4))
ggplot(bgc_mean_coverage, aes(reorder(BCG,-Mean), Mean )) + geom_bar(stat = "identity", fill="red4") + theme_classic() + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) + ylab ("Mean coverage")  + xlab ("Biosynthetic genes clusters")

select_coverage=droplevels(bgc_mean_coverage$BCG[bgc_mean_coverage$Mean>0.05])


#Collapse RPKM to pathways

bgc_rpkm=select(bgc, contains(".RPKM"))
colnames(bgc_rpkm)=gsub(".RPKM","",colnames(bgc_rpkm))
cc_rd$id2=rownames(cc_rd)
row.names(cc_rd)=cc_rd$IDs_MGS
cc_rd$row.names=NULL

bgc_rpkm2=merge(cc_rd,t(bgc_rpkm), by="row.names")
rownames(bgc_rpkm2)=bgc_rpkm2$id2
bgc_rpkm2$id2=NULL
bgc_rpkm2$Row.names=NULL
bgc_rpkm2$IDs_MGS=NULL
bgc_rpkm2$PF_RD=NULL
bgc_rpkm3=subset(bgc_rpkm2, row.names(bgc_rpkm2) %in% row.names(my_phenos))

#Subset MGC based on coverage

bgc_rpkm4=bgc_rpkm3[,colnames(bgc_rpkm3)%in%select_coverage]
bgc_rpkm4=as.data.frame(t(bgc_rpkm4))
bgc_rpkm4$type=str_split_fixed(row.names(bgc_rpkm4),"--", 3)[,2]
bgc_rpkm4= bgc_rpkm4 %>% group_by(type) %>% summarise_all(sum)
bgc_rpkm4=as.data.frame(bgc_rpkm4)
row.names(bgc_rpkm4)=make.names(bgc_rpkm4$type)
bgc_rpkm4$type=NULL
bgc_rpkm4=as.data.frame(t(bgc_rpkm4))


#Transform 
mgc_fil_trans=transform_and_filter_taxa(bgc_rpkm4, samples_row = T, method = "clr",missing_filter = 20)

```

3. Analyis MGC IBD vs Controls

```{r}

#Test differences in IBD

row.names(cc_rd)=cc_rd$PID
cc_rd$PID=NULL
#cc_rd$IDs_MGS=NULL

my_phenos=merge(my_phenos,cc_rd, by="row.names")
row.names(my_phenos)=my_phenos$Row.names
my_phenos$Row.names=NULL

case_control_mgc=merge(my_phenos,mgc_fil_trans, by="row.names")

rownames(case_control_mgc)=case_control_mgc$Row.names
case_control_mgc$Row.names=NULL
case_control_mgc$LC.COLUMN=NULL
case_control_mgc$Amount_sample_gram=NULL
case_control_mgc$metabolon_Month_in_freezer=NULL
#colnames(case_control_mgc)=make.names(colnames(case_control_mgc))




flag=1
for ( i in 1:11) {
  my_pheno=colnames(case_control_mgc)[i]
  if (my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="PF_RD" | my_pheno=="host_BMI" | my_pheno=="host_Sex" | my_pheno=="host_Age") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 12:ncol(case_control_mgc)){
      my_trait=colnames(case_control_mgc)[a]
      my_uni_test=case_control_mgc[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","PF_RD",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","clinical_BowelMovementADayDef", "host_BMI","PF_RD",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      #my_uni_test[,8]=as.numeric(as.character(my_uni_test[,8]))
      my_samples=nrow(my_uni_test)
      my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
        my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
        my_univariate_results=my_lm_pheno
        flag=5
      }
    }
  }
}

associations_ibd_mgc2=my_univariate_results
write.table(associations_ibd_mgc2, "~/Desktop/Metabolomics_v2/2.Input/BGC_with_IBD_NEW_collapsed.txt", sep = "\t", quote = F)

```


4. Test MGC vs metabolites in IBD and controls cohorts
----
  
  
```{r}

#####################################################
### TEST within IBD BGC - Metabolite associations ###
#####################################################
regressors=case_control_quantitative[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","ibd_IBD")]
regressors_IBD=subset(regressors,regressors$ibd_IBD=="IBD")
regressors_IBD$ibd_IBD=NULL
mtb=case_control_quantitative[,c(14:ncol(case_control_quantitative))]

ibd_test_mgc=merge(regressors_IBD,mgc_fil_trans, by="row.names")
rownames(ibd_test_mgc)=ibd_test_mgc$Row.names
ibd_test_mgc$Row.names=NULL

ibd_test_mgc=merge(ibd_test_mgc,mtb, by="row.names")
rownames(ibd_test_mgc)=ibd_test_mgc$Row.names
ibd_test_mgc$Row.names=NULL


#ibd_test_bgc=read.table("~/Desktop/Metabolomics_for_cluster/2.Input/intermed/ibd_BGC.txt", sep = "\t", header = T, row.names = 1)

flag=1
for ( i in 1:141) {
  my_pheno=colnames(ibd_test_mgc)[i]
  if (my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="LC.COLUMN" | my_pheno=="host_BMI"| my_pheno=="Amount_sample_gram"| my_pheno=="metabolon_Month_in_freezer" | my_pheno=="host_Sex" | my_pheno=="host_Age") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 142:ncol(ibd_test_mgc)){
      my_trait=colnames(ibd_test_mgc)[a]
      my_uni_test=ibd_test_mgc[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","clinical_BowelMovementADayDef", "host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      #my_uni_test[,8]=as.numeric(as.character(my_uni_test[,8]))
      my_samples=nrow(my_uni_test)
      my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
        my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
        my_univariate_results=my_lm_pheno
        flag=5
      }
    }
  }
}

associations_mgc_wIBD=my_univariate_results


associations_bgc_wIBD$FDR=p.adjust(associations_bgc_wIBD$`Pr(>|t|)`,method = "BH")
associations_bgc_wIBD$Bonferroni=p.adjust(associations_bgc_wIBD$`Pr(>|t|)`,method = "bonferroni")

write.table(associations_bgc_wIBD, "~/Desktop/Metabolomics_v2/3.Preliminary_results/7.Biosynthetic_gene_clusters_and_enzyme_comissions/within_IBD_associations_BGC_quantitative.txt", sep = "\t", quote = F)


##########################################################
### TEST within Controls BGC - Metabolite associations ###
##########################################################

regressors_CNT=subset(regressors,regressors$ibd_IBD!="IBD")
regressors_CNT$ibd_IBD=NULL

cnt_test_mgc=merge(regressors_CNT,mgc_fil_trans, by="row.names")
rownames(cnt_test_mgc)=cnt_test_mgc$Row.names
cnt_test_mgc$Row.names=NULL

cnt_test_mgc=merge(cnt_test_mgc,mtb, by="row.names")
rownames(cnt_test_mgc)=cnt_test_mgc$Row.names
cnt_test_mgc$Row.names=NULL

#cnt_test_bgc=read.table("~/Desktop/Metabolomics_for_cluster/2.Input/intermed/controls_BGC.txt", sep = "\t", header = T, row.names = 1)

flag=1
for ( i in 1:141) {
  my_pheno=colnames(cnt_test_mgc)[i]
  if (my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="LC.COLUMN" | my_pheno=="host_BMI"| my_pheno=="Amount_sample_gram"| my_pheno=="metabolon_Month_in_freezer" | my_pheno=="host_Sex" | my_pheno=="host_Age") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 142:ncol(cnt_test_mgc)){
      my_trait=colnames(cnt_test_mgc)[a]
      my_uni_test=cnt_test_mgc[,c("host_Sex","host_Age","clinical_BowelMovementADayDef","host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","clinical_BowelMovementADayDef", "host_BMI","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      #my_uni_test[,8]=as.numeric(as.character(my_uni_test[,8]))
      my_samples=nrow(my_uni_test)
      my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
        my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
        my_univariate_results=my_lm_pheno
        flag=5
      }
    }
  }
}

associations_mgc_CNT=my_univariate_results
#associations_controls=subset(associations_controls, associations_controls$phenotype!="BioMK_ChromograninA")
associations_bgc_CNT$FDR=p.adjust(associations_bgc_CNT$`Pr(>|t|)`,method = "BH")
associations_bgc_CNT$Bonferroni=p.adjust(associations_bgc_CNT$`Pr(>|t|)`,method = "bonferroni")

write.table(associations_bgc_CNT, "~/Desktop/Metabolomics_v2/3.Preliminary_results/7.Biosynthetic_gene_clusters_and_enzyme_comissions/Controls_associations_BGC_quantitative.txt", sep = "\t", quote = F)



################################################################
### TEST within IBD BGC - Metabolite associations prevalence ###
################################################################

regressors=ibd_test_phenos2[,c(1,2,3,4,6,7,8,34)]
mtb=ibd_test_phenos2[,c(345:ncol(ibd_test_phenos2))]

#colnames(bcg_fil_trans2)=make.names(colnames(bcg_fil_trans2))

ibd_test_phenos2=merge(regressors,bcg_fil_trans2, by="row.names")
rownames(ibd_test_phenos2)=ibd_test_phenos2$Row.names
ibd_test_phenos2$Row.names=NULL

ibd_test_phenos2=merge(ibd_test_phenos2,mtb, by="row.names")
rownames(ibd_test_phenos2)=ibd_test_phenos2$Row.names
ibd_test_phenos2$Row.names=NULL


#ibd_test_phenos2=read.table("~/Desktop/Metabolomics_for_cluster/2.Input/intermed/controls_BGC_prev2.txt", sep = "\t", header = T, row.names = 1)

ibd_test_phenos2$clinical_BowelMovementADayDef.1=NULL

flag=1
for ( i in 1:323) {
  my_pheno=colnames(ibd_test_phenos2)[i]
  if (my_pheno=="run_day_cat" | my_pheno=="clinical_BowelMovementADayDef" | my_pheno=="LC.COLUMN"| my_pheno=="Amount_sample_gram"| my_pheno=="metabolon_Month_in_freezer" | my_pheno=="host_Sex" | my_pheno=="host_Age" | my_pheno=="host_BMI") {
    print (paste("Skipping phenotype:",my_pheno))
  } else {
    print (paste("Testing:",my_pheno))
    for (a in 324:ncol(ibd_test_phenos2)){
      my_trait=colnames(ibd_test_phenos2)[a]
      my_uni_test=ibd_test_phenos2[,c("host_Sex","host_Age","run_day_cat","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","clinical_BowelMovementADayDef", "host_BMI",my_pheno,my_trait)]
      my_preds=c("host_Sex","host_Age","host_BMI", "run_day_cat","LC.COLUMN","Amount_sample_gram","metabolon_Month_in_freezer","clinical_BowelMovementADayDef",my_pheno)
      my_f=as.formula(paste(my_trait, paste(my_preds, collapse = " + "), sep = " ~ "))
      my_uni_test[,10]=as.numeric(as.character(my_uni_test[,10]))
      my_samples=nrow(my_uni_test)
      #my_lm=summary(lm(my_f,data = my_uni_test ))
      my_lm=summary(glm(my_f, family = binomial(link="logit"), data =my_uni_test))
      my_lm_coef=as.data.frame(my_lm$coefficients)
      my_lm_pheno=try(my_lm_coef[grep(my_pheno, rownames(my_lm_coef)), ])
      my_lm_pheno$metabolite=my_trait
      my_lm_pheno$phenotype=my_pheno
      my_lm_pheno$factor=rownames(my_lm_pheno)
      if (flag!=1){
        my_univariate_results=rbind(my_univariate_results,my_lm_pheno)
      }else{
        my_univariate_results=my_lm_pheno
        flag=5
      }
    }
  }
}

associations_BGC_ibd_prev=my_univariate_results
associations_BGC_ibd_prev$FDR=p.adjust(associations_BGC_ibd_prev$`Pr(>|z|)`,method = "BH")
associations_BGC_ibd_prev$Bonferroni=p.adjust(associations_BGC_ibd_prev$`Pr(>|z|)`,method = "bonferroni")
write.table(associations_BGC_ibd_prev, "~/Desktop/Metabolomics_v2/3.Preliminary_results/7.Biosynthetic_gene_clusters_and_enzyme_comissions/IBD_associations_BGC_prevalence.txt", sep = "\t", quote = F)

```








