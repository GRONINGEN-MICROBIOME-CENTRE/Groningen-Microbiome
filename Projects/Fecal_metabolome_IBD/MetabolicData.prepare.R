# =============================================================================================
#  The first script is for explore metabolic data (based on Arnau)
# =============================================================================================
library(corrplot)
library(plotly)
library(pheatmap)
library(NbClust)
library(reshape2)
library(ggrepel)
library(dplyr)
library(eulerr)
library(vegan)
library(ape)
library(ggridges)
library (UpSetR)
library(ggplot2)
library(psych)
library(ggord)
library(mice)
library(caret)
library(factoextra)
library(Rtsne)
library(heatmaply)
library(ggpubr)
library(ggbiplot)
library(coin)
library(DT)
library(pvclust)
library(compositions)
library(RColorBrewer)

#Inverse rank transformation function. 
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
`%ni%` <- Negate(`%in%`)

#transformed metabolites (raw data transformed to adjust distrubution, medians = 1)
all_metabolites <- read.delim("./InputFiles/all_metabolites2.txt", row.names=1)
#raw values (AUC of peaks)
all_metabolites_raw <- read.delim("./InputFiles/all_metabolites_raw2.txt", row.names=1)
#Metabolites annotation
annot <- read.delim("./InputFiles/info_metabolites_v2.txt", row.names=1)
annot$SUPER.PATHWAY.1=NULL
#Phenotypes
phenos_ibd=read.delim("./InputFiles/Merged_IBD_phenos_clean_v2.txt",row.names = 1)
phenos_lld=read.delim("./InputFiles/Merged_LLD_phenos_clean_v2.txt", row.names=1)
cc_pheno=read.delim("./InputFiles/mini_pheno.txt", row.names=1)
# Storage
#translate id's 
IDs_LLD <- read.delim("./InputFiles/IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("./InputFiles/IDs_IBD.txt", header=FALSE)

# scale and impute
metabolite_test=as.data.frame(t(all_metabolites_raw))
#scale date, mean =0 , sd =1
scaled.test <- scale(metabolite_test)
#apply log10 transformation + scale
metabolite_test_log=log10(metabolite_test)
scaled_log_metabolites <- scale(metabolite_test_log)
#Impute na's


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
rownames(all_new_ID_raw)=all_new_ID_raw$V5
all_new_ID_raw$Row.names=NULL
all_new_ID_raw$V5=NULL
row.names(all_new_ID)=all_new_ID$V5
all_new_ID$V5=NULL
all_new_ID$Row.names=NULL

remove_pouch_stoma=rownames(phenos_ibd)[phenos_ibd$ibd_CurrentStomaOrPouch=="yes"]
#Keep data with stoma
all_new_ID_with_stoma=all_new_ID
all_new_ID_raw_with_stoma=all_new_ID_raw
#Remove
all_new_ID=all_new_ID[row.names(all_new_ID)%ni%remove_pouch_stoma,]
all_new_ID_raw=all_new_ID_raw[row.names(all_new_ID_raw)%ni%remove_pouch_stoma,]

#Subset disease cohort
CD=rownames(phenos_ibd)[phenos_ibd$ibd_Diagnosis=="CD"]
UC=rownames(phenos_ibd)[phenos_ibd$ibd_Diagnosis=="UC"]
CT=rownames(phenos_lld)
all_new_ID_raw_cd=all_new_ID_raw[rownames(all_new_ID_raw)%in%CD,]
all_new_ID_raw_uc=all_new_ID_raw[rownames(all_new_ID_raw)%in%UC,]
all_new_ID_raw_ct=all_new_ID_raw[rownames(all_new_ID_raw)%in%CT,]
all_new_ID_cd=all_new_ID[rownames(all_new_ID)%in%CD,]
all_new_ID_uc=all_new_ID[rownames(all_new_ID)%in%UC,]
all_new_ID_ct=all_new_ID[rownames(all_new_ID)%in%CT,]
##Choose method for metabolite transformation
#Tranform to inv-rank values
norm_metab=as.data.frame(apply(all_new_ID,2,invrank))
norm_metab_CD=as.data.frame(apply(all_new_ID_cd,2,invrank))
norm_metab_UC=as.data.frame(apply(all_new_ID_uc,2,invrank))
norm_metab_CT=as.data.frame(apply(all_new_ID_ct,2,invrank))


# summary
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

#Calculate presence
summary_short=summary_met[,c("Metabolite", "Non_NAs")]
summary_short2=summary_met[,c("Metabolite","Non_NAs","Non_NAs_CT","Non_NAs_CD", "Non_NAs_UC")]
summary_short2$prevalence=((as.numeric(as.character(summary_short2$Non_NAs))/nrow(all_new_ID))*100)
summary_short2$prevalence_CT=((as.numeric(as.character(summary_short2$Non_NAs_CT))/nrow(all_new_ID_ct))*100)
summary_short2$prevalence_CD=((as.numeric(as.character(summary_short2$Non_NAs_CD))/nrow(all_new_ID_cd))*100)
summary_short2$prevalence_UC=((as.numeric(as.character(summary_short2$Non_NAs_UC))/nrow(all_new_ID_raw_uc))*100)
summary_short2$Non_NAs=NULL
summary_short2$Non_NAs_CT=NULL
summary_short2$Non_NAs_CD=NULL
summary_short2$Non_NAs_UC=NULL
summary_short2[,-1] <-round(summary_short2[,-1],0)
summary_short3=melt(summary_short2)

to_remove_90=subset(summary_short2,prevalence_CT<10 | prevalence_CD<10 | prevalence_UC<10)
keepCT=subset(summary_short2,prevalence_CT>10)
keepCT=factor(keepCT$Metabolite)
keepCD=subset(summary_short2,prevalence_CD>10)
keepCD=factor(keepCD$Metabolite)
keepUC=subset(summary_short2,prevalence_UC>10)
keepUC=factor(keepUC$Metabolite)

to_remove=subset(summary_short2,prevalence_CT<10 & prevalence_CD<10 & prevalence_UC<10)
to_remove1=to_remove
rownames(to_remove1)=to_remove1$Metabolite
to_remove1=merge(annot,to_remove1,by="row.names")
to_remove2=factor(to_remove$Metabolite)
datatable(to_remove1)

#Filter/remove selected metabolites
norm_CD_filt=subset(norm_metab_CD,select = names(norm_metab_CD) %ni% to_remove2)
norm_CT_filt=subset(norm_metab_CT,select = names(norm_metab_CT) %ni% to_remove2)
norm_UC_filt=subset(norm_metab_UC,select = names(norm_metab_UC) %ni% to_remove2)
norm_filt=subset(norm_metab,select = names(norm_metab) %ni% to_remove2)
annot_v2=subset(annot,rownames(annot) %in% names(norm_filt))
#Merge phenotyes
CD_f_t=merge(cc_pheno,norm_CD_filt,by="row.names")
CT_f_t=merge(cc_pheno,norm_CT_filt,by="row.names")
UC_f_t=merge(cc_pheno,norm_UC_filt,by="row.names")
all_f_t=merge(cc_pheno,norm_filt,by="row.names")

# =============================================================================================
#  Next script is for preparing metabolic data for mmQTL
# =============================================================================================
# prepare covariates
Covariate_CD=phenos_ibd[,c("host_Age_sampling","host_BMI","host_Sex")]
Covariate_CD=Covariate_CD[rownames(Covariate_CD) %in% rownames(norm_CD_filt),]
Covariate_UC=phenos_ibd[,c("host_Age_sampling","host_BMI","host_Sex")]
Covariate_UC=Covariate_UC[rownames(Covariate_UC) %in% rownames(norm_UC_filt),]
Covariate_CT=phenos_lld[,c("host_age","host_BMI","host_Sex")]
Covariate_CT=Covariate_CT[rownames(Covariate_CT) %in% rownames(norm_CT_filt),]
Covariate_CD$host_Sex=as.character(Covariate_CD$host_Sex)
Covariate_UC$host_Sex=as.character(Covariate_UC$host_Sex)
Covariate_CD$host_Sex[Covariate_CD$host_Sex=="female"]=1
Covariate_UC$host_Sex[Covariate_UC$host_Sex=="female"]=1
Covariate_CD$host_Sex[Covariate_CD$host_Sex=="male"]=2
Covariate_UC$host_Sex[Covariate_UC$host_Sex=="male"]=2
Covariate_CD$host_Sex=as.numeric(Covariate_CD$host_Sex)
Covariate_UC$host_Sex=as.numeric(Covariate_UC$host_Sex)

colnames(Covariate_CT)=c("host_Age_sampling","host_BMI","host_Sex")
Covariate_all=rbind(Covariate_CD,Covariate_UC)
Covariate_all=rbind(Covariate_all,Covariate_CT)
  
# imputate phenotype using median value
Covariate_CD=apply(Covariate_CD,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
Covariate_CD=as.data.frame(Covariate_CD)
Covariate_UC=apply(Covariate_UC,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
Covariate_UC=as.data.frame(Covariate_UC)
Covariate_CT=apply(Covariate_CT,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
Covariate_CT=as.data.frame(Covariate_CT)
Covariate_all=apply(Covariate_all,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
Covariate_all=as.data.frame(Covariate_all)

# re-order samples
Covariate_CD=Covariate_CD[order(rownames(Covariate_CD)),]
Covariate_UC=Covariate_UC[order(rownames(Covariate_UC)),]
Covariate_CT=Covariate_CT[order(rownames(Covariate_CT)),]
Covariate_all=Covariate_all[order(rownames(Covariate_all)),]
norm_CD_filt=norm_CD_filt[order(rownames(norm_CD_filt)),]
norm_UC_filt=norm_UC_filt[order(rownames(norm_UC_filt)),]
norm_CT_filt=norm_CT_filt[order(rownames(norm_CT_filt)),]
norm_filt=norm_filt[rownames(norm_filt) %in% rownames(Covariate_all),]
norm_filt=norm_filt[order(rownames(norm_filt)),]

# correct for age, gender and BMI
norm_CD_filt_correct = apply(norm_CD_filt,2,function(x){
  
  x = resid(lm(x ~ .,data = Covariate_CD))

  return(x)
})
norm_CD_filt_correct=as.data.frame(norm_CD_filt_correct)
norm_UC_filt_correct = apply(norm_UC_filt,2,function(x){
  
  x = resid(lm(x ~ .,data = Covariate_UC))
  
  return(x)
})
norm_UC_filt_correct=as.data.frame(norm_UC_filt_correct)
norm_CT_filt_correct = apply(norm_CT_filt,2,function(x){
  
  x = resid(lm(x ~ .,data = Covariate_CT))
  
  return(x)
})
norm_CT_filt_correct=as.data.frame(norm_CT_filt_correct)
norm_filt_correct = apply(norm_filt,2,function(x){
  
  x = resid(lm(x ~ .,data = Covariate_all))
  
  return(x)
})
norm_filt_correct=as.data.frame(norm_filt_correct)

# write out phenotype files
norm_CD_filt_correct = as.data.frame(t(norm_CD_filt_correct))
norm_CD_filt_correct = cbind(rownames(norm_CD_filt_correct),norm_CD_filt_correct)
colnames(norm_CD_filt_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(norm_CD_filt_correct),
                   Gene = rownames(norm_CD_filt_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(norm_CD_filt_correct, file = "CD_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "CD_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
write.table(colnames(norm_CD_filt_correct),file = "CD.list.txt",sep="\t",row.names=F,quote = F)

norm_UC_filt_correct = as.data.frame(t(norm_UC_filt_correct))
norm_UC_filt_correct = cbind(rownames(norm_UC_filt_correct),norm_UC_filt_correct)
colnames(norm_UC_filt_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(norm_UC_filt_correct),
                   Gene = rownames(norm_UC_filt_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(norm_UC_filt_correct, file = "UC_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "UC_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)
write.table(colnames(norm_UC_filt_correct),file = "UC.list.txt",sep="\t",row.names=F,quote = F)

norm_CT_filt_correct = as.data.frame(t(norm_CT_filt_correct))
norm_CT_filt_correct = cbind(rownames(norm_CT_filt_correct),norm_CT_filt_correct)
colnames(norm_CT_filt_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(norm_CT_filt_correct),
                   Gene = rownames(norm_CT_filt_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(norm_CT_filt_correct, file = "CT_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "CT_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)

norm_filt_correct = as.data.frame(t(norm_filt_correct))
norm_filt_correct = cbind(rownames(norm_filt_correct),norm_filt_correct)
colnames(norm_filt_correct)[1] = "-"
annot = data.frame(platform = "RDP",
                   HT12v3.ArrayAddress = rownames(norm_filt_correct),
                   Gene = rownames(norm_filt_correct),
                   Chr = 4,
                   ChrStart = 1000,
                   ChrEnd = 1000)
colnames(annot)[2] = "HT12v3-ArrayAddress"
write.table(norm_filt_correct, file = "All_numeric.metabolic.txt",sep="\t",row.names = F,quote = F)
write.table(annot,file = "All_numeric.metabolic.txt.annot",sep="\t",row.names=F,quote = F)

