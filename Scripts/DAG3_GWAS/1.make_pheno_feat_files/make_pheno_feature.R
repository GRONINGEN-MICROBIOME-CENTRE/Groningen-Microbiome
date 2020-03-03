###################################
### build microbiome features-phenotypes file
### version 2.0
### date: 26-17-2020
### Author:EALM
###################################
##New

############----------------Packages and formulas---------------################
library(data.table)
library(tidyverse)
library(gridExtra)

### formulas

#####WersmaLab function for normalizing microbiome
normalize = function(data,samples.in.rows = F,to.abundance = T,transform = "log",move.log = "min"){
  if (samples.in.rows == F){
    data = t(data)
    data[is.na(data)] = 0
  } else {
    data= as.matrix(data)
  }
  if(to.abundance == TRUE){
    data = sweep(data,1,rowSums(data),"/")
    data[is.nan(data)] = 0
  }
  if(!is.element(transform,c("log","asinsqrt","none"))) stop("Wrong transformation! Use 'log', 'asinsqrt' or 'none'")
  if (transform == "log"){
    if(move.log == "min"){
      min = as.numeric(min(data[data>0],na.rm = T))
      data = log(data+min)
      if(any(data<0,na.rm =T)) {data = data - min(data)}
    } else if(is.numeric(move.log)){
      data = log(data + move.log)
      if(any(data<0,na.rm =T)) {data = data - min(data)}
    } else{
      stop("please define the value to prevent -Inf from log transformation. can be numeric or 'min'")
    }
  } else if (transform == "asinsqrt" & to.abundance ==T) {
    data =asin(sqrt(data))
  } else if (transform == "asinsqrt" & any(data>1)) {
    stop("using asinsqrt without abundance transformation is not.allowed")
  } else {
    data = data
  }
  return(t(data))
} 

###############--------------------- main ----------------######################
### 1. Getting the appropiate microbiome data

## microbiome data previously fitered
#setwd("/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/")
inDF <- read.table('DAG3_metaphlan_bacteria_archaea_filtered.txt',sep='\t',header=T)
colnames(inDF)
nrow(inDF)
## we read from this file;
## A total of 316 microbiome features from abundance information only

### 2. Selecting the population for GWAS analysis
## we take only the individuals from DAG3-UGLI population
#### samples with genetic information
dag3.pairing<-fread("pairing_DAG3_UGLI")
nrow(dag3.pairing)
### retireve QC information: samples that passed microbiome QC?
dagqc.samples<-fread("DAG3_QC_final.csv")
###get medication and BMI data (for now stored in gearshift later deleted)
phe.data<-fread("bmi.med.data",data.table = F )
head(phe.data)
### retrieve pathways form Metacyc (previously processed by ranko)
metacyc<-fread("DAG3_humann_Metacyc_filtered.txt",data.table = F )
nrow(metacyc)
head(metacyc[,c(1:5)])
ncol(metacyc)
head(metacyc[,c(335:339)])


### 3. merge information from all data

### add information about samples with genetic info
inDF$genetic_info<-ifelse(inDF$ID %in% dag3.pairing$sample,1,0)
nrow(inDF[inDF$genetic_info=="1",])
nrow(inDF)

### add information about Europeans
inDF$EUR<-"NA"
inDF$EUR<-ifelse((inDF$ID %in% dag3.pairing$sample & dag3.pairing$EUR[match(inDF$ID,dag3.pairing$sample)]=="1") ,"1","0")
nrow(inDF[inDF$EUR=="1",])
nrow(inDF[inDF$EUR=="0",])

### add information about microbiome QC
inDF$passed_mb_QC<-ifelse((inDF$ID %in% dagqc.samples$sample  &  dagqc.samples$QC.status[match(inDF$ID,dagqc.samples$sample)]=="OK") ,1,0)
nrow(inDF[inDF$passed_mb_QC=="1",])
nrow(inDF[inDF$passed_mb_QC=="0",])
inDF$QC.status<-dagqc.samples$QC.status[match(inDF$ID,dagqc.samples$sample)]
table(inDF$QC.status)
nrow(inDF[is.na(inDF$QC.status),])

### rename database
all.data<-inDF

### merge antropometric and medical data
colnames(phe.data)[1]<-"ID"
alldata2<-left_join(all.data,phe.data,by="ID")
nrow(alldata2)
summary(alldata2$ANTHRO.BMI)

### there are at least 21 samples with NA information in BMI... not sure why. List has been passed to Ranko
## look for the ones that are not using antibiotics
table(alldata2$MED.MEDS.Antibacterials_ATC_J01)

### add pathway information
alldata2<-left_join(alldata2,metacyc,by="ID")
nrow(alldata2)
ncol(alldata2)

###combine QC data
colnames(dagqc.samples)[1]<-"ID"
alldata2<-left_join(alldata2,dagqc.samples,by="ID")
nrow(alldata2)
ncol(alldata2)


########### select the population for analysis ############
## in gearshift
#setwd("/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/features_and_phenotypes")
#alldata2<-fread("DAG3_gen_BMI_merged_untransformed_v2.txt",header = T,data.table = F)
nrow(alldata2[which(alldata2$genetic_info=="1" &
              alldata2$passed_mb_QC=="1"&
                alldata2$EUR=="1" &
                is.na(alldata2$ANTHRO.BMI)==F),])

alldata3<-alldata2[which(alldata2$genetic_info=="1" &
                           alldata2$passed_mb_QC=="1"&
                           alldata2$EUR=="1" &
                           is.na(alldata2$ANTHRO.BMI)==F),]
nrow(alldata3)
ncol(alldata3)
## 7738 individuals will be used for gwas


### 4. 
################----------- normalize data ----------------#####################

#source(normalize.functions.R)
colsToNor <- c()
colsToNor <- c(colsToNor,grep('__',colnames(alldata3)))
colsToNor <- c(colsToNor,grep('PWY',colnames(alldata3)))
##################
#####using wersmaLabfunction to Log_normalize microbiome

merged_log = normalize(alldata3[,colsToNor],to.abundance = T,samples.in.rows = T)
merged.frame.log<-data.frame(t(merged_log))
merged.frame.log<-data.frame(merged.frame.log,alldata3[,-colsToNor])


## plot
log1<-ggplot(data=merged.frame.log,aes(merged.frame.log$k__Bacteria.p__Actinobacteria))+geom_density()+xlab("bac1")
log2<-ggplot(data=merged.frame.log,aes(merged.frame.log$k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Rikenellaceae.g__Alistipes.s__Alistipes_senegalensis))+geom_density()+xlab("bac2")
log3<-ggplot(data=merged.frame.log,aes(merged.frame.log$k__Bacteria.p__Firmicutes.c__Negativicutes.o__Selenomonadales.f__Veillonellaceae.g__Veillonella.s__Veillonella_parvula))+geom_density()+xlab("bac3")

pat1<-ggplot(data=merged.frame.log,aes(merged.frame.log$PWY.3841..folate.transformations.II))+geom_density()+xlab("pw1")
pat2<-ggplot(data=merged.frame.log,aes(merged.frame.log$PWY.5136..fatty.acid..beta..oxidation.II..peroxisome.))+geom_density()+xlab("pw2")
pat3<-ggplot(data=merged.frame.log,aes(merged.frame.log$PWY.5138..unsaturated..even.numbered.fatty.acid..beta..oxidation))+geom_density()+xlab("pw3")

colnames(merged.frame.log)[400:450]

tiff("normalized.pwys.1.tiff",  
     width = 2000, height = 3000, 
     units = "px", res = 300, compression = "lzw")
grid.arrange( ncol=2, 
              arrangeGrob(log1,log2,log3,top="log-transformed taxa"),
              arrangeGrob(pat1,pat2,pat3,top="log-trans. pwys"))

dev.off()

#######--------write the normalized data for further analysis----####
write.table(merged.frame.log,file = paste0('DAG3_gen_BMI_merged_log_v2.txt'),sep='\t',row.names = F,quote = F)

### 5. list features
####list mb features to analyze as outcome###
col_taxa<-grep('__',colnames(alldata3))
col_pwy<-grep('PWY',colnames(alldata3))

write.table(colnames(merged.frame.log[,col_taxa]),"all.taxa.v2.list",quote=F,row.names = F,col.names = F)
write.table(colnames(merged.frame.log[,col_pwy]),"all.pwy.v2.list",quote=F,row.names = F,col.names = F)

