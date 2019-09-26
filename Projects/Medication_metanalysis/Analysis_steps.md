Medication Project 
===================

Creator: Arnau Vich

Year: 2017

1.Raw pathaways, first, we filter the stratified pathways, keeping the information for the overall pathway. 
------------------------------------------------------------------------------------------------------------
```{bash}
less $input_humann2.tsv | head -1 >> $input_unstrat.tsv
less $input_humann2.tsv | grep -v "|" >> $input_unstrat.tsv
## Remove unmapped | unaligned pathways
```

2.Remove duplicate sample id's in the second row + remove path description, just keeping the metacyc id (split by ":")
-------------------------------------------------------------------------------------------------------------------------


3.Clean headers (terminal/bash)
--------------------------------

```{bash}

less IBD_humann2_unstrat.txt | sed -e s/_kneaddata_merged_Abundance//g > tmp1.txt 
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_humann2_unstrat_clean.txt
rm tmp1.txt


less LLD_humann2_unstrat.txt | sed -e s/_kneaddata_merged_Abundance//g > tmp1.txt 
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_LLD.txt tmp1.txt | tr " " "\t" > LLD_humann2_unstrat_clean.txt
rm tmp1.txt


less MIBS_humann2_unstrat.tsv | sed -e s/_kneaddata_merged_Abundance//g >  MIBS_humann2_unstrat_clean.txt

less MIBS_taxonomy_metaphlan2_092017.tsv | head -1 | sed -e s/_metaphlan//g | tr " " "\t" >  MIBS_taxonomy_clean.txt
less MIBS_taxonomy_metaphlan2_092017.tsv | awk -F "\t" '{if(NR>2){print $0}}' >>  MIBS_taxonomy_clean.txt

less IBD_taxonomy_metaphlan2_082017.txt | head -1 | sed -e s/_metaphlan//g  | tr " " "\t" >  tmp1.txt
less IBD_taxonomy_metaphlan2_082017.txt | awk -F "\t" '{if(NR>2){print $0}}' >>  tmp1.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_IBD.txt tmp1.txt | tr " " "\t" > IBD_taxonomy_unstrat_clean.txt
rm tmp1.txt

less LLD_taxonomy_metaphlan2_092017.txt | head -1 | sed -e s/_metaphlan//g  | tr " " "\t" >  tmp1.txt
less LLD_taxonomy_metaphlan2_092017.txt | awk -F "\t" '{if(NR>2){print $0}}' >>  tmp1.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ../rename_LLD.txt tmp1.txt | tr " " "\t" > LLD_taxonomy_unstrat_clean.txt
rm tmp1.txt

```


4.Metadata and summary statistics in R
----------------------------------------
```{R}
 #Filter metadata script: https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/Filter_metadata.R
 #Filter Min RD > 10 M reads (remove 30 samples)
 phenos=read.table("../phenotypes.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
 filter_metadata(phenos, 3, 10000000,100)

 ##Filter stomas / pouches
phenos=read.table("../filtered_metadata.txt", header = T, sep = "\t",check.names = F, row.names = 1)
stomas=read.table("../remove_stoma.txt")
remove=stomas$V1
final_phenos=phenos[!row.names(phenos)%in%remove,]

## Summary statistics
# Script: https://github.com/WeersmaLabIBD/General_Tools/blob/master/Metadata_summary.R

cat_table$name=row.names(final_phenos)
cat_table$cat=final_phenos$cohort
cat_table=unique(cat_table)
rownames(cat_table)=cat_table$name
cat_table$name=NULL
summary_statistics_metadata(final_phenos,cat_table)
summary_statistics_metadata(final_phenos)
samples_to_keep=row.names(final_phenos)

```



5.Filtering in R 
------------------------------------

```{r}
#Use script: https://github.com/WeersmaLabIBD/Microbiome/blob/master/Tools/normalize_data.R
#Import pathway data
path_IBD=read.table("../1.1.Cleaned_tables/IBD_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_IBD=path_IBD[,colnames(path_IBD)%in%samples_to_keep]
path_IBD2= t(path_IBD)
path_IBD2[is.na(path_IBD2)] = 0
path_IBD2 = sweep(path_IBD2,1,rowSums(path_IBD2),"/")
path_IBD2[is.nan(path_IBD2)] = 0

## Output is a filtered table, summary statistics, graph filtering statistics 

filtering_taxonomy(path_IBD2,0.00000001,10)
#colnames(path_IBD) = sub("_kneaddata_merged_Abundance","",colnames(path_IBD))

path_MIBS=read.table("../1.1.Cleaned_tables/MIBS_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_MIBS=path_MIBS[,colnames(path_MIBS)%in%samples_to_keep]
path_MIBS2= t(path_MIBS)
path_MIBS2[is.na(path_MIBS2)] = 0
path_MIBS2 = sweep(path_MIBS2,1,rowSums(path_MIBS2),"/")
path_MIBS2[is.nan(path_MIBS2)] = 0
filtering_taxonomy(path_MIBS2,0.00000001,10)

path_LLD=read.table("../1.1.Cleaned_tables/LLD_humann2_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
path_LLD=path_LLD[,colnames(path_LLD)%in%samples_to_keep]
path_LLD2= t(path_LLD)
path_LLD2[is.na(path_LLD2)] = 0
path_LLD2 = sweep(path_LLD2,1,rowSums(path_LLD2),"/")
path_LLD2[is.nan(path_LLD2)] = 0
filtering_taxonomy(path_LLD2,0.00000001,10)


## Create join table
path_tmp=merge(path_IBD, path_MIBS, by="row.names", all = T)
rownames(path_tmp)=path_tmp$Row.names
path_tmp$Row.names=NULL
path_all=merge(path_tmp, path_LLD, by="row.names", all = T)
rownames(path_all)=path_all$Row.names
path_all$Row.names=NULL
path_all=path_all[,colnames(path_all)%in%samples_to_keep]
path_all= t(path_LLD)
path_all[is.na(path_all)] = 0
path_all = sweep(path_all,1,rowSums(path_all),"/")
path_all[is.nan(path_all)] = 0
filtering_taxonomy(path_all,0.00000001,10)




## Now for taxonomies 

tax_IBD=read.table("../1.1.Cleaned_tables/IBD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_IBD=tax_IBD[,colnames(tax_IBD)%in%samples_to_keep]
tax_IBD=t(tax_IBD)
filtering_taxonomy(tax_IBD,0.01,10)

tax_MIBS=read.table("../1.1.Cleaned_tables/MIBS_taxonomy_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_MIBS=tax_MIBS[,colnames(tax_MIBS)%in%samples_to_keep]
tax_MIBS=t(tax_MIBS)
filtering_taxonomy(tax_MIBS,0.01,10)

tax_LLD=read.table("../1.1.Cleaned_tables/LLD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
tax_LLD=tax_LLD[,colnames(tax_LLD)%in%samples_to_keep]
tax_LLD=t(tax_LLD)
filtering_taxonomy(tax_LLD,0.01,10)


tax_tmp=merge(t(tax_IBD), t(tax_MIBS), by="row.names", all = T)
#View(tax_tmp)
rownames(tax_tmp)=tax_tmp$Row.names
tax_tmp$Row.names=NULL
tax_all=merge(tax_tmp, t(tax_LLD), by="row.names", all = T)
rownames(tax_all)=tax_all$Row.names
tax_all$Row.names=NULL
tax_all=t(tax_all)
filtering_taxonomy(tax_all,0.01,10)
```


6.Normalize and merge with phenos
-----------------------------

```{R}

path = "./"
out.file<-""
file.names <- dir(path, pattern ="path.txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=F)
  #Run log transformation keeping 0's on pathways: Log transformation keeping 0's. 
  file=normalize(file, samples.in.rows = T, to.abundance = F, transform = "log", move.log = "min")
  file=t(file)
  out.file <- merge(final_phenos, file, by="row.names")
  rownames(out.file)=out.file$Row.names
  out.file$Row.names=NULL
  write.table(out.file, paste(file.names[i], "_pheno.txt", sep = ""),  quote=F, row.names=T, col.names=T, sep="\t")
}


path = "./"
out.file<-""
file.names <- dir(path, pattern ="taxonomy.txt")
for(i in 1:length(file.names)){
  file <- read.table(file.names[i],header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names=F)
  #Normalization step!
  file=file/100
  file= asin(sqrt(file))
  out.file <- merge(final_phenos, file, by="row.names")
  rownames(out.file)=out.file$Row.names
  out.file$Row.names=NULL
  write.table(out.file, paste(file.names[i], "_pheno.txt", sep = ""),  quote=F, row.names=T, col.names=T, sep="\t")
}
```

7.Run linear models
----------------

**Loop testing 2 things: 
1) Associations per drug with Age, Sex, BMI, Read Depth as co-variates
2) Associations considering all drugs together Age, Sex, BMI, Read Depth as co-variates**

**Example for IBD associations**
```{R}
# Read files
IBD=read.table("MIBS_filtered_taxonomy.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/MIBS_filtered_taxonomy_pheno.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/LLD_filtered_taxonomy_pheno.txt", sep="\t", header = T, row.names = 1)
# Pathways
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/IBD_filtered_path_pheno.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/MIBS_filtered_path_pheno.txt", sep="\t", header = T, row.names = 1)
#IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/6.Input_files/LLD_filtered_path_pheno.txt", sep="\t", header = T, row.names = 1)

#Start flags and duplicate datasets with different names
myOutliers=TRUE
IBD2=IBD
#IBD2$cohort=NULL
IBD3=IBD2
#Create an empty matrix, number of rows = number of bacteria to test, NUMBER 48 indicates the first column containing taxa information. 
#Number of columns will depend on the number of test and the information that you want to keep in the output files
results_IBD=matrix(,ncol = 16, nrow = length(colnames(IBD2)[49:ncol(IBD2)]))
# x= will indicate the number of row to write the results in each iteration
x=0
#d=4 <- Not used anymore
#Change 48 for the first column containing taxa
#Loop detecting outliers
for (i in 49:ncol(IBD2)){
  x=x+1
  #Get initial statistics of number of participants and non-0 values per taxa
  results_IBD[x,1]=length(IBD2[,i])
  results_IBD[x,2]=sum(IBD2[,i]!=0)
  #Test outliers, but if all are marked as outliers keep the original values... If you want to perform strict analysis you may want to remove them or pre-filter the table
  while (any(myOutliers==T)){
    if (sum(is.na(IBD2[,i]))==length(IBD2[,i])){
        IBD2[,i]=IBD3[,i]
        break
    }
    outliers=grubbs.test(IBD2[,i])
    myOutliers = outlier( IBD2[,i], logical = TRUE )
    # Threshold for outlier identification
    if (outliers$p.value < 0.05){
      #Transform relative abundances to NA's in case of outliers
      IBD2[,i][myOutliers] <- NA
    } else {
      break
    }
  }
}
#Change 5:47 for the phenotypes to test (here 41 meds)
#Loop per medication
for (a in 5:48){
  #Report the number of users and non-users
  results_IBD[,3]=sum(IBD2[,a]=="User")
  results_IBD[,4]=sum(IBD2[,a]!="User")
  results_IBD[,5]=0
  results_IBD[,6]=0
  results_IBD[,7]=0
  results_IBD[,8]=0
  results_IBD[,9]=0
  results_IBD[,10]=0
  results_IBD[,11]=0
  results_IBD[,12]=0
  results_IBD[,13]=0
  results_IBD[,14]=0
  results_IBD[,15]=0
  # Keep info about of the drug name
  drug=colnames(IBD2[a])
  results_IBD[,16]=drug
  #results_IBD[,12]=drug
  z=0
  # Calculate frequency table of medication use per sex. 
    # If there's not medication users, skip analyses. 
  if (colnames(IBD2[a]) !="cohort"){
   if (sum(IBD2[,a]=="User")>=4 ){
  #Loop per taxa / feature
  #d=d+1 <- Not used anymore
    for (b in 49:ncol(IBD2)){
    z=z+1
    #Create a subset including Taxa, age, bmi, read-depth, sex, 1 medication
    temp=IBD2[,c(b,1:4,a,18)] 
    # Use complete cases -> remove participant containing NAs
    df2<-temp[complete.cases(temp),]
    # Create frequency table medication ~ sex
    #test=table(df2[,6],df2$Sex )
    # If there's no users, skip analysis, print NA
    if (sum(subset(df2,df2[,6]=="User")[,1]!=0)>2){
    # Perform an strict test, testing 1 taking into account age, bmi, read-depth, sex and all medication    
    ##CHECK HERE AS WELL: If variable less than 2 non-zero values will give an error      
    #LLD
        #v=lm(IBD2[,b]~Age+BMI+PFReads+Sex+ACE_inhibitor+alpha_blockers+angII_receptor_antagonist+anti_androgen_oral_contraceptive+anti_epileptics+anti_histamine+antibiotics_merged+benzodiazepine_derivatives_related+beta_blockers+beta_sympathomimetic_inhaler+bisphosphonates+ca_channel_blocker+calcium+cohort+ferrum+folic_acid+insulin+IUD_that_includes_hormones+K_saving_diuretic+laxatives+melatonine+mesalazines+metformin+methylphenidate+NSAID+opiat+oral_anti_diabetics+oral_contraceptive+oral_steroid+other_antidepressant+paracetamol+platelet_aggregation_inhibitor+PPI+SSRI_antidepressant+statin+steroid_inhaler+steroid_nose_spray+thiazide_diuretic+thyrax+tricyclic_antidepressant+triptans+vitamin_D+vitamin_K_antagonist, data = IBD2)
    #IBD
        #v=lm(IBD2[,b]~Age+BMI+PFReads+Sex+ACE_inhibitor+alpha_blockers+angII_receptor_antagonist+anti_androgen_oral_contraceptive+anti_epileptics+anti_histamine+antibiotics_merged+benzodiazepine_derivatives_related+beta_blockers+beta_sympathomimetic_inhaler+bisphosphonates+ca_channel_blocker+calcium+cohort+ferrum+folic_acid+insulin+IUD_that_includes_hormones+K_saving_diuretic+laxatives+melatonine+mesalazines+metformin+methylphenidate+NSAID+opiat+oral_anti_diabetics+oral_contraceptive+oral_steroid+other_antidepressant+paracetamol+platelet_aggregation_inhibitor+PPI+Ranitidine+SSRI_antidepressant+statin+steroid_inhaler+steroid_nose_spray+thiazide_diuretic+thyrax+tricyclic_antidepressant+triptans+vitamin_D+vitamin_K_antagonist, data = IBD2)
    #MIBS
        v=lm(IBD2[,b]~Age+BMI+PFReads+Sex+ACE_inhibitor+alpha_blockers+angII_receptor_antagonist+anti_androgen_oral_contraceptive+anti_epileptics+anti_histamine+antibiotics_merged+benzodiazepine_derivatives_related+beta_blockers+beta_sympathomimetic_inhaler+bisphosphonates+ca_channel_blocker+calcium+cohort+laxatives+mesalazines+metformin+NSAID+opiat+oral_anti_diabetics+oral_contraceptive+oral_steroid+other_antidepressant+paracetamol+platelet_aggregation_inhibitor+PPI+SSRI_antidepressant+statin+steroid_inhaler+steroid_nose_spray+thiazide_diuretic+thyrax+tricyclic_antidepressant+triptans+vitamin_D+vitamin_K_antagonist, data = IBD2)
        vv=summary(v)
        bla=as.data.frame(vv$coefficients)
        selection_drug=try(bla[grep(drug, rownames(bla)), ])
        #Check numnber of users >1
        results_IBD[z,13]=selection_drug[nrow(selection_drug),4]
        results_IBD[z,14]=selection_drug[nrow(selection_drug),1]
        results_IBD[z,15]=selection_drug[nrow(selection_drug),2]

        # Perform an less strict test, taxa ~ 1 Medication, Age, BMI, sequenced reads, sex
        y=lm(df2[,1]~df2[,6]+Age+BMI+PFReads+Sex+cohort, data = df2)
        yy=summary(y)
        # Extract pvalues, coefficients and SE
        results_IBD[z,5]=yy$coefficients[2,1]
        results_IBD[z,6]=yy$coefficients[2,2]
        results_IBD[z,7]=yy$coefficients[2,4]
      } 
      else {
        results_IBD[z,5]=NA
        results_IBD[z,6]=NA
        results_IBD[z,7]=NA
        results_IBD[z,8]=NA
        results_IBD[z,9]=NA
        results_IBD[z,10]=NA
        results_IBD[z,11]=NA
        results_IBD[z,12]=NA
        results_IBD[z,13]=NA
        results_IBD[z,14]=NA
        results_IBD[z,15]=NA
      }
    }
   }
 }
  # Adjust the less strict test using fdr method and p.adjust function
  results_IBD[,8]=p.adjust( results_IBD[,7], method = "fdr")
  #define column names
  colnames(results_IBD)=c("N","non-zeros", "Users", "Non-users", "Coef", "StdError","pvalue", "qvalue","pval_drug_interact","coeff_drug_interact","SE.drug_interact", "pval_drug_sex_test", "pval_correcting_all", "coeff_correcting_all", "SE.correcting_all", "drug")
  #define rownames 
  rownames(results_IBD)=colnames(IBD2)[49:ncol(IBD2)]
  ## CHANGE PREFIX OF THE RESULT FILE 
  #assign(paste('MIBS',drug , sep = '_'), results_IBD)
  write.table(results_IBD, file=paste('MIBS',drug , sep = '_'), sep="\t", quote=F)
}


```

8.Meta-analysis and Heterogeneity mesure
--------------------------------------------

```{R}
#library (meta)
#setwd("../paths/")
#drug_list=c("laxatives")
drug_list=c("ACE_inhibitor","alpha_blockers","angII_receptor_antagonist","anti_androgen_oral_contraceptive","anti_epileptics","anti_histamine","antibiotics_merged","benzodiazepine_derivatives_related","beta_blockers","beta_sympathomimetic_inhaler","bisphosphonates","ca_channel_blocker","calcium","laxatives","metformin","NSAID","opiat","oral_anti_diabetics","oral_contraceptiveZ","oral_steroid","other_antidepressant","paracetamol","platelet_aggregation_inhibitor","PPI","SSRI_antidepressant","statin","steroid_inhaler","steroid_nose_spray","thiazide_diuretic","thyrax","tricyclic_antidepressant","triptans","vitamin_D","vitamin_K_antagonist")
#drug_list=c("ACE_inhibitor","angII_receptor_antagonist","anti_epileptics","anti_histamine","antibiotics_merged","benzodiazepine_derivatives_related","beta_blockers","beta_sympathomimetic_inhaler","ca_channel_blocker","calcium","laxatives","metformin","NSAID","opiat","oral_contraceptiveZ","oral_steroid","other_antidepressant","paracetamol","platelet_aggregation_inhibitor","PPI","SSRI_antidepressant","statin","steroid_inhaler","steroid_nose_spray","thiazide_diuretic","thyrax","tricyclic_antidepressant","triptans","vitamin_D","vitamin_K_antagonist")
path="./"
flag=1
#Import and merge files
for (a in drug_list){
	file.names = dir(path,pattern=a)
	my_IBD=grep("IBD",file.names, value = T )
	my_MIBS=grep("MIBS",file.names, value = T )
	my_LLD=grep("LLD",file.names, value = T )
	IBD=read.table(my_IBD, sep="\t", header = T, row.names = 1)
	LLD=read.table(my_LLD, sep="\t", header = T, row.names = 1)
	MIBS=read.table(my_MIBS, sep="\t", header = T, row.names = 1)
	IBD_MIBS=merge(IBD,MIBS, by = "row.names")
	rownames(IBD_MIBS)=IBD_MIBS$Row.names
	IBD_MIBS$Row.names=NULL
	all=merge(IBD_MIBS,LLD, by = "row.names")
 
	selection=all
	#CHANGE column names to add cohort of origin. 
	#colnames(selection)=c("Taxa","Factor.IBD","Coef.IBD","N.IBD","N0.IBD","Pval.IBD","Qval.IBD","SE.IBD","NonUsers.IBD","Users.IBD","Neff.IBD","Ref.IBD","Alt.IBD","Factor.MIBS","Coef.MIBS","N.MIBS","N0.MIBS","Pval.MIBS","Qval.MIBS","SE.MIBS","NonUsers.MIBS","Users.MIBS","Neff.MIBS","Ref.MIBS","Alt.MIBS","Factor.LLD","Coef.LLD","N.LLD","N0.LLD","Pval.LLD","Qval.LLD","SE.LLD","NonUsers.LLD","Users.LLD","Neff.LLD","Ref.LLD","Alt.LLD")  
	colnames(selection)=c("Taxa", "N.IBD","non-zeros.IBD", "Users.IBD", "Non-users.IBD", "Coef.IBD", "SE.IBD","Pval.IBD", "Qval.IBD","pval.drug.interact.IBD", "coeff.drug.interact.IBD","SE.drug.interact.IBD" ,"pval.sex.drug.IBD", "pval.correcting.all.IBD", "coeff.correcting.all.IBD","SE.correcting_all.IBD","drug.IBD","N.MIBS","non-zeros.MIBS", "Users.MIBS", "Non-users.MIBS", "Coef.MIBS", "SE.MIBS","Pval.MIBS", "Qval.MIBS","pval.drug.interact.MIBS", "coeff.drug.interact.MIBS","SE.drug.interact.MIBS" ,"pval.sex.drug.MIBS", "pval.correcting.all.MIBS", "coeff.correcting.all.MIBS","SE.correcting_all.MIBS","drug.MIBS","N.LLD","non-zeros.LLD", "Users.LLD", "Non-users.LLD", "Coef.LLD", "SE.LLD","Pval.LLD", "Qval.LLD","pval.drug.interact.LLD", "coeff.drug.interact.LLD","SE.drug.interact.LLD" ,"pval.sex.drug.LLD", "pval.correcting.all.LLD", "coeff.correcting.all.LLD","SE.correcting_all.LLD","drug.LLD")
	if (sum(selection$SE.correcting_all.IBD,na.rm = T)==0){
	#Calculate Inverse variance
	#51
	selection$inverse_var.mibs=1/selection$SE.MIBS^2
	#52
	selection$inverse_var.lld=1/selection$SE.LLD^2
#	selection$inverse_var.mibs=1/selection$SE.correcting_all.MIBS^2
#	selection$inverse_var.lld=1/selection$SE.correcting_all.LLD^2

	#Calculate SE  #53
	selection$se=sqrt(1/(selection$inverse_var.mibs+selection$inverse_var.lld))

	#Calculate Beta #54
	selection$beta=(selection$inverse_var.mibs*selection$Coef.MIBS+selection$inverse_var.lld*selection$Coef.LLD)/(selection$inverse_var.mibs+selection$inverse_var.lld)
#	selection$beta=(selection$inverse_var.mibs*selection$coeff.correcting.all.MIBS+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.mibs+selection$inverse_var.lld)
	
	#Calculate Z-score #55
	selection$Z=selection$beta/selection$se
	
	#Calculate meta p-value #56
	selection$P=2*pnorm(-abs(selection$Z))
	#Remove data only avalable in one cohort
	selection=selection[complete.cases(selection$P), ]

	#Adjust pvalue with FDR #57
	selection$FDR=p.adjust(selection$P,method = "fdr")
	
	#Create empty columns
	selection$Het.Q="NS" #58
	selection$Het.I2="NS" #59
	selection$Het.Pval="NS" #60
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			#Select coefficients
			TE=c( selection[i,22], selection[i,38])
#			TE=c( selection[i,31], selection[i,47])
			#Select Standart error
			SE=c( selection[i,23], selection[i,39])
#			SE=c( selection[i,32], selection[i,48])
			het=metagen(TE,SE)
			#Change the number here depending of the number of columns in your data, should match with column Het.I2
			selection[i,59]=het$I2
			#Match Het.Q column number
			selection[i,58]=het$Q 
			#Calculate p-value from Q calculation // Het.Pval
			selection[i,60]=pchisq(het$Q,df=2,lower.tail=F)
		}
	}
	} else if (sum(selection$SE.correcting_all.MIBS,na.rm = T)==0) {
	#Calculate Inverse variance
	#50
	selection$inverse_var.ibd=1/selection$SE.IBD^2
	#52
	selection$inverse_var.lld=1/selection$SE.LLD^2

#	selection$inverse_var.ibd=1/selection$SE.correcting_all.IBD^2
#	selection$inverse_var.lld=1/selection$SE.correcting_all.LLD^2

	#Calculate SE  #53
	selection$se=sqrt(1/(selection$inverse_var.ibd+selection$inverse_var.lld))

	#Calculate Beta #54
	selection$beta=(selection$inverse_var.ibd*selection$Coef.IBD+selection$inverse_var.lld*selection$Coef.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.lld)
	
	#Calculate Z-score #55
	selection$Z=selection$beta/selection$se
	
	#Calculate meta p-value #56
	selection$P=2*pnorm(-abs(selection$Z))
	#Remove data only avalable in one cohort
	selection=selection[complete.cases(selection$P), ]

	#Adjust pvalue with FDR #57
	selection$FDR=p.adjust(selection$P,method = "fdr")
	
	#Create empty columns
	selection$Het.Q="NS" #58
	selection$Het.I2="NS" #59
	selection$Het.Pval="NS" #60
	
	
	#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			#Select coefficients
			TE=c( selection[i,6], selection[i,38])
#			TE=c( selection[i,15], selection[i,47])
			#Select Standart error
			SE=c( selection[i,7], selection[i,39])
#			SE=c( selection[i,16], selection[i,48])
			het=metagen(TE,SE)
			#Change the number here depending of the number of columns in your data, should match with column Het.I2
			selection[i,59]=het$I2
			#Match Het.Q column number
			selection[i,58]=het$Q 
			#Calculate p-value from Q calculation // Het.Pval
			selection[i,60]=pchisq(het$Q,df=2,lower.tail=F)
		}
	}
	} else if (sum(selection$SE.correcting_all.LLD,na.rm = T)==0) {

	#Calculate Inverse variance
	#50
	selection$inverse_var.ibd=1/selection$SE.IBD^2
	#51
	selection$inverse_var.mibs=1/selection$SE.MIBS^2
#	selection$inverse_var.ibd=1/selection$SE.correcting_all.IBD^2
#	selection$inverse_var.mibs=1/selection$SE.correcting_all.MIBS^2
	#Calculate SE  #53
	selection$se=sqrt(1/(selection$inverse_var.ibd+selection$inverse_var.mibs))

	#Calculate Beta #54
	selection$beta=(selection$inverse_var.ibd*selection$Coef.IBD+selection$inverse_var.mibs*selection$Coef.MIBS)/(selection$inverse_var.ibd+selection$inverse_var.mibs)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.mibs*selection$coeff.correcting.all.MIBS)/(selection$inverse_var.ibd+selection$inverse_var.mibs)
	#Calculate Z-score #55
	selection$Z=selection$beta/selection$se
	#Calculate meta p-value #56
	selection$P=2*pnorm(-abs(selection$Z))
	#Remove data only avalable in one cohort
	selection=selection[complete.cases(selection$P), ]
	#Adjust pvalue with FDR #57
	selection$FDR=p.adjust(selection$P,method = "fdr")
	#Create empty columns
	selection$Het.Q="NS" #58
	selection$Het.I2="NS" #59
	selection$Het.Pval="NS" #60
	#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			#Select coefficients
			TE=c( selection[i,6], selection[i,22])
#			TE=c( selection[i,15], selection[i,31])
			#Select Standart error
			SE=c( selection[i,7], selection[i,23])
#			SE=c( selection[i,16], selection[i,32])
			het=metagen(TE,SE)
			#Change the number here depending of the number of columns in your data, should match with column Het.I2
			selection[i,59]=het$I2
			#Match Het.Q column number
			selection[i,58]=het$Q 
			#Calculate p-value from Q calculation // Het.Pval
			selection[i,60]=pchisq(het$Q,df=2,lower.tail=F)
		}
	} 
	}else {
	#Calculate Inverse variance
	#50
	selection$inverse_var.ibd=1/selection$SE.IBD^2
	#51
	selection$inverse_var.mibs=1/selection$SE.MIBS^2
	#52
	selection$inverse_var.lld=1/selection$SE.LLD^2

#	selection$inverse_var.ibd=1/selection$SE.correcting_all.IBD^2
#	selection$inverse_var.mibs=1/selection$SE.correcting_all.MIBS^2
#	selection$inverse_var.lld=1/selection$SE.correcting_all.LLD^2

	#Calculate SE  #53
	selection$se=sqrt(1/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld))

	#Calculate Beta #54
	selection$beta=(selection$inverse_var.ibd*selection$Coef.IBD+selection$inverse_var.mibs*selection$Coef.MIBS+selection$inverse_var.lld*selection$Coef.LLD)/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld)
#	selection$beta=(selection$inverse_var.ibd*selection$coeff.correcting.all.IBD+selection$inverse_var.mibs*selection$coeff.correcting.all.MIBS+selection$inverse_var.lld*selection$coeff.correcting.all.LLD)/(selection$inverse_var.ibd+selection$inverse_var.mibs+selection$inverse_var.lld)
	
	#Calculate Z-score #55
	selection$Z=selection$beta/selection$se
	
	#Calculate meta p-value #56
	selection$P=2*pnorm(-abs(selection$Z))
	#Remove data only avalable in one cohort
	selection=selection[complete.cases(selection$P), ]

	#Adjust pvalue with FDR #57
	selection$FDR=p.adjust(selection$P,method = "fdr")
	
	#Create empty columns
	selection$Het.Q="NS" #58
	selection$Het.I2="NS" #59
	selection$Het.Pval="NS" #60
	
	
	#Heterogeneity using Cochran's Q-test for meta-FDR < 0.1
	for (i in 1:length(rownames(selection))){
		if (selection$FDR[i]<0.1){
			#Select coefficients
			TE=c( selection[i,6], selection[i,22], selection[i,38])
#			TE=c( selection[i,15], selection[i,31], selection[i,47])
			#Select Standart error
			SE=c( selection[i,7], selection[i,23], selection[i,39])
#			SE=c( selection[i,16], selection[i,32], selection[i,48])
			het=metagen(TE,SE)
			#Change the number here depending of the number of columns in your data, should match with column Het.I2
			selection[i,59]=het$I2
			#Match Het.Q column number
			selection[i,58]=het$Q 
			#Calculate p-value from Q calculation // Het.Pval
			selection[i,60]=pchisq(het$Q,df=2,lower.tail=F)
		}
	}
	}
	selection$qval.correcting.all.IBD=p.adjust( selection$pval.correcting.all.IBD, method = "fdr")
	selection$qval.correcting.all.MIBS=p.adjust( selection$pval.correcting.all.MIBS, method = "fdr")
	selection$qval.correcting.all.LLD=p.adjust( selection$pval.correcting.all.LLD, method = "fdr")
	selection$pval.drug.interact.IBD=NULL
	selection$coeff.drug.interact.IBD=NULL
	selection$SE.drug.interact.IBD=NULL
	selection$pval.sex.drug.IBD=NULL
	selection$pval.drug.interact.MIBS=NULL
	selection$coeff.drug.interact.MIBS=NULL
	selection$SE.drug.interact.MIBS=NULL
	selection$pval.sex.drug.MIBS=NULL
	selection$pval.drug.interact.LLD=NULL
	selection$coeff.drug.interact.LLD=NULL
	selection$SE.drug.interact.LLD=NULL
	selection$pval.sex.drug.LLD=NULL
write.table(selection, file=paste(a, "_meta.txt", sep=""), sep="\t", quote=F)
}
```


9.Diversity measuraments
--------------------------

Calculate Shannon Index and create boxplots per cohorts

```{r}
library(vegan)
library(ggplot2)
library(RColorBrewer)


IBD=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/IBD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
LLD=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/LLD_taxonomy_unstrat_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
MIBS=read.table("~/Desktop/PPI_v2/00.With_new_data/1.1.Cleaned_tables/MIBS_taxonomy_clean.txt", header = T, sep = "\t",check.names = F, as.is = T, row.names = 1)
     
phenos=read.table("~/Desktop/PPI_v2/00.With_new_data/filtered_metadata.txt", header = T, sep = "\t",check.names = F, row.names = 1, stringsAsFactors = T)
stomas=read.table("~/Desktop/PPI_v2/00.With_new_data/remove_stoma.txt")
remove=stomas$V1
final_phenos=phenos[!row.names(phenos)%in%remove,]
samples_to_keep=row.names(final_phenos)

tax_IBD=IBD[,colnames(IBD)%in%samples_to_keep]
tax_MIBS=MIBS[,colnames(MIBS)%in%samples_to_keep]
tax_LLD=LLD[,colnames(LLD)%in%samples_to_keep]
tax_IBD2 = as.data.frame(t(tax_IBD[!duplicated(tax_IBD, fromLast = T), ]))
colnames(tax_IBD2) = gsub("[|]",".",colnames(tax_IBD2))
tax_MIBS2 = as.data.frame(t(tax_MIBS[!duplicated(tax_MIBS, fromLast = T), ]))
colnames(tax_MIBS2) = gsub("[|]",".",colnames(tax_MIBS2))
tax_LLD2 = as.data.frame(t(tax_LLD[!duplicated(tax_LLD, fromLast = T), ]))
colnames(tax_LLD2) = gsub("[|]",".",colnames(tax_LLD2))

filtering_taxonomy = function(x){
  x = rev(x)
  result = c()
  for(i in 1:length(x)){
    rmstr = grep("t__",x[i])
    if (length(rmstr) >0) next
    check_prev = grep(x[i], result[length(result)])
    if (length(check_prev) == 0) result = append(result,x[i])
  }
  
  rev(result)
}

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_MIBS2))
tax_MIBS3= tax_MIBS2[,colnames_taxa2_filtered]

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_IBD2))
tax_IBD3= tax_IBD2[,colnames_taxa2_filtered]

colnames_taxa2_filtered = filtering_taxonomy(colnames(tax_LLD2))
tax_LLD3= tax_LLD2[,colnames_taxa2_filtered]

coded_phenos=final_phenos
coded_phenos[coded_phenos=="User"]=1
coded_phenos[coded_phenos=="Non_user"]=0
coded_phenos2=coded_phenos
coded_phenos2$cohort=NULL
coded_phenos3=coded_phenos2[5:46]
coded_phenos3$total_meds=0
coded_phenos3=data.matrix(coded_phenos3)
coded_phenos3=as.data.frame(coded_phenos3)
coded_phenos3$total_meds <- rowSums(coded_phenos3)

MIBS_SI <- as.data.frame(diversity(tax_MIBS3, index="shannon"))
LLD_SI <- as.data.frame(diversity(tax_LLD3, index="shannon"))
IBD_SI <- as.data.frame(diversity(tax_IBD3, index="shannon"))

coded_phenos4=coded_phenos3
coded_phenos4[1:42]=NULL



IBD_SI2=merge(coded_phenos4,IBD_SI, by="row.names")
MIBS_SI2=merge(coded_phenos4,MIBS_SI, by="row.names")
LLD_SI2=merge(coded_phenos4,LLD_SI, by="row.names")

colnames(MIBS_SI2)=c("SID","Meds","Div.")
colnames(IBD_SI2)=c("SID","Meds","Div.")
colnames(LLD_SI2)=c("SID","Meds","Div.")
ggplot( MIBS_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")
ggplot( LLD_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")
ggplot( IBD_SI2, aes(x=as.factor(Meds),y=Div.)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index")

library(ggplot2)
library(reshape2)

#Create figure shannon number medication x cohort
IBD_SI2$cohort="IBD"
MIBS_SI2$cohort="MIBS"
LLD_SI2$cohort="LLD"
tmpSI=rbind(IBD_SI2,MIBS_SI2)
SI2=rbind(tmpSI,LLD_SI2)
SI2$Meds2=SI2$Meds
SI2[SI2$Meds>5,]$Meds2="+5"
ggplot (SI2, aes (Meds2, Div., fill=cohort)) + geom_boxplot(outlier.size = NA, alpha=0.7, colour="grey36") + geom_jitter(shape=16, alpha=0.2, position=position_jitter(0.2)) + theme_bw() + geom_smooth(method = "lm", color="black", aes(group=1), size=0.8, fill="purple") + facet_grid(~cohort)

fit = lm (Div. ~ Meds, data= IBD_SI2)
summary(fit)
p=0.88
estimated=-0.001029 
r2=-0.002168


fit = lm (Div. ~ Meds, data= LLD_SI2)
summary(fit)
p=0.77
estimated=-0.0016 
r2=-0.00080

fit = lm (Div. ~ Meds, data= MIBS_SI2)
summary(fit)
p=0.059
estimated=-0.016 
r2=-0.0089


# Test shannon x cohort x medication
row.names(SI2)=SI2$SID
SI2$SID=NULL
tmpSI=merge(SI2, coded_phenos3, by="row.names")
tmpSI3=melt(tmpSI, id.vars = c("Row.names", "Meds", "Div.", "Meds2", "cohort"))
ggplot(tmpSI3, aes(as.factor(value), Div., fill=as.factor(value))) + geom_boxplot(outlier.size = 0.2) + facet_grid(cohort~variable) +  scale_fill_manual(values=c ("grey36", "grey90")) + xlab("Users yes/no") + ylab("Shannon Index") + theme_bw()


tmpLLD=tmpSI[tmpSI$cohort=="LLD", ]
tmpMIBS=tmpSI[tmpSI$cohort=="MIBS", ]
tmpIBD=tmpSI[tmpSI$cohort=="IBD", ]

shannon_medi_lld=matrix(nrow = ncol(tmpLLD)-5, ncol = 6)
cnt=1
for (i in 6:ncol(tmpLLD)){
  #mean_users
  shannon_medi_lld[cnt,2]=mean(tmpLLD$Div.[tmpLLD[,i]==1])
  #sd_users
  shannon_medi_lld[cnt,3]=sd(tmpLLD$Div.[tmpLLD[,i]==1])
  #mean_n_users
  shannon_medi_lld[cnt,4]=mean(tmpLLD$Div.[tmpLLD[,i]==0])
  #sd_n_users
  shannon_medi_lld[cnt,5]=sd(tmpLLD$Div.[tmpLLD[,i]==0])
  test=wilcox.test(tmpLLD$Div.[tmpLLD[,i]==0],tmpLLD$Div.[tmpLLD[,i]==1] )
  shannon_medi_lld[cnt,6]=test$p.value
  cnt=cnt+1
}
shannon_medi_lld=as.data.frame(shannon_medi_lld)
shannon_medi_lld[,1]=colnames(tmpLLD)[6:ncol(tmpLLD)]
colnames(shannon_medi_lld)=c("Meds", "MeanUsers","SDUsers", "MeanNonUsers", "SDNonUsers", "pvalue")
shannon_medi_lld$qval=0
shannon_medi_lld$qval=p.adjust(shannon_medi_lld$pvalue, method = "fdr")

shannon_medi_IBD=matrix(nrow = ncol(tmpIBD)-5, ncol = 6)
cnt=1
for (i in 6:ncol(tmpIBD)){
  #mean_users
  shannon_medi_IBD[cnt,2]=mean(tmpIBD$Div.[tmpIBD[,i]==1])
  #sd_users
  shannon_medi_IBD[cnt,3]=sd(tmpIBD$Div.[tmpIBD[,i]==1])
  #mean_n_users
  shannon_medi_IBD[cnt,4]=mean(tmpIBD$Div.[tmpIBD[,i]==0])
  #sd_n_users
  shannon_medi_IBD[cnt,5]=sd(tmpIBD$Div.[tmpIBD[,i]==0])
  test=wilcox.test(tmpIBD$Div.[tmpIBD[,i]==0],tmpIBD$Div.[tmpIBD[,i]==1] )
  shannon_medi_IBD[cnt,6]=test$p.value
  cnt=cnt+1
}
shannon_medi_IBD=as.data.frame(shannon_medi_IBD)
shannon_medi_IBD[,1]=colnames(tmpIBD)[6:ncol(tmpIBD)]
colnames(shannon_medi_IBD)=c("Meds", "MeanUsers","SDUsers", "MeanNonUsers", "SDNonUsers", "pvalue")
shannon_medi_IBD$qval=0
shannon_medi_IBD$qval=p.adjust(shannon_medi_IBD$pvalue, method = "fdr")


shannon_medi_MIBS=matrix(nrow = ncol(tmpMIBS)-5, ncol = 6)
cnt=1
for (i in 6:ncol(tmpMIBS)){
  #mean_users
  if (length(tmpMIBS$Div.[tmpMIBS[,i]==1])==0){
    cnt=cnt+1
    next
    
  }
  shannon_medi_MIBS[cnt,2]=mean(tmpMIBS$Div.[tmpMIBS[,i]==1])
  #sd_users
  shannon_medi_MIBS[cnt,3]=sd(tmpMIBS$Div.[tmpMIBS[,i]==1])
  #mean_n_users
  shannon_medi_MIBS[cnt,4]=mean(tmpMIBS$Div.[tmpMIBS[,i]==0])
  #sd_n_users
  shannon_medi_MIBS[cnt,5]=sd(tmpMIBS$Div.[tmpMIBS[,i]==0])
  test=wilcox.test(tmpMIBS$Div.[tmpMIBS[,i]==0],tmpMIBS$Div.[tmpMIBS[,i]==1] )
  shannon_medi_MIBS[cnt,6]=test$p.value
  cnt=cnt+1
}
shannon_medi_MIBS=as.data.frame(shannon_medi_MIBS)
shannon_medi_MIBS[,1]=colnames(tmpMIBS)[6:ncol(tmpMIBS)]
colnames(shannon_medi_MIBS)=c("Meds", "MeanUsers","SDUsers", "MeanNonUsers", "SDNonUsers", "pvalue")
shannon_medi_MIBS$qval=0
shannon_medi_MIBS$qval=p.adjust(shannon_medi_MIBS$pvalue, method = "fdr")


```


10.Calculate beta-diversity, create pcoa plots and MANOVA (adonis) analyses
---


```{R}
LLD_SI2$cohort="1_LLD"
MIBS_SI2$cohort="2_MIBS"
IBD_SI2$cohort="3_IBD"
b1=rbind(LLD_SI2,MIBS_SI2)
b2=rbind(b1,IBD_SI2)

ggplot( b2, aes(x=as.factor(Meds),y=Div., fill=cohort)) + geom_boxplot() + theme_classic() + xlab("Number of different medication") + ylab("Shannon Index") + ylim(0,3) + facet_grid(. ~ cohort) + scale_fill_manual(breaks=b2$cohort, values=my_col)


beta<- vegdist(tax_IBD3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_IBD=as.data.frame(mypcoa)

beta<- vegdist(tax_MIBS3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_MIBS=as.data.frame(mypcoa)

beta<- vegdist(tax_LLD3, method="bray")
mypcoa=cmdscale(beta, k = 5)
beta_LLD=as.data.frame(mypcoa)

b1=rbind(beta_LLD,beta_MIBS)
b3=rbind(b1,beta_IBD)

rownames(b2)=b2$SID
b2$SID=NULL
b4=merge(b2,b3,by="row.names")

ggplot(b4,aes(V3,V4))  +  geom_point (size=2, aes(colour=as.factor(Meds))) + theme_classic() + labs(x="PCoA1", y="PCoA3") + facet_grid(. ~ cohort) + scale_color_manual("Total Medication",values = blues_fun(14))

adonis_IBD <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_IBD3),]
ad <- adonis(formula = tax_IBD3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_IBD[1,1]=aov_table[1,1]
adonis_IBD[1,2]=aov_table[1,5]
adonis_IBD[1,3]=aov_table[1,6]
colnames(adonis_IBD) = c("Df", "R2", "Pr(>F)")


adonis_MIBS <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_MIBS3),]
ad <- adonis(formula = tax_MIBS3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_MIBS[1,1]=aov_table[1,1]
adonis_MIBS[1,2]=aov_table[1,5]
adonis_MIBS[1,3]=aov_table[1,6]
colnames(adonis_MIBS) = c("Df", "R2", "Pr(>F)")


adonis_LLD <- matrix(ncol = 3, nrow=1) 
b5_pheno <- b2[rownames(tax_LLD3),]
ad <- adonis(formula = tax_LLD3 ~ b5_pheno[,1] , data = b5_pheno, permutations = 10000, method = "bray")
aov_table <- ad$aov.tab
adonis_LLD[1,1]=aov_table[1,1]
adonis_LLD[1,2]=aov_table[1,5]
adonis_LLD[1,3]=aov_table[1,6]
colnames(adonis_LLD) = c("Df", "R2", "Pr(>F)")

```


11.Figures
--------------------------

```{R}
#Figures 
#Figure1
xx <- read.delim(".counts_results.txt")
ggplot(xx,aes(Meds,value,fill=Cohort)) + geom_bar(stat="identity", position=position_dodge()) + theme_minimal() + scale_fill_manual(values=c('red','black','gray46','gray72' )) + facet_grid(Test~Feature, drop=T, scales = "free") + coord_flip() + scale_x_discrete(limits = rev(levels(xx$Meds)))



#Figure2
lld <- read.delim("./fig2/lld.txt", row.names=1, check.names = F)
mibs<- read.delim("./fig2/mibs.txt", row.names=1, check.names = F)
ibd<- read.delim("./fig2/ibd.txt", row.names=1, check.names = F)

rel_lld=sweep(lld, 2, colSums(lld), '/')
rel_ibd=sweep(ibd, 2, colSums(ibd), '/')
rel_mibs=sweep(mibs, 2, colSums(mibs), '/')

pheno_ibd=read.delim("./IBD_phenos.txt", row.names=1, check.names = F)
pheno_lld=read.delim("./LLD_phenos.txt", row.names=1, check.names = F)
pheno_mibs=read.delim("./MIBS_phenos.txt", row.names=1, check.names = F)

lld2=merge(pheno_lld, as.data.frame(t(rel_lld)), by = "row.names")
ibd2=merge(pheno_ibd, as.data.frame(t(rel_ibd)), by = "row.names")
mibs2=merge(pheno_mibs, as.data.frame(t(rel_mibs)), by = "row.names")

ibd3=as.data.frame(t(ibd2[,c(1,19,38,49:ncol(ibd2))]))
mibs3=as.data.frame(t(mibs2[,c(1,19,38,49:ncol(mibs2))]))
lld3=as.data.frame(t(lld2[,c(1,19,38,49:ncol(lld2))]))


x1=merge(ibd3, mibs3, by = "row.names", all.x = T)
row.names(x1)=x1$Row.names
x1$Row.names=NULL
x2=merge(x1, lld3, by = "row.names", all.x = T)

x3=as.data.frame(t(x2))

colnames(x3)=unlist(x3[1,])
x3=x3[-1,]

x4=melt(x3,id.vars = c("Row.names","cohort","PPI"))

x4$variable2=gsub(".*s__","",x4$variable)
x4[x4$variable=="PWY0-1297:_superpathway_of_purine_deoxyribonucleosides_degradation|unclassified",]$variable2="unclassified"
ggplot(x4,aes(variable2,as.numeric(as.character(value2)),fill=PPI)) + geom_boxplot(outlier.size = 0.2)+facet_grid(cohort2~.) + theme_minimal()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Relative contribution to Purine deoxyribonucleosides degradation (PWY0-1297)") + xlab("")

```



12.Test stratified pathways (bacterial contribution to each pathways) and gene familes (pathways associated to drug categories in the previous steps)
----------------

```{bash}
#Obtain stratified pathways (pathways + contributing species)
#rename sample id's
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ./rename_IBD.txt ./IBD_humann2_pathways_uniref90_082017.txt | tr " " "\t" > IBD_renamed.txt
awk 'NR==FNR{d[$1]=$2;next}FNR==1{for(i=1;i<=NF;i++)$i=d[$i]?d[$i]:$i}7' ./rename_LLD.txt ./LLD_humann2_Uniref90_092017.txt | tr " " "\t" > LLD_renamed.txt

#subset samples used in the analysis
awk -f extract.awk selected.txt LLD_renamed.txt | tr " " "\t" > LLD_strat_samples.txt
awk -f extract.awk selected.txt MIBS_humann2_uniref90_092017.tsv | tr " " "\t" > MIBS_strat_samples.txt 
awk -f extract.awk selected.txt IBD_renamed.txt | tr " " "\t" > IBD_strat_samples.txt

#select pathways associated to each drug (e.g PPI)
less MIBS_strat_samples.txt | head -1 > PPI_MIBS2.txt
less IBD_strat_samples.txt | head -1 > PPI_IBD2.txt
less LLD_strat_samples.txt | head -1 > PPI_LLD2.txt
for i in `cat PPI.txt` ; do less MIBS_strat_samples.txt | grep ^$i >> PPI_MIBS2.txt ; done
for i in `cat PPI.txt` ; do less IBD_strat_samples.txt | grep ^$i >> PPI_IBD2.txt ; done
for i in `cat PPI.txt` ; do less LLD_strat_samples.txt | grep ^$i >> PPI_LLD2.txt ; done

#Re-normalize pathways, calculating the relative contribution of each species to the pathways 
humann2_renorm_table -i PPI_IBD2.txt -o PPI_IBD_norm.txt -u relab -m levelwise
humann2_renorm_table -i PPI_MIBS2.txt -o PPI_MIBS_norm.txt -u relab -m levelwise
humann2_renorm_table -i PPI_LLD2.txt -o PPI_LLD_norm.txt -u relab -m levelwise
```

```{bash}
#Obtain gene families
humann2_unpack_pathways --input-pathways all_pahtways.tsv --output all_gene_families.tsv

#Extract gene families related to the drug-associated pathways
for i in `cat oral_steroids.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_oral_steroids.txt ; done
for i in `cat SSRI.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_SSRI.txt ; done
for i in `cat laxatives.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_laxatives.txt ; done
for i in `cat metformin.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_metformin.txt ; done
for i in `cat PPI.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_PPI.txt ; done
for i in `cat AB.txt` ; do less all_gene_families.tsv | grep ^$i | awk -F "|" '{print $3}' | awk '{print $1}'>> gf_AB.txt ; done

#Remove non alphanumerical characters from the gene families names (optional)
less MIBSx.txt | tr -d '(),:#' | tr -d "/" | tr -d \' | awk -F "\t" '{if(NR!=2){print $0}}'> MIBS.txt

#Filter gene families absent in more than 90% of the samples (repeat per cohort)
less MIBS.txt | awk -F "\t" '{ zeroes=0; for( i=1; i<NF; i++ ) { if( $i == 0 ) {zeroes++} } if( zeroes / (NF-1) < 0.9 ) { print $0 } }' >> MIBS_filt.txt

#Subset per cohort and medication
head -1 MIBS_filt.txt > ./new2/MIBS_filt_metformin.txt
awk  'FNR==NR {arr[$0];next} $1 in arr' met_mibs.txt MIBS_filt.txt > ./new2/MIBS_filt_metformin.txt
head -1 MIBS_filt.txt > ./new2/MIBS_filt_laxatives.txt
awk 'FNR==NR {arr[$0];next} $1 in arr' lax_mibs.txt MIBS_filt.txt > ./new2/MIBS_filt_laxatives.txt
head -1 MIBS_filt.txt > ./new2/MIBS_filt_ppi.txt
awk  'FNR==NR {arr[$0];next} $1 in arr' ppi_mibs.txt MIBS_filt.txt > ./new2/MIBS_filt_ppi.txt
head -1 MIBS_filt.txt > ./new2/MIBS_filt_ab.txt
awk  'FNR==NR {arr[$0];next} $1 in arr' ab_mibs.txt MIBS_filt.txt > ./new2/MIBS_filt_ab.txt

#Normalize data

for i in *txt; do humann2_renorm_table -i $i -o $i_norm.txt -u relab -m levelwise; done

```

```{r}
#Repeat per drug and gene familes / stratified pathways
library (coin)
#Load phenotypes
IBD_pheno=read.table("./IBD_phenos.txt",header = T,row.names = 1)
#Load normalized data
IBD=read.table("./ibd_filt_ab.txt", header = T, check.names = F)
#Load raw data (for statistics, e.g. median, non-zeros, etc.)
IBDn=read.table("./IBD_filt_ab.txt", header = T, check.names = F)
IBD=unique(IBD)
IBDn=unique(IBDn)
rownames(IBD)=IBD$Gene_Family
rownames(IBDn)=IBDn$Gene_Family
IBD$Gene_Family=NULL
IBDn$Gene_Family=NULL
IBD_n=as.data.frame(t(IBDn))
test_n=merge(IBD_pheno, IBD_n, by="row.names")
#Log transform
IBD=as.data.frame(t(normalize(IBD,samples.in.rows = F,to.abundance = F, transform = "log", move.log = "min")))
test=merge(IBD_pheno, IBD, by="row.names")

x=0
my_res=matrix(,ncol = 10, nrow = length(colnames(test)[49:ncol(test)]))

for (i in 49:ncol(test)){
	x=x+1
	my_res[x,2]=median(test_n[,i][test_n$antibiotics_merged=="User"])
  	my_res[x,3]=median(test_n[,i][test_n$antibiotics_merged!="User"])
  	my_res[,4]=sum(test_n$antibiotics_merged=="User")
  	my_res[,5]=sum(test_n$antibiotics_merged!="User")
  	my_res[x,6]=sum(test_n[,i][test_n$antibiotics_merged=="User"]!=0)
  	my_res[x,7]=sum(test_n[,i][test_n$antibiotics_merged!="User"]!=0)
  	#Wilcoxon effect size
	my_res[x,8]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE))
	my_res[x,9]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE)) / sqrt (nrow(test))
	asoc=wilcox.test(test[,i]~test$antibiotics_merged)
	my_res[x,1]=colnames(test)[i]
	my_res[x,10]=asoc$p.value
}
my_res2=as.data.frame(my_res[my_res[,4]!="NaN", ])
colnames(my_res2)=c("Path","Median_Users","Median_Non_Users", "Number_users", "Number_non-users", "Non_zeros_User", "Non_zeros_Non_Users", "Wilcoxon Z (Users vs Non-users)","Wilcoxon effect size","p.value")
my_res2$q.value=p.adjust(as.numeric(as.character(my_res2$p.value)),method = "fdr")
my_res2$cohort="IBD"
res_IBD=my_res2


IBD_pheno=read.table("./LLD_phenos.txt",header = T,row.names = 1)
IBD=read.table("./lld_filt_ab.txt", header = T, check.names = F)
IBDn=read.table("./LLD_filt_ab.txt", header = T, check.names = F)
IBD=unique(IBD)
IBDn=unique(IBDn)
rownames(IBD)=IBD$Gene_Family
rownames(IBDn)=IBDn$Gene_Family
IBD$Gene_Family=NULL
IBDn$Gene_Family=NULL
IBD_n=as.data.frame(t(IBDn))
test_n=merge(IBD_pheno, IBD_n, by="row.names")
IBD=as.data.frame(t(normalize(IBD,samples.in.rows = F,to.abundance = F, transform = "log", move.log = "min")))
test=merge(IBD_pheno, IBD, by="row.names")

x=0
my_res=matrix(,ncol = 10, nrow = length(colnames(test)[49:ncol(test)]))

for (i in 49:ncol(test)){
	x=x+1
	my_res[x,2]=median(test_n[,i][test_n$antibiotics_merged=="User"])
  	my_res[x,3]=median(test_n[,i][test_n$antibiotics_merged!="User"])
  	my_res[,4]=sum(test_n$antibiotics_merged=="User")
  	my_res[,5]=sum(test_n$antibiotics_merged!="User")
  	my_res[x,6]=sum(test_n[,i][test_n$antibiotics_merged=="User"]!=0)
  	my_res[x,7]=sum(test_n[,i][test_n$antibiotics_merged!="User"]!=0)
  	#Wilcoxon effect size
	my_res[x,8]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE))
	my_res[x,9]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE)) / sqrt (nrow(test))
	asoc=wilcox.test(test[,i]~test$antibiotics_merged)
	my_res[x,1]=colnames(test)[i]
	my_res[x,10]=asoc$p.value
}
my_res2=as.data.frame(my_res[my_res[,4]!="NaN", ])
colnames(my_res2)=c("Path","Median_Users","Median_Non_Users", "Number_users", "Number_non-users", "Non_zeros_User", "Non_zeros_Non_Users", "Wilcoxon Z (Users vs Non-users)","Wilcoxon effect size","p.value")
my_res2$q.value=p.adjust(as.numeric(as.character(my_res2$p.value)),method = "fdr")
my_res2$cohort="LLD"
res_LLD=my_res2



IBD_pheno=read.table("./MIBS_phenos.txt",header = T,row.names = 1)
IBD=read.table("./mibs_filt_ab.txt", header = T, check.names = F)
IBDn=read.table("./MIBS_filt_ab.txt", header = T, check.names = F)
IBD=unique(IBD)
IBDn=unique(IBDn)
rownames(IBD)=IBD$Gene_Family
rownames(IBDn)=IBDn$Gene_Family
IBD$Gene_Family=NULL
IBDn$Gene_Family=NULL
IBD_n=as.data.frame(t(IBDn))
test_n=merge(IBD_pheno, IBD_n, by="row.names")
IBD=as.data.frame(t(normalize(IBD,samples.in.rows = F,to.abundance = F, transform = "log", move.log = "min")))
test=merge(IBD_pheno, IBD, by="row.names")

x=0
my_res=matrix(,ncol = 10, nrow = length(colnames(test)[50:ncol(test)]))

for (i in 50:ncol(test)){
	x=x+1
	my_res[x,2]=median(test_n[,i][test_n$antibiotics_merged=="User"])
  	my_res[x,3]=median(test_n[,i][test_n$antibiotics_merged!="User"])
  	my_res[,4]=sum(test_n$antibiotics_merged=="User")
  	my_res[,5]=sum(test_n$antibiotics_merged!="User")
  	my_res[x,6]=sum(test_n[,i][test_n$antibiotics_merged=="User"]!=0)
  	my_res[x,7]=sum(test_n[,i][test_n$antibiotics_merged!="User"]!=0)
  	#Wilcoxon effect size
	my_res[x,8]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE))
	my_res[x,9]=statistic(wilcox_test(test[,i]~test$antibiotics_merged, conf.int=TRUE)) / sqrt (nrow(test))
	asoc=wilcox.test(test[,i]~test$antibiotics_merged)
	my_res[x,1]=colnames(test)[i]
	my_res[x,10]=asoc$p.value
}
my_res2=as.data.frame(my_res[my_res[,4]!="NaN", ])
colnames(my_res2)=c("Path","Median_Users","Median_Non_Users", "Number_users", "Number_non-users", "Non_zeros_User", "Non_zeros_Non_Users", "Wilcoxon Z (Users vs Non-users)","Wilcoxon effect size","p.value")
my_res2$q.value=p.adjust(as.numeric(as.character(my_res2$p.value)),method = "fdr")
my_res2$cohort="MIBS"
res_MIBS=my_res2

rownames(res_MIBS)=res_MIBS$Path
rownames(res_LLD)=res_LLD$Path
rownames(res_IBD)=res_IBD$Path

res_MIBS$Path=NULL
res_LLD$Path=NULL
res_IBD$Path=NULL

res_1=merge(res_MIBS,res_LLD, by="row.names", all=T)
rownames(res_1)=res_1$Row.names
res_1$Row.names=NULL
res_2=merge(res_1,res_IBD, by="row.names", all=T)
res_3=res_2
res_3$q.value.x=p.adjust(as.numeric(as.character(res_3$p.value.x)),method = "fdr")
res_3$q.value.y=p.adjust(as.numeric(as.character(res_3$p.value.y)),method = "fdr")
res_3$q.value=p.adjust(as.numeric(as.character(res_3$p.value)),method = "fdr")
```
