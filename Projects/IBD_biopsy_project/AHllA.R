# HAllA analysis

# ======================================================================================================================================
# halla files for gearshift run
# ======================================================================================================================================

# genus
genus=read.table("OutputTable/CLR.genus.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)
genus=genus[rownames(genus) %in% rownames(covariate_bac),]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(genus),]

str(covariate_bac)
covariate_bac=covariate_bac[,c("Diagnosis","Inflammation","Location_rough","age_at_biopsy","sex","BMI","smoking_DB",
                               "Aminosalicylates","Thiopurines",
                               "Steroids","Methotrexaat","resec_part_colon","resec_ileocec","resec_part_small")]
covariate_bac$CD.vs.control=NA
covariate_bac$CD.vs.control[covariate_bac$Diagnosis=="CD"]=1
covariate_bac$CD.vs.control[covariate_bac$Diagnosis=="Control"]=0

covariate_bac$UC.vs.control=NA
covariate_bac$UC.vs.control[covariate_bac$Diagnosis=="UC"]=1
covariate_bac$UC.vs.control[covariate_bac$Diagnosis=="Control"]=0

covariate_bac$CD.vs.UC=NA
covariate_bac$CD.vs.UC[covariate_bac$Diagnosis=="CD"]=1
covariate_bac$CD.vs.UC[covariate_bac$Diagnosis=="UC"]=0

covariate_bac$IBD.vs.Control=NA
covariate_bac$IBD.vs.Control[covariate_bac$Diagnosis!="Control"]=1
covariate_bac$IBD.vs.Control[covariate_bac$Diagnosis=="Control"]=0

covariate_bac$Diagnosis=NULL
covariate_bac$Location_rough=as.character(covariate_bac$Location_rough)
covariate_bac$Location_rough[covariate_bac$Location_rough=="colon"]=1
covariate_bac$Location_rough[covariate_bac$Location_rough=="ileum"]=0
covariate_bac$Location_rough=as.numeric(covariate_bac$Location_rough)
str(covariate_bac)

covariate_bac=as.data.frame(t(covariate_bac))
genus=as.data.frame(t(genus))

covariate_bac=covariate_bac[,order(colnames(covariate_bac))]
genus=genus[,order(colnames(genus))]

write.table(covariate_bac,file = "AHllA/Covariate.genus.txt",sep = "\t",row.names = T,quote = F)
write.table(genus,file = "AHllA/CLR.genus.txt",sep = "\t",row.names = T,quote = F)

# ======================================================================================================================================
# halla comparisons between 1000IBD, HMP2
# ======================================================================================================================================
umcg=read.table("AHllA/all_associations.txt",sep = "\t",header = T)
hmp=read.table("AHllA/HMP2/all_associations.txt",sep = "\t",header = T)

group_umcg=umcg[umcg$Y_features=="CD.vs.UC",]
group_hmp=hmp[hmp$Y_features=="CD.vs.UC",]

group_umcg=group_umcg[order(group_umcg$q.values),]
group_umcg=group_umcg[1:10,]
group_hmp=group_hmp[group_hmp$X_features %in% group_umcg$X_features,]
















