# HAllA analysis

# ======================================================================================================================================
# halla files 
# ======================================================================================================================================

# genus
genus=read.table("OutputTable/CLR.genus.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F)
covariate_bac=read.table("OutputTable/Covariate_bac.organized.txt",sep = "\t",header = T,stringsAsFactors = T,row.names = 1)

covariate_bac=as.data.frame(t(covariate_bac))
genus=as.data.frame(t(genus))

covariate_bac=covariate_bac[,order(colnames(covariate_bac))]
genus=genus[,order(colnames(genus))]

write.table(covariate_bac,file = "AHllA/Covariate.genus.txt",sep = "\t",row.names = T,quote = F)
write.table(genus,file = "AHllA/CLR.genus.txt",sep = "\t",row.names = T,quote = F)
