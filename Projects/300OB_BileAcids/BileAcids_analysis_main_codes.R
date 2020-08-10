#####main codes for 300OB cohort bile acids analysis (for questions, please contact Lianmin Chen, lianminchen@yeah.net)
#1.QTL mapping (GWAS)
#2.microbiome wide association (MWAS)
#3.LD clumping


###1.QTL mapping (Genome wide association, GWAS)
#QTL mapping is performed by calculating Spearman correlation between SNP dosage and bile acid entities.
#detailed describtion of the JAVA based QTL mapping tool is here: https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline


###2.micorbiome bile acid entities associations (MWAS)
#pairwise Spearman correlation is used to check association between micorbial features and bile acid entities (use 1000x permutation to calculate FDR)
cor_spearman_permutation=function(microbiome, metabolites,n){
  microbiome=data.frame(microbiome,check.names = F)
  metabolites=data.frame(metabolites,check.rows = F)
  microbiome=microbiome[intersect(row.names(microbiome),row.names(metabolites)),]
  metabolites=metabolites[row.names(microbiome),]
  pvals = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  cors = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  number = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  for(i in 1:ncol(microbiome)){
    for(j in 1:ncol(metabolites)){
      a=microbiome[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,i]
      b=metabolites[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,j]
      try(cor<-cor.test(a,b,method = "spearman"),silent = T)
      try(pvals[j,i]<-cor$p.value,silent = T)
      try(cors[j,i]<-cor$estimate,silent = T)
      number[j,i] = length(a)
    }
  }
  pvals[is.na(pvals)]=1
  cors[is.na(cors)]=0
  colnames(cors) = colnames(microbiome)
  rownames(cors) = colnames(metabolites)
  ind = which(pvals <=1,arr.ind = T)
  association = data.frame(metabolites = colnames(metabolites)[ind[,1]],microbiome = colnames(microbiome)[ind[,2]],Coefficient_r = cors[ind],P_value = pvals[ind],Q_value = NA,number = number[ind],stringsAsFactors = F)
  association$name=paste(association[,1],association[,2],sep = "_with_")
  per_p=matrix(NA,nrow = nrow(association),ncol = n)
  colnames(per_p)=1:n
  per_cor=matrix(NA,nrow = nrow(association),ncol = n)
  colnames(per_cor)=1:n
  row.names(per_cor)=association$name
  per_number=matrix(NA,nrow = nrow(association),ncol = n)
  colnames(per_number)=1:n
  row.names(per_number)=association$name
  for(k in 1:n){
    print(k)
    per_sp=matrix(sample(as.matrix(microbiome)),nrow = nrow(microbiome),ncol = ncol(microbiome))
    per_path=matrix(sample(as.matrix(metabolites)),nrow = nrow(metabolites),ncol = ncol(metabolites))
    pvals=matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
    cors=matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
    number=matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
    for(i in 1:ncol(per_sp)){
      for(j in 1:ncol(per_path)){
        a=per_sp[is.na(per_sp[,i])==F&is.na(per_path[,j])==F,i]
        b=per_path[is.na(per_sp[,i])==F&is.na(per_path[,j])==F,j]
        try(cor<-cor.test(a,b,method = "spearman"),silent = T)
        try(pvals[j,i]<-cor$p.value,silent = T)
        try(cors[j,i]<-cor$estimate,silent = T)
        try(number[j,i]<-length(a),silent = T)
      }
    }
    pvals[is.na(pvals)]=1
    cors[is.na(cors)]=0
    ind = which(pvals<=1,arr.ind = T)
    per_p[,k]=pvals[ind]
    per_cor[,k]=cors[ind]
    per_number[,k]=number[ind]
  }
  association=association[order(association[,4],decreasing = F),]
  for(i in 1:nrow(association)){
    association[,5][i]=(length(which(per_p<=association[,4][i]))/n)/i
  }
  return(association)
}

###3.select independent QTLs (LD clumping (PLINK clumping method),500kb,r2=0.05) for variance prediction by using 5-fold CV
library(devtools)
library(TwoSampleMR)
trait=colnames(ba) #names of 39 bile acid entities
cv1=read.table("man_cv1.txt",header = T,sep = "\t",fill = T) #QTLs mapping results (GWAS results)
cv1=cv1[which(cv1$PValue<1e-5),] #only run clumping on QTLs with P<1e-5, i.e., pass genome wide suggestive line
result=NULL
for(i in 1:length(trait)){
  tmp=cv1[which(cv1$ProbeName==trait[i]),]
  tmp=tmp[,c(2,3,4,5,1)] #make it in a good order to suit the LD clumping function
  colnames(tmp)=c("SNP","chr_name","chrom_start","id.exposure","pval.exposure")
  clump=clump_data(tmp, clump_kb = 500, clump_r2 = 0.05, clump_p1 = 1,clump_p2 = 1)
  result=rbind(result,clump)
}
