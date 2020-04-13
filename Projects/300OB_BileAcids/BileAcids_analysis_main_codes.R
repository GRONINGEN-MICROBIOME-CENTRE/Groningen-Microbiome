#####main codes for 300OB cohort bile acids analysis (for questions, please contact Lianmin Chen, lianminchen@yeah.net)
#1.QTL mapping (GWAS)
#2.microbiome wide association (MWAS)
#3.LD clumping
#4.cross validation based variance estimation


###QTL mapping (Genome wide association, GWAS)
#QTL mapping is performed by calculating Spearman correlation between SNP dosage and bile acid entities.
#detailed describtion of the JAVA based QTL mapping tool is here: https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline


###micorbiome bile acid entities associations (MWAS)
#pairwise Spearman correlation is used to check association between micorbial features and bile acid entities


###select independent QTLs (LD clumping (PLINK clumping method),500kb,r2=0.05) for variance prediction by using 5-fold CV
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
