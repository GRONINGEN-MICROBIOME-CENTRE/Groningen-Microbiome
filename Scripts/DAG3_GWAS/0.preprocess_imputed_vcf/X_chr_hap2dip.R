

#################form vcf#####################
library(data.table)
library(tidyverse)
library(optparse)

###corrected loop nrow(0) no loop
###corrected loop form n to ncol (not nrow)


#### note: input is an output from this command
# bcftools +fixploidy X.DAG3_ID_imputed_info_filtered.vcf.gz  -- -s samples_sex -p ploidy.file > X.DAG3.dip.vcf
## which already converted haploid GT into diploid
## Also, the header of vcf file was cutted with these commands
# grep -n "#CHROM" X.DAG3.dip.vcf |cut -d":" -f1
# head -n 114 X.DAG3.dip.vcf >header.vcf

## Then we got the number of total snps
# wc -l  X.DAG3.dip.vcf
### 551872-114=551758 total snps

#### extracting the last snp
###tail -n 10 X.DAG3.dip.vcf |less -S -x25
##154929412
#### extracting the first snp
#less -S -x25 X.DAG3.dip.vcf |head -n 116
### first snp position 2699968 
##snp difference is 154929412-2699968 =152229444

###calculon
wkdir<-"/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/X_proc"

##gearshift
#wkdir<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf/X_proc"



################################################################################


## first we have to cut the file into chunks ...of 100K
##lets call the total number of SNPs
header.file<-file.path(wkdir,"header.vcf")
input.xfile<-file.path(wkdir,"X.DAG3.dip.vcf.gz")
sex.file<-file.path(wkdir,"samples_sex")

samples.sex<-fread(sex.file,data.table = F,header=T)

numvar=154929412
m<-0
lim<-0
while (lim < numvar){
  lim<-ifelse(1000000*(m+1)+2699968>numvar,numvar,1000000*(m+1)+2699968)
  #### correct the nfrom
  nfrom<-(1000000*(m)+2699968)
  cutted.x.file<-file.path(wkdir,paste0("X_part",m+1,".vcf"))
  
  ###cut with bcftools
  bcfcut.call<-paste0("ml BCFtools; ",
                      " bcftools view -r X:",nfrom,"-",lim, " -o ",cutted.x.file," ",input.xfile)
  system(bcfcut.call)
  
  ##read the cutted file
  vcfx<-fread(cutted.x.file, skip="#CHROM",data.table=F )
  #n<-10
  if (nrow(vcfx)!=0){
    
    for (n in 10:ncol(vcfx)){
      ###condition for males only
      dagID<-colnames(vcfx)[n]
      if(samples.sex[samples.sex$ID==dagID,2]=="M"){
        indivgen<-as.vector(vcfx[,c(n)])
        splitted<-strsplit(indivgen, ":")
        splitted2<-lapply(splitted,FUN= function(x) {
          if(x[1]=="1/1"){
            x[2]<-2
          }
          if(x[1]=="0/0"){
            x[2]<-0
          }
          paste0(c(x[1],x[2],x[3]),collapse=":")
        } )
        indivgen2<-unlist(splitted2)
        vcfx[,c(n)]<-indivgen2
      }
    }
    file.temp.out<-file.path(wkdir,paste0("dip_X_part",m+1,"temp.vcf"))
    write.table(vcfx, file.temp.out,quote = F,sep = '\t',row.names = F)
    
    file.part.out<-file.path(wkdir,paste0("dip_X_part",m+1,".vcf"))
    
    ### Add vcf header bcftools
    add.header.call<-paste0("cat ",header.file," ", file.temp.out , " > ",file.part.out)
    system(add.header.call)
    ### index vcf file
    index.call<-paste0("ml BCFtools; ",
                       "bgzip -f"," ",file.part.out," ; ",
                       "bcftools  index ", file.part.out,".gz ;",
                       " echo ",file.part.out,".gz ", " >> ",wkdir,"/merge.list" )
    system(index.call)
    
  }
  
  m=m+1
  #  cat header.vcf DSdip.X.part.test.vcf >frank.vcf
  #  bgzip frank2.vcf 
  #  bcftools index frank2.vcf.gz 
  
}


merge.call<-paste0("ml BCFtools; ",
                   " bcftools concat -a -f ",wkdir,"/merge.list", " -o ", wkdir,"/X_DS_dip_merged.vcf")
system(merge.call)

#  bcftools concat -a  frank.vcf.gz frank2.vcf.gz -o merged_frank.vcf



### done