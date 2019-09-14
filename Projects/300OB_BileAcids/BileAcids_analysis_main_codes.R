####main codes for 300OB cohort bile acids analysis (for questions, please contact Lianmin Chen, lianminchen@yeah.net)
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


###estimate inter-individual variations of bile acid explained by predictors (use micorbiome as an example, i.e., included LASSO to selecte features) by using 5-fold cross validation (CV) 
predictor=microbiome #a matrix with subjects in row while microbial features in column
responder=ba #a matrix with subjects in row while bile acid entities in column
library(foreach)
library(glmnet)
nfolds=5
set.seed(666)
folds = split(sample(nrow(predictor)),1:nfolds)
check=NULL
threshold=c(0.05,0.01,0.005,0.001,0.0005,0.0001,0.00005,0.00001)
for(p in threshold){
  results=matrix(NA,nrow = ncol(responder),ncol = 9)
  colnames(results)=c("meta","r2_cv","sd_r2_cv","cv_1","cv_2","cv_3","cv_4","cv_5")
  for(i in 1:ncol(responder)){
    results[i,1]=colnames(responder)[i]
    test.pred=NULL#
    for(n in 1:5){
      predictor_p=predictor[-unlist(folds[n],use.names=FALSE),]
      responder_p=responder[-unlist(folds[n],use.names=FALSE),i]
      predictor_v=predictor[unlist(folds[n],use.names=FALSE),]
      responder_v=responder[unlist(folds[n],use.names=FALSE),i]
      pick=NULL
      for(j in 1:ncol(predictor_p)){
        try(cor<-cor.test(responder_p,predictor_p[,j],method = "spearman"),silent = T)
        try(pick<-c(pick,cor$p.value),silent = T)
      }
      if(length(which(pick<p))>0){
        tmp_pre=predictor_p[,which(pick<p)]
        tmp_val=predictor_v[,which(pick<p)] 
        if(length(which(pick<p))==1){
          Data.train = cbind(tmp_pre,responder_p)
          Data.test  = cbind(tmp_val,responder_v)
          lmf = lm(responder_p~.,data.frame(Data.train))
          tmp_val=data.frame(tmp_val)
          colnames(tmp_val)="tmp_pre"
          test.pred=c(test.pred,(1-sum((predict(lmf,tmp_val)-responder_v)^2)/sum((responder_v-mean(responder_v))^2))*100)
        }else{
          try(cv<-glmnet::cv.glmnet(as.matrix(tmp_pre), responder_p, alpha=1,type.measure='mse'),silent = T)#a LASSO step is included for to select microbial features (for genetic feature selection, we applied LD clumping) or lets say reduce colinearity
          try(tmp_pre<-tmp_pre[,which(coef(cv, s='lambda.min')!=0)],silent = T)
          try(tmp_val<-tmp_val[,which(coef(cv, s='lambda.min')!=0)],silent = T)
          Data.train = cbind(tmp_pre,responder_p)
          Data.test  = cbind(tmp_val,responder_v)
          if(length(which(coef(cv, s='lambda.min')!=0))==1){
            lmf = lm(responder_p~.,data.frame(Data.train))
            tmp_val=data.frame(tmp_val)
            colnames(tmp_val)="tmp_pre"
            test.pred=c(test.pred,(1-sum((predict(lmf,tmp_val)-responder_v)^2)/sum((responder_v-mean(responder_v))^2))*100)
          }else{
            lmf = lm(responder_p~.,data.frame(Data.train,check.names = F))
            test.pred=c(test.pred,(1-sum((predict(lmf,tmp_val)-responder_v)^2)/sum((responder_v-mean(responder_v))^2))*100)
          }
    }
      }else{test.pred=c(test.pred, NA)}
    }
    if(length(which(pick<p))>0){
      results[i,2]=mean(test.pred,na.rm = T)
      res=responder[,i]
      results[i,3]=sd(test.pred,na.rm = T)
      results[i,4:8]=test.pred
    }
  }
  check=cbind(check,results[,2])
}
row.names(check)=results[,1]
colnames(check)=threshold
colnames(check)=paste(colnames(check),"_p")
check=cbind(row.names(check),check)
colnames(check)[1]="meta"
write.table(check,file = "Vm_ba_all_p_threshold.txt",sep = "\t",quote = F,row.names = F)
