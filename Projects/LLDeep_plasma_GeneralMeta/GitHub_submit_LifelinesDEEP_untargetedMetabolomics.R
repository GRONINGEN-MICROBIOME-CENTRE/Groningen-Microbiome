####main codes for Lifelines plasma untargeted metabolomics (Metabolome platform) study
###for questions, please contact LianminChen (lianminchen@yeah.net)

############################################################  statistics  #############################################################
###diet, snp and microbiome estimate metabolite variance by using linear model with lasso for feature selection 
#input 'asso_diet','asso_microbiome','asso_genetics' are metabolite associations to dietary, genetic and microbial features. 'data' and 'data_rep' are data frames that contain metabolites, genetics, microbiome and dietary habits from all participants
lasso_variance=function(tmp, data,data_rep){
  if(length(tmp)==0){
    fit=c(0,0,1,0,1)
  }else{
    if(length(tmp)==1){
      tmp=data[,c(meta[i],tmp)]
      fit=summary(lm(tmp[,1]~.,data.frame(data.frame(tmp[,2]))))
      tmp_rep=data_rep[,colnames(tmp)]
      fit_rep=summary(lm(tmp_rep[,1]~.,data.frame(data.frame(tmp_rep[,2]))))
      if(length(unique(tmp_rep[,2]))==1){
        fit=c(1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,1)
      }else{
        fit=c(1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,pf(fit_rep$fstatistic[1],fit_rep$fstatistic[2],fit_rep$fstatistic[3],lower.tail = F))
      }
    }else{
      if(length(tmp)>1){
        tmp=data[,c(meta[i],tmp)]
        cv=glmnet::cv.glmnet(as.matrix(tmp[,2:ncol(tmp)]), tmp[,1], alpha=1, nfolds=10, type.measure='mse')
        if(length(which(as.vector(coef(cv, s='lambda.min'))!=0))==1){
          fit=c(0,0,1,0,1)
        }else{
          tmp=tmp[,which(as.vector(coef(cv, s='lambda.min'))!=0)]
          fit=summary(lm(tmp[,1]~.,data.frame(tmp[,2:ncol(tmp)])))
          tmp_rep=data_rep[,colnames(tmp)]
          fit_rep=summary(lm(tmp_rep[,1]~.,data.frame(tmp_rep[,2:ncol(tmp_rep)])))
          fit=c(ncol(tmp)-1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,pf(fit_rep$fstatistic[1],fit_rep$fstatistic[2],fit_rep$fstatistic[3],lower.tail = F))
        }
      }
    }
  }
  return(fit)
}
result=data.frame(meta=NA,n_diet=NA,r2.adj_diet=NA,r2.adj_diet_p=NA,r2.adj_diet_rep=NA,r2.adj_diet_rep_p=NA,n_microbe=NA,r2.adj_microbe=NA,r2.adj_microbe_p=NA,r2.adj_microbe_rep=NA,r2.adj_microbe_rep_p=NA,n_snp=NA,r2.adj_snp=NA,r2.adj_snp_p=NA,r2.adj_snp_rep=NA,r2.adj_snp_rep_p=NA)
meta=colnames(data)[1:1183]
for(i in 1:length(meta)){
  tmp=c(meta[i],lasso_variance(tmp=unique(as.character(asso_diet$diet[which(asso_diet$meta==meta[i])])),data = data,data_rep = data_fup),
        lasso_variance(tmp=unique(as.character(asso_microbiome$microbiome[which(asso_microbiome$meta==meta[i])])),data = data,data_rep = data_fup),
        variance_snp(tmp=unique(as.character(asso_genetics$SNPName[which(asso_genetics$ProbeName==meta[i])])),data = data,data_rep = data_fup))
  result=rbind(result,tmp)
}

###metabolite predict dietary quality score
#input 'tmp' is a data frame taht contains all metabolites associated with dietary quality score,'data' and 'data_rep' are data frames that contain metabolites and dietary quality score from all participants
lasso_variance_diet_score=function(tmp, data,data_rep){
  if(length(tmp)==0){
    fit=c(0,0,1,0,1)
  }else{
    if(length(tmp)==1){
      tmp=data[,c(meta[i],tmp)]
      fit=summary(lm(tmp[,1]~.,data.frame(data.frame(tmp[,2]))))
      tmp_rep=data_rep[,colnames(tmp)]
      fit_rep=summary(lm(tmp_rep[,1]~.,data.frame(data.frame(tmp_rep[,2]))))
      if(length(unique(tmp_rep[,2]))==1){
        fit=c(1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,1)
      }else{
        fit=c(1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,pf(fit_rep$fstatistic[1],fit_rep$fstatistic[2],fit_rep$fstatistic[3],lower.tail = F))
      }
    }else{
      if(length(tmp)>1){
        tmp=data[,c(meta[i],tmp)]
        cv=glmnet::cv.glmnet(as.matrix(tmp[,2:ncol(tmp)]), tmp[,1], alpha=1, nfolds=10, type.measure='mse')
        if(length(which(as.vector(coef(cv, s='lambda.min'))!=0))==1){
          fit=c(0,0,1,0,1)
        }else{
          tmp=tmp[,which(as.vector(coef(cv, s='lambda.min'))!=0)]
          fit=summary(lm(tmp[,1]~.,data.frame(tmp[,2:ncol(tmp)])))
          model=lm(tmp[,1]~.,data.frame(tmp[,2:ncol(tmp)]))
          tmp_rep=data_rep[,colnames(tmp)]
          pre_value=predict(model,tmp_rep)
          fit_rep=summary(lm(tmp_rep[,1]~.,data.frame(tmp_rep[,2:ncol(tmp_rep)])))
          fit=c(ncol(tmp)-1,fit$adj.r.squared,pf(fit$fstatistic[1],fit$fstatistic[2],fit$fstatistic[3],lower.tail = F),fit_rep$adj.r.squared,pf(fit_rep$fstatistic[1],fit_rep$fstatistic[2],fit_rep$fstatistic[3],lower.tail = F))
        }
      }
    }
  }
  return(list(fit,pre_value))
}

###mQTL mapping (Genome wide association for metabolites, mGWAS)
#QTL mapping is performed by calculating Spearman correlation between SNP dosage and plasma metabolites
#detailed describtion of the JAVA based QTL mapping tool is here: https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline

###bi-directional Mendelian randomization analysis
#'mr_37microbial_abundance_45metabolites.txt' contains associations between genetics, micorbiome and metabolites
library(TwoSampleMR)
mr=read.delim("../result/mr_37microbial_abundance_45metabolites.txt",header = T,sep = "\t")
exposure=unique(as.character(mr$ProbeName))
outcome=unique(as.character(mr$meta))
result_mr=NULL
for(i in 1:length(exposure)){
  tmp=mr[which(mr$ProbeName==as.character(exposure[i])),]
  for(j in 1:length(unique(as.character(tmp$meta)))){
    tmp_specific=tmp[which(tmp$meta==unique(as.character(tmp$meta))[j]),]
    #discovery
    data.exposure=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$IncludedDatasetsCorrelationCoefficient,se=rep(1/sqrt(933-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.exposure=format_data(data.exposure, type="exposure")
    data.outcome=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$outcome_baseline_r,se=rep(1/sqrt(933-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.outcome=format_data(data.outcome, type="outcome")
    data=harmonise_data(exposure_dat =data.exposure,outcome_dat = data.outcome)
    data$mr_keep='TRUE'
    data$mr_keep=as.logical(data$mr_keep)
    result.dis=mr(data)
    result.dis$outcome=unique(as.character(tmp$meta))[j]
    result.dis$exposure=as.character(exposure)[i]
    result.dis$heterogeneity=NA
    result.dis$pleiotropy=NA
    result.dis$heterogeneity[c(1,3)]=mr_heterogeneity(data)[,8]
    result.dis$pleiotropy[1]=mr_pleiotropy_test(data)[1,7]
    colnames(result.dis)[7:11]=paste("discovery_",colnames(result.dis)[7:11],sep = "")
    #replication
    data.exposure=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$LLD_fup_rep_r,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.exposure=format_data(data.exposure, type="exposure")
    data.outcome=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$outcome_fup_r,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.outcome=format_data(data.outcome, type="outcome")
    data=harmonise_data(exposure_dat =data.exposure,outcome_dat = data.outcome)
    data$mr_keep='TRUE'
    data$mr_keep=as.logical(data$mr_keep)
    result.rep=mr(data)
    colnames(result.rep)[7:9]=paste("rep_",colnames(result.rep)[7:9],sep = "")
    result=cbind(result.dis[,c(4,3,5:11)],result.rep[,7:9])
    result_mr=rbind(result_mr,result)
  }
}

###mediation analysis
library(mediation)
#mediation
colnames(data)=c("X","Y","M")
model.m=lm(M~X,data)
model.y=lm(Y~X+M,data)
summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
#inverse mediation
colnames(data)=c("X","M","Y")
model.m=lm(M~X,data)
model.y=lm(Y~X+M,data)
summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))