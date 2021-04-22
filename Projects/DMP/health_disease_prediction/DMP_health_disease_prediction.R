# ==================================================
# By: Weersma Group, UMCG (2020)
# DMP project, training and testing of
# models for prediction of health and diseases
# ==================================================

## libraries and functions
library(glmnet)
library(pROC)
library(vegan)

## function do_clr_extrnalWeighting
#
# function performs CLR transformation 
#  using external weighting
# ===========================================
do_clr_externalWeighting = function(interest_matrix, core_matrix){
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  
  # estimate weighting parameter
  gm_mean = function(x, na.rm=TRUE){
    exp(sum(log(x), na.rm=na.rm) / length(x))
  }
  Gmean_core = apply(core_matrix, 1, gm_mean)
  
  #do transformation
  data_prepared = cbind(Gmean_core,interest_matrix)
  data_transformed = t(apply(data_prepared,1,function(x){
    log(x / x[1])[-1]
  }))
  colnames(data_transformed) = colnames(data_transformed)
  rownames(data_transformed) = rownames(data_transformed)
  data_transformed
}

# ===========================
# ===========================
# MAIN
# ===========================
# ===========================

# load all necessary data
covar = read.table("Association.analysis.input/covariates.txt")
taxa = read.table("Association.analysis.input/taxa.txt")
shannon.div = diversity(taxa[,grep("[.]s__",colnames(taxa))],index = "shannon")
pathways = read.table("Association.analysis.input/pathways_metacyc.txt")
load("Association.analysis.input/phenos.RData")
dis2 = read.table("diseases.final.txt",as.is = T)

# select diseases
pheno_disease = pheno[,dis2[,1]]

# transform data using CLR
taxa_transformed = do_clr_externalWeighting(taxa,taxa[,grep("[.]s__",colnames(taxa))])
taxa_transformed = taxa_transformed[,colSums(taxa>0)>nrow(taxa) * 0.05]
pathways_transformed = do_clr_externalWeighting(pathways,pathways)
pathways_transformed = pathways_transformed[,colSums(pathways>0)>nrow(pathways) * 0.05]

## split data into training and test sets
set.seed(12348)
data.pred = data.frame(covar[,1:2],BMI = pheno[,1],shannon = shannon.div,taxa_transformed,pathways_transformed)
train.set = sample(1:nrow(data.pred),size = round(0.9 * nrow(data.pred)))
train.x = as.matrix(data.pred[train.set,])
test.x = as.matrix(data.pred[-train.set,])

# clinical parameters (age, sex, BMI)
clin.train.x = train.x;clin.train.x[,4:ncol(clin.train.x)] = 0
clin.test.x = test.x;clin.test.x[,4:ncol(clin.test.x)] = 0

# microbiome factors (taxa and pathways)
microb.train.x = train.x;microb.train.x[,1:3] = 0
microb.test.x = test.x;microb.test.x[,1:3] = 0

# outcome
train.y = pheno_disease[train.set,]
test.y = pheno_disease[-train.set,]

## use cross-validated glmnet for training and parameter selection
# =================================================================

# (prediction of health == no diseases)
# ==============================================
model.health = cv.glmnet(x=train.x,
                         y = train.y[,"MED.DISEASES.None.No.Diseases"],
                         nfold = 10,
                         alpha = 0.5,
                         penalty.factor = c(0,0,0,rep(1,ncol(train.x)-3)),
                         family = "binomial")

# use models to predict results on training and test sets
train.prediction.full.min = predict(model.health,newx = train.x,s= "lambda.min")
train.prediction.full.1se = predict(model.health,newx = train.x,s= "lambda.1se")
  
test.prediction.full.min = predict(model.health,newx = test.x,s= "lambda.min")
test.prediction.full.1se = predict(model.health,newx = test.x,s= "lambda.1se")

train.prediction.clin.min = predict(model.health,newx = clin.train.x,s= "lambda.min")
train.prediction.clin.1se = predict(model.health,newx = clin.train.x,s= "lambda.1se")

test.prediction.clin.min = predict(model.health,newx = clin.test.x,s= "lambda.min")
test.prediction.clin.1se = predict(model.health,newx = clin.test.x,s= "lambda.1se")

train.prediction.microb.min = predict(model.health,newx = microb.train.x,s= "lambda.min")
train.prediction.microb.1se = predict(model.health,newx = microb.train.x,s= "lambda.1se")

test.prediction.microb.min = predict(model.health,newx = microb.test.x,s= "lambda.min")
test.prediction.microb.1se = predict(model.health,newx = microb.test.x,s= "lambda.1se")

# create result tables
# > training set
predictions.train = cbind(train.prediction.microb.min,
                    train.prediction.microb.1se,
                    train.prediction.clin.min,
                    train.prediction.clin.1se,
                    train.prediction.full.min,
                    train.prediction.full.1se)
# > test set
predictions.test = cbind(test.prediction.microb.min,
                    test.prediction.microb.1se,
                    test.prediction.clin.min,
                    test.prediction.clin.1se,
                    test.prediction.full.min,
                    test.prediction.full.1se
                    )

# calculate AUCs for no-disease (healthy) phenotype
aucs.train = apply(predictions.train,2,
      function(x){auc(train.y[,"MED.DISEASES.None.No.Diseases"],x)})
aucs.test = apply(predictions.test,2,
                   function(x){auc(test.y[,"MED.DISEASES.None.No.Diseases"],x)})

names(aucs.train) = c("microbes.min","microbes.1se","clin.min","clin.1se","full.min","full.1se")
names(aucs.test) = c("microbes.min","microbes.1se","clin.min","clin.1se","full.min","full.1se")

##calculate AUCs for prediction of diseases, parallel implementation
library(doSNOW)
cl = makeCluster(8)
registerDoSNOW(cl)
disease_pred_models = foreach(i = 1:36) %dopar% {
  complete.x = train.x[!is.na(train.y[,i]),]
  complete.y = train.y[!is.na(train.y[,i]),]
  cv1 = glmnet::cv.glmnet(complete.x,complete.y[,i],alpha = 0.5,nfolds = 5, family = "binomial")
  
}

disease_preds_microbes_train = do.call(cbind,lapply(disease_pred_models,function(x){s =train.x;s[,1:3] = 0;predict(x,newx = s,s = "lambda.min")[,1]}))
disease_preds_clin_train = do.call(cbind,lapply(disease_pred_models,function(x){s =train.x;s[,4:608] = 0;predict(x,newx = s,s = "lambda.min")[,1]}))
disease_preds_full_train = do.call(cbind,lapply(disease_pred_models,function(x){s =train.x;predict(x,newx = s,s = "lambda.min")[,1]}))
disease_preds_microbes_test = do.call(cbind,lapply(disease_pred_models,function(x){s =test.x;s[,1:3] = 0;predict(x,newx = s,s = "lambda.min")[,1]}))
disease_preds_clin_test = do.call(cbind,lapply(disease_pred_models,function(x){s =test.x;s[,4:608] = 0;predict(x,newx = s,s = "lambda.min")[,1]}))
disease_preds_full_test = do.call(cbind,lapply(disease_pred_models,function(x){s =test.x;predict(x,newx = s,s = "lambda.min")[,1]}))

disease_aucs = foreach(i = 1:36,.combine = rbind)%do%{
  data.frame(train.microbes = auc(train.y[,i],disease_preds_microbes_train[,i]),
             train.clin = auc(train.y[,i],disease_preds_clin_train[,i]),
             train.full = auc(train.y[,i],disease_preds_full_train[,i]),
             test.microbes = auc(test.y[,i],disease_preds_microbes_test[,i]),
             test.clin = auc(test.y[,i],disease_preds_clin_test[,i]),
             test.full = auc(test.y[,i],disease_preds_full_test[,i])
  )
             
};rownames(disease_aucs) = colnames(pheno_disease)[1:36]

aucs_byHealth_lambdaMin = foreach(i = 1:36,.combine = rbind)%do%{
  data.frame(train.microbes = auc(train.y[,i],predictions.train[,1]),
             train.clin = auc(train.y[,i],predictions.train[,3]),
             train.full = auc(train.y[,i],predictions.train[,5]),
             test.microbes = auc(test.y[,i],predictions.test[,1]),
             test.clin = auc(test.y[,i],predictions.test[,3]),
             test.full = auc(test.y[,i],predictions.test[,5])
  )
  
};rownames(aucs_byHealth_lambdaMin) = colnames(pheno_disease)[1:36]

aucs_byHealth_lambda1se = foreach(i = 1:36,.combine = rbind)%do%{
  data.frame(train.microbes = auc(train.y[,i],predictions.train[,2]),
             train.clin = auc(train.y[,i],predictions.train[,4]),
             train.full = auc(train.y[,i],predictions.train[,6]),
             test.microbes = auc(test.y[,i],predictions.test[,2]),
             test.clin = auc(test.y[,i],predictions.test[,4]),
             test.full = auc(test.y[,i],predictions.test[,6])
  )
  
};rownames(aucs_byHealth_lambda1se) = colnames(pheno_disease)[1:36]
