# load data 

metaphlan =read.table("alex_stuff/clean/data/microbe/metaphlan_raw.txt",header=T,as.is = T,sep="\t")
pheno = read.table("alex_stuff/clean/data/phenotypes/allPhenotypes_clear.txt",header=T,as.is = T,sep="\t")

library(vegan)

metaphlan = metaphlan[,colSums(metaphlan)>0]

metaphlan_filtered = metaphlan[,-grep("t__",colnames(metaphlan))]
library(vegan)

# wilcoxon test / 
wilcox_pvals = c()
for (i in 1:ncol(metaphlan_filtered)){
  w1 = wilcox.test(metaphlan_filtered[,1] ~ pheno$antrop_gender.F1M2)
  wilcox_pvals = c(wilcox_pvals,w1$p.value) 
}
wilcox_results = data.frame(taxon = colnames(metaphlan_filtered),P = wilcox_pvals) 

#log transformation
metaphlan_filtered_addedLittle = metaphlan_filtered + min(metaphlan_filtered[metaphlan_filtered>0])/2

metaphlan_log = log(metaphlan_filtered_addedLittle,base = 10)


# external function
do_clr_externalWeighting = function(interest_matrix, core_matrix) {
  if(any(interest_matrix==0)) interest_matrix = interest_matrix + min(interest_matrix[interest_matrix>0])/2
  if(any(core_matrix==0)) core_matrix = core_matrix + min(core_matrix[core_matrix>0])/2
  #estimate weighting parameter
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

metaphlan_clr = do_clr_externalWeighting(metaphlan_filtered, 
                                         metaphlan_filtered[,grep("s__",colnames(metaphlan_filtered))])




#linear regression with abudance
results = data.frame()

for(i in 1:10){
  for(j in 3:10) {
  lm1 = lm(metaphlan_clr[,i] ~ pheno$antrop_age+pheno$antrop_gender.F1M2 + pheno[,j])
  summary1 = summary(lm1)
  summary1$coefficients[4,1:4]
  oneReport = data.frame(bac = colnames(metaphlan_clr)[i],pheno = colnames(pheno)[j],
             coef = summary1$coefficients[4,1],
             SE = summary1$coefficients[4,2],
             T = summary1$coefficients[4,3],
             P = summary1$coefficients[4,4]
             )
  results = rbind(results,oneReport)
  }
}
results$FDR = p.adjust(results$P,method = "bonferroni")



# alpha diversity

metaphlan_species = metaphlan_filtered[,grep("s__",colnames(metaphlan_filtered))]

library(foreach)

#foreach vesrions
alpha_divs = foreach(div = c("shannon","simpson","invsimpson"),.combine = cbind)%do%
  { diversity(metaphlan_species,index = div)
    
  }
  
#for
alpha_divs_for = matrix(nrow = nrow(metaphlan_species),ncol = 3)
rownames(alpha_divs_for) = rownames(metaphlan_species)
colnames(alpha_divs_for) = c("shannon","simpson","invsimpson")

for(i in c("shannon","simpson","invsimpson")){
  alpha_divs_for[,i] = diversity(metaphlan_species,index = i)
  
}

#inverse rank function
invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}


alpha_transformed = apply(alpha_divs_for,2,invrank)

results_alpha = data.frame()
for(i in 3:3){
  for(j in 3:10) {
    lm1 = lm(alpha_transformed[,i] ~ pheno[,j]+ pheno$antrop_age+pheno$antrop_gender.F1M2)
    summary1 = summary(lm1)
    summary1$coefficients[4,1:4]
    oneReport = data.frame(bac = colnames(alpha_transformed)[i],pheno = colnames(pheno)[j],
                           coef = summary1$coefficients[4,1],
                           SE = summary1$coefficients[4,2],
                           T = summary1$coefficients[4,3],
                           P = summary1$coefficients[4,4]
    )
    results_alpha = rbind(results_alpha,oneReport)
  }
}
results_alpha$FDR = p.adjust(results$P,method = "bonferroni")


#beta diversity
metaphlan_species = metaphlan_species[
  setdiff(
  rownames(metaphlan_species),c("LLDeep_0071", "LLDeep_0076","LLDeep_0593")),]

pheno = pheno[setdiff(
  rownames(metaphlan_species),c("LLDeep_0071", "LLDeep_0076","LLDeep_0593")),]

betadist = vegdist(metaphlan_species,method = "bray")

ad1 = adonis(betadist ~ pheno$antrop_age+pheno$antrop_gender.F1M2 + pheno$antrop_BMI)
ad2 = adonis(betadist ~ pheno$antrop_gender.F1M2+ pheno$antrop_age+ pheno$antrop_BMI)

adX = adonis(betadist ~ pheno$antrop_age+pheno$antrop_gender.F1M2)
adX = adonis(betadist ~ pheno$antrop_gender.F1M2+pheno$antrop_age)


library(doMC)
registerDoMC(10)
ad1 = adonis(betadist ~ pheno$antrop_age+pheno$antrop_gender.F1M2 + pheno$antrop_BMI,parallel = T)


ad1$aov.tab

#PCoA plots
meta_f = metaphlan_species[-which(rowSums(metaphlan_species) ==0),]
pheno_f = pheno[-which(rowSums(metaphlan_species) ==0),]

betadist = vegdist(meta_f)


pcoa1 = capscale(betadist ~ 1)
pcoa2 = capscale(meta_f ~ 1, distance = "bray")

envfit(pcoa2 ~ pheno_f$antrop_gender.F1M2 + pheno_f$antrop_age)

pcoa2constrained = capscale(meta_f ~ pheno_f$antrop_age+pheno_f$antrop_gender.F1M2, distance = "bray")


#NMDS
nmds1 = metaMDS(meta_f)
plot(envfit(nmds1 ~ pheno_f$antrop_age+pheno_f$antrop_gender.F1M2))






