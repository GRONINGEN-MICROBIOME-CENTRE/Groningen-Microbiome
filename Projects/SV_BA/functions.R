## BA Categories
primary.group <- c('CA','GCA','TCA','CDCA','GCDCA', 'TCDCA')
second.group  <- c('DCA','GDCA','TDCA','LCA','GLCA','TLCA') # 'TLCA_3S','GLCA_3S'
urso.group    <- c('UDCA', 'GUDCA', 'TUDCA')

ca.group      <- c('CA', "GCA", 'TCA', 'DCA','GDCA', 'TDCA')
cdca.group    <- c('CDCA','GCDCA', 'TCDCA',"LCA", 'GLCA', 'TLCA') # 'TLCA_3S','GLCA_3S'
udca.group    <- c('UDCA', 'GUDCA',"TUDCA")
unconjugated  <- c('CA','CDCA','DCA','LCA')
conjugated    <- c('TCA','GCA','TCDCA','GCDCA','TDCA', 'GDCA', 'TLCA','GLCA') # 'TLCA_3S','GLCA_3S'
ca.deconjugated    <- c('CA','TCA', 'GCA')
ca.dehydroxylation <- c('DCA','TDCA', 'GDCA')
#ca.deconjugated    <- c('DCA','TDCA', 'GDCA')
#ca.dehydroxylation <- c('CA','TCA', 'GCA')
glycine.group <- c('GCA', 'GCDCA','GDCA', 'GLCA') # 'GLCA_3S'
taurine.group <- c('TCA','TCDCA','TDCA','TLCA') # 'TLCA_3S'

primary.group_p <- paste(primary.group,'_p',sep = '')
second.group_p  <- paste(second.group, '_p',sep = '')
urso.group_p    <- paste(urso.group, '_p',sep = '')

concentration.group <- c(primary.group,second.group,urso.group)
proportion.group    <- c(primary.group_p,second.group_p,urso.group_p)

## color panels
mycolor2_blue_yellow <- c("#08B2E3","#FED766")
mycolor2_green_blue  <- c("#2EC4B6","#235789")
mycolor2_blue_red    <- c("#08B2E3","#ee6352")

mycolor3 <- c("#2EC4B6","#235789", "grey50") # c("#4f94cd","#ff4040", "grey50")
mycolor3_red_blue_yellow <- c("#ee6352","#08B2E3","#FED766")
mycolor3_green_red_yellow<- c("#2EC4B6","#ee6352","#FED766")

## libraty packages
library(tidyverse)
library(ggpubr)
library(GGally)
library(ggmosaic)
library(cowplot)
library(ggsci)
library(ggthemes)
library(gplots)
library(VennDiagram)
library(circlize)
library(viridis)
library(wesanderson)
library(ggalluvial)

library(reshape2)
library(Hmisc)
library(vegan)
library(microbiome)
library(mediation)

## Convert SV name
changeSVname<-function(SVrawid){
  testname     <- SVrawid
  species_name <- as.character(taxonomy$organism[match(str_replace_all(testname, "\\..*",""), taxonomy$X)])
  region       <- str_replace_all(testname, ".*\\:","") 
  region_list  <- str_split(region,";") 
  
  region_suf   <- NULL
  i<-1
  for (i in c(1:length(region_list))){
    if(length(region_list[[i]]) <= 2){
      region_suf[i] <- paste(":",region[i],sep = "")
    }else{
      region_suf[i] <- paste(":",region_list[[i]][1]," and ",length(region_list[[i]])-1," segments", sep = "")
    }
    i <- i+1
  }
  paste(species_name,region_suf,sep = "")
}


## Calculate SV size
calcSVSize <- function(SVrawid){
  testname <- SVrawid
  region   <- str_replace_all(testname, ".*\\:","") 
  region_list <- str_split(region,";") 
  
  sv_size  <- NULL
  for (i in c(1:length(region_list))) {
    frag_df    <- str_split(region_list[[i]], "_") %>% unlist %>% as.numeric  %>% matrix(byrow = T, ncol = 2) %>% as.data.frame 
    sv_size[i] <- sum(as.numeric(frag_df$V2)-as.numeric(frag_df$V1))
  }
  sv_size
}



shared_sv_dis<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    #inVec<-inList[[3]]
    Vec_i<-Vec_i
    insvdf<-data.frame(inVec,Vec_i) %>% na.omit
    sv_dis<- vegdist(t(insvdf), method = "canberra")
    #length(na.omit(Vec_i+inVec))
    return(sv_dis)
  }
  
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-inList[[i]]
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
}

## Calculate SE
se <- function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}


## Summary statistics
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=TRUE,conf.interval=.95, .drop=TRUE) {
  require(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N            = length2(xx[[col]], na.rm=na.rm),
                     nonZero_N    = sum(xx[[col]]!=0,na.rm = na.rm),
                     nonZero_rate = sum(xx[[col]]!=0,na.rm = na.rm)/length2(xx[[col]], na.rm=na.rm),
                     mean         = mean   (xx[[col]], na.rm=na.rm),
                     sd           = sd     (xx[[col]], na.rm=na.rm),
                     mindata      = min(xx[[col]], na.rm=na.rm),
                     maxdata      = max(xx[[col]], na.rm=na.rm),
                     quant1_  = quantile(xx[[col]], na.rm=na.rm)[2] ,
                     quant2_  = quantile(xx[[col]], na.rm=na.rm)[4] ,
                     quantmid = quantile(xx[[col]], na.rm=na.rm)[3] ,
                     mea      = mean   (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  
  ciMult   <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  datac$xd <- datac$mea-datac$se
  datac$xu<- datac$mea+datac$se
  return(datac)
}

## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}

## pie chart
my_pie<-function(in1,item,mycol = c("#F24D4D","#4D7EAE")){
  require(tibble)
  require(showtext)
  require(Cairo)
  require(ggsci)
  require(tibble)
  require(scales)
  require(ggrepel)
  require(forcats)
  require(scatterpie)
  
  showtext_auto()
  
  in1.tibble            <- as_tibble(subset(in1, in1$items%in%item))
  in1.tibble$categories <- fct_reorder(in1.tibble$categories, in1.tibble$value)
  in1.tibble            <- in1.tibble[order(in1.tibble$value, decreasing = TRUE), ]
  piepercent            <- round(100*in1.tibble$value/sum(in1.tibble$value), 2)
  my_labels             <- tibble(x.breaks = seq(1.3, 1.3, length.out = length(piepercent)),
                                  y.breaks = cumsum(in1.tibble$value) - in1.tibble$value/2,
                                  labels = paste(in1.tibble$categories, "\n",in1.tibble$value,", ",piepercent, "%", sep = ""),
                                  categories = in1.tibble$categories)
  
  pdfName    <- paste(item, ".pie.pdf", sep = "")
  ggplot(in1.tibble, aes(x = 1, y = value, fill = categories)) +
    ggtitle(paste(item)) +
    geom_bar(stat="identity", color='white') + 
    coord_polar(theta='y') + 
    theme(legend.position = "None",
          axis.ticks=element_blank(),  # the axis ticks
          axis.title=element_blank(),  # the axis labels
          axis.text.y=element_blank(), # the 0.75, 1.00, 1.25 labels.
          axis.text.x = element_blank(),
          panel.background = element_blank(),
          panel.grid = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    #    scale_fill_brewer(palette = "Set3", direction = -1)+
    scale_fill_manual(values=mycol) + # set fill color manually
    geom_text_repel(data = my_labels, 
                    aes(x = x.breaks, y = y.breaks, label = labels),
                    size = 7,        # label text size
                    show.legend = FALSE,
                    inherit.aes = FALSE,
                    arrow = arrow(length=unit(0.01, "npc")),
                    force = 1,
                    segment.color = 'grey50'
    )
}


## Batch adonis
my_adonis_terms<-function(inDist, inDf, covDf=NULL, covar=NULL){
  require(vegan)
  ## test
  #inDist <- lld_msv_dist_std[[1]]
  #inDf   <- lld_ba
  #covDf  <- lld_covar
  #ovar  <- covar1
  ## test
  
  covDf<- covDf[,covar]
  
  inDist_rmna <- dist_rmna(inDist)
  covDf_rmna   <- na.omit(covDf)
  
  my_adonis<-function(inMat, inVec, covdf = NULL){
    ## test
    #inMat<-inDist_rmna
    #inVec<-inDf[,1]
    #covdf<-covDf_rmna
    ## test
    inMat_covdf_inter <- intersect(rownames(inMat),rownames(covdf))
    inMat_covdf_inVec_inter <- intersect(rownames(inDf)[!is.na(inVec)], inMat_covdf_inter)
    
    inMat_rmna <- inMat[match(inMat_covdf_inVec_inter,rownames(inMat)),
                        match(inMat_covdf_inVec_inter,colnames(inMat))]
    inVec_rmna <- inVec[match(inMat_covdf_inVec_inter, rownames(inDf))]
    covdf_rmna <- covdf[match(inMat_covdf_inVec_inter, rownames(covdf)),]
    
    in_cov <- data.frame(inVec_rmna,covdf_rmna)
    
    sample_size<-nrow(in_cov)
    
    adonis_res <- adonis(as.dist(inMat_rmna)~.,data = in_cov,parallel = 4)
    
    return(list(adonis_res$aov.tab[1,5],adonis_res$aov.tab[1,6], sample_size))
  }
  
  adonis_res_batch <- apply(inDf, 2, my_adonis,inMat = inDist_rmna, covdf = covDf_rmna)
  adonis_res_batch.unlist <- matrix(unlist(adonis_res_batch), ncol = 3, byrow = T)
  colnames(adonis_res_batch.unlist)<-c("R2", "P", "N")
  
  return(adonis_res_batch.unlist)
}


my_adonis_terms_adjAbun<-function(distList, inDf, covDf, covar,info){
  #distList<-lld_msv_dist_std
  #inDf<-lld_ba
  #covDf<-lld_covar
  #covar<-covar1
  #info<-info
  
  
  res_table<-NULL
  for (i in 1:length(distList)) { # 
    #i<-2
    cat(paste(i,"\n")) 
    
    n<-names(distList)[i] %>% 
      str_replace_all("msv_","") %>% 
      match(., info$organism)
    
    inCovar<-c(covar, info$organism[n])

    
    adonis_res<-my_adonis_terms(distList[[i]], inDf, covDf = covDf, covar = inCovar)
    
    adonis_res_df<-data.frame(Species = rep(info$Short_name[n], dim(adonis_res)[1]),
                              BA = colnames(inDf),
                              as.data.frame(adonis_res))
    
    res_table<-rbind(res_table, adonis_res_df)
  }
  
  res_table$Species<-as.character(res_table$Species)
  res_table$BA<-as.character(res_table$BA)
  
  res_table$fdr<-p.adjust(res_table$P,method = 'fdr')
  res_table$bonferroni<-p.adjust(res_table$P,method = 'bonferroni')
  
  # R2
  r2_mat<-matrix(data  = res_table$R2,
                 nrow  = length(unique(res_table$Species)),
                 ncol  = length(unique(res_table$BA)),
                 byrow = T)
  
  rownames(r2_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(r2_mat)<-colnames(inDf)
  
  # P
  p_mat<-matrix(data  = res_table$P,
                nrow  = length(unique(res_table$Species)),
                ncol  = length(unique(res_table$BA)),
                byrow = T)
  
  rownames(p_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(p_mat)<-colnames(inDf)
  
  # fdr
  fdr_mat<-matrix(data  = res_table$fdr,
                  nrow  = length(unique(res_table$Species)),
                  ncol  = length(unique(res_table$BA)),
                  byrow = T)
  
  rownames(fdr_mat)<-names(distList)[1:length(distList)] %>% str_replace_all("msv_", "") %>% match(.,info$organism) %>% info$Short_name[.]
  colnames(fdr_mat)<-colnames(inDf)
  
  return(list(table = res_table,
              r2 = r2_mat,
              P = p_mat,
              FDR = fdr_mat))
}

my_meta_p <- function(inVec, study_name, n_col, p_col) {
  require(metap)
  #require(metafor)
  
  #inVec<-cbind_adonis_res.table[20,]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #_col<-c(9,14)
  
  study_n <- inVec[n_col] %>% as.numeric
  study_p <- inVec[p_col] %>% as.numeric
  
  inDf <- data.frame(study_name, study_n, study_p)
  
  Wi=sqrt(study_n)
  #Convert p-values to Z-scores
  Zi=qnorm(1-(study_p/2))
  # Meta-zscore  
  Z=(sum(Zi*Wi)/sqrt(sum(Wi^2)))
  # Convert Z-score to p-value
  MetaP=2*pnorm(-abs(Z))
  
  #Cochran Q-test (test heterogeneity of your meta-analysis)
  
  WeightedZ= sum(sqrt(study_n)*Zi)
  totalSample=sum(study_n)
  
  #Calculate expected Z
  expected_Z= sqrt(study_n)*WeightedZ/totalSample
  het_Sum= sum((Zi-expected_Z)*(Zi-expected_Z))
  
  #Get p-value of heterogeneity test!
  my_pvalue_co=pchisq(het_Sum, lower.tail = F, df=length(study_p)-1)
  
  return(
    list(
      Meta.p = MetaP,
      Meta.hetero = het_Sum,
      Meta.hetero.p = my_pvalue_co
    )
  )
}
## Batch meta-analysis
my_batch_meta_p <- function(inDf,study_name,n_col,p_col,row_var_col = 1,col_var_col = 2) {
  #inDf<-cbind_adonis_res.table[c(1:10),]
  #study_name<-c("LLD", "300OB")
  #n_col<-c(10,15)
  #p_col<-c(9,14)
  #row_var_col<-1
  #col_var_col<-2
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_p,
    study_name = study_name,
    n_col = n_col,
    p_col = p_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 3, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.p',  "Meta.hetero", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  rowName <- unique(inDf[, row_var_col])
  colName <- unique(inDf[, col_var_col])
  
  N_row <- length(rowName)
  N_col <- length(colName)
  
  
  batch_res.p <- matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
  colnames(batch_res.p) <- colName
  rownames(batch_res.p) <- rowName
  batch_res.p[is.na(batch_res.p)] <- 0
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  # fdr matrix
  batch_res.fdr <- matrix(
    data  = batch_res_edge$Meta.fdr.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  colnames(batch_res.fdr) <- colName
  rownames(batch_res.fdr) <- rowName
  batch_res.fdr[is.na(batch_res.fdr)] <- 1
  
  # bonferroni matrix
  batch_res.bon <- matrix(
    data  = batch_res_edge$Meta.bonferroni.p,
    nrow  = N_row,
    ncol  = N_col,
    byrow = T
  )
  
  colnames(batch_res.bon) <- colName
  rownames(batch_res.bon) <- rowName
  batch_res.bon[is.na(batch_res.bon)] <- 1
  
  return(
    list(
      table      = batch_res_edge,
      p          = batch_res.p,
      fdr        = batch_res.fdr,
      bonferroni = batch_res.bon
    )
  )
}


## The Normal Quantile Transformation
qtrans<-function(x){
  k<-!is.na(x)
  k<-which(x!="-999")
  ran<-rank(as.numeric(x[k]))
  y<-qnorm((1:length(k)-0.5)/length(k))
  x[k]<-y[ran]
  x
}



## linear model
lm_btw_mats<-function(y_mat,x_mat,cov_mat,covar){
  require(reshape2)
  require(R.utils)
  
  ## test block
  #y_mat<- y_mat_i #lld_intri[,c(1:3)]
  #x_mat<- x_mat #lld_vsv[,c(1:3)]
  #cov_mat<- cov_mat# lld_covar
  #covar<- covar_i #covar
  ## test block
  
  my_lm<-function(y,x){
   # y<-y_mat[,1]
  #  x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Y"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"X"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    try(lm_res <- summary(lm(Y~.,data = lm_input)), silent = T)
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 10, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2,1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}




## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, se_col) {
  require(meta)
  require(metafor)
  
  ## test start
  #inVec<-vsv_crass_lm_res.edge[100,]
  #study_name<-c("LLD", "300OB", "IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  
  ## test end
  
  study_beta <- inVec[beta_col] %>% as.numeric
  study_se <- inVec[se_col] %>% as.numeric
  
  #study_beta<-c(0.7, 0.65)
  #study_se<-c(0.07, 0.08)
  #stydy_n<-c(1000, 300)
  
  
  inDf <- data.frame(study_name, study_beta, study_se)
  
  m.hksj <- metagen(
    study_beta,
    study_se,
    data = inDf,
    studlab = study_name
  )
  
  m.hksj.res <- summary(m.hksj)
  
  return(
    list(
      Meta.beta = m.hksj.res$random$TE,
      Meta.se = m.hksj.res$random$seTE,
      Meta.p = m.hksj.res$random$p,
      Meta.I2 = m.hksj.res$I2$TE,
      Meta.hetero.p = m.hksj$pval.Q
    )
  )
}



## Batch meta-analysis
my_batch_meta_lm <- function(inDf,study_name,beta_col,se_col) {
  ## test start
  #inDf<-vsv_crass_lm_res.edge[c(1:10),]
  #study_name<-c("LLD", "300OB","IBD")
  #beta_col<-c(3,10,17)
  #se_col<-c(4,11,18)
  ## test end
  
  batch_res <- apply(
    inDf,
    1,
    my_meta_lm,
    study_name = study_name,
    beta_col = beta_col,
    se_col = se_col
  )
  
  batch_res.unlist <- matrix(unlist(batch_res), ncol = 5, byrow = T)
  colnames(batch_res.unlist) <-
    c('Meta.beta', 'Meta.SE', "Meta.p", "Meta.I2", "Meta.hetero.p")
  batch_res_edge <- cbind(inDf, batch_res.unlist)
  
  # add p adjust
  batch_res_edge <- data.frame(
    as.data.frame(batch_res_edge),
    Meta.fdr.p = p.adjust(batch_res_edge$Meta.p, method = "fdr"),
    Meta.bonferroni.p = p.adjust(batch_res_edge$Meta.p, method = "bonferroni")
  )
  
  return(batch_res_edge)
}


lm_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
  #y_mat<-ob_contin[,c(1:10)]
  #x_mat<-ob_vsv[,c(1:100)]
  #cov_mat<-ob_basic
  #covar<-covar
  #abun<-ob_abun_clr
  #info<-info
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun, info, abun_col){
  #y_mat<-vsgv_lld
  #x_mat<-lld_exp[,c(1:10)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun_clr
  #info<-info
  #abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    i<-1
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lm_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat, info, abun_col){
  #  y_mat<-lld_ba[,c(1:10)]
  #  x_mat<-vsgv_lld[,c(1:100)]
  #  cov_mat<-lld_basic
  #  covar<-covar
  #  abun<-lld_abun_clr
  #  pc_mat<-lld_msv_pc_cum0.6
  #  info<-info
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
#    i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))] 
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lm_btw_mats_adjAbunPCs2<-function(y_mat, x_mat, cov_mat, covar, abun, pc_mat, info, abun_col){
#    y_mat<-vsgv_lld[,c(1:10)]
#    x_mat<-lld_exp[,c(1:10)]
#    cov_mat<-lld_basic
#    covar<-covar
#    abun<-lld_abun_clr
#    pc_mat<-lld_msv_pc_cum0.6
#    info<-info
#    abun_col<-9
  
  cov_mat<-cbind(cov_mat, abun,pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
#    i<-2
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]%>% as.matrix
    colnames(y_mat_i)<-colnames(y_mat)[grep(spe_name,colnames(y_mat))]
    
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"Yes"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lm_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x_i$Adjust_Abundance<-"No"
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

## Get density of 2-demision dataset
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}


## linear model
lr_btw_mats<-function(y_mat,x_mat,cov_mat,covar, direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  
  require(reshape2)
  require(glmnet)
  require(R.utils)
  
  ## test block
  #y_mat <- lld_disea[,c(1:5)]
  #x_mat <- lld_vsv[,c(1:3)]
  #cov_mat <- lld_covar
  #covar   <- covar #covar
  #direction<-c(1,1)
  ## test block
  
  my_lm<-function(y,x){
    #y<-y_mat[,2]
    #x<-x_mat[,1]
    
    beta    <- NA
    se      <- NA
    p.value <- NA
    N       <- NA
    
    y_uniq_N        <- NA
    x_uniq_N        <- NA
    y_non_zero_N    <-NA
    x_non_zero_N    <-NA
    y_non_zero_rate <-NA
    x_non_zero_rate <-NA
    y_0_N<-NA
    y_1_N<-NA
    
    lm_input<-data.frame(Y = y, X = x,cov_mat[,covar]) %>% sapply(as.numeric) %>% na.omit
    N        <- nrow(lm_input)
    
    y_uniq_N   <- length(unique(lm_input[,"Y"]))
    x_uniq_N   <- length(unique(lm_input[,"X"]))
    
    y_non_zero_N <- sum(!isZero(lm_input[,"Gender"]))
    x_non_zero_N <- sum(!isZero(lm_input[,"Gender"]))
    
    y_non_zero_rate<-y_non_zero_N/N
    x_non_zero_rate<-x_non_zero_N/N
    
    y_0_N<-sum(lm_input[,"Y"]==0)
    y_1_N<-sum(lm_input[,"Y"]==1)
    
    lm_input_tmp <- apply(lm_input[,c(2:ncol(lm_input))], 2, qtrans) %>% as.data.frame
    lm_input<-data.frame(Y = lm_input[,1],lm_input_tmp)
    
    try(lm_res <- summary(glm(Y~., data = lm_input, family = 'binomial')))
    
    indv<-'X'
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta,
                    se = se,
                    p.value = p.value,
                    N = N,
                    y_uniq_N     = y_uniq_N,
                    x_uniq_N     = x_uniq_N,
                    y_non_zero_N = y_non_zero_N,
                    x_non_zero_N = x_non_zero_N,
                    y_non_zero_rate = y_non_zero_rate,
                    x_non_zero_rate= x_non_zero_rate,
                    y_0_N = y_0_N,
                    y_1_N = y_1_N)),
        silent = T)
    
  }
  
  
  y_x<-sapply( 
    as.data.frame(y_mat),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(x_mat)
    )
  )
  y_x.unlist <- matrix(unlist(y_x), ncol = 12, byrow = T)
  
  # beta matrix
  y_x.beta<-matrix(y_x.unlist[,1],ncol = ncol(y_mat), byrow = F)
  colnames(y_x.beta)<-colnames(y_mat)
  rownames(y_x.beta)<-colnames(x_mat)
  y_x_edge_beta  <- melt(y_x.beta)
  
  y_x_edge<-data.frame(y_x_edge_beta[,c(2:1)],
                       as.data.frame(y_x.unlist),
                       fdr.p = p.adjust(y_x.unlist[,3], method = "fdr"),
                       bonferroni.p = p.adjust(y_x.unlist[,3], method = "bonferroni"))
  
  
  colnames(y_x_edge)<-c("Y", "X", "Beta","SE", "p","N","y_uniq_N","x_uniq_N", "y_non_zero_N", "x_non_zero_N", "y_non_zero_rate","x_non_zero_rate","y_0_N","y_1_N","fdr.p","bonferroni.p")
  
  return(y_x_edge)
}


lr_btw_mats_adjAbun<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
  #y_mat<-lld_intri[,c(1:10)]
  #x_mat<-lld_vsv[,c(1:100)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun
  #info<-info
  #direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lr_btw_mats_adjAbun2<-function(y_mat, x_mat, cov_mat, covar, abun,info, abun_col){
  #y_mat<-lld_intri[,c(1:10)]
  #x_mat<-lld_vsv[,c(1:100)]
  #cov_mat<-lld_basic
  #covar<-covar
  #abun<-lld_abun
  #info<-info
  #direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}


lr_btw_mats_adjAbunPCs<-function(y_mat, x_mat, cov_mat, covar, abun,pc_mat,info, abun_col){
 # y_mat<-lld_intri[,c(1:10)]
#  x_mat<-lld_vsv[,c(1:100)]
# cov_mat<-lld_basic
#  covar<-covar
#  abun<-lld_abun_clr
#  pc_mat<-
#  info<-info
#  direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun, pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    x_mat_i<-x_mat[,grep(spe_name,colnames(x_mat))]
    if(dim(x_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,info[i,abun_col])
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat, x_mat_i, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}

lr_btw_mats_adjAbunPCs2<-function(y_mat, x_mat, cov_mat, covar, abun,pc_mat,info, abun_col){
  # y_mat<-lld_intri[,c(1:10)]
  #  x_mat<-lld_vsv[,c(1:100)]
  # cov_mat<-lld_basic
  #  covar<-covar
  #  abun<-lld_abun_clr
  #  pc_mat<-
  #  info<-info
  #  direction<-c(1,1)
  
  cov_mat<-cbind(cov_mat, abun, pc_mat)
  
  y_x.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-19
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% str_replace_all("\\]","\\\\\\]")
    y_mat_i<-y_mat[,grep(spe_name,colnames(y_mat))]%>% as.matrix
    colnames(y_mat_i)<-colnames(y_mat)[grep(spe_name,colnames(y_mat))]
    
    if(dim(y_mat_i)[2]>0){
      if(info[i,abun_col]%in%colnames(abun)){
        covar_i<-c(covar,colnames(cov_mat)[grep(spe_name, colnames(cov_mat))])
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }else{
        covar_i<-covar
        y_x_i <- lr_btw_mats(y_mat_i, x_mat, cov_mat, covar_i)
        y_x.edge<-rbind(y_x.edge,y_x_i)
      }
    }
  }
  
  
  y_x.edge<-as.data.frame(y_x.edge)
  y_x.edge$fdr.p<-p.adjust(y_x.edge$p,method = 'fdr')
  y_x.edge$bonferroni.p<-p.adjust(y_x.edge$p,method = 'bonferroni')
  
  return(y_x.edge)
}



## Mediation analysis for linear model
my_lm_mediation<-function(input.inv,input.med, input.dv,  covDf){
  #input.inv<-lld_exp$melatonine
  #input.med<-lld_vsv$`Faecalibacterium cf. prausnitzii KLE1255:1373_1377`
  #input.dv <-lld_ba$CA_dehydro_deconju_ratio
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Faecalibacterium_prausnitzii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2]) # input.inv
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna)
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
    }
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
  }
  
  return(res_list)
}

## Bidirectional mediation analysis for linear model
lm_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-exp_vsv_ba_df[1,]
  #indvDf<-lld_exp
  #dvDf1<-lld_vsv
  #dvDf2<-lld_ba
  #covDf<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  # test block
  
  if(is.na(inVec[4])){
    covar<-covar
  }else{
    covar<-c(covar, colnames(covDf)[grep(inVec[4],colnames(covDf))])
  }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lm_mediation(indv, dv1, dv2, covDf[,covar])
  dir2_res <- my_lm_mediation(indv, dv2, dv1, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}


## Mediation analysis for logistic model
my_lr_mediation_1<-function(inpu.inv, input.dv, input.med, covDf){
  #input.inv<-indv
  #input.med<-dv2
  #input.dv <-dv1
  #covDf<- covDf[,covar]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  med<-input.df.rmna$input.med
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.med<-med
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=lm(input.dv~.,input.df.rmna[,-3])
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=glm(input.med~., input.df.rmna[,-2], family = "binomial")
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=lm(input.dv~.,input.df.rmna) 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}


my_lr_mediation_2<-function(inpu.inv, input.dv, input.med, covDf){
  #input.inv<-lld_exp$cereals
  #input.med<-lld_ba$GDCA
  #input.dv <-lld_dsv$`[Eubacterium] hallii DSM 3353:2969_2983`
  #covDf<-lld_covar[,c('Gender','Age','BMI','Reads_number','s__Eubacterium_hallii')]
  
  input.df<-data.frame(input.inv,input.dv,input.med,covDf)
  input.df.rmna<-na.omit(input.df)
  dv<-input.df.rmna$input.dv
  input.df.rmna<-apply(input.df.rmna, 2, qtrans) %>% as.data.frame
  input.df.rmna$input.dv<-dv
  
  if(length(table(input.df.rmna$input.inv))>1 &
     length(table(input.df.rmna$input.dv))>1 &
     length(table(input.df.rmna$input.med))>1
  ){
    if(length(table(input.df.rmna$input.inv)) < 3 &
       sort(table(input.df.rmna$input.inv),decreasing = T)[2] <= 10){
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=NA,    inv_med.p=NA,
                       med_dv.beta=NA,     med_dv.p=NA,
                       inv_dv.beta=NA,     inv_dv.p=NA,
                       ACME.beta=NA,       ACME.p=NA,
                       ADE.beta=NA,        AED.p = NA,
                       Total.effect=NA,    Total.effet.p=NA,
                       Prop.mediated = NA, Prop.mediated.p = NA)
      
    }else{
      fit.totaleffect=glm(input.dv~.,input.df.rmna[,-3], family = "binomial")
      fit.totaleffect.res<-summary(fit.totaleffect)
      
      fit.mediator=lm(input.med~., input.df.rmna[,-2])
      fit.mediator.res<-summary(fit.mediator)
      
      fit.dv=glm(input.dv~.,input.df.rmna, family = "binomial") 
      fit.dv.res<-summary(fit.dv)
      
      results <- mediate(fit.mediator, fit.dv, covariates = colnames(covDf),
                         treat='input.inv', mediator='input.med', boot=T)
      res<-summary(results)
      
      res_list <- list(N = nrow(input.df.rmna),
                       inv_med.beta=fit.mediator.res$coefficients[2,1],   inv_med.p=fit.mediator.res$coefficients[2,4],
                       med_dv.beta=fit.dv.res$coefficients[3,1],          med_dv.p=fit.dv.res$coefficients[3,4],
                       inv_dv.beta=fit.totaleffect.res$coefficients[2,1], inv_dv.p=fit.totaleffect.res$coefficients[2,4],
                       ACME.beta=res$d0,                                  ACME.p=res$d0.p,
                       ADE.beta=res$z0,                                   AED.p = res$z0.p,
                       Total.effect=res$tau.coef,                         Total.effet.p=res$tau.p,
                       Prop.mediated = res$n0,                            Prop.mediated.p = res$n0.p)
      
    }
    
  }else{
    res_list <- list(N = nrow(input.df.rmna),
                     inv_med.beta=NA,    inv_med.p=NA,
                     med_dv.beta=NA,     med_dv.p=NA,
                     inv_dv.beta=NA,     inv_dv.p=NA,
                     ACME.beta=NA,       ACME.p=NA,
                     ADE.beta=NA,        AED.p = NA,
                     Total.effect=NA,    Total.effet.p=NA,
                     Prop.mediated = NA, Prop.mediated.p = NA)
    
  }
  
  return(res_list)
}





lr_bimediation<-function(inVec, indvDf, dvDf1, dvDf2, covDf, covar ){
  ## test block
  #inVec<-exp_dsv_ba_df[1,]
  #indvDf<-lld_exp
  #dvDf1<-lld_dsv
  #dvDf2<-lld_ba
  #covDf<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  # test block
  
  if(is.na(inVec[4])){
    covar<-covar
  }else{
    covar<-c(covar, inVec[4])
  }
  
  indv<-indvDf[,match(inVec[1], colnames(indvDf))]
  dv1<-dvDf1[,match(inVec[2], colnames(dvDf1))]
  dv2<-dvDf2[,match(inVec[3], colnames(dvDf2))]
  
  dir1_res <- my_lr_mediation_1(indv, dv2, dv1, covDf[,covar])
  dir2_res <- my_lr_mediation_2(indv, dv1, dv2, covDf[,covar])
  
  names(dir1_res)<-paste("dir1.",names(dir1_res),sep = "")
  names(dir2_res)<-paste("dir2.",names(dir2_res),sep = "")
  
  MediationDirection<-"none"
  if(!is.na(dir1_res$dir1.Prop.mediated.p) & !is.na(dir2_res$dir2.Prop.mediated.p)){
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "both"}
    if( dir1_res$dir1.Prop.mediated.p<0.05 &  dir2_res$dir2.Prop.mediated.p>0.05){MediationDirection <- "indv_dv1_dv2"}
    if( dir1_res$dir1.Prop.mediated.p>0.05 &  dir2_res$dir2.Prop.mediated.p<0.05){MediationDirection <- "indv_dv2_dv1"}
  }
  
  bires<-list(indv=inVec[1], dv1=inVec[2], dv2=inVec[3],MediationDirection = MediationDirection)
  
  res<-c(bires,dir1_res,dir2_res)
  
  return(res)
}


get_PCs<-function(inDist, eig.cutoff = 0.6){
  #inDist <- all_msv_dist_std[[1]]
  #eig.cutoff <- 0.6
  #eig.cutoff<-as.numeric(eig.cutoff)
  
  pcoa_res <- cmdscale(inDist, k=30, eig = T)
  
  total_eig <- 0
  i <- 0
  
  while (total_eig<eig.cutoff & i<10) {
    i<-i+1
    total_eig <- total_eig + pcoa_res$eig[i]/100
  }
  
  if(i<=1){
    pcoa<-data.frame(X1 = pcoa_res$points[, 1])
  }else{
    pcoa <- data.frame(pcoa_res$points)[,c(1:i)]
  }
  
  return(list(PCoA= pcoa,PC_num = i,Total_eig = total_eig))
  
}


my_uniq<-function(x){
  length(unique(x))
}



## Recalculate gene position
reCalcuPos<-function(spe, spe_scaf){
  spe_gene<-spe
  
  spe_gene<-data.frame(Taxonomy_id = str_replace_all(spe$contig_id, "\\..*", ""),
                       Contig_id = str_replace_all(spe$contig_id, ".*\\.", ""),
                       spe)
  Base_length<-NULL
  Base_length[1]<-0
  con_id<-spe_gene[1,2]
  for (i in 1:nrow(spe_gene)) {
    if(spe_gene[i,7] > spe_gene[i, 8]){
      tmp<-spe_gene[i,7]
      spe_gene[i,7]<-spe_gene[i, 8]
      spe_gene[i,8]<-tmp
    }
    
    if(i > 1){
      if( spe_gene[i,2]==con_id){
        Base_length[i]<-Base_length[i-1]
      }else{
        Base_length[i]<-spe_scaf$Length[match(con_id,spe_scaf$Contig)]+Base_length[i-1]
        con_id<-spe_gene[i,2]
      }
    }
  }
  
  spe_gene<-data.frame(spe_gene,
                       Base_length = Base_length,
                       Total_start = spe_gene$start + Base_length,
                       Total_stop  = spe_gene$stop  + Base_length)
  
  
}

