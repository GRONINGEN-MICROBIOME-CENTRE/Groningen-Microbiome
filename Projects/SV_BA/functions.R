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

library(ggthemes)
library(GGally)
library(ggpubr)
library(ggord)
library(ggrepel)
library(gplots)
library(ggalluvial)
library(circlize)
library(VennDiagram)
library(networkD3)
library(ggmosaic)
library(gg.gap)
library(rbokeh)
library(dplyr)
library(cowplot)
library(ggsci)
library(viridis)
library(wesanderson)
library(RColorBrewer)
library(scales)

library(MASS)
library(ade4)
library(vegan)
library(phytools)
library(reshape2)
library(meta)
library(metafor)
library(mediation)
library(corrplot)
library(ggcorrplot2)
library(microbiome)

library(R.utils)
library(Cairo)
library(extrafont)
library(gridExtra)
library(showtext)
library(beepr)
library(future.apply)

## General setting
showtext_auto()

## Calculate SE
se <- function(x) {x<-na.omit(x);sqrt(var(x)/length(x))}

## Get density of 2-demision dataset
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
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




get_residual_lm<-function(y, x){
  #y<-vsgv_i[,1]
  #x<-abun_i
  
  xy<-data.frame(y,x)
  xy_norm <- apply(xy, 2, qtrans) %>% as.data.frame
  
  lm_res <- lm(y~x,data = xy_norm)
  xy_norm$y.residuals<-rep(NA, nrow(xy_norm))
  xy_norm$y.residuals[match(names(lm_res$residuals),row.names(xy_norm))]<-lm_res$residuals
  return(xy_norm$y.residuals)
}

get_residual_lr<-function(y, x){
  #y<-dsgv_i[,1]
  #x<-abun_i
  
  xy<-data.frame(y,x)
  xy_norm <- xy
  xy_norm$x<-qtrans(x)
  
  lr_res <- glm(y~x, data = xy_norm, family = 'binomial')
  
  xy_norm$y.residuals<-rep(NA, nrow(xy_norm))
  xy_norm$y.residuals[match(names(lr_res$residuals),row.names(xy_norm))]<-lr_res$residuals
  return(xy_norm$y.residuals)
}


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

## Prepare panphlan profile
prepData <- function(inDf,getLevel="S",isHMP2=F,minBac=95,doDropZero=T,doRescale=T) {
  if (isHMP2) {
    taxRaw <- inDf
  } else {
    taxRaw <- as.data.frame(t(inDf))
  }
  if (doDropZero) {
    taxF <- taxRaw[rowSums(taxRaw) > 0,]
  } else {
    taxF <- taxRaw
  }
  
  #  taxF <- taxF[taxRaw$k__Bacteria >= minBac,]
  # take only bacteria
  taxF <- taxF[,grep('k__Bacteria',colnames(taxF))]
  # take only species
  taxF <- taxF[,grep('t__',colnames(taxF),invert = T)]
  taxF.spec <- taxF[,grep('s__',colnames(taxF))]
  # take only genera
  taxF <- taxF[,grep('s__',colnames(taxF),invert = T)]
  taxF.gen <- taxF[,grep('g__',colnames(taxF))]
  # take only families
  taxF <- taxF[,grep('g__',colnames(taxF),invert = T)]
  taxF.family <- taxF[,grep('f__',colnames(taxF))]
  # take only orders
  taxF <- taxF[,grep('f__',colnames(taxF),invert = T)]
  taxF.order <- taxF[,grep('o__',colnames(taxF))]
  # take only classes
  taxF <- taxF[,grep('o__',colnames(taxF),invert = T)]
  taxF.class <- taxF[,grep('c__',colnames(taxF))]
  # take only phyla
  taxF <- taxF[,grep('c__',colnames(taxF),invert = T)]
  taxF.phyla <- taxF[,grep('p__',colnames(taxF))]
  # rescale
  if (getLevel=="S") {taxF.final <- taxF.spec
  } else if (getLevel=="G") {taxF.final <- taxF.gen
  } else if (getLevel=="F") {taxF.final <- taxF.family
  } else if (getLevel=="O") {taxF.final <- taxF.order
  } else if (getLevel=="C") {taxF.final <- taxF.class
  } else if (getLevel=="P") {taxF.final <- taxF.phyla}
  if (doRescale) {
    taxF.final.rescaled <- taxF.final/rowSums(taxF.final)
  } else {
    taxF.final.rescaled <- taxF.final
  }
  # return
  colnames(taxF.final.rescaled) <- gsub('\\|','.',colnames(taxF.final.rescaled))
  taxF.final.rescaled
}

## shorten the metaphlan2 taxa name
shortenNames <- function(inDF,sep = "\\|",direction = 2) {
  if (direction == 2) {
    for (c in grep('__', rownames(inDF))) {
      cn <- rownames(inDF)[c]
      cns <- unlist(strsplit(cn, sep))
      cnsf <- cns[length(cns)]
      rownames(inDF)[c] <- cnsf
    }
    
    
  } else{
    for (c in grep('__', colnames(inDF))) {
      cn <- colnames(inDF)[c]
      cns <- unlist(strsplit(cn, sep))
      cnsf <- cns[length(cns)]
      colnames(inDF)[c] <- cnsf
    }
  }
  inDF
}

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

## Get average distance matrix of multi distance matrices
mergeDist<-function(inDistList){
  require(matrixStats)
  tmp<-as.matrix(do.call(cbind, inDistList))
  tmp<-array(tmp, dim=c(dim(inDistList[[1]]), length(inDistList)))
  tmp<-colMeans(aperm(tmp, c(3, 1, 2)), na.rm = TRUE) # can be adjust
  row.names(tmp)<-row.names(inDistList[[1]])
  colnames(tmp)<-colnames(inDistList[[1]])
  tmp[is.nan(as.matrix(tmp))]<-NA
  tmp<-tmp/max(tmp,na.rm = T)
  return(tmp)
}

## Remove samples with NA in distance matrix
dist_rmna<-function(inDist){
  while(sort(colSums(is.na(inDist)),decreasing = T)[1] > 0){
    rmid<-names(sort(colSums(is.na(inDist)),decreasing = T)[1])
    inDist<-inDist[-match(rmid, rownames(inDist)),-match(rmid, colnames(inDist))]
  }
  return(inDist)
}

## convert distance matrix to GCTA grm file
make_grm<-function(inMat,marker_n){
  grm.bin<-data.frame(matrix(NA, nrow = nrow(inMat)*(nrow(inMat)+1)/2, ncol = 4))
  grm.bin[1,]<-c(1,1,marker_n,inMat[1,1])
  for (i in c(2:nrow(inMat))) {
    grm_start<-i*(i-1)/2+1
    grm_end  <-i*(i+1)/2
    grm.bin[c(grm_start:grm_end),1]<-i
    grm.bin[c(grm_start:grm_end),2]<-c(1:i)
    grm.bin[c(grm_start:grm_end),3]<-marker_n[i, 1:i]
    grm.bin[c(grm_start:grm_end),4]<-inMat[i,1:i]
  }
  return(grm.bin)
}

shared_marker_n<-function(inDf){
  #inDf<-all_sv
  inList <- split(as.matrix(inDf), seq(nrow(inDf)))
  
  shared_n_func<-function(inVec,Vec_i){
    shared_n<- sum(!is.na(inVec)&Vec_i)
      #length(na.omit(Vec_i+inVec))
    return(shared_n)
  }
  
  marker_n_mat<-matrix(NA, nrow = nrow(inDf), ncol = nrow(inDf))
  for (i in 1:nrow(inDf)) { #nrow(inDf)
    #i<-2
    Vec_i<-!is.na(inList[[i]])
    shared_n_i<-sapply(inList, shared_n_func,Vec_i = Vec_i)
    marker_n_mat[i,]<-shared_n_i
  }
  rownames(marker_n_mat)<-rownames(inDf)
  colnames(marker_n_mat)<-rownames(inDf)
  
  return(marker_n_mat)
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


## linear model
lm_btw_mats<-function(mat0,mat1,mat2,covar, direction = c(1,1),y_mat = 0){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  # y_mat: 0 means mat0 is y in linear model; 1 mean mat1 is y in linear model
  
  require(reshape2)
  
  ## test block
  #mat0<-lld_tmao[,c(1:3)]
  #mat1<-lld_vsv[,c(1:3)]
  #mat2<-lld_basic
  #covar<-covar1
  #direction<-c(1,1)
  #y_mat <- 0
  ## test block
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  
  col_var<-mat0
  row_var<-mat1
  
  my_lm<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,1]
    
    beta    <- NA
    p.value <- NA
    N       <- NA
    uniq_N  <- NA
    
    lm_input<-data.frame(phen = a,
                         mbio = b,
                         mat2[,covar])
    #
    #lm_input$mbio[!is.na(lm_input$mbio)]<-0
    #
    
    lm_input <- sapply(lm_input, as.numeric)
    
    lm_input <- na.omit(lm_input)
    N        <- nrow(lm_input)
    uniq_N   <- length(unique(lm_input[,1]))
    lm_input <- apply(lm_input, 2, qtrans) %>% as.data.frame
    
    
    if(y_mat==0){
      try(lm_res <- summary(lm(phen~.,data = lm_input)), silent = T)
      indv<-'mbio'
    }else{
      try(lm_res <- summary(lm(mbio~.,data = lm_input)), silent = T)
      indv<-'phen'
    }
    
    try(beta    <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),1],silent = T)
    try(se      <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),2], silent = T)
    try(p.value <- lm_res$coefficients[match(indv,rownames(lm_res$coefficients)),4],silent = T)
    
    
    try(return(list(beta = beta, se = se, p.value = p.value, N = N, uniq_N = uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_lm(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 5, byrow = T)
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # se matrix
  col_var_row_var.se<-matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.se)<-colnames(col_var)
  rownames(col_var_row_var.se)<-colnames(row_var)
  #col_var_row_var.se[is.na(col_var_row_var.se)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,3],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  
  
  # N matrix
  col_var_row_var.N <- matrix(col_var_row_var.unlist[,4],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.N)<-colnames(col_var)
  rownames(col_var_row_var.N)<-colnames(row_var)
  col_var_row_var.N[is.na(col_var_row_var.N)]<-0
  
  # uniq_N matrix
  col_var_row_var.uN <- matrix(col_var_row_var.unlist[,5],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.uN)<-colnames(col_var)
  rownames(col_var_row_var.uN)<-colnames(row_var)
  col_var_row_var.uN[is.na(col_var_row_var.uN)]<-0
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_se     <- melt(col_var_row_var.se)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  col_var_row_var_N          <- melt(col_var_row_var.N)
  col_var_row_var_uN         <- melt(col_var_row_var.uN)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,
                              col_var_row_var_edge_se,
                              col_var_row_var_edge_p,
                              col_var_row_var_N,
                              col_var_row_var_uN)[,-c(4,5,7,8,10,11,13,14)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "Beta","SE", "p","N","uniq_N")
  
  
  # add p adjust
  col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
                                   fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
                                   bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$fdr.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var_edge$bonferroni.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  
  colnames(col_var_row_var.bon)<-colnames(col_var)
  rownames(col_var_row_var.bon)<-colnames(row_var)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  return(list(table      = col_var_row_var_edge,
              beta       = col_var_row_var.beta,
              se         = col_var_row_var.se,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon,
              N          = col_var_row_var.N,
              uniq_N     = col_var_row_var.uN))
}

lm_btw_mats_adjAbun<-function(mat0, mat1, mat2, covar, abun,info,res_order,
                              direction=c(1,1),y_mat =0){
  #mat0<-lld_ba
  #mat1<-lld_vsv
  #mat2<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  #abun<-lld_abun_clr
  #info<-info
  #res_order<-paste(lld_vsv_ba_lm_res$table$Taxa,lld_vsv_ba_lm_res$table$Phenotype)
  #direction<-c(1,1)
  #y_mat <- 0
  
  col_var_row_var.edge<-NULL
  for (i in c(1:nrow(info))){
    #  i<-16
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% 
      str_replace_all("\\]","\\\\\\]")
    mat1_i<-mat1[,grep(spe_name,colnames(mat1))]
    if(dim(mat1_i)[2]>0){
      if(info$Metaphlan3_name[i]%in%colnames(abun)){
        covar_i<-c(covar,info$Metaphlan3_name[i])
        col_var_row_var_i <- lm_btw_mats(mat0, mat1_i, mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }else{
        covar_i<-covar
        col_var_row_var_i <- lm_btw_mats(mat0, mat1_i, mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }
    }
  }
  col_var_row_var.edge$fdr.p<-p.adjust(col_var_row_var.edge$p,method = 'fdr')
  col_var_row_var.edge$bonferroni.p<-p.adjust(col_var_row_var.edge$p,method = 'bonferroni')
  
  col_var_row_var.edge<-col_var_row_var.edge[match(res_order,
                                                          paste(col_var_row_var.edge$Taxa,col_var_row_var.edge$Phenotype)),]
  
  # beta matrix
  col_var_row_var.r<-matrix(data  = col_var_row_var.edge$Beta,
                            nrow  = ncol(mat1),
                            ncol  = ncol(mat0),
                            byrow = F)
  colnames(col_var_row_var.r)<-colnames(mat0)
  rownames(col_var_row_var.r)<-colnames(mat1)
  col_var_row_var.r[is.na(col_var_row_var.r)]<-0
  
  
  # p matrix
  col_var_row_var.p<-matrix(data  = col_var_row_var.edge$p,
                            nrow  = ncol(mat1),
                            ncol  = ncol(mat0),
                            byrow = F)
  colnames(col_var_row_var.p)<-colnames(mat0)
  rownames(col_var_row_var.p)<-colnames(mat1)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var.edge$fdr.p,
                              nrow  = ncol(mat1),
                              ncol  = ncol(mat0),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(mat0)
  rownames(col_var_row_var.fdr)<-colnames(mat1)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var.edge$bonferroni.p,
                              nrow  = ncol(mat1),
                              ncol  = ncol(mat0),
                              byrow = F)
  colnames(col_var_row_var.bon)<-colnames(mat0)
  rownames(col_var_row_var.bon)<-colnames(mat1)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  
  return(list(table      = col_var_row_var.edge,
              beta       = col_var_row_var.r,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon))
}




lm_btw_mats_adjAbun2<-function(mat0, mat1, mat2, covar, abun,info,
                              direction=c(1,1),y_mat =0){
  #mat0<-lld_vsv
  #mat1<-lld_exp
  #mat2<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  #abun<-lld_abun_clr
  #info<-info
  #direction<-c(1,1)
  #y_mat <- 0
  
  col_var_row_var.edge<-NULL
  for (i in c(1:nrow(info))){
    #i<-1
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% 
      str_replace_all("\\]","\\\\\\]")
    mat0_i<-mat0[,grep(spe_name,colnames(mat0))]
    
    if(dim(mat0_i)[2]>0){
      if(info$Metaphlan3_name[i]%in%colnames(abun)){
        covar_i<-c(covar,info$Metaphlan3_name[i])
        col_var_row_var_i <- lm_btw_mats(mat0_i,mat1 , mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }else{
        covar_i<-covar
        col_var_row_var_i <- lm_btw_mats(mat0_i, mat1, mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }
    }
  }
  col_var_row_var.edge$fdr.p<-p.adjust(col_var_row_var.edge$p,method = 'fdr')
  col_var_row_var.edge$bonferroni.p<-p.adjust(col_var_row_var.edge$p,method = 'bonferroni')
  
  return(list(table      = col_var_row_var.edge))
}

## logistic model
lr_btw_mats<-function(mat0,mat1,mat2,covar, direction = c(1,1),y_mat = 0){
  # mat0: binary variables, mat1: continuous variables, mat2: covar
  # direction: 1 means samples in row and variables in column; 2 means samples in column and variables in row
  # y_mat: 0 means mat0 is y in linear model; 1 mean mat1 is y in linear model
  
  require(reshape2)
  require(glmnet)
  
  ## test block
  #mat0<-mat0_i
  #mat1<-mat1
  #mat2<-mat2
  #covar<-covar_i
  #direction<-c(1,1)
  #y_mat <- 0
  ## test block
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  
  col_var<-mat0
  row_var<-mat1
  
  my_lr<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,1]
    
    beta    <- NA
    p.value <- NA
    N       <- NA
    uniq_N  <- NA
    
    lr_input<-data.frame(dv = a,
                         indv = b,
                         mat2[,covar])
    lr_input <- sapply(lr_input, as.numeric)
    
    lr_input <- na.omit(lr_input)
    N        <- nrow(lr_input)
    uniq_N   <- length(unique(lr_input[,1]))
    lr_input[,c(2:ncol(lr_input))] <- apply(lr_input[,c(2:ncol(lr_input))], 2, qtrans)
    lr_input<-as.data.frame(lr_input)
    
    
    if(y_mat==0){
      try(lr_res <- summary(glm(dv~., data = lr_input, family = 'binomial')), silent = T)
      indv<-'mbio'
    }else{
      try(lr_res <- summary(glm(dv~., data = lr_input, family = 'binomial')), silent = T)
      indv<-'phen'
    }
    
    #try(beta   <-lr_res$coefficients[2,1],silent = T)
    #try(p.value<-lr_res$coefficients[2,4],silent = T)
    try(beta<-lr_res$coefficients[match("indv",rownames(lr_res$coefficients)),1],silent = T)
    try(p.value<-lr_res$coefficients[match("indv",rownames(lr_res$coefficients)),4],silent = T)

    
    try(return(list(beta = beta, p.value = p.value, N = N, uniq_N = uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_lr(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 4, byrow = T)
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # N matrix
  col_var_row_var.N <- matrix(col_var_row_var.unlist[,3],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.N)<-colnames(col_var)
  rownames(col_var_row_var.N)<-colnames(row_var)
  col_var_row_var.N[is.na(col_var_row_var.N)]<-0
  
  # uniq_N matrix
  col_var_row_var.uN <- matrix(col_var_row_var.unlist[,4],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.uN)<-colnames(col_var)
  rownames(col_var_row_var.uN)<-colnames(row_var)
  col_var_row_var.uN[is.na(col_var_row_var.uN)]<-0
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  col_var_row_var_N          <- melt(col_var_row_var.N)
  col_var_row_var_uN         <- melt(col_var_row_var.uN)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,
                              col_var_row_var_edge_p,
                              col_var_row_var_N,
                              col_var_row_var_uN)[,-c(4,5,7,8,10,11)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "Beta", "p","N","uniq_N")
  
  
  # add p adjust
  col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
                                   fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
                                   bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$fdr.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var_edge$bonferroni.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  
  colnames(col_var_row_var.bon)<-colnames(col_var)
  rownames(col_var_row_var.bon)<-colnames(row_var)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  return(list(table      = col_var_row_var_edge,
              beta       = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon,
              N          = col_var_row_var.N,
              uniq_N     = col_var_row_var.uN))
}

lr_btw_mats_adjAbun2<-function(mat0, mat1, mat2, covar, abun,info,
                               direction=c(1,1),y_mat =0){
  #mat0<-lld_dsv[,c(1:10)]
  #mat1<-lld_exp[,c(1:10)]
  #mat2<-lld_covar
  #covar<-c('Gender','Age','BMI','Reads_number')
  #abun<-lld_abun_clr
  #info<-info
  #direction<-c(1,1)
  #y_mat <- 0
  
  col_var_row_var.edge<-NULL
  for (i in c(1:nrow(info))){
   # i<-18
    cat(paste(i,info$organism[i],"\n"))
    spe_name<-str_replace_all(info$organism[i],"\\[","\\\\\\[") %>% 
      str_replace_all("\\]","\\\\\\]")
    mat0_i<-mat0[,grep(spe_name,colnames(mat0))]
    
    if(dim(mat0_i)[2]>0){
      if(info$Metaphlan3_name[i]%in%colnames(abun)){
        covar_i<-c(covar,info$Metaphlan3_name[i])
        col_var_row_var_i <- lr_btw_mats(mat0_i,mat1 , mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }else{
        covar_i<-covar
        col_var_row_var_i <- lr_btw_mats(mat0_i, mat1, mat2, covar_i, y_mat = 0)
        col_var_row_var.edge<-rbind(col_var_row_var.edge,col_var_row_var_i$table)
      }
    }
  }
  col_var_row_var.edge$fdr.p<-p.adjust(col_var_row_var.edge$p,method = 'fdr')
  col_var_row_var.edge$bonferroni.p<-p.adjust(col_var_row_var.edge$p,method = 'bonferroni')
  
  return(list(table      = col_var_row_var.edge))
}


## kruskal-wallis test
kw_btw_mats<-function(mat0,mat1,direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  require(reshape2)
  require(rstatix)
  
  #mat0<-all_ba[,c(1:10)]
  #mat1<-all_msv_cluster_sub[,c(1:10)]
  #mat2<-covar
  #direction<-c(1,1)
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  col_var<-mat0
  row_var<-mat1
  
  my_kw<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,1]
    
    effSize  <- NA
    p.value  <- NA
    N        <- NA
    N_minor  <- NA
    uniq_N   <- NA
    
    kw_input <-data.frame(phen = a,
                         mbio = b)
    kw_input<-na.omit(kw_input)
    
    res.kw <- kruskal_test(phen ~ mbio, data = kw_input)
    res.kw.eff<-kruskal_effsize(phen ~ mbio, data = kw_input)
    
    effSize   <- res.kw.eff$effsize
    p.value   <- res.kw$p
    N         <- nrow(kw_input)
    N_minor   <- table(kw_input$mbio)[(length(table(kw_input$mbio))+1-rank(table(kw_input$mbio)))==2]
    uniq_N    <- length(unique(kw_input$phen))
    
    try(return(list(effSize = effSize, p.value = p.value,N=N,N_minor=N_minor,uniq_N=uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_kw(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 5, byrow = T)
  
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,col_var_row_var.unlist)[,-c(3)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "effectSize", "p","N","N_minor","uniq_N")
  
  
  # add p adjust
  col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
                                   fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
                                   bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$fdr.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  # bonferroni matrix
  col_var_row_var.bon<-matrix(data  = col_var_row_var_edge$bonferroni.p,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  
  colnames(col_var_row_var.bon)<-colnames(col_var)
  rownames(col_var_row_var.bon)<-colnames(row_var)
  col_var_row_var.bon[is.na(col_var_row_var.bon)]<-1
  
  return(list(table      = col_var_row_var_edge,
              effectSize = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr,
              bonferroni = col_var_row_var.bon))
  
}

## Permutation kruskal-wallis test
permKW_btw_mats<-function(mat0,mat1,direction = c(1,1)){
  # mat0: phenotypes, mat1: microbiome, mat2: covar
  require(reshape2)
  require(rstatix)
  require(coin)
  
  #mat0<-all_ba[,c(1:10)]
  #mat1<-all_msv_cluster_sub[,c(1:10)]
  #direction<-c(1,1)
  
  if(direction[1]==2){mat0 =t(mat0)}
  if(direction[2]==2){mat1 =t(mat1)}
  col_var<-mat0
  row_var<-mat1
  
  my_kw<-function(a,b){
    #a<-col_var[,1]
    #b<-row_var[,3]
    
    effSize  <- NA
    p.value  <- NA
    perm.fdr <- NA
    N        <- NA
    N_minor  <- NA
    uniq_N   <- NA
    
    kw_input <-data.frame(phen = a,
                          mbio = b)
    kw_input<-na.omit(kw_input)
    
    res.kw <- kruskal_test(phen ~ as.factor(mbio), data = kw_input)
    res.kw.eff<-kruskal_effsize(phen ~ as.factor(mbio), data = kw_input)
    res.kw.perm <- oneway_test(phen ~ as.factor(mbio), data=kw_input)
    
    effSize   <- res.kw.eff$effsize
    p.value   <- pvalue(res.kw)
    perm.fdr  <- pvalue(res.kw.perm)
    N         <- nrow(kw_input)
    N_minor   <- table(kw_input$mbio)[(length(table(kw_input$mbio))+1-rank(table(kw_input$mbio)))==2]
    uniq_N    <- length(unique(kw_input$phen))
    
    try(return(list(effSize = effSize,
                    p.value = p.value,
                    perm.fdr=perm.fdr , 
                    N=N,
                    N_minor=N_minor,
                    uniq_N=uniq_N)), silent = T)
  }
  
  
  col_var_row_var<-sapply(
    as.data.frame(col_var),
    function(x) Map(function(a,b) my_kw(a,b),
                    list(x),
                    as.data.frame(row_var)
    )
  )
  col_var_row_var.unlist <- matrix(unlist(col_var_row_var), ncol = 6, byrow = T)
  
  
  # beta matrix
  col_var_row_var.beta<-matrix(col_var_row_var.unlist[,1],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.beta)<-colnames(col_var)
  rownames(col_var_row_var.beta)<-colnames(row_var)
  col_var_row_var.beta[is.na(col_var_row_var.beta)]<-0
  
  # p matrix
  col_var_row_var.p <- matrix(col_var_row_var.unlist[,2],ncol = ncol(col_var), byrow = F)
  colnames(col_var_row_var.p)<-colnames(col_var)
  rownames(col_var_row_var.p)<-colnames(row_var)
  col_var_row_var.p[is.na(col_var_row_var.p)]<-1
  
  # convert matrix to edge list
  col_var_row_var_edge_p     <- melt(col_var_row_var.p)
  col_var_row_var_edge_beta  <- melt(col_var_row_var.beta)
  
  col_var_row_var_edge<-cbind(col_var_row_var_edge_beta,col_var_row_var.unlist)[,-c(3)]
  colnames(col_var_row_var_edge)<-c("Taxa", "Phenotype", "effectSize", "p","perm.fdr","N","N_minor","uniq_N")
  
  
  # add p adjust
  #col_var_row_var_edge<-data.frame(as.data.frame(col_var_row_var_edge),
  #                                 fdr.p = p.adjust(col_var_row_var_edge$p, method = "fdr"),
  #                                 bonferroni.p = p.adjust(col_var_row_var_edge$p, method = "bonferroni"))
  
  # fdr matrix
  col_var_row_var.fdr<-matrix(data  = col_var_row_var_edge$perm.fdr,
                              nrow  = nrow(col_var_row_var.beta),
                              ncol  = ncol(col_var_row_var.beta),
                              byrow = F)
  colnames(col_var_row_var.fdr)<-colnames(col_var)
  rownames(col_var_row_var.fdr)<-colnames(row_var)
  col_var_row_var.fdr[is.na(col_var_row_var.fdr)]<-1
  
  
  return(list(table      = col_var_row_var_edge,
              effectSize = col_var_row_var.beta,
              p          = col_var_row_var.p,
              fdr        = col_var_row_var.fdr))
  
}

## Meta-analysis
my_meta_lm <- function(inVec, study_name, beta_col, se_col) {
  require(meta)
  require(metafor)
  
  #inVec<-cbind_vsv_ba_lm_edge[1,]
  #study_name<-c("LLD", "300OB")
  #beta_col<-c(3,10)
  #se_col<-c(4,11)
  
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
my_batch_meta_lm <- function(inDf,study_name,beta_col,se_col,row_var_col = 1,col_var_col = 2) {
    #inDf<-cbind_vsv_ba_lm_edge[c(1:10),]
    #study_name<-c("LLD", "300OB")
    #beta_col<-c(3,10)
    #se_col<-c(4,11)
    #row_var_col<-1
    #col_var_col<-2
    
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
    
    rowName <- unique(inDf[, row_var_col])
    colName <- unique(inDf[, col_var_col])
    
    N_row <- length(rowName)
    N_col <- length(colName)
    
    batch_res.beta <-
      matrix(batch_res.unlist[, 1], ncol = N_col, byrow = F)
    colnames(batch_res.beta) <- colName
    rownames(batch_res.beta) <- rowName
    batch_res.beta[is.na(batch_res.beta)] <- 0
    
    batch_res.p <- matrix(batch_res.unlist[, 3], ncol = N_col, byrow = F)
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
      byrow = F
    )
    colnames(batch_res.fdr) <- colName
    rownames(batch_res.fdr) <- rowName
    batch_res.fdr[is.na(batch_res.fdr)] <- 1
    
    # bonferroni matrix
    batch_res.bon <- matrix(
      data  = batch_res_edge$Meta.bonferroni.p,
      nrow  = N_row,
      ncol  = N_col,
      byrow = F
    )
    
    colnames(batch_res.bon) <- colName
    rownames(batch_res.bon) <- rowName
    batch_res.bon[is.na(batch_res.bon)] <- 1
    
    return(
      list(
        table      = batch_res_edge,
        beta       = batch_res.beta,
        p          = batch_res.p,
        fdr        = batch_res.fdr,
        bonferroni = batch_res.bon
      )
    )
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



## Mediation analysis for linear model
my_lm_mediation<-function(input.inv, input.dv, input.med, covDf){
  #input.inv<-lld_exp$melatonine
  #input.dv <-lld_ba$CA_dehydro_deconju_ratio
  #input.med<-lld_vsv$`Faecalibacterium cf. prausnitzii KLE1255:1373_1377`
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
    covar<-c(covar, inVec[4])
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
    
    if(info$Abundance_available[n] == "Yes"){
      inCovar<-c(covar, info$Metaphlan3_name[n])
    }else{
      inCovar<-covar
    }
    
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


## Clustering analysis
my_cluster<-function(inDist, ps.cutoff=0.8,my_seed = 1){
  require(NbClust)
  require(fpc)
  require(tsne)
  require(ggsci)
  
  # test
  #inDist<-all_msv_dist_std[[5]]
  # test
  
  clu_n<-prediction.strength(as.dist(inDist), Gmin=2, Gmax=10, M=50,
                             clustermethod=claraCBI ,usepam=T,diss = T,
                             classification="centroid",cutoff=ps.cutoff,
                             distances=T,count=F)
  clu<-claraCBI(as.dist(inDist), clu_n$optimalk, usepam=T,diss=T)
  clu_df<-as.data.frame(clu$partition)
  
  rownames(clu_df)<-rownames(inDist)
  
  pcoa_res<-cmdscale(inDist, k=5, eig = T)
  pcoa <- data.frame(pcoa_res$points)
  
  set.seed(my_seed)
  tsne_res<-Rtsne::Rtsne(inDist, is_distance = TRUE,
                         perplexity = as.integer((nrow(inDist)-1)/3),
                         theta = 0, pca = T,
                         eta=as.integer(nrow(inDist)/12))
  
  tsne_df <- tsne_res$Y %>%
    data.frame() %>%
    setNames(c("X", "Y"))
  tsne_df <- data.frame(Cluster = as.factor(clu$partition),tsne_df)
  
  return(list(pcoa=pcoa, pcoa_res = pcoa_res, clu_n=clu_n,tsne_df=tsne_df))
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


