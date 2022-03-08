####main codes for Lifelines plasma untargeted metabolomics study (GeneralMetabolomics platform) 
###for questions, please contact LianminChen (lianminchen@yeah.net)


############################################################  statistics  #############################################################
###distance matrix‒based variance estimation
#input 'meta'is a matrix with metabolites abundance and 'data' is a matrix with phenotypes (snp,microbial abundance and dietary habits)
library(vegan)
set.seed(666)
bloodmeta=scale(meta)
bloodmeta_distance = as.matrix(vegdist(bloodmeta, method="euclidean"))
results=NULL
for(i in colnames(data)){
  adonis=adonis(bloodmeta_distance_fup~data[,i], permutations = 1000)
  results=rbind(results,adonis$aov.tab[1,])
}

###metabolite associations to genetics, micorbiome and diet
#associations to genetics was done by using QTL mapping pipeline: https://github.com/molgenis/systemsgenetics/blob/master/eqtl-mapping-pipeline
#associations to diet and microbiome was done by Spearman correlation
cor_spearman=function(microbiome/diet,metabolites,name){
  microbiome=microbiome[intersect(row.names(microbiome),row.names(metabolites)),]
  metabolites=metabolites[row.names(microbiome),]
  pvals = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  cors = matrix(NA,ncol = ncol(microbiome),nrow = ncol(metabolites))
  for(i in 1:ncol(microbiome)){
    for(j in 1:ncol(metabolites)){
      a=microbiome[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,i]
      b=metabolites[is.na(microbiome[,i])==F&is.na(metabolites[,j])==F,j]
      if(length(a)>2&length(b)>2){
        cor = cor.test(a,b,method = "spearman")
        pvals[j,i] = cor$p.value
        cors[j,i] = cor$estimate
      }
    }
  }
  qvals = matrix(p.adjust(pvals,method = "BH"),ncol = ncol(pvals))
  colnames(qvals) = colnames(microbiome)
  rownames(qvals) = colnames(metabolites)
  colnames(cors) = colnames(microbiome)
  rownames(cors) = colnames(metabolites)
  ind = which(pvals <=1,arr.ind = T)
  association = data.frame(metabolites = colnames(metabolites)[ind[,1]],microbiome = colnames(microbiome)[ind[,2]],Cor = cors[ind],Pval = pvals[ind],Qval = qvals[ind],stringsAsFactors = F)
  association$name=paste(association$metabolites,association$microbiome,sep = "_with_")
  colnames(association)[1:2]=name
  qvals[is.na(qvals)]=1
  cors=cors[row.names(qvals),colnames(qvals)]
  pvals[pvals<0.05]="."
  pvals[qvals<0.05]="*"
  pvals[pvals>0.05]=NA
  qvals[qvals<0.05]="*"
  qvals[qvals>0.05]=NA
  return(list(association,cors,qvals,pvals))
}

###estimating variance of individual metabolites
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
ElasticNet_variance=function(tmp, data,data_rep){
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
        cv=glmnet::cv.glmnet(as.matrix(tmp[,2:ncol(tmp)]), tmp[,1], alpha=0.5, nfolds=10, type.measure='mse')
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

###Lifelines diet quality score prediction
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

###interaction analysis
#the dataframe 'tmp' contains variables
summary(lm(meta~diet+microbe+diet*microbe,tmp))$coefficients[4,4]

###mediation analysis
library(mediation)
#mediation
colnames(data)=c("X","Y","M")
model.m=lm(M~X,data)
model.y=lm(Y~X+M+X*M,data)
summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
#inverse mediation
colnames(data)=c("X","M","Y")
model.m=lm(M~X,data)
model.y=lm(Y~X+M+X*M,data)
summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))


############################################################  plots  #############################################################
###cirplot
cirplot=function(edge, plot_name){
  colnames(edge)=c("from","to","fre")
  from=data.frame(from=unique(as.character(edge$from)))
  for(i in 1:nrow(from)){
    from$fre[i]=sum(edge$fre[which(edge$from==as.character(from$from[i]))])
  }
  from=from[order(from$fre,decreasing = T),]
  from=as.character(from$from)
  to=data.frame(to=unique(as.character(edge$to)))
  for(i in 1:nrow(to)){
    to$fre[i]=sum(edge$fre[which(edge$to==as.character(to$to[i]))])
  }
  to=to[order(to$fre,decreasing = T),]
  to=as.character(to$to)
  color=colorRampPalette(c('darkblue', 'lightblue','yellow', 'red','pink', 'purple'))(n=length(unique(edge$from)))
  grid.col =  c(color[1:length(from)],rep('grey',length(to)))
  names(grid.col)=c(from,to)
  library(circlize)
  library(wesanderson)
  pdf(plot_name, width = 30, height = 30,useDingbats = F)
  circos.par(start.degree = 0)
  chordDiagram(edge,annotationTrack = "grid",
               grid.col=grid.col,
               order = c(rev(from),rev(to)),
               big.gap = 5,
               preAllocateTracks = list(track.margin = c(0, uh(100, "mm")), 
                                        track.height = max(strwidth(unlist(dimnames(edge)))))
  )
  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
  }, bg.border = NA)
  dev.off()
  circos.clear()
}

###sankey plot
p=ggplot(net,aes(axis1 = net$diet, axis2 = net$microbe_abb, axis3 = net$meta,y= net$fre))+
  scale_x_discrete(limits = c("Diet", "Microbiome", "Meta")) +
  geom_alluvium(aes(fill = net$meta),cex=0.1)+
  geom_stratum() +theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+geom_text(stat = "stratum",cex=0.8,aes(label = after_stat(stratum))) 

#manhattan plot
library(lattice)
manhattan.plot<-function(chr, pos, pvalue, 
                         sig.level=NA, annotate=NULL, ann.default=list(),
                         should.thin=T, thin.pos.places=2, thin.logp.places=2, 
                         xlab="Chromosome", ylab=expression(-log[10](p-value)),
                         col=c("gray","darkgray"), panel.extra=NULL, pch=20, cex=0.8,...) {  
  if (length(chr)==0) stop("chromosome vector is empty")
  if (length(pos)==0) stop("position vector is empty")
  if (length(pvalue)==0) stop("pvalue vector is empty")
  #make sure we have an ordered factor
  if(!is.ordered(chr)) {
    chr <- ordered(chr)
  } else {
    chr <- chr[,drop=T]
  }
  #make sure positions are in kbp
  if (any(pos>1e6)) pos<-pos/1e6;
  #calculate absolute genomic position
  #from relative chromosomal positions
  posmin <- tapply(pos,chr, min);
  posmax <- tapply(pos,chr, max);
  posshift <- head(c(0,cumsum(posmax)),-1);
  names(posshift) <- levels(chr)
  genpos <- pos + posshift[chr];
  getGenPos<-function(cchr, cpos) {
    p<-posshift[as.character(cchr)]+cpos
    return(p)
  }
  #parse annotations
  grp <- NULL
  ann.settings <- list()
  label.default<-list(x="peak",y="peak",adj=3, pos=3, offset=0.5, 
                      col=NULL, fontface=NULL, fontsize=10, show=F)
  parse.label<-function(rawval, groupname) {
    r<-list(text=groupname)
    if(is.logical(rawval)) {
      if(!rawval) {r$show <- F}
    } else if (is.character(rawval) || is.expression(rawval)) {
      if(nchar(rawval)>=1) {
        r$text <- rawval
      }
    } else if (is.list(rawval)) {
      r <- modifyList(r, rawval)
    }
    return(r)
  }
  if(!is.null(annotate)) {
    if (is.list(annotate)) {
      grp <- annotate[[1]]
    } else {
      grp <- annotate
    } 
    if (!is.factor(grp)) {
      grp <- factor(grp)
    }
  } else {
    grp <- factor(rep(1, times=length(pvalue)))
  }
  ann.settings<-vector("list", length(levels(grp)))
  ann.settings[[1]]<-list(pch=pch, col=col, cex=cex, fill=col, label=label.default)
  if (length(ann.settings)>1) { 
    lcols<-trellis.par.get("superpose.symbol")$col 
    lfills<-trellis.par.get("superpose.symbol")$fill
    for(i in 2:length(levels(grp))) {
      ann.settings[[i]]<-list(pch=pch, 
                              col=lcols[(i-2) %% length(lcols) +1 ], 
                              fill=lfills[(i-2) %% length(lfills) +1 ], 
                              cex=cex, label=label.default);
      ann.settings[[i]]$label$show <- T
    }
    names(ann.settings)<-levels(grp)
  }
  for(i in 1:length(ann.settings)) {
    if (i>1) {ann.settings[[i]] <- modifyList(ann.settings[[i]], ann.default)}
    ann.settings[[i]]$label <- modifyList(ann.settings[[i]]$label, 
                                          parse.label(ann.settings[[i]]$label, levels(grp)[i]))
  }
  if(is.list(annotate) && length(annotate)>1) {
    user.cols <- 2:length(annotate)
    ann.cols <- c()
    if(!is.null(names(annotate[-1])) && all(names(annotate[-1])!="")) {
      ann.cols<-match(names(annotate)[-1], names(ann.settings))
    } else {
      ann.cols<-user.cols-1
    }
    for(i in seq_along(user.cols)) {
      if(!is.null(annotate[[user.cols[i]]]$label)) {
        annotate[[user.cols[i]]]$label<-parse.label(annotate[[user.cols[i]]]$label, 
                                                    levels(grp)[ann.cols[i]])
      }
      ann.settings[[ann.cols[i]]]<-modifyList(ann.settings[[ann.cols[i]]], 
                                              annotate[[user.cols[i]]])
    }
  }
  rm(annotate)
  #reduce number of points plotted
  if(should.thin) {
    thinned <- unique(data.frame(
      logp=round(-log10(pvalue),thin.logp.places), 
      pos=round(genpos,thin.pos.places), 
      chr=chr,
      grp=grp)
    )
    logp <- thinned$logp
    genpos <- thinned$pos
    chr <- thinned$chr
    grp <- thinned$grp
    rm(thinned)
  } else {
    logp <- -log10(pvalue)
  }
  rm(pos, pvalue)
  gc()
  #custom axis to print chromosome names
  axis.chr <- function(side,...) {
    if(side=="bottom") {
      panel.axis(side=side, outside=T,
                 at=((posmax+posmin)/2+posshift),
                 labels=levels(chr), 
                 ticks=F, rot=0,
                 check.overlap=F
      )
    } else if (side=="top" || side=="right") {
      panel.axis(side=side, draw.labels=F, ticks=F);
    }
    else {
      axis.default(side=side,...);
    }
  }
  #make sure the y-lim covers the range (plus a bit more to look nice)
  prepanel.chr<-function(x,y,...) { 
    A<-list();
    maxy<-ceiling(max(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)))+1.5;
    miny<-ceiling(min(y, ifelse(!is.na(sig.level), -log10(sig.level), 0)));
    A$ylim=c(miny,maxy);
    A;
  }
  xyplot(logp~genpos, chr=chr, groups=grp,
         axis=axis.chr, ann.settings=ann.settings, 
         #scal y axis
         prepanel=prepanel.chr, scales=list(axs="i"),
         panel=function(x, y, ..., getgenpos) {
           if(!is.na(sig.level)) {
             #add significance line (if requested)
             panel.abline(h=-log10(sig.level), lty=2);
           }
           panel.superpose(x, y, ..., getgenpos=getgenpos);
           if(!is.null(panel.extra)) {
             panel.extra(x,y, getgenpos, ...)
           }
         },
         panel.groups = function(x,y,..., subscripts, group.number) {
           A<-list(...)
           #allow for different annotation settings
           gs <- ann.settings[[group.number]]
           A$col.symbol <- gs$col[(as.numeric(chr[subscripts])-1) %% length(gs$col) + 1]    
           A$cex <- gs$cex[(as.numeric(chr[subscripts])-1) %% length(gs$cex) + 1]
           A$pch <- gs$pch[(as.numeric(chr[subscripts])-1) %% length(gs$pch) + 1]
           A$fill <- gs$fill[(as.numeric(chr[subscripts])-1) %% length(gs$fill) + 1]
           A$x <- x
           A$y <- y
           do.call("panel.xyplot", A)
           #draw labels (if requested)
           if(gs$label$show) {
             gt<-gs$label
             names(gt)[which(names(gt)=="text")]<-"labels"
             gt$show<-NULL
             if(is.character(gt$x) | is.character(gt$y)) {
               peak = which.max(y)
               center = mean(range(x))
               if (is.character(gt$x)) {
                 if(gt$x=="peak") {gt$x<-x[peak]}
                 if(gt$x=="center") {gt$x<-center}
               }
               if (is.character(gt$y)) {
                 if(gt$y=="peak") {gt$y<-y[peak]}
               }
             }
             if(is.list(gt$x)) {
               gt$x<-A$getgenpos(gt$x[[1]],gt$x[[2]])
             }
             do.call("panel.text", gt)
           }
         },
         xlab=xlab, ylab=ylab, 
         panel.extra=panel.extra, getgenpos=getGenPos, ...
  );
}
