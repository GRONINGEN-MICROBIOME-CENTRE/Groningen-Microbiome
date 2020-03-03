# ============================================================================================
#       This script is to anlayze outfiles from StrainPhlan
# ============================================================================================
require(ggtree)
require(tidyverse)
library(optparse)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(randomcoloR)
library(crayon)
library(cowplot)
library(ggsignif)
library(foreach)

InputFolder="/Users/hushixian/Desktop/PhD/Liver project/Analysis/StrainPhlanInput/"
filenames <- list.files(InputFolder,  full.names=F)
metadata = read.delim ("Cas_TxL_annotation.csv",sep = ",",stringsAsFactors = T)
Color<-(distinctColorPalette(length(unique(metadata$Subject))))
names(Color)=unique(metadata$Subject)
metadata=merge(Color,metadata,by="Subject")
metadata=metadata[,c(3,1,2,4,5)]

coupling=read.table("Coupling.file.txt",header = T)
coupling$MGS=gsub("\\.","-",coupling$MGS)
coupling=coupling[coupling$sample %in% metadata$SampleID,]
coupling=coupling[order(coupling$sample),]
metadata=metadata[order(metadata$SampleID),]
metadata$SampleID=coupling$MGS
metadata$Subject=as.factor(metadata$Subject)

plotNames=list()
count = 1
for(i in filenames){
  
  taxa=i
  
  if (file.exists(paste0(InputFolder,'/',taxa,'/',"RAxML_bestTree.",taxa,".tree",sep = ""))) {
    
    BestTree=read.tree(paste0(InputFolder,'/',taxa,'/',"RAxML_bestTree.",taxa,".tree",sep = ""))
    meta=metadata
    meta=meta[meta$SampleID %in% BestTree$tip.label,]
    
    culledtree <- BestTree
    TreePLot=ggtree(culledtree)
    TreePLot <- TreePLot %<+% meta
    
    plotNames[[count]]=TreePLot+geom_tippoint( size = 4, aes( color = Subject,shape=TimePrePost ),na.rm = T ) + 
      aes( branch.length = 'length' ) +
      theme_tree2() + theme(legend.position="right")+scale_color_manual(values = Color)+ggtitle(taxa)
    
    cat(yellow(taxa,"plot is done","\n"))
    count = count + 1
    
  }
  
}
pdf("StrainPhlan.TreePlot.pdf",width = 30,height = 30)
cowplot::plot_grid(plotlist=plotNames,ncol = 5,scale = 0.7)
dev.off()






plotNames=list()
count = 1

InputFolder="/Users/hushixian/Desktop/PhD/Liver project/Analysis/StrainPhlanInput/"
filenames <- list.files(InputFolder,  full.names=F)

for(taxa in filenames){
  
  if (file.exists(paste0(InputFolder,'/',taxa,'/',"RAxML_bestTree.",taxa,".tree",sep = ""))) {
    
    Dmat=read.csv(paste0(InputFolder,'/',taxa,'/',taxa,"_dmat_Rready.csv"))
    metadata = read.delim ("Cas_TxL_annotation.csv",sep = ",",stringsAsFactors = T)
    colnames(Dmat)=lapply(colnames(Dmat),function(x){
      strsplit(x,"\\.")[[1]][1]})
    Dmat$name=rownames(Dmat)
    Dmat$name=lapply(Dmat$name,function(x){
      strsplit(x,"\\.")[[1]][1]})
    Dmat=Dmat[!duplicated(Dmat$name),]
    rownames(Dmat)=lapply(rownames(Dmat),function(x){
      strsplit(x,"\\.")[[1]][1]})
    meta=metadata[metadata$SampleID %in% colnames(Dmat),]
    Dmat$name=NULL
    Dmat=Dmat[,!duplicated(colnames(Dmat)),drop=F]
    meta$Subject=as.factor(meta$Subject)
    
    evolution= foreach(i=1:length(unique(meta$Subject)),.combine = rbind) %do%  {
      
      subject=as.character(unique(meta$Subject)[i])
      if(length(unique(meta[meta$Subject==subject,]$TimePrePost))>1){
        tmp=meta[meta$Subject==subject,]
        if(length(tmp$TimePrePost[tmp$TimePrePost=="Post"])>1){
          pre=as.character(tmp$SampleID[which(tmp$TimePrePost=="Pre")])
          post=as.character(tmp$SampleID[which(tmp$TimePrePost=="Post")])
          pre_post=c()
          for(n in 1:length(post)){
            pre_post=append(pre_post,Dmat[rownames(Dmat)==pre,colnames(Dmat)==post[n]])
          }
          
          pre_post=data.frame(Dissimilarity=pre_post,Group="Pre.vs.Post" )
          post_post=Dmat[rownames(Dmat) %in% post,colnames(Dmat) %in% post,drop=F]
          post_post=(post_post[lower.tri(post_post)])
          post_post=data.frame(Dissimilarity=post_post,Group="Post.vs.Post")
          
          compare=rbind(pre_post,post_post)
          compare$Subject=subject
          
        }else{compare=NULL}
      }else{compare=NULL}
      return.string=compare
    }
    
    if(!is.null(nrow(evolution))){
      plotNames[[count]]=ggplot(evolution, aes(x=Group, y=Dissimilarity, fill=Group)) + 
        geom_boxplot(width=0.1)+scale_fill_npg()+theme_bw()+
        geom_jitter()+
        geom_signif(comparisons = list(c("Pre.vs.Post","Post.vs.Post")), 
                    map_signif_level=function(p)sprintf("p = %.2g", p),
                    test = "wilcox.test",vjust=-0.2)+xlab("")
      count=count+1
      cat(yellow(taxa,"plot is done","\n"))
    }

  }
  
}

cowplot::plot_grid(plotlist=plotNames,ncol = 2)
