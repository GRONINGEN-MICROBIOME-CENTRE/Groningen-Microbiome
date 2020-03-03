###################################
### merge GWAS and plot
### version 1.0
### date: 01-29-2020
### Author:EALM
###################################
##New
## 06-02-2020
## changed axes in qqplot (corrected the X-Y swap)
## corrected code to make lambdas.list for all microbes (changed from 1 to n, and 23 for mb  in line 127)
## added shortnames to lambdas plot
## added chucnk to plot all lambdas by 20
## added chunk to plot together the gwas plot by 10
## corrected n reaassignment in the qqplot, that caused the loop to fail
## added if loop to plot only manhattan for gwas <0.00000001
## 26-02-2020
## added argument parser to work together with a sbatch launcher
## changed chromosome reading as numeric for character, to take also chromosome X


####Packages
library(data.table)
library(tidyverse)
library(gridExtra)
library(ggpubr)
library(optparse)

#####paths & variables
##gearshift (does not plot yet)
#mblist<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/heritable.list"
#input<-"/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/outputs/"

##calculon
#input<-"/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output"
#opt$input<-"/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2"
#mblist<-"/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/heritable.list"
#opt$mblist<-"/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2/list.to.graph"

######################## --------arguments -----------##########
option_list = list(
  
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="working directory directory where the folders of each feature are located", metavar="character"),
  
  make_option(c("-m", "--mblist"), type="character", default=NULL,
              help="list file containing the names of the features evaluated, with NO HEADER", metavar="character")
  
)

opt_parser  <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

input<-opt$input
mblist<-opt$mblist
######-----------main-----------------------#############
mblist<-fread(mblist,data.table = F,header=F)
lambdas.list<-data.frame(name=NA,lambda=NA,short.name=NA)
img.list<-list()
#n=6
for (n in (1:nrow(mblist))) {
  
  
  mb<-mblist[n,"V1"]
  chrfile.path <- file.path(input,paste0(mb,"_gwas"))
  chr.files <- list.files(chrfile.path, pattern = ".chr", full.names = TRUE)

  ## make single file with all chromosoemes
  chr.dat.list <- lapply(chr.files, fread, data.table=FALSE)
  chr.dat.list <- lapply(chr.dat.list, function(x){x[,c("id","chr", "pos","AF.alt", "AC.alt", "num",  "beta","SE", "pval")]})
  #chr.dat.list <- lapply(chr.dat.list, function(x){head(x[,c("id","chr", "pos","AF.alt", "AC.alt", "num",  "beta","SE", "pval")],1000)})#test
  chr.dat.list <- lapply(chr.dat.list, function(df){mutate_at(df, .vars = c("chr"), as.character)})
  chr.dat.list <- lapply(chr.dat.list, function(df){mutate_at(df, .vars = c("pval"), as.numeric)})
  whole.genome <- bind_rows(chr.dat.list)
  rm(chr.dat.list)
  nrow(whole.genome)
  ###-----------Manhatan plot---------------
  whole.genome<-whole.genome[order(whole.genome$chr),]
  
  nCHR <- length(unique(whole.genome$chr))
  whole.genome$BPcum <- NA
  s <- 0
  nbp <- c()

  for (i in unique(sort(as.factor(whole.genome$chr)))){
    nbp[i] <- max(whole.genome[whole.genome$chr == i,]$pos)
    whole.genome[whole.genome$chr == i,"BPcum"] <- whole.genome[whole.genome$chr == i,"pos"] + s
    s <- s + nbp[i]
  }
  
  axis.set <- whole.genome %>% 
    group_by(chr) %>% 
    summarize(center = (max(BPcum) + min(BPcum)) / 2)
  ylim <- abs(floor(log10(min(whole.genome$pval)))) + 10
  sig <- 5e-8
  
  #### make short names for plot 
  mb.names<-strsplit(mb,"\\.")
  names<-lapply(mb.names, function(x){paste0(x[length(x)-1],".",x[length(x)])})
  short.name<-as.character(unlist(names))
  ### make manhatan plot if there is any GWAS significant signal
  if (any(whole.genome$pval<0.00000001)){
    manhplot <- ggplot(whole.genome, aes(x = BPcum, y = -log10(pval), 
                                         color = as.factor(chr), size = -log10(pval))) +
      geom_point(alpha = 0.75) +
      geom_hline(yintercept = -log10(sig), color = "red", linetype = "dashed") + 
      scale_x_continuous(label = axis.set$chr, breaks = axis.set$center) +
      scale_y_continuous(expand = c(0,0), limits = c(0, ylim)) +
      scale_color_manual(values = rep(c("#276FBF", "#183059"), nCHR)) +
      scale_size_continuous(range = c(0.5,3)) +
      labs(x = NULL, 
           y = "-log10(p)") + 
      ggtitle(short.name)+
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90, size = 8, vjust = 0.5)
      )  
    
    mhn.plot.file<-file.path(chrfile.path,paste0("mhn.plot.tiff"))
    tiff(mhn.plot.file,  
         width = 1500, height = 1500, 
         units = "px", res = 300, compression = "lzw")
    plot(manhplot)
    dev.off()
    img.list[[n]]<- manhplot #save the plot in a list to plot together later
    
  }
  
      
  ## ----QQplot------------------------------------------------------
  g <- nrow(whole.genome)
  # expected p-values
  exp.pval <- (1:g) / g
  # observed p-values
  obs.pval <- sort(whole.genome$pval)
  qq.plot.file<-file.path(chrfile.path,paste0("qq.plot.tiff"))
  tiff(qq.plot.file,  
       width = 1500, height = 1500, 
       units = "px", res = 300, compression = "lzw")
  
  plot(-log10(exp.pval), -log10(obs.pval), xlab="-log10(expected P)", ylab="-log10(observed P)", cex=2/3)+
    abline(0, 1, col="red", lty=2)
  
  dev.off()
  
  ###----calculate lambda-------
  whole.genome$chisq<-qchisq(1-whole.genome$pval,1)
  lambda.table<-whole.genome %>% 
    group_by(chr) %>% 
    summarise(lambda=median(chisq)/qchisq(0.5,1))
  summary.row<-c(mb,median( whole.genome$chisq)/qchisq(0.5,1))
  lambda.table<-data.frame(rbind(lambda.table,summary.row))
  lambda.file<-file.path(chrfile.path,paste0("lambda.table"))
  write.table(lambda.table,lambda.file,quote = F,row.names = F,sep = '\t')
  ##lambda is 1.03568. to my point of view there is not a lot of inflation, but... would this change when using it  with all the chromosomes?
  lambdas.list[n,]<-lambda.table[lambda.table$chr==mb,]### changed from 1 to n and 23 for lambda.table$chr==mb
  lambdas.list[n,"short.name"]<-short.name
  ####save merged gwas results
  gwas.file<-file.path(chrfile.path,paste0("whole.genome.gwas"))
  write.table(whole.genome,gwas.file,quote = F,row.names = F,sep = '\t')
  
}

###----------make the lambda table for all microbes in list----------#####
lambdas.list.file<-file.path(input,paste0("lambdas.list"))
write.table(lambdas.list,lambdas.list.file,quote = F,row.names = F,sep = '\t')
#lambdas.list<-fread(lambdas.list.file,data.table=F,header=T)

  ### plot lambdas in barplot
  lambdas.list$lambda<-as.numeric(lambdas.list$lambda)
  lambda.bar<-ggplot(lambdas.list,aes(x=lambdas.list$short.name,y=lambda))+
    geom_bar(stat="identity",width=0.98)+
    geom_hline(yintercept = 1, color = "red", linetype = "dashed")+
    geom_hline(yintercept = 1.1, color = "blue", linetype = "dotted")+
    geom_hline(yintercept = 0.9, color = "blue", linetype = "dotted")+
    scale_y_continuous(limits=c(0, 1.5),breaks = c(seq(0,1.5,0.1)))+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90, size = 8))+
    xlab("Microbial feature")+
    ylab("Lambda")
  
  lambda.plot.file<-file.path(input,"lambda.bar.plot.tiff")
  tiff(lambda.plot.file,  
       width = 6000, height = 2000, 
       units = "px", res = 300, compression = "lzw")
  plot(lambda.bar)
  dev.off()

  ### plot in frequency histogram
  lambda.his<-ggplot(lambdas.list,aes(x=lambda))+
    geom_histogram(binwidth = 0.005)+
    geom_vline(xintercept = 1, color = "red", linetype = "dashed")+
    geom_vline(xintercept = 1.1, color = "blue", linetype = "dotted")+
    geom_vline(xintercept = 0.9, color = "blue", linetype = "dotted")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, size = 8))+
    xlab("Lambda")+
    ylab("Frequency")
  lambda.plot.file<-file.path(input,"lambda.his.plot.tiff")
  tiff(lambda.plot.file,  
       width = 6000, height = 2000, 
       units = "px", res = 300, compression = "lzw")
  plot(lambda.his)
  dev.off()
  
###plot every 10 gwas
m=0
glimit=0
while(glimit<length(img.list)){
  glimit<-ifelse(10*(m+1)>length(img.list),length(img.list),10*(m+1))
  plist<-img.list[(10*(m)+1):glimit]
  ### join gwas
  image.file.plt.file<-file.path(input,paste0((m+1),".mhn.all.tiff"))
  tiff(image.file.plt.file, width = 6000, height = 3000, units = "px", res = 300, compression = "lzw")
  do.call("grid.arrange", c(plist, ncol=floor(sqrt(10))))
  dev.off()  
  m=m+1
}







