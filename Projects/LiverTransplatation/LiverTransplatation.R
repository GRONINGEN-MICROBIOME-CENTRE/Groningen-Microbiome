# ================================== Liver ==================================

library(vegan)
library(ggsci)
library(ggplot2)
library(crayon)
library(foreach)
library(pheatmap)
library(RColorBrewer)
library(randomcoloR)
library(ggsignif)
library(heatmaply)
library(ggpubr)
library(Rtsne)
source(file = "Microbiome.function.R")

# import data
taxa=read.table("./MicrobialData/Metaphlan.merged.txt",sep = "\t",row.names = 1,header = T,stringsAsFactors = F,fill = T)
taxa=as.data.frame(t(taxa))
taxa[is.na(taxa)]=0
taxa=SelectLevel(taxa, level = "species")
PhenoClass=read.table("./Pheno.latest/Phenotype.category.txt",header = T,sep = "\t")

# merge phenotype data, extract only LT and control, rename all smaples
pheno_pre=read.table("./Pheno.latest/Pheno_pre.txt",sep = "\t",header = T,stringsAsFactors = F,fill =T, check.names = F)
pheno_post=read.table("./Pheno.latest/Pheno_post.txt",sep = "\t",header = T,stringsAsFactors = F,fill =T, check.names = F)

pheno_pre=pheno_pre[,colnames(pheno_pre) %in% colnames(pheno_post)]
pheno_post=pheno_post[,colnames(pheno_post) %in% colnames(pheno_pre)]
pheno_pre$`LT-status`="LT_pre"
pheno_post$`LT-status`="LT_post"
pheno=rbind(pheno_pre,pheno_post)

coupling=read.table("Coupling.file.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
pheno=pheno[pheno$sample %in% coupling$sample,]
coupling=coupling[coupling$sample %in% pheno$sample,]
taxa_LT=taxa[rownames(taxa) %in% coupling$MGS,]
taxa_LT=taxa_LT[order(match(rownames(taxa_LT),coupling$MGS)),]
rownames(taxa_LT)=coupling$sample

# remove problmatic samples (QC)
taxa_LT=taxa_LT[rownames(taxa_LT)!="TxL2_A11_177",]
taxa_LT=taxa_LT[rownames(taxa_LT)!="TxL1_F02_14",]
taxa_LT=taxa_LT[rownames(taxa_LT)!="TxL2_E06_141",]
taxa_LT=taxa_LT[rownames(taxa_LT)!="TxL5_A05_510",]
taxa_LT=taxa_LT[rownames(taxa_LT)!="TxL5_F06_523",]
taxa_HC=taxa[grep("LL",rownames(taxa)),]
taxa_combine=rbind(taxa_HC,taxa_LT)
pheno=pheno[pheno$sample %in% rownames(taxa_LT),]

# clean taxa, present rate >10%, relative abundance >0.01%
taxa_combine=Filter(taxa_combine,present = 0.1,relative = 0.0001)
taxa_LT=Filter(taxa_LT,present = 0.1,relative = 0.0001)
taxa_HC=Filter(taxa_HC,present = 0.1,relative = 0.0001)

# alpha diversity
shannon_index=data.frame(row.names = rownames(taxa_combine),shannon=diversity(taxa_combine,"shannon"))
shannon_index=merge(shannon_index,pheno,by.x="row.names",by.y="sample",all=T)
shannon_index$`LT-status`[is.na(shannon_index$`LT-status`)]="Population_DAG3"
shannon_index$new_group[is.na(shannon_index$`new_group`)]="Population_DAG3"
shannon_index$`LT-status`=factor(shannon_index$`LT-status`,levels = c("LT_pre","LT_post","Population_DAG3"))
ggplot(shannon_index, aes(x=`LT-status`, y=shannon, fill=`LT-status`)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+scale_fill_npg()+theme_bw()+
  geom_signif(comparisons = list(c("Population_DAG3", "LT_post"),c("Population_DAG3","LT_pre")), 
              map_signif_level=function(p)sprintf("p = %.2g", p),
              test = "wilcox.test",vjust=-0.2,y_position = c(4.5,4.8))+
  xlab("")
ggsave("./Plot/Status.alpha.pdf")

shannon_index$new_group=factor(shannon_index$new_group,levels = c("Pre-LT","Peri-operative phase","Early post-LT (< 1 yr.)",
                                                                  "1-5 yrs. post-LT","5-10 yrs. post-LT","Long-term post-LT (>10 yrs.)","Population_DAG3"))
ggplot(shannon_index, aes(x=new_group, y=shannon, fill=new_group)) + 
  geom_violin(trim=FALSE)+
  geom_boxplot(width=0.1, fill="white")+scale_fill_npg()+theme_bw()+
  geom_signif(comparisons = list(c("Population_DAG3", "Pre-LT"),c("Population_DAG3","Long-term post-LT (>10 yrs.)"),
                                 c("Pre-LT","Peri-operative phase"),c("Pre-LT","Early post-LT (< 1 yr.)"),
                                 c("Pre-LT","Long-term post-LT (>10 yrs.)"),c("Pre-LT","1-5 yrs. post-LT")), 
              map_signif_level=function(p)sprintf("p = %.2g", p),
              test = "wilcox.test",vjust=-0.2,y_position = c(5.2,4.4,4.4,4.6,5.0,4.8),size=0.1,textsize=2)+
  xlab("")+theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave("./Plot/Group.alpha.pdf")

shannon_index$`years since tx`=as.numeric(shannon_index$`years since tx`)
ggplot(shannon_index[!is.na(shannon_index$`years since tx`),], aes(x = `years since tx`, y = shannon)) + 
  geom_point(aes(color=`years since tx`)) +scale_color_gradient(low="#00A087B2", high="#91D1C2B2")+
  geom_smooth()
ggsave("./Plot/YearSinceLT.alpha.pdf")

# overall beta-diversity
beta_diversity=vegdist(taxa_combine,method = "bray")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates
pcoas=merge(pcoas,pheno,by.x="Sample",by.y="sample",all=T)
pcoas$`LT-status`[is.na(pcoas$`LT-status`)]="Population_DAG3"
pcoas$new_group[is.na(pcoas$new_group)]="Population_DAG3"
ggplot (pcoas, aes(PCoA1,PCoA2,color=`LT-status`)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_lancet()+theme_bw()+scale_fill_lancet()+
  stat_ellipse(aes(fill=`LT-status`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))
ggsave("./Plot/Species.beta.diversity.pdf")

pcoas$new_group=factor(pcoas$new_group,levels = c("Pre-LT","Peri-operative phase","Early post-LT (< 1 yr.)",
                                                  "1-5 yrs. post-LT","5-10 yrs. post-LT","Long-term post-LT (>10 yrs.)","Population_DAG3"))
ggplot (pcoas, aes(PCoA1,PCoA2,color=new_group)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_npg()+theme_bw()+scale_fill_npg()+
  stat_ellipse(aes(fill=new_group),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))
ggsave("./Plot/Species.beta.diversity.newGroup.pdf")

# Tsne plot
set.seed(9)  
tsne_model_1 = Rtsne(as.matrix(beta_diversity), check_duplicates=FALSE, pca=TRUE, perplexity=30, theta=0.5, dims=2)
d_tsne_1 = as.data.frame(tsne_model_1$Y) 
rownames(d_tsne_1)=rownames(taxa_combine)
d_tsne_1=merge(d_tsne_1,pheno,by.x="row.names",by.y="sample",all=T)
d_tsne_1$`LT-status`[is.na(d_tsne_1$`LT-status`)]="Population_DAG3"
d_tsne_1$new_group[is.na(d_tsne_1$new_group)]="Population_DAG3"
ggplot (d_tsne_1, aes(V1,V2,color=`LT-status`)) + 
  geom_point() + theme_bw() + 
  scale_color_lancet()+theme_bw()+scale_fill_lancet()+
  stat_ellipse(aes(fill=`LT-status`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))

# beta diversity on pre_LT
pheno_pre=read.table("./Pheno.latest/Pheno_pre.txt",sep = "\t",header = T,stringsAsFactors = F,fill =T, check.names = F)
pheno_pre=pheno_pre[pheno_pre$sample %in% rownames(taxa_LT),]
pheno_pre=pheno_pre[order(pheno_pre$sample),]
taxa_pre=taxa_LT[rownames(taxa_LT) %in% pheno_pre$sample,]
taxa_pre=taxa_pre[order(rownames(taxa_pre)),]
beta_diversity=vegdist(taxa_pre,method = "bray")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates

pcoas=merge(pcoas,pheno_pre,by.x="Sample",by.y="sample",all=F)
ggplot (pcoas, aes(PCoA1,PCoA2,color=disease)) + 
  geom_point() + theme_bw()+
  stat_ellipse(aes(fill=disease),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))


# beta diversity on post_LT
taxa_LT=taxa_LT[order(rownames(taxa_LT)),]
pheno=pheno[order(pheno$sample),]
beta_diversity=vegdist(taxa_LT[pheno$`LT-status`=="LT_post",],method = "bray")
PCoAList=do_PCoA(beta_diversity)
pcoas=PCoAList$Coordinates
pcoas=merge(pcoas,pheno,by.x="Sample",by.y="sample",all=F)
p1=ggplot (pcoas, aes(PCoA1,PCoA2,color=`re-tx`)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_npg()+theme_bw()+scale_fill_npg()+
  stat_ellipse(aes(fill=`re-tx`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))+
  theme(legend.position = c(0.85, 0.8),legend.background = element_rect(fill=NA))
p2=ggplot (pcoas, aes(PCoA1,PCoA2,color=dum_antibiotics)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_npg()+theme_bw()+scale_fill_npg()+
  stat_ellipse(aes(fill=`dum_antibiotics`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))+
  theme(legend.position = c(0.85, 0.8),legend.background = element_rect(fill=NA))
p3=ggplot (pcoas, aes(PCoA1,PCoA2,color=immunesuppresent)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_npg()+theme_bw()+scale_fill_npg()+
  stat_ellipse(aes(fill=`immunesuppresent`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))+
  theme(legend.position = c(0.85, 0.8),legend.background = element_rect(fill=NA))
p4=ggplot (pcoas, aes(PCoA1,PCoA2,color=ppi)) + 
  geom_point() + theme_bw() + 
  xlab(label = paste("PCoA1",PCoAList$Variance1,sep = " ")) + ylab(label = paste("PCoA2",PCoAList$Variance2,sep = " ")) + 
  scale_color_npg()+theme_bw()+scale_fill_npg()+
  stat_ellipse(aes(fill=`ppi`),level=0.8,geom="polygon",alpha=0.2)+scale_linetype_manual(values=c(5,5,5))+
  theme(legend.position = c(0.85, 0.8),legend.background = element_rect(fill=NA))
ggarrange(p1,p2,p3,p4)
ggsave("./Plot/Post_LT.beta.diversity.pdf",width = 10,height = 8)

# variance explanation, overall LT samples
variable_overall=pheno[,c(which(colnames(pheno)=="LT-status"):ncol(pheno))]
rownames(variable_overall)=pheno$sample
variable_overall=variable_overall[order(rownames(variable_overall)),]
variable_overall=variable_overall[,colSums(!is.na(variable_overall))>nrow(variable_overall)*0.8]
pheno=pheno[order(pheno$sample),]
taxa_LT=taxa_LT[order(rownames(taxa_LT)),]
adonis_overall=matrix(nrow = ncol(variable_overall),ncol = 4)
adonis_overall=as.data.frame(adonis_overall)
colnames(adonis_overall)=c("Factor","F_model","R2","Pvalue")
for(n in 1:ncol(variable_overall)){
  
  variable.sub.sub=variable_overall[,n,drop=F]
  variable.sub.sub=na.omit(variable.sub.sub)
  taxa.sub=taxa_LT[rownames(taxa_LT) %in% rownames(variable.sub.sub),]
  taxa.sub=taxa.sub[order(rownames(taxa.sub)),]
  variable.sub.sub=variable.sub.sub[order(rownames(variable.sub.sub)),,drop=F]
  if(length(unique(variable.sub.sub[,1])) >1 & nrow(variable.sub.sub) >20){
    ad=adonis(taxa.sub ~ variable.sub.sub[,1],permutations = 1000,method = "bray")
    ad_table=ad$aov.tab
    
    adonis_overall$Factor[n]=colnames(variable.sub.sub)[1]
    adonis_overall$F_model[n]=ad_table$F.Model[1]
    adonis_overall$R2[n]=ad_table$R2[1]
    adonis_overall$Pvalue[n]=ad_table$`Pr(>F)`[1]
    
  }else{
    adonis_overall$Factor[n]=colnames(variable.sub.sub)[1]
    adonis_overall$F_model[n]="No power"
    adonis_overall$R2[n]="No power"
    adonis_overall$Pvalue[n]="No power"
  }
  cat(green("Assess variable","-----",colnames(variable_overall)[n],"\n"))
  
}
adonis_overall$Sig=NA
adonis_overall$Sig[adonis_overall$Pvalue<0.05]="Significant"
adonis_overall$Sig[adonis_overall$Pvalue>=0.05]="Non-Significant"
adonis_overall$R2=as.numeric(adonis_overall$R2)
adonis_overall$Pvalue=as.numeric(adonis_overall$Pvalue)
adonis_overall=merge(adonis_overall,PhenoClass,by.x = "Factor",by.y = "Phenotype",all=F)
adonis_overall=adonis_overall[order(adonis_overall$R2,decreasing = F),]
adonis_overall$Factor=factor(adonis_overall$Factor,levels = adonis_overall$Factor)
ggplot(data=adonis_overall, aes(x=Factor, y=R2,fill=Sig)) +
  geom_bar(stat="identity")+
  coord_flip()+theme_bw()+facet_grid(Category ~ .,scales="free", drop = TRUE)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 5),
        axis.text.y = element_text(hjust = 1,vjust = 1,size = 6))+ylab("")+
  scale_fill_lancet()+theme(strip.text.y = element_text(size=5))+
  theme(strip.background =element_rect(fill="#91D1C2B2"))+guides(fill=F)
ggsave("./Plot/Adonis.overall.pdf",width = 5,height = 8)

# variance explanation, differenct timepoints of LT
variable=pheno[,c(which(colnames(pheno)=="re-tx"):ncol(pheno))]
rownames(variable)=pheno$sample
variable=variable[order(rownames(variable)),]
pheno=pheno[order(pheno$sample),]
taxa_LT=taxa_LT[order(rownames(taxa_LT)),]
taxa_LT=taxa_LT[,colnames(taxa_LT) %in% colnames(taxa_combine)]

pheno$new_group=factor(pheno$new_group,levels = c("Pre-LT","Peri-operative phase","Early post-LT (<1y)",
                                                  "1-5y post-LT","5-10y post-LT","Long-term post-LT (>10y)"))
variable=variable[,colSums(!is.na(variable))>nrow(variable)*0.8]

adonis= foreach(i=1:6,.combine = rbind) %do%  {
  
  time=levels(pheno$new_group)[i]
  variable.sub=variable[which(pheno$new_group==time),]
  
  tmp=matrix(ncol = 6, nrow=ncol(variable))
  tmp=as.data.frame(tmp)
  colnames(tmp)=c("Time","Factor","DF","F_model","R2","Pvalue")
  for(n in 1:ncol(variable.sub)){

    variable.sub.sub=variable.sub[!is.na(variable.sub[,n]),]
    taxa.sub=taxa_LT[rownames(taxa_LT) %in% rownames(variable.sub.sub),]
    taxa.sub=taxa.sub[order(rownames(taxa.sub)),]
    variable.sub.sub=variable.sub.sub[order(rownames(variable.sub.sub)),]
    if(length(unique(variable.sub.sub[,n])) >1 & nrow(variable.sub.sub) >20){
    ad=adonis(taxa.sub ~ variable.sub.sub[,n],permutations = 1000,method = "bray")
    ad_table=ad$aov.tab
    
    tmp$Time[n]=time
    tmp$Factor[n]=colnames(variable.sub.sub)[n]
    tmp$DF[n]=ad_table$Df[1]
    tmp$F_model[n]=ad_table$F.Model[1]
    tmp$R2[n]=ad_table$R2[1]
    tmp$Pvalue[n]=ad_table$`Pr(>F)`[1]
    
    }else{
      tmp$Time[n]=time
      tmp$Factor[n]=colnames(variable.sub.sub)[n]
      tmp$DF[n]="No power"
      tmp$F_model[n]="No power"
      tmp$R2[n]="No power"
      tmp$Pvalue[n]="No power"
    }
    cat(yellow(levels(pheno$new_group)[i],"---",colnames(variable)[n],"\n"))
    
  }
  
  return.string=tmp
}

adonis$Time=factor(adonis$Time,levels = c("Pre-LT","Peri-operative phase","Early post-LT (<1y)",
                                          "1-5y post-LT","5-10y post-LT","Long-term post-LT (>10y)"))
adonis=adonis[adonis$Time!="Peri-operative phase",]
adonis$R2=as.numeric(adonis$R2)
adonis$Pvalue=as.numeric(adonis$Pvalue)
adonis$R2[is.na(adonis$R2)]=0
adonis$Pvalue[is.na(adonis$Pvalue)]=1
adonis$Sig=NA
adonis$Sig[adonis$Pvalue<0.05]="Significant"
adonis$Sig[adonis$Pvalue>=0.05]="Non-Significant"
adonis=merge(adonis,PhenoClass,by.x = "Factor",by.y = "Phenotype",all=F)
adonis=adonis[order(adonis$R2,decreasing = F),]
adonis$Factor=factor(adonis$Factor,levels = unique(adonis$Factor))
ggplot(data=adonis, aes(x=Factor, y=R2,fill=Sig)) +
  geom_bar(stat="identity")+
  coord_flip()+theme_bw()+facet_grid(Category ~ Time,scales="free", drop = TRUE)+
  theme(axis.text.x = element_text(angle = 45,hjust = 1,vjust = 1,size = 5),
        axis.text.y = element_text(hjust = 1,vjust = 1,size = 4))+
  scale_fill_lancet()+guides(fill=F)+ylab("")+ylim(0,0.07)+xlab("")+
  theme(strip.text = element_text(size=5))+
  theme(strip.background.y =element_rect(fill="#91D1C2B2"),strip.background.x =element_rect(fill="#8491B4B2"))
ggsave("./Plot/Adonis.timepoints.pdf",width = 8,height = 5)


taxa_LT=taxa_LT[order(rownames(taxa_LT)),]
pheno=pheno[order(pheno$sample),]
pre=(CompositionTable(taxa_LT[pheno$new_group=="Pre-LT",],n=ncol(taxa_LT)))
time1=(CompositionTable(taxa_LT[pheno$new_group=="Peri-operative phase",],n=ncol(taxa_LT)))
time2=(CompositionTable(taxa_LT[pheno$new_group=="Early post-LT (<1y)",],n=ncol(taxa_LT)))
time3=(CompositionTable(taxa_LT[pheno$new_group=="1-5y post-LT",],n=ncol(taxa_LT)))
time4=(CompositionTable(taxa_LT[pheno$new_group=="5-10y post-LT",],n=ncol(taxa_LT)))
time5=(CompositionTable(taxa_LT[pheno$new_group=="Long-term post-LT (>10y)",],n=ncol(taxa_LT)))
control=CompositionTable(taxa_HC,n=ncol(taxa_HC))
pre=AverageTable(pre)
time1=AverageTable(time1)
time2=AverageTable(time2)
time3=AverageTable(time3)
time4=AverageTable(time4)
time5=AverageTable(time5)
control=AverageTable(control)
colnames(pre)[1]="Pre-LT"
colnames(time1)[1]="Peri-operative phase"
colnames(time2)[1]="Early post-LT (<1y)"
colnames(time3)[1]="1-5y post-LT"
colnames(time4)[1]="5-10y post-LT"
colnames(time5)[1]="Long-term post-LT (>10y)"
colnames(control)[1]="Population_DGA3"
top=Reduce(union, list(pre$Taxa[1:20],time1$Taxa[1:20],time3$Taxa[1:20],time4$Taxa[1:20],time5$Taxa[1:20],control$Taxa[1:20])) 
heatmap_data=Reduce(function(x, y) merge(x, y, all=TRUE,by="Taxa"), list(pre[rownames(pre) %in% top,,drop=F],
                                                               time2[rownames(time2) %in% top,,drop=F],
                                                               time3[rownames(time3) %in% top,,drop=F],
                                                               time4[rownames(time4) %in% top,,drop=F],
                                                               time5[rownames(time5) %in% top,,drop=F]))
control=control[rownames(control) %in% top,,drop=F]
heatmap_data=merge(heatmap_data,control,by="Taxa",all=T)
rownames(heatmap_data)=heatmap_data$Taxa
heatmap_data$Taxa=NULL
heatmap_data[is.na(heatmap_data)]=0
pdf("./Plot/Abundant.bacteria.timepoints.pdf",width = 10)
pheatmap(heatmap_data)
dev.off()


pre=as.data.frame(t(taxa_LT[pheno$new_group=="Pre-LT",colnames(taxa_LT) %in% top]))
time1=as.data.frame(t(taxa_LT[pheno$new_group=="Peri-operative phase",colnames(taxa_LT) %in% top]))
time2=as.data.frame(t(taxa_LT[pheno$new_group=="Early post-LT (<1y)",colnames(taxa_LT) %in% top]))
time3=as.data.frame(t(taxa_LT[pheno$new_group=="1-5y post-LT",colnames(taxa_LT) %in% top]))
time4=as.data.frame(t(taxa_LT[pheno$new_group=="5-10y post-LT",colnames(taxa_LT) %in% top]))
time5=as.data.frame(t(taxa_LT[pheno$new_group=="Long-term post-LT (>10y)",colnames(taxa_LT) %in% top]))
control=as.data.frame(t(taxa_HC[,colnames(taxa_HC) %in% top]))
annotation=data.frame(sample=(colnames(control)),group="Population_DGA3")
annotation2=data.frame(sample=pheno$sample,group=pheno$new_group)
annotation=rbind(annotation,annotation2)
rownames(annotation)=annotation$sample
annotation$sample=NULL
pre$taxa=rownames(pre)
time1$taxa=rownames(time1)
time2$taxa=rownames(time2)
time3$taxa=rownames(time3)
time4$taxa=rownames(time4)
time5$taxa=rownames(time5)
control$taxa=rownames(control)
heatmap_data=Reduce(function(x, y) merge(x, y, all=TRUE,by="taxa"), list(pre,time1,time2,time3,time4,time5,control))
heatmap_data[is.na(heatmap_data)]=0
rownames(heatmap_data)=heatmap_data$taxa
heatmap_data$taxa=NULL
pheatmap(heatmap_data, annotation_col = annotation)

# =================================================
# Comparison
# =================================================

invrank= function(x) {qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))}
taxa_LT=taxa_LT[order(rownames(taxa_LT)),]

covariate=read.table("./Pheno.latest/Covariates.maaslin.txt",sep = "\t",header = T,stringsAsFactors = F)
covariate[covariate==""]=NA
rownames(covariate)=covariate$sample

# long term vs. pre
covariate_tmp=covariate[!is.na(covariate$Maaslin.longterm),]
rownames(covariate_tmp)=covariate_tmp$sample
covariate_tmp=covariate_tmp[,c("Maaslin.longterm","age","gender","BMI","ppi","dum_antibiotics","dum_laxative"),drop=F]

covariate_tmp=apply(covariate_tmp,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
covariate_tmp=as.data.frame(covariate_tmp)

taxa_tmp=taxa_LT[rownames(taxa_LT) %in% rownames(covariate_tmp),]
taxa_tmp[taxa_tmp==0]=NA
taxa_tmp=apply(taxa_tmp,2,function(x){invrank(x)})
taxa_tmp = apply(taxa_tmp,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_tmp[!is.na(x),2:7,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100

  return(x)
})

longterm = foreach(i=1:ncol(taxa_tmp),.combine = rbind) %do%  {
    x=as.data.frame(taxa_tmp[,i,drop=F])
    x=na.omit(x)
    y=as.data.frame(covariate_tmp[rownames(covariate_tmp) %in% rownames(x),1,drop=F])
    x=x[order(rownames(x)),,drop=F]
    y=y[order(rownames(y)),,drop=F]
    mm=cor.test(x[,1],y[,1],method = "spearman")
    cat(yellow(colnames(taxa_tmp)[i],"\n"))
  return.string = data.frame(Bac = colnames(taxa_tmp)[i], Cor=mm$estimate,Pvalue=mm$p.value,Group="LongTerm")
}
longterm$FDR=p.adjust(longterm$Pvalue)
longterm=longterm[order(longterm$FDR),]

#  1 year vs. pre
covariate_tmp=covariate[!is.na(covariate$Maaslin.1year),]
rownames(covariate_tmp)=covariate_tmp$sample
covariate_tmp=covariate_tmp[,c("Maaslin.1year","age","gender","BMI","ppi","dum_antibiotics","dum_laxative"),drop=F]

covariate_tmp=apply(covariate_tmp,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
covariate_tmp=as.data.frame(covariate_tmp)

taxa_tmp=taxa_LT[rownames(taxa_LT) %in% rownames(covariate_tmp),]
taxa_tmp[taxa_tmp==0]=NA
taxa_tmp=apply(taxa_tmp,2,function(x){invrank(x)})
taxa_tmp = apply(taxa_tmp,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_tmp[!is.na(x),2:7,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  
  return(x)
})

oneyear = foreach(i=1:ncol(taxa_tmp),.combine = rbind) %do%  {
  x=as.data.frame(taxa_tmp[,i,drop=F])
  x=na.omit(x)
  y=as.data.frame(covariate_tmp[rownames(covariate_tmp) %in% rownames(x),1,drop=F])
  x=x[order(rownames(x)),,drop=F]
  y=y[order(rownames(y)),,drop=F]
  mm=cor.test(x[,1],y[,1],method = "spearman")
  cat(yellow(colnames(taxa_tmp)[i],"\n"))
  return.string = data.frame(Bac = colnames(taxa_tmp)[i], Cor=mm$estimate,Pvalue=mm$p.value,Group="OneYear")
}
oneyear$FDR=p.adjust(oneyear$Pvalue)
oneyear=oneyear[order(oneyear$FDR),]

#  1=5-10 year vs. pre
covariate_tmp=covariate[!is.na(covariate$Maaslin.5.10year),]
rownames(covariate_tmp)=covariate_tmp$sample
covariate_tmp=covariate_tmp[,c("Maaslin.5.10year","age","gender","BMI","ppi","dum_antibiotics","dum_laxative"),drop=F]

covariate_tmp=apply(covariate_tmp,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
covariate_tmp=as.data.frame(covariate_tmp)

taxa_tmp=taxa_LT[rownames(taxa_LT) %in% rownames(covariate_tmp),]
taxa_tmp[taxa_tmp==0]=NA
taxa_tmp=apply(taxa_tmp,2,function(x){invrank(x)})
taxa_tmp = apply(taxa_tmp,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_tmp[!is.na(x),2:7,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  
  return(x)
})

fivetotenyear = foreach(i=1:ncol(taxa_tmp),.combine = rbind) %do%  {
  x=as.data.frame(taxa_tmp[,i,drop=F])
  x=na.omit(x)
  y=as.data.frame(covariate_tmp[rownames(covariate_tmp) %in% rownames(x),1,drop=F])
  x=x[order(rownames(x)),,drop=F]
  y=y[order(rownames(y)),,drop=F]
  mm=cor.test(x[,1],y[,1],method = "spearman")
  cat(yellow(colnames(taxa_tmp)[i],"\n"))
  return.string = data.frame(Bac = colnames(taxa_tmp)[i], Cor=mm$estimate,Pvalue=mm$p.value,Group="FiveToTenYear")
}
fivetotenyear$FDR=p.adjust(fivetotenyear$Pvalue)
fivetotenyear=fivetotenyear[order(fivetotenyear$FDR),]

#  1-5 year vs. pre
covariate_tmp=covariate[!is.na(covariate$Maaslin.1.5year),]
rownames(covariate_tmp)=covariate_tmp$sample
covariate_tmp=covariate_tmp[,c("Maaslin.1.5year","age","gender","BMI","ppi","dum_antibiotics","dum_laxative"),drop=F]

covariate_tmp=apply(covariate_tmp,2,function(x){
  
  x[is.na(x)]=median(x[!is.na(x)])
  return(x)
  
})
covariate_tmp=as.data.frame(covariate_tmp)

taxa_tmp=taxa_LT[rownames(taxa_LT) %in% rownames(covariate_tmp),]
taxa_tmp[taxa_tmp==0]=NA
taxa_tmp=apply(taxa_tmp,2,function(x){invrank(x)})
taxa_tmp = apply(taxa_tmp,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_tmp[!is.na(x),2:7,drop = FALSE]
  
  x.resid = resid(lm(x.subset ~ .,data = covariate.subset))
  x[!is.na(x)] = x.resid+100
  
  return(x)
})

onetofiveyear = foreach(i=1:ncol(taxa_tmp),.combine = rbind) %do%  {
  x=as.data.frame(taxa_tmp[,i,drop=F])
  x=na.omit(x)
  y=as.data.frame(covariate_tmp[rownames(covariate_tmp) %in% rownames(x),1,drop=F])
  x=x[order(rownames(x)),,drop=F]
  y=y[order(rownames(y)),,drop=F]
  mm=cor.test(x[,1],y[,1],method = "spearman")
  cat(yellow(colnames(taxa_tmp)[i],"\n"))
  return.string = data.frame(Bac = colnames(taxa_tmp)[i], Cor=mm$estimate,Pvalue=mm$p.value,Group="OneToFiveYear")
}
onetofiveyear$FDR=p.adjust(onetofiveyear$Pvalue)
onetofiveyear=onetofiveyear[order(onetofiveyear$FDR),]

# plot most changeble bacteria
Bacteria_change=Reduce(union,list(oneyear$Bac[oneyear$FDR<2],onetofiveyear$Bac[onetofiveyear$FDR<2],fivetotenyear$Bac[fivetotenyear$FDR<2],longterm$Bac[longterm$FDR<2]))
palette <- distinctColorPalette(length(Bacteria_change))
Bacteria_change=taxa_LT[,colnames(taxa_LT) %in% Bacteria_change]

Bacteria_change=Bacteria_change[order(rownames(Bacteria_change)),]
pheno=pheno[order(pheno$sample),]
pre=(CompositionTable(Bacteria_change[pheno$new_group=="Pre-LT",],n=ncol(Bacteria_change)))
time1=(CompositionTable(Bacteria_change[pheno$new_group=="Peri-operative phase",],n=ncol(Bacteria_change)))
time2=(CompositionTable(Bacteria_change[pheno$new_group=="Early post-LT (< 1 yr.)",],n=ncol(Bacteria_change)))
time3=(CompositionTable(Bacteria_change[pheno$new_group=="1-5 yrs. post-LT",],n=ncol(Bacteria_change)))
time4=(CompositionTable(Bacteria_change[pheno$new_group=="5-10 yrs. post-LT",],n=ncol(Bacteria_change)))
time5=(CompositionTable(Bacteria_change[pheno$new_group=="Long-term post-LT (>10 yrs.)",],n=ncol(Bacteria_change)))

pre$group="Pre-LT"
time1$group="Peri-operative phase"
time2$group="Early post-LT (<1y)"
time3$group="1-5y post-LT"
time4$group="5-10y post-LT"
time5$group="Long-term post-LT (>10y)"

change_plot=rbind(pre,time1,time2,time3,time4,time5)
change_plot=merge(change_plot,pheno[,c("sample","new_group","years since tx"),drop=F],by.x = "ID",by.y = "sample",all=T)
change_plot$`years since tx`=as.numeric(change_plot$`years since tx`)
change_plot$new_group=factor(change_plot$new_group,levels = c("Pre-LT","Peri-operative phase","Early post-LT (< 1 yr.)",
                                                                  "1-5 yrs. post-LT","5-10 yrs. post-LT","Long-term post-LT (>10 yrs.)"))
change_plot$group[change_plot$new_group=="Pre-LT"]=-1
change_plot$group[change_plot$new_group=="Peri-operative phase"]=1
change_plot$group[change_plot$new_group=="Early post-LT (< 1 yr.)"]=2
change_plot$group[change_plot$new_group=="1-5 yrs. post-LT"]=3
change_plot$group[change_plot$new_group=="5-10 yrs. post-LT"]=4
change_plot$group[change_plot$new_group=="Long-term post-LT (>10 yrs.)"]=5
change_plot$group=as.numeric(change_plot$group)
ggplot(change_plot, aes(x = group, y = Relative)) + 
  geom_point(aes(color=Level)) +
  geom_smooth(aes(color=Level))+
  facet_wrap(Level ~ .,ncol=3,scales = "free")+guides(color=F)+scale_color_manual(values = palette)
ggsave("./Plot/Bacteria_change.pdf",width = 10,height = 100,limitsize = FALSE)
