###load data
annotation=read.delim("./key_lld_1183meta_annotation.txt",header = T,quote = "",sep = "\t",check.names = F)
annotation$super_class=gsub("\"","",annotation$super_class)
meta=read.delim("./data_1442samples_LLD_baseline_1183plasma_metabolites.txt",header = T,sep = "\t")
meta=meta[order(row.names(meta)),]
meta_fup=read.delim("./data_342samples_LLD_fup_1183plasma_metabolites.txt",header = T,sep = "\t")
meta_fup=meta_fup[intersect(row.names(meta),row.names(meta_fup)),]
row.names(meta_fup)=paste(row.names(meta_fup),"_F",sep = "")
meta_fup=meta_fup[order(row.names(meta_fup)),]
meta_rep=read.delim("./data_314samples_239LLD_80GoNL_1183plasma_metabolites.txt",header = T,row.names = 1,sep = "\t")
meta_rep=meta_rep[order(row.names(meta_rep)),]
pheno=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_207pheno_log_imp_1135samples.txt",header = T,sep = "\t")
pheno_fup=read.delim("~/Documents/UMCG/manuscript/LLD_fup_paper/1_LLD_5years_data/phenotype/data_pheno_LLD_base_fup_338pairs_62pheno_18med_17disease_log_imput_min_10participants.txt",header = T,row.names = 1,sep = "\t")
pheno_rep=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_1135samples_345fup_mgs_key_raw_reads.txt",header = T,row.names = 1,sep = "\t")
pheno_rep=pheno_rep[which(row.names(pheno_rep)%in%row.names(pheno_fup)),]
row.names(pheno_rep)=paste(row.names(pheno_rep),"_F",sep = "")
pheno_fup=merge(pheno_rep[,5:6],pheno_fup,by="row.names",all = F)
row.names(pheno_fup)=pheno_fup[,1]
colnames(pheno_fup)[3]="totalreads"
pheno_fup=pheno_fup[,-c(1,2,4)]
pheno_rep=read.delim("./LLD_1440samples_age_sex_bmi_smk.txt",header = T,row.names = 1,sep = "\t")
pheno_rep=pheno_rep[row.names(meta_rep),]
library(microbiome)
sp=read.delim("~/Documents/work.UMCG/research-projects/microbiome_data/merged_microbime_data/LLD_1135base_338fup_558species.txt",header = T,sep = "\t",check.names = F)
sp=sp/rowSums(sp)
sp=sp[,colSums(sp)>0]#467
clr=data.frame(t(abundances(t(sp),transform = "clr")),check.names = F)
sp=sp[rowSums(sp)>0,colSums(sp>0)>nrow(sp)*0.1]#156
sp=clr[,colnames(sp)]
path=read.delim("~/Documents/UMCG/research-projects/microbiome_data/merged_microbime_data/LLD_1135base_338fup_467path_abundance.txt",quote = "",header = T,sep = "\t",check.names = F)
path=path[,colSums(path)>0]#467
clr=data.frame(t(abundances(t(path),transform = "clr")),check.names = F)
path=path[rowSums(path)>0,colSums(path>0)>nrow(path)*0.1]#343
path=clr[,colnames(path)]
dsv=read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_SVs/LLD_base_fup_SVs/01.cleanData/20200615_LLD12_deletionStructuralVariation_1471samples.tsv",header = T,sep = "\t",check.names = F)
dsv=dsv[,colSums(is.na(dsv)==F)>nrow(dsv)*0.1]
dsv=merge(sp[,1:2],dsv,by="row.names",all = T)
row.names(dsv)=dsv[,1]
dsv=dsv[,-c(1:3)]
vsv=read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_SVs/LLD_base_fup_SVs/01.cleanData/20200615_LLD12_variableStructuralVariation_1471samples.tsv",header = T,sep = "\t",check.names = F)
vsv=vsv[,colSums(is.na(vsv)==F)>nrow(vsv)*0.1]
vsv=merge(sp[,1:2],vsv,by="row.names",all = T)
row.names(vsv)=vsv[,1]
vsv=vsv[,-c(1:3)]
id_base=intersect(row.names(meta),row.names(pheno))[grep("LLDeep_0186",intersect(row.names(meta),row.names(pheno)),invert = T)]
id_fup=intersect(intersect(row.names(meta_fup),paste(id_base,"_F",sep = "")),row.names(sp))
id_rep=row.names(meta_rep)[1:237]
id_gonl=row.names(meta_rep)[238:314]
meta=meta[id_base,]
meta_fup=meta_fup[id_fup,]
sp_fup=sp[id_fup,]
sp=sp[id_base,]
path_fup=path[id_fup,]
path=path[id_base,]
dsv_fup=dsv[intersect(row.names(dsv),id_fup),]
dsv=dsv[intersect(row.names(dsv),id_base),]
vsv_fup=vsv[intersect(row.names(vsv),id_fup),]
vsv=vsv[intersect(row.names(vsv),id_base),]
pheno=pheno[row.names(meta),]
pheno_fup=pheno_fup[id_fup,]
bgc=data.frame(t(read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_BGC/LLD_base_BGC/LLD_baseline_1135samples_all_RPKMs_norm.tsv",header = T,sep = "\t",row.names = 1,check.names = F)),check.names = F)
name=data.frame(t(read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_BGC/LLD_base_BGC/LLD_CGR+HMP+Clos_RPKM.txt",header = T,sep = "\t",row.names = 1,check.names = F)),check.names = F)
colnames(bgc)=colnames(name)
name=NA
bgc=bgc[,colSums(bgc>0)>nrow(bgc)*0.1]
bgc=bgc[row.names(pheno),]
bgc_fup=read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_BGC/LLD_base_BGC/LLD_fup_338samples_all_RPKMs_norm.tsv",header = T,sep = "\t",row.names = 1,check.names = F)
row.names(bgc_fup)=paste(row.names(bgc_fup),"_F",sep = "")
bgc_fup=bgc_fup[row.names(sp_fup),colnames(bgc)]
diet_score=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_1447samples_diet_score.txt",header = T,sep = "\t")
diet_score=na.omit(diet_score)
diet_score_rep=diet_score[intersect(row.names(diet_score),row.names(meta_rep)),]
diet_score=diet_score[intersect(row.names(diet_score),row.names(meta)),]

############################################analysis
#diet summary
diet=pheno[,48:125]
result=data.frame(matrix(NA,nrow = ncol(diet),ncol = 7))
colnames(result)=c("diet","mean_female","sd_female","mean_male","sd_male","mean_all","sd_all")
for(i in 1:ncol(diet)){
  result$diet[i]=as.character(colnames(diet)[i])
  if(i %in% c(26:30, 39:42)){
    result$mean_female[i]=length(which(diet[which(pheno$antrop_gender.F1M2==1),i]==1))
    result$sd_female[i]=length(which(diet[which(pheno$antrop_gender.F1M2==1),i]==0))
    result$mean_male[i]=length(which(diet[which(pheno$antrop_gender.F1M2==2),i]==1))
    result$sd_male[i]=length(which(diet[which(pheno$antrop_gender.F1M2==2),i]==0))
    result$mean_all[i]=length(which(diet[which(pheno$antrop_gender.F1M2!=0),i]==1))
    result$sd_all[i]=length(which(diet[which(pheno$antrop_gender.F1M2!=0),i]==0))
  }else{
    result$mean_female[i]=mean(diet[which(pheno$antrop_gender.F1M2==1),i],na.rm = T)
    result$sd_female[i]=sd(diet[which(pheno$antrop_gender.F1M2==1),i],na.rm = T)
    result$mean_male[i]=mean(diet[which(pheno$antrop_gender.F1M2==2),i],na.rm = T)
    result$sd_male[i]=sd(diet[which(pheno$antrop_gender.F1M2==2),i],na.rm = T)
    result$mean_all[i]=mean(diet[,i],na.rm = T)
    result$sd_all[i]=sd(diet[,i],na.rm = T)
  }
}
write.table(result,file = "../result/summary_78diet_1054samples.txt",quote = F,sep = "\t",row.names = F)
#stability
result=data.frame(matrix(NA,nrow = ncol(meta),ncol = 8))
colnames(result)=c("meta","mean_base","sd_base","mean_fup","sd_fup","P_wilcoxon","r_spearman","P_spearman")
for(i in 1:ncol(meta)){
  result$meta[i]=colnames(meta)[i]
  result$mean_base[i]=mean(meta[gsub("_F","",row.names(meta_fup)),i])
  result$sd_base[i]=sd(meta[gsub("_F","",row.names(meta_fup)),i])
  result$mean_fup[i]=mean(meta_fup[,i])
  result$sd_fup[i]=sd(meta_fup[,i])
  result$P_wilcoxon[i]=wilcox.test(meta[gsub("_F","",row.names(meta_fup)),i],meta_fup[,i],paired = T)$p.value
  result$r_spearman[i]=cor.test(meta[gsub("_F","",row.names(meta_fup)),i],meta_fup[,i],method = "spearman")$estimate
  result$P_spearman[i]=cor.test(meta[gsub("_F","",row.names(meta_fup)),i],meta_fup[,i],method = "spearman")$p.value
}
result=merge(result,annotation[,c(1,11,12,13)],by="meta")
position=data.frame(cbind(unique(as.character(result$super_class)),unique(as.character(result$super_class))))
for(i in 1:nrow(position)){
  position$mean[i]=mean(result$r_spearman[which(result$super_class==as.character(position$X1[i]))])
  position$median[i]=median(result$r_spearman[which(result$super_class==as.character(position$X1[i]))])
}
position=position[order(position$median,decreasing = F),]
position$col=rep(c("blue4","gray","brown3"),6)
result$col=ifelse(result$super_class%in%position$X1[c(1,4,7,10,13,16)],"blue4",ifelse(result$super_class%in%position$X1[c(2,5,8,11,14,17)],"gray","brown3")    )
p4=ggplot(result, aes(x=super_class, y=r_spearman,fill=col))+coord_flip()+geom_boxplot(width=0.2,outlier.size = 0.0001,cex=0.05,alpha=0.7,position = position_dodge(0.9))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =0.1,alpha=0.1,binwidth=0.015)+guides(colour=FALSE)+theme_bw()+scale_fill_manual(values=c("blue4","gray","brown3"))+theme(axis.text.x = element_text(angle =0,hjust = 1,vjust = 0.5),text = element_text(size=6),axis.title = element_text(size = 6),panel.grid=element_blank(),legend.key.size =unit(0.4,"line"))+labs(x=NULL, y="r (Spearman)")+scale_x_discrete(limits = position$X1)+theme(axis.text.x = element_text(angle = 0,hjust = 0.5,vjust = 1))
result$order=as.numeric(gsub("meta_","",result$meta))
result=result[order(result$order,decreasing = F),]

#############pca
#1054
library(vegan)
euclidean_distance = as.matrix(vegdist(scale(meta), method="euclidean"))
euclidean.cmd=cmdscale(euclidean_distance,k=100,eig = T) 
pc=matrix(NA,nrow = 100,ncol = 2)
colnames(pc)=c("pc","%")
for(i in 1:nrow(pc)){
  pc[i,1]=i
  pc[i,2]=round(euclidean.cmd$eig[i]/sum(euclidean.cmd$eig)*100,digits = 2)
}
sum(pc[,2])#72.66%
sum(pc[1:10,2])#38.93%
euclidean.cmd=euclidean.cmd$points
colnames(euclidean.cmd)=1:ncol(euclidean.cmd)
for(i in 1:ncol(euclidean.cmd)){
  colnames(euclidean.cmd)[i]=paste("PC",i,sep = "")
}
euclidean.cmd=data.frame(euclidean.cmd)
euclidean.cmd$subset=ifelse(row.names(euclidean.cmd)%in%gsub("_F","",id_fup),"with","without")
write.table(euclidean.cmd,file = "../result/result_1054sample_1183meta_100pc.txt",quote = F,sep = "\t")

#pheno~100pc
euclidean.cmd=read.delim("../result/result_1054sample_1183meta_100pc.txt",header = T,sep = "\t")
result=cor_spearman(pheno[,c(3,2,5,43,48:208)], euclidean.cmd[,1:100],c("PCs","phenotype"))
result=result[[1]]
result=result[which(result$Qval<0.05),]
write.table(result,"../result/result_asso_100pc_phenotype.txt",quote = F,sep = "\t",row.names = F)
result=cor_spearman(meta, euclidean.cmd[,1:2],c("PCs","meta"))
result=result[[1]]
result=result[which(result$Qval<0.05),]
result=result[which(abs(result$Cor)>0.2),]
result=result[which(result$meta%in%result$meta[duplicated(result$meta)]),]
result=merge(result[,1:5],annotation,by.x ="meta",by.y = "row.names",all = F)
variance=read.delim("../result/table_variation_lasso_meta_snp_diet_microbiome.txt",header = T,sep = "\t",row.names = 1)
variance=variance[which(row.names(variance)%in%result$meta),]
variance=merge(variance,result[,c(1:3,16,24)],by.x = "row.names",by.y ="meta",all = F)
library(ggplot2)
library(gridExtra)
library(ggExtra)
p1=ggplot (euclidean.cmd, aes(PC1,PC2)) + geom_point(alpha=1,size = 0.05)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+ theme(legend.position = "none")+stat_ellipse(type = "norm",linetype = 2)+lims(x=c(-30,30),y = c(-20, 20))
p2=ggplot (euclidean.cmd, aes(PC1,PC2,colour=subset)) + geom_point(alpha=1,size = 0.05)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+ theme(legend.position = "none")+stat_ellipse(geom = "polygon", alpha = 0.2,color=NA,aes(fill = subset))+scale_color_manual(values=c("brown2","gray"))+scale_fill_manual(values=c("brown2","gray"))+lims(x=c(-30,30),y = c(-20, 20))
p3=ggplot (euclidean.cmd, aes(PC1,PC2,color=meta$meta_360)) + geom_point(alpha=1,size = 1)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+scale_color_gradient(low="blue2",high="brown2")+lims(x=c(-30,30),y = c(-20, 20))
p4=ggplot (euclidean.cmd, aes(PC1,PC2,color=meta$meta_360)) + geom_point(alpha=1,size = 1)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+ theme(legend.position = "none")+scale_color_gradient(low="blue2",high="brown2")+lims(x=c(-30,30),y = c(-20, 20))
p5=ggMarginal(p4,colour = "gray95",fill="gray80",type="densigram",size = 10,cex=0,binwidth = 1)
p6=ggplot (euclidean.cmd, aes(PC1,PC2,color=meta$meta_365)) + geom_point(alpha=1,size = 1)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+scale_color_gradient(low="blue2",high="brown2")+lims(x=c(-30,30),y = c(-20, 20))
p7=ggplot (euclidean.cmd, aes(PC1,PC2,color=meta$meta_365)) + geom_point(alpha=1,size = 1)+xlab(label = "PC1 11.5%") + ylab(label = "PC2 5.6%")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank(),axis.line = element_line(colour = "black"),panel.border = element_blank(),axis.text = element_text(size = 5),axis.title = element_text(size = 6))+ theme(legend.position = "none")+scale_color_gradient(low="blue2",high="brown2")+lims(x=c(-30,30),y = c(-20, 20))
p8=ggMarginal(p7,colour = "gray95",fill="gray80",type="densigram",size = 10,cex=0,binwidth = 1)
pdf(file="../figure/figure_pca_1054samples_pca_raw_1183meta.pdf",useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,ncol=4,nrow = 4)
dev.off()

############################################################# mQTLs ###############################################################
##baseline mQTLs
#adj age, sex, smk, oral
colnames(pheno)[c(3,2,43,152)]
covariate=pheno[,c(3,2,43,152)]
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi
colnames(pheno)[c(3,2,43,152,5)]
covariate=pheno[,c(3,2,43,152,5)]
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet
colnames(pheno)[c(3,2,43,152,5,48:125)]
covariate=pheno[,c(3,2,43,152,5,48:125)]
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet, medication
colnames(pheno)[c(3,2,43,152,5,48:151,153:169)]
covariate=pheno[,c(3,2,43,152,5,48:151,153:169)]
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_med_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet, medication, disease
colnames(pheno)[c(3,2,43,152,5,48:151,153:208)]
covariate=pheno[,c(3,2,43,152,5,48:151,153:208)]
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_med_disease_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet, medication, disease, sp
covariate=cbind(pheno[,c(3,2,43,152,5,48:151,153:208)],sp[row.names(pheno),])
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_med_disease_sp_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet, medication, disease, sp, path
covariate=cbind(pheno[,c(3,2,43,152,5,48:151,153:208)],sp[row.names(pheno),],path[row.names(pheno),])
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_med_disease_sp_path_qtl.txt",sep = "\t",quote = F,row.names = F)
#adj age, sex, smk, oral, bmi, diet, medication, disease, sp, path, gutSMASH
write.table(bgc,file = "../data/data_1054samples_LLD_baseline_3075gutSMASH.txt",sep = "\t",quote = F)
bgc=read.delim("../data/data_1054samples_LLD_baseline_3075gutSMASH.txt",sep = "\t",header = T)
covariate=cbind(pheno[,c(3,2,43,152,5,48:151,153:208)],sp[row.names(pheno),],path[row.names(pheno),],bgc[row.names(pheno),])
save=linear_adj_meta(meta,covariate)
write.table(save,file = "./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_bmi_diet_med_disease_sp_path_gutSMASH_qtl.txt",sep = "\t",quote = F,row.names = F)
#replication adj age, sex and smk
covariate=pheno_rep[,c(1,2,4)]
save=linear_adj_meta(meta_rep,covariate)
write.table(save,file = "./mQTL/data_314samples_LLD_baseline_replication_1183plasma_metabolites_adj_age_sex_smk_qtl.txt",sep = "\t",quote = F,row.names = F)
#replication adj age, sex, smk and oral
covariate=pheno_fup[,c(3,2,27,37)]
save=linear_adj_meta(meta_fup,covariate)
write.table(save,file = "./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",sep = "\t",quote = F,row.names = F)
##follow-up delta
#delta adj age and sex
covariate=pheno[,c(3,2)]
row.names(covariate)=paste(row.names(covariate),"_F",sep = "")
save=linear_adj_delta(meta_fup-meta[gsub("_F","",id_fup),],covariate[row.names(meta_fup-meta[gsub("_F","",id_fup),]),])
colnames(save)=gsub("_F","",colnames(save))
write.table(save,file = "./mQTL/data_311samples_LLD_baseline_fup_1183plasma_metabolites_adj_baseline_age_sex_qtl.txt",sep = "\t",quote = F,row.names = F)

##metaQTLs replication
meta_dis=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_adj=data.frame(t(read.delim("./mQTL/data_314samples_LLD_baseline_replication_1183plasma_metabolites_adj_age_sex_smk_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
geno_lld=t(read.delim("../data/lld_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno_lld)=gsub("1_LLDeep","LLDeep",row.names(geno_lld))
geno_gonl=t(read.delim("../data/gonl_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
key=read.delim("../data/gte_gonl.txt",header = F,sep = "\t",check.names = F)
geno_gonl=merge(key,geno_gonl,by.x = "V1",by.y = "row.names",all = F)
row.names(geno_gonl)=geno_gonl[,2]
geno_gonl=geno_gonl[,-c(1:2)]
geno_gonl_type=t(read.delim("../data/gonl_snp_list.genotypes.txt",row.names = 1,header = T,sep = "\t",check.names = F))
geno_gonl_type=merge(key,geno_gonl_type,by.x = "V1",by.y = "row.names",all = F)
row.names(geno_gonl_type)=geno_gonl_type[,2]
geno_gonl_type=geno_gonl_type[,-c(1:2)]
qtl=read.delim("../result/result_73independent_qtl_63meta_57snp_1054samples_all_correction_modles.txt",header = T,sep = "\t")
#base discovery check
meta=meta_dis
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#933
meta=meta[row.names(snp),]
for(i in 1:nrow(qtl)){
  cor=cor.test(snp[,as.character(qtl$SNPName[i])],meta[,as.character(qtl$ProbeName[i])],method = "spearman")
  qtl$LLD_dis_P[i]=cor$p.value
  qtl$LLD_dis_r[i]=cor$estimate
}
#fup rep
meta=meta_fup
row.names(meta)=gsub("_F","",row.names(meta))
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#281
meta=meta[row.names(snp),]
for(i in 1:nrow(qtl)){
  cor=cor.test(snp[,as.character(qtl$SNPName[i])],meta[,as.character(qtl$ProbeName[i])],method = "spearman")
  qtl$LLD_fup_rep_P[i]=cor$p.value
  qtl$LLD_fup_rep_r[i]=cor$estimate
}
#LLD2 rep
meta=meta_adj
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#237
meta=meta[row.names(snp),]
for(i in 1:nrow(qtl)){
  cor=cor.test(snp[,as.character(qtl$SNPName[i])],meta[,as.character(qtl$ProbeName[i])],method = "spearman")
  qtl$LLD2_rep_P[i]=cor$p.value
  qtl$LLD2_rep_r[i]=cor$estimate
}
#gonl rep
library(stringr)
meta=meta_adj
snp=geno_gonl[intersect(row.names(meta),row.names(geno_gonl)),]#77
snp_type=geno_gonl_type[row.names(snp),]
meta=meta[row.names(snp),]
for(i in 1:nrow(qtl)){
  cor=cor.test(snp[,as.character(qtl$SNPName[i])],meta[,as.character(qtl$ProbeName[i])],method = "spearman")
  qtl$gonl_rep_P[i]=cor$p.value
  qtl$gonl_rep_r[i]=cor$estimate
  qtl$gonl_assessed_allele[i]=str_split_fixed(as.character(snp_type[which(snp[,as.character(qtl$SNPName[i])]==max(snp[,as.character(qtl$SNPName[i])],na.rm = T))[1],as.character(qtl$SNPName[i])]),"/",2)[1]
}
qtl$direction=ifelse(qtl$IncludedDatasetsCorrelationCoefficient*qtl$LLD_dis_r>0,1,-1)
qtl$LLD_fup_rep_r=qtl$LLD_fup_rep_r*qtl$direction
qtl$LLD2_rep_r=qtl$LLD2_rep_r*qtl$direction
qtl$direction_gonl=ifelse(qtl$AlleleAssessed==qtl$gonl_assessed_allele,1,-1)
qtl$gonl_rep_r=qtl$gonl_rep_r*qtl$direction_gonl
pdf("../figure/figure_metaQTL_replication.pdf",width=9, height=10, useDingbats = F)
par(mfrow=c(3,3),mgp=c(1.5,0.5,0))
plot(qtl$IncludedDatasetsCorrelationCoefficient,qtl$LLD_fup_rep_r,pch = 19,cex=0.5,col=ifelse(qtl$LLD_fup_rep_P>0.05,"gray","black"),xlim = c(-1,1),ylim = c(-1,1),xlab = "metaQTL effect size in LLD",ylab = "metaQTL effect size in LLD fup")
abline(a=0,b=1,h=0,v=0,lty=2,col="gray70")
plot(qtl$IncludedDatasetsCorrelationCoefficient,qtl$LLD2_rep_r,pch = 19,cex=0.5,col=ifelse(qtl$LLD2_rep_P>0.05,"gray","black"),xlim = c(-1,1),ylim = c(-1,1),xlab = "metaQTL effect size in LLD",ylab = "metaQTL effect size in LLD subset")
abline(a=0,b=1,h=0,v=0,lty=2,col="gray70")
plot(qtl$IncludedDatasetsCorrelationCoefficient,qtl$gonl_rep_r,pch = 19,cex=0.5,col=ifelse(qtl$gonl_rep_P>0.05,"gray","black"),xlim = c(-1,1),ylim = c(-1,1),xlab = "metaQTL effect size in LLD",ylab = "metaQTL effect size in GoNL")
abline(a=0,b=1,h=0,v=0,lty=2,col="gray70")
dev.off()
#replication rate 100% in at least one cohort
nrow(qtl[which(qtl$IncludedDatasetsCorrelationCoefficient*qtl$LLD_fup_rep_r>0&qtl$LLD_fup_rep_P<0.05|qtl$IncludedDatasetsCorrelationCoefficient*qtl$LLD2_rep_r>0&qtl$LLD2_rep_P<0.05|qtl$IncludedDatasetsCorrelationCoefficient*qtl$gonl_rep_r>0&qtl$gonl_rep_P<0.05),]) 
qtl=merge(qtl,annotation[,c(1,11:19)],by.x ="ProbeName",by.y = "meta",all = F)
qtl=qtl[,c(1,4:8,3,10,11,17:22,24,12:14,26:34)]
write.table(qtl,"../result/result_73independent_qtl_63meta_57snp_1054samples_all_correction_modles_replication.txt",quote = F,sep = "\t",row.names = F)
qtl=read.delim("../result/result_73independent_qtl_63meta_57snp_1054samples_all_correction_modles_replication.txt",header = T,sep = "\t")

################################################ Microbiome associations to metabolites ###########################################
#adj function
#adj age, sex, smk and oral from metabolites
meta_adj=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_delta_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_baseline_fup_1183plasma_metabolites_adj_baseline_age_sex_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
row.names(meta_delta_adj)=paste(row.names(meta_delta_adj),"_F",sep = "")
#adj age, sex, smk, ppi, antibiotics and laxatives from microbiome
covariate=pheno[,c(3,2,43,158,132,144)]
sp_adj=linear_adj_microbiome(sp,covariate)
path_adj=linear_adj_microbiome(path,covariate)
vsv_adj=linear_adj_microbiome(vsv,covariate)
bgc_adj=linear_adj_microbiome(bgc,covariate)
write.table(sp_adj,file = "../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(path_adj,file = "../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(vsv_adj,file = "../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(bgc_adj,file = "../data/data_1054samples_LLD_baseline_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
covariate=pheno_fup[,c(3,2,27,41)]
sp_fup_adj=linear_adj_microbiome(sp_fup,covariate)
path_fup_adj=linear_adj_microbiome(path_fup,covariate)
vsv_fup_adj=linear_adj_microbiome(vsv_fup,covariate)
bgc_fup_adj=linear_adj_microbiome(bgc_fup,covariate)
write.table(sp_fup_adj,file = "../data/data_311samples_LLD_fup_156sp_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(path_fup_adj,file = "../data/data_311samples_LLD_fup_343path_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(vsv_fup_adj,file = "../data/data_311samples_LLD_fup_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(bgc_fup_adj,file = "../data/data_311samples_LLD_fup_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
#adj baseline age and sex from delta abundance
covariate=pheno[,c(3,2)]
row.names(covariate)=paste(row.names(covariate),"_F",sep = "")
sp_delta_adj=linear_adj_microbiome(sp_fup-sp[gsub("_F","",row.names(sp_fup)),],covariate[row.names(sp_fup-sp[gsub("_F","",row.names(sp_fup)),]),])
path_delta_adj=linear_adj_microbiome(path_fup-path[gsub("_F","",row.names(path_fup)),],covariate[row.names(path_fup-path[gsub("_F","",row.names(path_fup)),]),])
vsv_delta_adj=linear_adj_microbiome(vsv_fup-vsv[gsub("_F","",row.names(vsv_fup)),],covariate[row.names(vsv_fup-vsv[gsub("_F","",row.names(vsv_fup)),]),])
bgc_delta_adj=linear_adj_microbiome(bgc_fup-bgc[gsub("_F","",row.names(bgc_fup)),],covariate[row.names(bgc_fup-bgc[gsub("_F","",row.names(bgc_fup)),]),])
write.table(sp_delta_adj,file = "../data/data_311samples_LLD_delta_156sp_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(path_delta_adj,file = "../data/data_311samples_LLD_delta_343path_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(vsv_delta_adj,file = "../data/data_311samples_LLD_delta_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
write.table(bgc_delta_adj,file = "../data/data_311samples_LLD_delta_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",quote = F,sep = "\t")
##cross-sectional and delta association
#bgc
bgc_adj=read.delim("../data/data_1054samples_LLD_baseline_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
bgc_fup_adj=read.delim("../data/data_311samples_LLD_fup_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
bgc_delta_adj=read.delim("../data/data_311samples_LLD_delta_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
result_dis=cor_spearman(meta_adj,bgc_adj,c("bgc","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,bgc_fup_adj,c("bgc","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta=cor_spearman(meta_delta_adj,bgc_delta_adj,c("bgc","meta"))[[1]]
colnames(result_delta)[3:5]=paste(colnames(result_delta)[3:5],"delta",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta[,3:6],by="name",all = T)
write.table(result,file = "../result/table_cor_meta_bgc_dis_rep_delta_associations_all.txt",quote = F,sep = "\t",row.names = F)
length(which(result$Qval_baseline<0.05))#55991
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0))#23662
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta>0))#4430
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=merge(result[,c(2:8,10,11)],annotation[,c(1,11:17)],by="meta",all = F)
result=result[which(result$Qval_baseline<0.05),]
write.table(result,file = "../result/table_cor_meta_bgc_dis_rep_delta_associations.txt",quote = F,sep = "\t",row.names = F)
#sp
result_dis=cor_spearman(meta_adj,sp_adj,c("species","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,sp_fup_adj,c("species","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta=cor_spearman(meta_delta_adj,sp_delta_adj,c("species","meta"))[[1]]
colnames(result_delta)[3:5]=paste(colnames(result_delta)[3:5],"delta",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta[,3:6],by="name",all = T)
write.table(result,file = "../result/table_cor_meta_sp_dis_rep_delta_associations_all.txt",quote = F,sep = "\t",row.names = F)
length(which(result$Qval_baseline<0.05))#3012
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0))#1373
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta>0))#234
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=merge(result[,c(2:8,10,11)],annotation[,c(1,11:17)],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_sp_dis_rep_delta_associations.txt",quote = F,sep = "\t",row.names = F)
#path
result_dis=cor_spearman(meta_adj,path_adj,c("pathway","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,path_fup_adj,c("pathway","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta=cor_spearman(meta_delta_adj,path_delta_adj,c("pathway","meta"))[[1]]
colnames(result_delta)[3:5]=paste(colnames(result_delta)[3:5],"delta",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta[,3:6],by="name",all = T)
write.table(result,file = "../result/table_cor_meta_path_dis_rep_delta_associations_all.txt",quote = F,sep = "\t",row.names = F)
length(which(result$Qval_baseline<0.05))#4965
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0))#2839
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta>0))#225
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=merge(result[,c(2:8,10,11)],annotation[,c(1,11:17)],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_path_dis_rep_delta_associations.txt",quote = F,sep = "\t",row.names = F)
#vsv
result_dis=cor_spearman(meta_adj,vsv_adj,c("vSVs","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,vsv_fup_adj,c("vSVs","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta=cor_spearman(meta_delta_adj,vsv_delta_adj,c("vSVs","meta"))[[1]]
colnames(result_delta)[3:5]=paste(colnames(result_delta)[3:5],"delta",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta[,3:6],by="name",all = T)
length(which(result$Qval_baseline<0.05))#
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0))#
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta>0))#
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=merge(result[,c(2:8,10,11)],annotation[,c(1,11:17)],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_vSV_dis_rep_delta_associations.txt",quote = F,sep = "\t",row.names = F)
#dsv
dsv_base=dsv[gsub("_F","",row.names(dsv_fup)),]
dsv_delta=dsv_fup
dsv_delta[is.na(dsv_delta)==F]=NA
for(i in 1:ncol(dsv_delta)){
  dsv_delta[,i][which(dsv_fup[,i]==0&dsv_base[,i]==0)]=0
  dsv_delta[,i][which(dsv_fup[,i]==0&dsv_base[,i]==1)]=1
  dsv_delta[,i][which(dsv_fup[,i]==1&dsv_base[,i]==0)]=2
  dsv_delta[,i][which(dsv_fup[,i]==1&dsv_base[,i]==1)]=3
}
dsv_delta_keep=dsv_delta
dsv_delta_keep[dsv_delta_keep==1]=NA
dsv_delta_keep[dsv_delta_keep==2]=NA
dsv_delta_change=dsv_delta
dsv_delta_change[dsv_delta_change==0]=NA
dsv_delta_change[dsv_delta_change==3]=NA
result_dis=cor_spearman(meta_adj,dsv,c("dSVs","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,dsv_fup,c("dSVs","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta_keep=cor_spearman(meta_delta_adj,dsv_delta_keep,c("dSVs","meta"))[[1]]
colnames(result_delta_keep)[3:5]=paste(colnames(result_delta_keep)[3:5],"delta_keep",sep = "_")
result_delta_change=cor_spearman(meta_delta_adj,dsv_delta_change,c("dSVs","meta"))[[1]]
colnames(result_delta_change)[3:5]=paste(colnames(result_delta_change)[3:5],"delta_change",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta_keep[,3:6],by="name",all = T)
result=merge(result,result_delta_change[,3:6],by="name",all = T)
length(which(result$Qval_baseline<0.05))#
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0))#
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta>0))#
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=merge(result[,c(2:8,10,11,13,14)],annotation[,c(1,11:17)],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_dSV_dis_rep_delta_associations.txt",quote = F,sep = "\t",row.names = F)

################################################ diet associations to metabolites ###########################################
#adj age, sex, smk and oral from metabolites
meta_adj=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
diet=pheno[,48:125]
result_dis=cor_spearman(meta_adj,diet,c("diet","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
write.table(result_dis,file = "../result/table_cor_meta_diet_associations.txt",quote = F,sep = "\t",row.names = F)
result=result_dis[which(result_dis$Qval_baseline<0.05),]
result=merge(annotation[,c(1,12:17)],result[,1:5],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_diet_associations_FDR.txt",quote = F,sep = "\t",row.names = F)
#adj age, sex and smk from diet score
covariate=pheno[,c(3,2,43)]
diet_score_adj=linear_adj_microbiome(diet_score,covariate)
write.table(diet_score_adj,file = "../data/data_1019samples_LLD_baseline_5diet_scores_adj_age_sex_smk.txt",quote = F,sep = "\t")
covariate=pheno_rep[,c(1,2,4)]
diet_score_rep_adj=linear_adj_microbiome(diet_score_rep,covariate)
write.table(diet_score_rep_adj,file = "../data/data_230samples_LLD_baseline_rep_5diet_scores_adj_age_sex_smk.txt",quote = F,sep = "\t")
###summary metabolites associations to diet, genetics and microbiome
annotation=read.delim("./key_lld_1183meta_annotation_type.txt",header = T,quote = "",sep = "\t",check.names = F)
annotation$name=gsub("\"","",annotation$name)
annotation$super_class=gsub("\"","",annotation$super_class)
annotation$class=gsub("\"","",annotation$class)
annotation$sub_class=gsub("\"","",annotation$sub_class)
annotation$direct_parent=gsub("\"","",annotation$direct_parent)
asso_diet=read.delim("../result/table_cor_meta_diet_associations_FDR.txt",header = T,sep = "\t",check.names = F,quote = "")
asso_diet=merge(asso_diet,annotation[,c(1,18,19)],by="meta",all = F)
asso_genetics=read.delim("../result/result_73independent_qtl_63meta_57snp_1054samples_all_correction_modles_replication.txt",header = T,sep = "\t",check.names = F,quote = "")
asso_microbiome=read.delim("../result/table_cor_meta_sp_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(asso_microbiome)[2]="microbiome"
tmp=read.delim("../result/table_cor_meta_path_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(tmp)[2]="microbiome"
asso_microbiome=rbind(asso_microbiome,tmp)
asso_microbiome=merge(asso_microbiome,annotation[,c(1,18,19,20)],by="meta",all = F)
length(intersect(as.character(asso_diet$meta),as.character(asso_genetics$ProbeName)))#49
length(intersect(intersect(as.character(asso_diet$meta),as.character(asso_genetics$ProbeName)),as.character(asso_microbiome$meta)))#6
length(intersect(intersect(as.character(asso_diet$meta),as.character(asso_genetics$ProbeName)),as.character(asso_microbiome$meta[which(asso_microbiome$Pval_delta<0.05)])))#6
intersect(intersect(as.character(asso_diet$meta),as.character(asso_genetics$ProbeName)),as.character(asso_microbiome$meta[which(asso_microbiome$Pval_delta<0.05)]))#2
##genetics of gut microbiome composition
meta_adj=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_delta_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_baseline_fup_1183plasma_metabolites_adj_baseline_age_sex_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
sp_adj=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
sp_fup_adj=read.delim("../data/data_311samples_LLD_fup_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
sp_delta_adj=read.delim("../data/data_311samples_LLD_delta_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
row.names(sp_delta_adj)=gsub("_F","",row.names(sp_delta_adj))
path_adj=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_fup_adj=read.delim("../data/data_311samples_LLD_fup_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_delta_adj=read.delim("../data/data_311samples_LLD_delta_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
row.names(path_delta_adj)=gsub("_F","",row.names(path_delta_adj))
diet=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_207pheno_log_imp_1135samples.txt",header = T,sep = "\t")
diet=diet[row.names(meta_adj),]
geno=t(read.delim("../data/lld_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno)=gsub("1_LLDeep","LLDeep",row.names(geno))
geno=geno[intersect(row.names(geno),row.names(meta_adj)),]
data=cbind(diet,meta_adj,sp_adj,path_adj)
data=merge(data,geno,by="row.names",all = T)
row.names(data)=data[,1]
data=data[,-1]
data_fup=cbind(meta_fup_adj,sp_fup_adj,path_fup_adj)
row.names(data_fup)=gsub("_F","",row.names(data_fup))
data_fup=merge(data_fup,geno[intersect(row.names(geno),row.names(data_fup)),],by="row.names",all = T)
row.names(data_fup)=data_fup[,1]
data_fup=data_fup[,-1]
data_delta=cbind(meta_delta_adj,sp_delta_adj,path_delta_adj)
data_delta=merge(data_delta,geno[intersect(row.names(geno),row.names(meta_delta_adj)),],by="row.names",all = T)
row.names(data_delta)=data_delta[,1]
data_delta=data_delta[,-1]
gene_micro=merge(asso_genetics,asso_microbiome[,1:9],by.x = "ProbeName",by.y = "meta",all = F)
gene_micro$Cor_baseline_g_m=NA
gene_micro$Pval_baseline_g_m=NA
gene_micro$Cor_fup_g_m=NA
gene_micro$Pval_fup_g_m=NA
for(i in 1:nrow(gene_micro)){
  gene_micro$Cor_baseline_g_m[i]=cor.test(data[,as.character(gene_micro$SNPName[i])],data[,as.character(gene_micro$microbiome[i])],method = "spearman")$estimate
  gene_micro$Pval_baseline_g_m[i]=cor.test(data[,as.character(gene_micro$SNPName[i])],data[,as.character(gene_micro$microbiome[i])],method = "spearman")$p.value
  gene_micro$Cor_fup_g_m[i]=cor.test(data_fup[,as.character(gene_micro$SNPName[i])],data_fup[,as.character(gene_micro$microbiome[i])],method = "spearman")$estimate
  gene_micro$Pval_fup_g_m[i]=cor.test(data_fup[,as.character(gene_micro$SNPName[i])],data_fup[,as.character(gene_micro$microbiome[i])],method = "spearman")$p.value
}
write.table(gene_micro,file = "../result/table_cor_genetics_on_microbiome.txt",quote = F,sep = "\t",row.names = F)
#ven plot
library(VennDiagram)
library(gridExtra)
diet=unique(as.character(asso_diet$meta))
genetics=unique(as.character(asso_genetics$ProbeName))
microbiome=unique(as.character(asso_microbiome$meta))
pdf("../figure/ven_genetics_diet_microbiome_overlap.pdf")
a=venn.diagram(list(Genetics=genetics, Microbiome=microbiome, Diet=diet),fill = c("brown3","blue4","gray90"),alpha = c(0.8,0.8,0.8), cex = 2,cat.fontface = 4,lty =1, filename = NULL)
grid.arrange(gTree(children = a),ncol=1)
dev.off()
tmp=merge(asso_diet[,c(1,8:11)],asso_genetics,by.x = "meta",by.y = "ProbeName",all = F)
tmp=merge(tmp,asso_microbiome[,1:9],by="meta",all = F)
write.table(tmp,file = "../result/table_cor_meta_diet_snp_microbiome_associations.txt",quote = F,sep = "\t",row.names = F)
#diet summary plot
data=data.frame(table(as.character(asso_diet$diet[which(asso_diet$Cor_baseline>0)])))
data$direction="pos"
tmp=data.frame(table(as.character(asso_diet$diet[which(asso_diet$Cor_baseline<0)])))
tmp$direction="neg"
data=rbind(data,tmp)
tmp=data
data=data.frame(table(as.character(asso_diet$diet)))
data=data[order(data$Freq,decreasing = T),]
data$id=1:nrow(data)
angle= 90-360*(data$id-0.5)/nrow(data)
data$hjust=ifelse(angle< -90, 1, 0)
data$angle=ifelse(angle < -90, angle+180, angle)
data$size=c(rep(1,17),rep(0.7,9),rep(0,48))
data=merge(tmp,data,by="Var1",all = F)
data$label=data$Var1
data$label=gsub("_log","",data$label)
data$label=gsub("_",".",data$label)
data$label[duplicated(as.character(data$label))]=NA
library(ggplot2)
library(viridis)
library(gridExtra)
#p=ggplot(data, aes(x=as.factor(id), y=Freq.x,fill=direction))+geom_bar(stat="identity", color="gray90")+theme_minimal()+theme(axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),plot.margin = unit(rep(-1,4), "cm"))+coord_polar(start = 0)+geom_text(data=data, aes(x=id, y=Freq.y+3, label=Var1, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=1, angle= data$angle, inherit.aes = FALSE)+scale_y_continuous(trans = 'log2')+scale_fill_viridis(discrete=TRUE)+annotate("text", x = max(data$id), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="brown3", size=2 , angle=0, fontface="bold", hjust=1)
#p=ggplot(data, aes(x=as.factor(id), y=Freq.x,fill=direction,alpha=0.2))+geom_bar(stat="identity", color="gray90",cex=0.1)+scale_fill_manual(values=c( "blue","brown3"))+theme_minimal()+theme(legend.position="none",axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),plot.margin = unit(rep(-1,4), "cm"))+coord_polar(start = 0)+scale_y_continuous(trans = 'log2')+geom_text(data=data, aes(x=id, y=Freq.y, label=label, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=0.7, angle= data$angle, inherit.aes = FALSE)+annotate("text", x = max(data$id), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="black", size=1 , angle=0, fontface="bold", hjust=1)
p=ggplot(data, aes(x=as.factor(id), y=Freq.x,fill=direction))+geom_bar(stat="identity", color="gray90",cex=0.1)+scale_fill_manual(values=c( "blue4","brown3"))+theme_minimal()+theme(legend.position="none",axis.text = element_blank(),axis.title = element_blank(),panel.grid = element_blank(),plot.margin = unit(rep(-1,4), "cm"))+coord_polar(start = 0)+geom_text(data=data, aes(x=id, y=Freq.y+3, label=label, hjust=hjust), color="black", fontface="bold",alpha=1, size=data$size, angle= data$angle, inherit.aes = FALSE)+annotate("text", x = max(data$id), y = c(0, 50, 100, 150, 200), label = c("0", "50", "100", "150", "200") , color="black", size=1, angle=0, fontface="bold", hjust=1)
pdf("../figure/figure_circul_barplot_meta_diet_asso.pdf",useDingbats = F)
grid.arrange(p,ncol=2,nrow = 2)
dev.off()

#adj age, sex, smk and oral from metabolites
meta_adj=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_delta_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_baseline_fup_1183plasma_metabolites_adj_baseline_age_sex_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_delta_adj=scale(meta_delta_adj)
sp_adj=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
sp_fup_adj=read.delim("../data/data_311samples_LLD_fup_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
sp_delta_adj=read.delim("../data/data_311samples_LLD_delta_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
row.names(sp_delta_adj)=gsub("_F","",row.names(sp_delta_adj))
path_adj=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_fup_adj=read.delim("../data/data_311samples_LLD_fup_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_delta_adj=read.delim("../data/data_311samples_LLD_delta_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
row.names(path_delta_adj)=gsub("_F","",row.names(path_delta_adj))
vsv_adj=read.delim("../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_fup_adj=read.delim("../data/data_311samples_LLD_fup_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_delta_adj=read.delim("../data/data_311samples_LLD_delta_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
row.names(vsv_delta_adj)=gsub("_F","",row.names(vsv_delta_adj))
diet=read.delim("../data/data_1054samples_LLD_baseline_78diets_adj_age_sex.txt",header = T,sep = "\t")
diet=diet[row.names(meta_adj),]
pheno=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_207pheno_log_imp_1135samples.txt",header = T,sep = "\t")
pheno=pheno[row.names(diet),]
pheno=pheno[,48:125]
diet[,c(26:30,39:42)]=pheno[,c(26:30,39:42)]
#geno=t(read.delim("../data/lld_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
geno=t(read.delim("../data/genotype_lld_generalmeta_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno)=gsub("1_LLDeep","LLDeep",row.names(geno))
geno=geno[intersect(row.names(geno),row.names(meta_adj)),]
data=cbind(diet,meta_adj,sp_adj,path_adj,vsv_adj)
data=merge(data,geno,by="row.names",all = T)
row.names(data)=data[,1]
data=data[,-1]
data_fup=cbind(meta_fup_adj,sp_fup_adj,path_fup_adj,vsv_fup_adj)
row.names(data_fup)=gsub("_F","",row.names(data_fup))
data_fup=merge(data_fup,geno[intersect(row.names(geno),row.names(data_fup)),],by="row.names",all = T)
row.names(data_fup)=data_fup[,1]
data_fup=data_fup[,-1]
data_delta=cbind(meta_delta_adj,sp_delta_adj,path_delta_adj,vsv_delta_adj)
data_delta=merge(data_delta,geno[intersect(row.names(geno),row.names(meta_delta_adj)),],by="row.names",all = T)
row.names(data_delta)=data_delta[,1]
data_delta=data_delta[,-1]
library(ggplot2)
library(gridExtra)
#diet
p1=ggplot(data,aes(x=coffee_log,y=meta_365))+geom_point(alpha=1,size=1.5,shape=21,fill="blue4",color="brown3")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Coffee intake", y="Quinic acid")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-3, 3))
p2=ggplot(data,aes(x=how_often_fish,y=meta_377))+geom_point(alpha=1,size=1.5,shape=21,fill="blue4",color="brown3")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="How often fish", y="2,6-Dimethoxy-4-propylphenol")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-3, 3))
p3=ggplot(data,aes(x=as.character(data$vegetarian),y=data$meta_284,colour=as.character(data$vegetarian)))+geom_violin(alpha=0.3,cex=0.4,colour="gray70")+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =0.5,alpha=0.8,binwidth=0.03,fill="white",stackratio=1.5)+scale_fill_manual(values=c("blue4","brown3"))+geom_boxplot(width=0.1,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Vegetarian", y="1-Methylhistidine")
#genetics
data_tmp1=data[is.na(data$rs1495741)==F,]
p4=ggplot(data_tmp1,aes(x=as.character(round(rs1495741,digits = 0)),y=meta_479,colour=as.character(round(rs1495741,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.01,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs1495741 (NAT2)", y="5-Acetylamino-6-formylamino-3-methyluracil")+lims(y = c(-0.7, 1.1))
data_tmp2=data[is.na(data$rs174564)==F,]
p5=ggplot(data_tmp2,aes(x=as.character(round(rs174564,digits = 0)),y=meta_1070,colour=as.character(round(rs174564,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.01,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs174564 (FADS2)", y="PA(16:0/22:5(4Z,7Z,10Z,13Z,16Z))")+lims(y = c(-0.6, 0.6))
data_tmp3=data[is.na(data$rs4149056)==F,]
p6=ggplot(data_tmp3,aes(x=as.character(round(rs4149056,digits = 0)),y=meta_1013,colour=as.character(round(rs4149056,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.01,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs4149056 (SLCO1B1)", y="5'-Carboxy-gamma-chromanol")+lims(y = c(-0.6, 0.6))
p14=ggplot(data_tmp3,aes(x=as.character(round(rs4149056,digits = 0)),y=meta_670,colour=as.character(round(rs4149056,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.01,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs4149056 (SLCO1B1)", y="3-O-Protocatechuoylceanothic acid")+lims(y = c(-0.6, 0.6))
data_tmp4=data[is.na(data$rs261332)==F,]
p7=ggplot(data_tmp4,aes(x=as.character(round(rs261332,digits = 0)),y=meta_1079,colour=as.character(round(rs261332,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.01,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs261332 (LIPC)", y="PC(15:0/18:4(6Z,9Z,12Z,15Z))")+lims(y = c(-0.6, 0.8))
data_tmp5=data_delta[is.na(data_delta$rs12331618)==F,]
p8=ggplot(data_tmp5,aes(x=as.character(round(rs12331618,digits = 0)),y=meta_628,colour=as.character(round(rs12331618,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=250,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs12331618 (KLKB1)", y="Arginyl-isoleucine changes")+lims(y = c(-10000,10000))
data_tmp6=data_delta[is.na(data_delta$rs12331618)==F,]
p9=ggplot(data_tmp6,aes(x=as.character(round(rs12331618,digits = 0)),y=meta_466,colour=as.character(round(rs12331618,digits = 0))))+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=150,fill="white",stackratio=1.2)+scale_fill_manual(values=c("blue4","brown3","gray60"))+geom_boxplot(width=0.5,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","gray60"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="rs12331618 (KLKB1)", y="N-Acetylgalactosamine changes")+lims(y = c(-6000, 5000))
#meta_628~rs12331618 r=0.46 p=6.6e-16; meta_466~rs12331618 r=0.44 p=1.6e-14     
#microbiome
p10=ggplot(data,aes(x=`k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Verrucomicrobiaceae|g__Akkermansia|s__Akkermansia_muciniphila`,y=meta_586))+geom_point(alpha=1,size=1.5,shape=21,fill="gray90",color="blue4")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Akkermansia muciniphila", y="Thiamine")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-2.2, 2.5))
p11=ggplot(data,aes(x=`PWY-5198: factor 420 biosynthesis`,y=meta_892))+geom_point(alpha=1,size=1.5,shape=21,fill="gray90",color="dodgerblue2")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Factor 420", y="2,3-Dehydrosilybin")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-0.5, 0.75))
p12=ggplot(data,aes(x=`Ruminococcus sp. SR1/5:300_305`,y=meta_457))+geom_point(alpha=1,size=1.5,shape=21,fill="gray90",color="brown3")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="vSV_Ruminococcus sp (300-305kb)", y="Tyrosol 4-sulfate")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-1.2, 2.1))
p13=ggplot(data,aes(x=`Ruminococcus sp. SR1/5:300_305`,y=meta_501))+geom_point(alpha=1,size=1.5,shape=21,fill="gray90",color="brown3")+geom_smooth(method = "lm",color="white",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="vSV_Ruminococcus sp (300-305kb)", y="(4-ethyl-2-methoxyphenyl)oxidanesulfonic acid")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-1.2, 3))
pdf(file="../figure/figure_example.pdf",useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p14,p6,p7,p8,p9,p10,p11,p12,p13,ncol=4,nrow = 4)
dev.off()
#MR correlation plots
data_tmp=data_fup
row.names(data_tmp)=paste(row.names(data_tmp),"_F",sep = "")
data_tmp=rbind(data_tmp,data[,colnames(data_tmp)])
data_tmp$cohort="baseline"
data_tmp$cohort[grep("_F",row.names(data_tmp))]="fup"
p1=ggplot(data_tmp,aes(x=data_tmp$`k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_rectale`,y=data_tmp$meta_63,color=data_tmp$cohort))+geom_point(alpha=1,size=0.1)+geom_smooth(method = "lm",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Eubacterium rectale abundance", y="p-Cresol level")+scale_fill_manual(values=c("blue4","brown3"))+scale_color_manual(values=c("blue4","brown3"))+lims(x = c(-2, 2),y = c(-3.5, 2.5))
p2=ggplot(data_tmp,aes(x=data_tmp$`k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_rectale`,y=data_tmp$meta_354,color=data_tmp$cohort))+geom_point(alpha=1,size=0.1)+geom_smooth(method = "lm",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Eubacterium rectale abundance", y="p-Cresol sulfate level")+scale_fill_manual(values=c("blue4","brown3"))+scale_color_manual(values=c("blue4","brown3"))+lims(x = c(-2, 2),y = c(-3.5, 2.5))
p3=ggplot(data_delta,aes(x=data_delta$`k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_rectale`,y=data_delta$meta_63))+geom_point(alpha=1,size=1.5,shape=21,fill="blue4",color="brown3")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Delat Eubacterium rectale abundance", y="Delta p-Cresol level")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))
p4=ggplot(data_delta,aes(x=data_delta$`k__Bacteria|p__Firmicutes|c__Clostridia|o__Clostridiales|f__Eubacteriaceae|g__Eubacterium|s__Eubacterium_rectale`,y=data_delta$meta_354))+geom_point(alpha=1,size=1.5,shape=21,fill="blue4",color="brown3")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Delat Eubacterium rectale abundance", y="Delta p-Cresol sulfate level")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))
#Mediation correlation plots
p5=ggplot(data,aes(x=as.character(data$white_wine),y=data$meta_126,colour=as.character(data$white_wine)))+geom_violin(alpha=0.3,cex=0.4,colour="gray70")+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =0.5,alpha=0.8,binwidth=0.015,fill="white",stackratio=1)+scale_fill_manual(values=c("blue4","brown3","blue4","brown3"))+geom_boxplot(width=0.1,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","blue4","brown3"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="White wine consumption", y="Pipecolic acid level")+lims(y = c(-0.5, 0.8))
p6=ggplot(data,aes(x=as.character(data$white_wine),y=data$`[Eubacterium] hallii DSM 3353:734_736;754_756`,colour=as.character(data$white_wine)))+geom_violin(alpha=0.3,cex=0.4,colour="gray70")+geom_dotplot(binaxis = "y", stackdir = "center", position = "dodge",dotsize =1,alpha=0.8,binwidth=0.05,fill="white",stackratio=1)+scale_fill_manual(values=c("blue4","brown3","blue4","brown3"))+geom_boxplot(width=0.1,outlier.size = 0.01,colour="black",cex=0.2,alpha=0.2)+scale_color_manual(values=c("blue4","brown3","blue4","brown3"))+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="White wine consumption", y="vSV_Eubacterium hallii:734_736;754_756")
#p5=ggplot(data,aes(x=data$white_wine,y=data$meta_126))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="White wine consumption", y="Pipecolic acid level")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))
#p6=ggplot(data,aes(x=data$white_wine,y=data$`[Eubacterium] hallii DSM 3353:734_736;754_756`))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="White wine consumption", y="vSV_Eubacterium hallii:734_736;754_756")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))
p7=ggplot(data_tmp,aes(x=data_tmp$`[Eubacterium] hallii DSM 3353:734_736;754_756`,y=data_tmp$meta_126,color=data_tmp$cohort))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99")+geom_smooth(method = "lm",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="vSV_Eubacterium hallii:734_736;754_756", y="Pipecolic acid level")+scale_fill_manual(values=c("blue4","brown3"))+scale_color_manual(values=c("blue4","brown3"))+lims(y = c(-1, 2))
p8=ggplot(data_delta,aes(x=data_delta$`[Eubacterium] hallii DSM 3353:734_736;754_756`,y=data_delta$meta_126))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Delat vSV_Eubacterium hallii:734_736;754_756", y="Delta Pipecolic acid level")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))
p9=ggplot(data,aes(x=data$coffee_log,y=data$meta_370))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Coffee consumption", y="trans-Ferulic acid")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-1, 1))
p10=ggplot(data,aes(x=data$coffee_log,y=data$`Ruminococcus sp. SR1/5:300_305`))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x="Coffee consumption", y="vSV_Ruminococcus sp:300_305")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(y = c(-3, 4))
p11=ggplot(data_tmp,aes(y=data_tmp$`Ruminococcus sp. SR1/5:300_305`,x=data_tmp$meta_370,color=data_tmp$cohort))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99")+geom_smooth(method = "lm",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(y="vSV_Ruminococcus sp:300_305", x="Pipecolic acid level")+scale_fill_manual(values=c("blue4","brown3"))+scale_color_manual(values=c("blue4","brown3"))+lims(x = c(-0.6, 0.6),y = c(-4, 5))
p12=ggplot(data_delta,aes(y=data_delta$`Ruminococcus sp. SR1/5:300_305`,x=data_delta$meta_370))+geom_point(alpha=1,size=1.5,shape=21,fill="gray99",color="blue4")+geom_smooth(method = "lm",color="brown3",size=0.5,linetype=2)+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(y="Delat vSV_Ruminococcus sp:300_305", x="Delta Pipecolic acid level")+scale_fill_manual(values=c("blue4","blue4"))+scale_color_manual(values=c("blue4","blue4"))+lims(x = c(-10000, 10000),y = c(-3, 3))
pdf(file="../figure/figure_example_MR_Mediation_correlation_plot.pdf",useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,ncol=4,nrow = 4)
dev.off()
#partial correlation
library("corpcor")
tmp=data[,c("white_wine","[Eubacterium] hallii DSM 3353:734_736;754_756","meta_126")]
tmp=na.omit(tmp)
cor=pcor(tmp,method = "spearman")$estimate
p=pcor(tmp,method = "spearman")$p.value
tmp=data[,c("coffee_log","Ruminococcus sp. SR1/5:300_305","meta_370")]
tmp=na.omit(tmp)
cor=pcor(tmp,method = "spearman")$estimate
p=pcor(tmp,method = "spearman")$p.value

#pheno
pheno_fup=pheno_fup[,c(5:26,46:62)]
pheno=pheno[,colnames(pheno_fup)]
pheno_base=pheno[gsub("_F","",row.names(pheno_fup)),]
pheno_delta=pheno_fup
for(i in 1:ncol(pheno)){
  if(i%in%c(1:22)){
    pheno_delta[,i]=pheno_fup[,i]-pheno_base[,i]
  }
  if(i%in%c(23:39)){
    pheno_delta[,i][which(pheno_fup[,i]==0&pheno_base[,i]==0)]=0
    pheno_delta[,i][which(pheno_fup[,i]==0&pheno_base[,i]==1)]=1
    pheno_delta[,i][which(pheno_fup[,i]==1&pheno_base[,i]==0)]=2
    pheno_delta[,i][which(pheno_fup[,i]==1&pheno_base[,i]==1)]=3
  }
}
pheno_delta_keep=pheno_delta
tmp=pheno_delta[,23:39]
tmp[tmp==1]=NA
tmp[tmp==2]=NA
pheno_delta_keep[,23:39]=tmp
pheno_delta_change=pheno_delta
tmp=pheno_delta[,23:39]
tmp[tmp==0]=NA
tmp[tmp==3]=NA
pheno_delta_keep[,23:39]=tmp
result_dis=cor_spearman(meta_adj,pheno,c("pheno","meta"))[[1]]
colnames(result_dis)[3:5]=paste(colnames(result_dis)[3:5],"baseline",sep = "_")
result_rep=cor_spearman(meta_fup_adj,pheno_fup,c("pheno","meta"))[[1]]
colnames(result_rep)[3:5]=paste(colnames(result_rep)[3:5],"fup",sep = "_")
result_delta_keep=cor_spearman(meta_delta_adj,pheno_delta_keep,c("pheno","meta"))[[1]]
colnames(result_delta_keep)[3:5]=paste(colnames(result_delta_keep)[3:5],"delta_keep",sep = "_")
result_delta_change=cor_spearman(meta_delta_adj,pheno_delta_change,c("pheno","meta"))[[1]]
colnames(result_delta_change)[3:5]=paste(colnames(result_delta_change)[3:5],"delta_change",sep = "_")
result_delta=cor_spearman(meta_delta_adj,pheno_delta,c("pheno","meta"))[[1]]
colnames(result_delta)[3:5]=paste(colnames(result_delta)[3:5],"delta",sep = "_")
result=merge(result_dis,result_rep[,3:6],by="name",all = T)
result=merge(result,result_delta[,3:6],by="name",all = T)
result=merge(result,result_delta_keep[,3:6],by="name",all = T)
result=merge(result,result_delta_change[,3:6],by="name",all = T)
tmp=result[result$Qval_baseline<0.05,]
tmp=data.frame(table(tmp$pheno))
length(which(result$Qval_baseline<0.05))#2585
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05))#1101
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_keep<0.05|result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_change<0.05))#367
length(which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_keep<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta_keep>0|result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_change<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta_change>0))#363
#result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Cor_baseline*result$Cor_fup>0),]
result=result[which(result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_keep<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta_keep>0|result$Qval_baseline<0.05&result$Pval_fup<0.05&result$Pval_delta_change<0.05&result$Cor_baseline*result$Cor_fup>0&result$Cor_baseline*result$Cor_delta_change>0),]
#result=result[which(result$Pval_delta<0.05),]

result=merge(result[,c(2:8,10,11,13,14,16,17)],annotation[,c(1,11:19)],by="meta",all = F)
write.table(result,file = "../result/table_cor_meta_path_dis_rep_delta_225associations.txt",quote = F,sep = "\t",row.names = F)

####two sample MR (genetics impact metabolites through microbiome)
library(TwoSampleMR)
library(MRInstruments)
#microbiomeQTL discovery and replication
source('./code/manhattan_function.R')
#microbiome qtl
geno_lld=t(read.delim("../data/lld_snp_list.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno_lld)=gsub("1_LLDeep","LLDeep",row.names(geno_lld))
sp_base=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_base=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_base=read.delim("../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
save=t(vsv_base[,1189:1782])
save=cbind(RowNames = rownames(save), save)
colnames(save)[1]=" "
write.table(save,file = "../data/data_1054samples_LLD_baseline_1782vsv_b3_adj_age_sex_smk_ppi_anti_lax_qtl.txt",sep = "\t",quote = F,row.names = F)
fake_annotation=data.frame(matrix(NA,nrow=nrow(save),ncol=6))
colnames(fake_annotation)=c("platform","id","name","chr","start","end")
fake_annotation$platform="metabo"
fake_annotation$id=row.names(save)
fake_annotation$name=row.names(save)
fake_annotation$chr=1
fake_annotation$start=1000
fake_annotation$end=1000
write.table(fake_annotation,file = "fake_annotation_vsv_b3.txt",quote = F,row.names = F,sep = "\t")
#all
manhattan=read.delim("../result/eQTL_adj_age_sex_smk_ppi_anti_lax_sp.txt",header = T,sep = "\t",fill = T)
all=manhattan
manhattan=read.delim("../result/eQTL_adj_age_sex_smk_ppi_anti_lax_path.txt",header = T,sep = "\t",fill = T)
all=rbind(all,manhattan)
manhattan=read.delim("../result/eQTL_adj_age_sex_smk_ppi_anti_lax_vsv1.txt",header = T,sep = "\t",fill = T)
all=rbind(all,manhattan)
manhattan=read.delim("../result/eQTL_adj_age_sex_smk_ppi_anti_lax_vsv2.txt",header = T,sep = "\t",fill = T)
all=rbind(all,manhattan)
manhattan=read.delim("../result/eQTL_adj_age_sex_smk_ppi_anti_lax_vsv3.txt",header = T,sep = "\t",fill = T)
all=rbind(all,manhattan)
all=all[which(all$PValue<1e-5),]
all=all[order(all$PValue,decreasing = F),]
manhattan=all
manhattan$name=paste(manhattan$SNPName,manhattan$ProbeName,sep="_with_")
all=all[,c(2:4,1,5)]
colnames(all)=c("SNP","CHR","BP","P","microbiome")
trait=unique(as.character(all$microbiome))
length(unique(trait))#166#2281
result=NULL
for(i in 1:length(trait)){
  tmp=all[which(all$microbiome==trait[i]),]
  tmp=tmp[,c(1:3,5,4)]
  colnames(tmp)=c("SNP","chr_name","chrom_start","id.exposure","pval.exposure")
  clump=clump_data(tmp, clump_kb = 500, clump_r2 = 0.05, clump_p1 = 1,clump_p2 = 1)
  result=rbind(result,clump)
}
result$name=paste(result$SNP,result$id.exposure,sep="_with_")
manhattan=manhattan[which(manhattan$name%in%result$name),]
length(unique(manhattan$SNPName))#167
write.table(manhattan,"../result/qtl_1054samples_snp_microbiome.txt",quote = F,sep = "\t",row.names = F)
write.table(unique(as.character(manhattan$SNPName)),"../result/result_snp_1054samples_snp_microbiome.txt",quote = F,sep = "\t",row.names = F)
#write.table(unique(as.character(manhattan$SNPName)),"../result/result_167snp_1054samples_snp_microbiome.txt",quote = F,sep = "\t",row.names = F)
#replication in lld fup
manhattan=read.delim("../result/qtl_1054samples_snp_microbiome.txt",header = T,sep = "\t")
geno_lld=t(read.delim("../data/lld_snp_list_microbiome.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno_lld)=gsub("1_LLDeep","LLDeep",row.names(geno_lld))
sp_base=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_base=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_base=read.delim("../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
microbiome_base=cbind(sp_base,path_base,vsv_base)
sp_fup=read.delim("../data/data_311samples_LLD_fup_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_fup=read.delim("../data/data_311samples_LLD_fup_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_fup=read.delim("../data/data_311samples_LLD_fup_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
microbiome_fup=cbind(sp_fup,path_fup,vsv_fup)
row.names(microbiome_fup)=gsub("_F","",row.names(microbiome_fup))
#base discovery check
meta=microbiome_base
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#933
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$ProbeName[i])],method = "spearman")
  manhattan$LLD_dis_P[i]=cor$p.value
  manhattan$LLD_dis_r[i]=cor$estimate
}
#fup rep
meta=microbiome_fup
row.names(meta)=gsub("_F","",row.names(meta))
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#281
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$ProbeName[i])],method = "spearman")
  manhattan$LLD_fup_rep_P[i]=cor$p.value
  manhattan$LLD_fup_rep_r[i]=cor$estimate
}
manhattan$direction=ifelse(manhattan$LLD_dis_r*manhattan$IncludedDatasetsCorrelationCoefficient>0,1,-1)
manhattan$LLD_fup_rep_r=manhattan$LLD_fup_rep_r*manhattan$direction
manhattan$replication=ifelse(manhattan$LLD_fup_rep_P<0.05&manhattan$LLD_fup_rep_r*manhattan$IncludedDatasetsCorrelationCoefficient>0,"yes","no")
manhattan=manhattan[,c(23,2:4,9,10,5,1,18,20,26:29)]
manhattan=manhattan[which(manhattan$replication=="yes"),]#3019
gene_name=read.delim("~/Documents/UMCG/research-projects/microbiome_data/database/NCBI37.3.gene.loc_MAGMA",header = F,sep = "")
manhattan$gene=NA
for(i in 1:nrow(manhattan)){
  tmp=gene_name[which(gene_name$V2==manhattan$SNPChr[i]&(manhattan$SNPChrPos[i]+100000)>gene_name$V3&(manhattan$SNPChrPos[i]-100000)<gene_name$V4),]
  if(ncol(tmp)>0){
    gene="NA_gene"
    for(j in 1:nrow(tmp)){
      gene=paste(gene,tmp$V6[j],sep = ", ")
    }
    manhattan$gene[i]=gene
  }else{
    manhattan$gene[i]=" "
  }
}
manhattan$gene=gsub("NA_gene, ","",manhattan$gene)
manhattan$gene_within=NA
for(i in 1:nrow(manhattan)){
  tmp=gene_name[which(gene_name$V2==manhattan$SNPChr[i]&(manhattan$SNPChrPos[i])>=gene_name$V3&(manhattan$SNPChrPos[i])<=gene_name$V4),]
  if(ncol(tmp)>0){
    gene="NA_gene"
    for(j in 1:nrow(tmp)){
      gene=paste(gene,tmp$V6[j],sep = ", ")
    }
    manhattan$gene_within[i]=gene
  }else{
    manhattan$gene_within[i]=" "
  }
}
manhattan$gene_within=gsub("NA_gene, ","",manhattan$gene_within)
manhattan$gene_within=gsub("NA, ","",manhattan$gene_within)
tmp=data.frame(table(as.character(manhattan$ProbeName)))
tmp=unique(as.character(tmp$Var1[which(tmp$Freq>=3)]))#389
manhattan=manhattan[which(manhattan$ProbeName%in%tmp),]
write.table(manhattan,"../result/result_1549independent_qtl_389microbe_1471snp_1054samples_replication.txt",quote = F,sep = "\t",row.names = F)
###microbiome association result
manhattan=read.delim("../result/result_1549independent_qtl_389microbe_1471snp_1054samples_replication.txt",header = T,sep = "\t")
asso_sp=read.delim("../result/table_cor_meta_sp_dis_rep_delta_associations.txt",header = T,sep = "\t",quote = "")
asso_sp=asso_sp[which(asso_sp$Qval_baseline<0.05&asso_sp$Pval_fup<0.05&asso_sp$Pval_delta<0.05&asso_sp$Cor_baseline*asso_sp$Cor_fup>0&asso_sp$Cor_baseline*asso_sp$Cor_delta>0),]
colnames(asso_sp)[2]="microbiome"
asso_path=read.delim("../result/table_cor_meta_path_dis_rep_delta_associations.txt",header = T,sep = "\t",quote = "")
asso_path=asso_path[which(asso_path$Qval_baseline<0.05&asso_path$Pval_fup<0.05&asso_path$Pval_delta<0.05&asso_path$Cor_baseline*asso_path$Cor_fup>0&asso_path$Cor_baseline*asso_path$Cor_delta>0),]
colnames(asso_path)[2]="microbiome"
#asso_vsv=read.delim("../result/table_cor_meta_vSV_dis_rep_delta_associations.txt",header = T,sep = "\t")
#colnames(asso_vsv)[2]="microbiome"
asso_microbiome=rbind(asso_sp,asso_path)
asso_microbiome=asso_microbiome[which(asso_microbiome$microbiome%in%intersect(as.character(asso_microbiome$microbiome),as.character(manhattan$ProbeName))),]
manhattan=manhattan[which(manhattan$ProbeName%in%intersect(as.character(asso_microbiome$microbiome),as.character(manhattan$ProbeName))),]
write.table(manhattan,file = "../result/result_microbeQTL_37microbe_45metabolite.txt",quote = F,sep = "\t",row.names = F)
length(unique(as.character(manhattan$ProbeName)))#37 microbial abundance
length(unique(as.character(asso_microbiome$meta)))#45 metabolites
###MR
geno_lld=t(read.delim("../data/lld_snp_list_microbiome.dosages.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno_lld)=gsub("1_LLDeep","LLDeep",row.names(geno_lld))
geno_lld=geno_lld[,unique(as.character(manhattan$SNPName))]
#meta_adj=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_adj=data.frame(t(read.delim("./mQTL/data_927samples_LLD_baseline_1183plasma_metabolites_adj_2g_all_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
row.names(meta_fup_adj)=gsub("_F","",row.names(meta_fup_adj))
exposure=unique(as.character(manhattan$ProbeName))
outcome=unique(as.character(asso_microbiome$meta))
mr=merge(manhattan[,2:ncol(manhattan)],asso_microbiome,by.x = "ProbeName",by.y = "microbiome",all = T)
for(i in 1:nrow(mr)){
  mr$outcome_baseline_r[i]=cor.test(geno_lld[intersect(row.names(geno_lld),row.names(meta_adj)),as.character(mr$SNPName[i])],meta_adj[intersect(row.names(geno_lld),row.names(meta_adj)),as.character(mr$meta[i])],method = "spearman")$estimate
  mr$outcome_baseline_P[i]=cor.test(geno_lld[intersect(row.names(geno_lld),row.names(meta_adj)),as.character(mr$SNPName[i])],meta_adj[intersect(row.names(geno_lld),row.names(meta_adj)),as.character(mr$meta[i])],method = "spearman")$p.value
  mr$outcome_fup_r[i]=cor.test(geno_lld[intersect(row.names(geno_lld),row.names(meta_fup_adj)),as.character(mr$SNPName[i])],meta_fup_adj[intersect(row.names(geno_lld),row.names(meta_fup_adj)),as.character(mr$meta[i])],method = "spearman")$estimate
  mr$outcome_fup_P[i]=cor.test(geno_lld[intersect(row.names(geno_lld),row.names(meta_fup_adj)),as.character(mr$SNPName[i])],meta_fup_adj[intersect(row.names(geno_lld),row.names(meta_fup_adj)),as.character(mr$meta[i])],method = "spearman")$p.value
}
mr$outcome_baseline_r=mr$outcome_baseline_r*mr$direction
mr$outcome_fup_r=mr$outcome_fup_r*mr$direction
write.table(mr,file = "../result/mr_37microbial_abundance_45metabolites_new.txt",quote = F,sep = "\t",row.names = F)
#read table
library(TwoSampleMR)
mr=read.delim("../result/mr_37microbial_abundance_45metabolites_new.txt",header = T,sep = "\t")
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
result_mr$name=paste(result_mr$exposure,result_mr$outcome,sep = "_with_")
write.table(result_mr,file = "../result/result_mr_microbiome_to_meta_new.txt",quote = F,sep = "\t",row.names = F)

##bi-directional MR for 10 significant links
#10 significant links
result_mr=read.delim("../result/result_mr_microbiome_to_meta_new.txt",header = T,sep = "\t")
tmp=result_mr[which(result_mr$method=="Inverse variance weighted"),]
tmp$discovery_pval_BH_adj=p.adjust(tmp$discovery_pval,method = "BH")
tmp=tmp[which(tmp$discovery_pval_BH_adj<0.05),]
tmp=tmp[which(tmp$rep_pval<0.05),]
tmp_wm=result_mr[which(result_mr$method=="Weighted median"),]
tmp_wm=tmp_wm[which(tmp_wm$name%in%tmp$name),]
colnames(tmp_wm)[7:9]=paste(colnames(tmp_wm)[7:9],"_WM",sep = "")
tmp_egger=result_mr[which(result_mr$method=="MR Egger"),]
tmp_egger=tmp_egger[which(tmp_egger$name%in%tmp$name),]
colnames(tmp_egger)[8:9]=paste(colnames(tmp_egger)[8:9],"_pval_Egger",sep = "")
tmp=merge(tmp,tmp_egger[,c(13,8,9)],by="name",all = F)
tmp=merge(tmp,tmp_wm[,c(13,7)],by="name",all = F)
tmp=tmp[,c(1,2,3,5:8,14,17,15,16,11:13)]
colnames(tmp)[c(5:8,12:14)]=paste(colnames(tmp)[c(5:8,12:14)],"IVW",sep = "_")
result_mr=tmp
result_micro_to_meta=tmp
write.table(tmp,file = "../result/result_mr_microbiome_to_meta_4links_new.txt",quote = F,sep = "\t",row.names = F)
#outcome > exposure qtl
qtl=read.delim("../result/eQTL_adj_2gpc_all_pheno.txt",header = T,sep = "\t")
qtl=qtl[which(qtl$PValue<1e-5),]
qtl=qtl[which(qtl$ProbeName%in%result_mr$outcome),]
qtl$name=paste(qtl$SNPName,qtl$ProbeName,sep="_with_")
all=qtl[,c(2:4,1,5)]
colnames(all)=c("SNP","CHR","BP","P","meta")
trait=unique(as.character(all$meta))
length(unique(trait))#9 4
result=NULL
for(i in 1:length(trait)){
  tmp=all[which(all$meta==trait[i]),]
  tmp=tmp[,c(1:3,5,4)]
  colnames(tmp)=c("SNP","chr_name","chrom_start","id.exposure","pval.exposure")
  clump=clump_data(tmp, clump_kb = 500, clump_r2 = 0.05, clump_p1 = 1,clump_p2 = 1)
  result=rbind(result,clump)
}
result$name=paste(result$SNP,result$id.exposure,sep="_with_")
qtl=qtl[which(qtl$name%in%result$name),]
length(unique(qtl$SNPName))#361 54
write.table(qtl,"../result/qtl_927samples_snp_meta_inverse_direction_MR.txt",quote = F,sep = "\t",row.names = F)
write.table(unique(as.character(qtl$SNPName)),"../result/result_snp_927samples_snp_meta_inverse_direction_MR.txt",quote = F,sep = "\t",row.names = F)
#replication in lld fup
manhattan=read.delim("../result/qtl_927samples_snp_meta_inverse_direction_MR.txt",header = T,sep = "\t")
geno_lld=t(read.delim("../data/lld_snp_list_inverse.dosages_new.txt",row.names = 1,header = T,sep = "\t",check.names = F))
row.names(geno_lld)=gsub("1_LLDeep","LLDeep",row.names(geno_lld))
meta_dis=data.frame(t(read.delim("./mQTL/data_927samples_LLD_baseline_1183plasma_metabolites_adj_2g_all_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta_fup_adj=data.frame(t(read.delim("./mQTL/data_311samples_LLD_fup_replication_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
row.names(meta_fup_adj)=gsub("_F","",row.names(meta_fup_adj))
sp_base=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_base=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_base=read.delim("../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
microbiome_base=cbind(sp_base,path_base,vsv_base)
sp_fup=read.delim("../data/data_311samples_LLD_fup_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path_fup=read.delim("../data/data_311samples_LLD_fup_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
vsv_fup=read.delim("../data/data_311samples_LLD_fup_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
microbiome_fup=cbind(sp_fup,path_fup,vsv_fup)
row.names(microbiome_fup)=gsub("_F","",row.names(microbiome_fup))
#base discovery check meta
meta=meta_dis
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#927
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$ProbeName[i])],method = "spearman")
  manhattan$LLD_dis_P[i]=cor$p.value
  manhattan$LLD_dis_r[i]=cor$estimate
}
#fup rep meta
meta=meta_fup_adj
row.names(meta)=gsub("_F","",row.names(meta))
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#281
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$ProbeName[i])],method = "spearman")
  manhattan$LLD_fup_rep_P[i]=cor$p.value
  manhattan$LLD_fup_rep_r[i]=cor$estimate
}
manhattan$direction=ifelse(manhattan$LLD_dis_r*manhattan$IncludedDatasetsCorrelationCoefficient>0,1,-1)
manhattan$LLD_fup_rep_r=manhattan$LLD_fup_rep_r*manhattan$direction
#manhattan$replication=ifelse(manhattan$LLD_fup_rep_P<0.05&manhattan$LLD_fup_rep_r*manhattan$IncludedDatasetsCorrelationCoefficient>0,"yes","no")
manhattan$replication=ifelse(manhattan$LLD_fup_rep_P<0.1&manhattan$LLD_fup_rep_r*manhattan$IncludedDatasetsCorrelationCoefficient>0,"yes","no")
manhattan=manhattan[,c(23,2:4,9,10,5,1,18,20,26:29)]
manhattan=manhattan[which(manhattan$replication=="yes"),]#42
manhattan=merge(result_mr[,2:3],manhattan,by.x="outcome",by.y = "ProbeName",all = F)

#base discovery check microbiome
meta=microbiome_base
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#933
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$exposure[i])],method = "spearman")
  manhattan$LLD_dis_P_microbiome[i]=cor$p.value
  manhattan$LLD_dis_r_microbiome[i]=cor$estimate
}
#fup rep meta
meta=microbiome_fup
row.names(meta)=gsub("_F","",row.names(meta))
snp=geno_lld[intersect(row.names(meta),row.names(geno_lld)),]#281
meta=meta[row.names(snp),]
for(i in 1:nrow(manhattan)){
  cor=cor.test(snp[,as.character(manhattan$SNPName[i])],meta[,as.character(manhattan$exposure[i])],method = "spearman")
  manhattan$LLD_fup_rep_P_microbiome[i]=cor$p.value
  manhattan$LLD_fup_rep_r_microbiome[i]=cor$estimate
}
manhattan$LLD_dis_r_microbiome=manhattan$LLD_dis_r_microbiome*manhattan$direction
manhattan$LLD_fup_rep_r_microbiome=manhattan$LLD_fup_rep_r_microbiome*manhattan$direction
write.table(manhattan,file = "../result/mr_inverse_10mr_link_new.txt",quote = F,sep = "\t",row.names = F)
#read table
library(TwoSampleMR)
mr=read.delim("../result/mr_inverse_10mr_link_new.txt",header = T,sep = "\t")
exposure=unique(as.character(mr$outcome))
outcome=unique(as.character(mr$exposure))
result_mr=NULL
for(i in 1:length(exposure)){
  tmp=mr[which(mr$outcome==as.character(exposure[i])),]
  for(j in 1:length(unique(as.character(tmp$exposure)))){
    tmp_specific=tmp[which(tmp$exposure==unique(as.character(tmp$exposure))[j]),]
    #discovery
    data.exposure=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$IncludedDatasetsCorrelationCoefficient,se=rep(1/sqrt(933-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.exposure=format_data(data.exposure, type="exposure")
    data.outcome=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$LLD_dis_r_microbiome,se=rep(1/sqrt(933-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.outcome=format_data(data.outcome, type="outcome")
    data=harmonise_data(exposure_dat =data.exposure,outcome_dat = data.outcome)
    data$mr_keep='TRUE'
    data$mr_keep=as.logical(data$mr_keep)
    result.dis=mr(data)
    result.dis$outcome=unique(as.character(tmp$exposure))[j]
    result.dis$exposure=as.character(exposure)[i]
    colnames(result.dis)[7:9]=paste("discovery_",colnames(result.dis)[7:9],sep = "")
    #replication
    data.exposure=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$LLD_fup_rep_r,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.exposure=format_data(data.exposure, type="exposure")
    data.outcome=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$LLD_fup_rep_r_microbiome,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.outcome=format_data(data.outcome, type="outcome")
    data=harmonise_data(exposure_dat =data.exposure,outcome_dat = data.outcome)
    data$mr_keep='TRUE'
    data$mr_keep=as.logical(data$mr_keep)
    result.rep=mr(data)
    colnames(result.rep)[7:9]=paste("rep_",colnames(result.rep)[7:9],sep = "")
    result=cbind(result.dis[,c(4,3,5:9)],result.rep[,7:9])
    result_mr=rbind(result_mr,result)
  }
}
result_mr$name=paste(result_mr$outcome,result_mr$exposure,sep = "_with_")
result_mr=result_mr[which(result_mr$method=="Inverse variance weighted"),]
result_mr=result_mr[,4:11]
colnames(result_mr)[1:7]=paste("inverse_MR_",colnames(result_mr)[1:7],sep = "")
colnames(result_mr)[1:7]=paste(colnames(result_mr)[1:7],"_IVW",sep = "")
result_meta_to_micro=result_mr
result_mr=merge(result_micro_to_meta,result_meta_to_micro,by="name",all = F)
write.table(result_mr,file = "../result/result_mr_microbial_meta_bidirectional_new.txt",quote = F,sep = "\t",row.names = F)

##leave one out
result_mr=read.delim("../result/result_mr_microbiome_to_meta_new.txt",header = T,sep = "\t")
tmp=result_mr[which(result_mr$method=="Inverse variance weighted"),]
tmp$discovery_pval_BH_adj=p.adjust(tmp$discovery_pval,method = "BH")
tmp=tmp[which(tmp$discovery_pval_BH_adj<0.05),]
tmp=tmp[which(tmp$rep_pval<0.05),]
result_mr=read.delim("../result/result_mr_microbial_meta_bidirectional_new.txt",header = T,sep = "\t")
#result_mr=result_mr[which(result_mr$discovery_pval_WM<0.05&result_mr$inverse_MR_discovery_pval_IVW>0.05&result_mr$inverse_MR_rep_pval_IVW>0.05),]
result_mr=tmp
mr=read.delim("../result/mr_37microbial_abundance_45metabolites_new.txt",header = T,sep = "\t")
mr$meta_name=mr$name
mr$name=paste(mr$ProbeName,mr$meta,sep = "_with_")
mr=mr[which(mr$name%in%result_mr$name),]
write.table(mr,file = "../result/result_mr_microbiome_to_meta_new.txt",quote = F,sep = "\t",row.names = F)
library(TwoSampleMR)
library(gridExtra)
exposure=unique(as.character(mr$ProbeName))
outcome=unique(as.character(mr$meta))
plot=list()
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
    plot=c(plot,mr_leaveoneout_plot(mr_leaveoneout(data)))
    #replication
    data.exposure=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$LLD_fup_rep_r,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.exposure=format_data(data.exposure, type="exposure")
    data.outcome=data.frame(SNP=as.character(tmp_specific$SNPName),beta=tmp_specific$outcome_fup_r,se=rep(1/sqrt(311-3),nrow(tmp_specific)),effect_allele = as.character(tmp_specific$AlleleAssessed))
    data.outcome=format_data(data.outcome, type="outcome")
    data=harmonise_data(exposure_dat =data.exposure,outcome_dat = data.outcome)
    data$mr_keep='TRUE'
    data$mr_keep=as.logical(data$mr_keep)
    plot=c(plot,mr_leaveoneout_plot(mr_leaveoneout(data)))
  }
}
pdf("../figure/figure_leave_one_out_7_MR_results_new.pdf",width=5, height=10, useDingbats = F)
do.call(grid.arrange, c(plot,nrow=7,ncol=2))
dev.off()


##MR plots
library(ggplot2)
library(gridExtra)
library(stringr)
result_mr=read.delim("../result/result_mr_microbiome_to_meta_new.txt",header = T,sep = "\t")
tmp=result_mr[which(result_mr$method=="Inverse variance weighted"),]
tmp$discovery_pval_BH_adj=p.adjust(tmp$discovery_pval,method = "BH")
tmp=tmp[which(tmp$discovery_pval_BH_adj<0.05),]
tmp=tmp[which(tmp$rep_pval<0.05),]
tmp$name=paste(tmp$exposure,tmp$outcome,sep = "_with_")
mr=read.delim("../result/mr_37microbial_abundance_45metabolites_new.txt",header = T,sep = "\t")
mr$name=paste(mr$ProbeName,mr$meta,sep = "_with_")
mr=mr[which(mr$name%in%tmp$name),]
write.table(mr,file = "../result/result_mr_microbiome_to_meta_4links_all_new.txt",quote = F,sep = "\t",row.names = F)
mr$se=1/sqrt(927-3)
mr$se_fup=1/sqrt(311-3)
mr$exposure=mr$IncludedDatasetsCorrelationCoefficient
mr$outcome=mr$outcome_baseline_r
data=data.frame(microbe=rep(as.character(mr$ProbeName),2),meta=rep(as.character(mr$meta),2),exposure=c(mr$IncludedDatasetsCorrelationCoefficient,mr$LLD_fup_rep_r),outcome=c(mr$outcome_baseline_r,mr$outcome_fup_r),se=c(mr$se,mr$se_fup),type=c(rep("dis",nrow(mr)),rep("rep",nrow(mr))))
data$microbe_abb[grep("s__",data$microbe)]=str_split_fixed(data$microbe[grep("s__",data$microbe)],"s__",2)[,2]
data$microbe_abb[grep(": ",data$microbe)]=str_split_fixed(data$microbe[grep(": ",data$microbe)],": ",2)[,1]
data$microbe=data$microbe_abb
exposure=unique(as.character(data$microbe))
outcome=unique(as.character(data$meta))
tmp1=data[which(as.character(data$microbe)==exposure[1]&as.character(data$meta)==outcome[1]),]
p1=ggplot(tmp1,aes(exposure,outcome,group=type,color=type))+geom_point(size=1,shape=21,fill="gray90")+geom_smooth(method = "lm",size=0.5,linetype=2,se=F)+geom_errorbar(aes(ymin = outcome-se,ymax = outcome+se),colour="gray60",alpha=1,width=0,size=0.1)+geom_errorbarh(aes(xmin = exposure-se,xmax = exposure+se),colour="gray60",alpha=1,height=0,size=0.1)+theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x=exposure[1], y=outcome[1])+scale_fill_manual(values=c("#34669a","#9ac43a"))+scale_color_manual(values=c("#34669a","#9ac43a"))
tmp2=data[which(as.character(data$microbe)==exposure[2]&as.character(data$meta)==outcome[2]),]
p2=ggplot(tmp2,aes(exposure,outcome,group=type,color=type))+geom_point(size=1,shape=21,fill="gray90")+geom_smooth(method = "lm",size=0.5,linetype=2,se=F)+geom_errorbar(aes(ymin = outcome-se,ymax = outcome+se),colour="gray60",alpha=1,width=0,size=0.1)+geom_errorbarh(aes(xmin = exposure-se,xmax = exposure+se),colour="gray60",alpha=1,height=0,size=0.1)+theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x=exposure[2], y=outcome[2])+scale_fill_manual(values=c("#34669a","#9ac43a"))+scale_color_manual(values=c("#34669a","#9ac43a"))
tmp3=data[which(as.character(data$microbe)==exposure[3]&as.character(data$meta)==outcome[3]),]
p3=ggplot(tmp3,aes(exposure,outcome,group=type,color=type))+geom_point(size=1,shape=21,fill="gray90")+geom_smooth(method = "lm",size=0.5,linetype=2,se=F)+geom_errorbar(aes(ymin = outcome-se,ymax = outcome+se),colour="gray60",alpha=1,width=0,size=0.1)+geom_errorbarh(aes(xmin = exposure-se,xmax = exposure+se),colour="gray60",alpha=1,height=0,size=0.1)+theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x=exposure[3], y=outcome[3])+scale_fill_manual(values=c("#34669a","#9ac43a"))+scale_color_manual(values=c("#34669a","#9ac43a"))
tmp4=data[which(as.character(data$microbe)==exposure[4]&as.character(data$meta)==outcome[4]),]
p4=ggplot(tmp4,aes(exposure,outcome,group=type,color=type))+geom_point(size=1,shape=21,fill="gray90")+geom_smooth(method = "lm",size=0.5,linetype=2,se=F)+geom_errorbar(aes(ymin = outcome-se,ymax = outcome+se),colour="gray60",alpha=1,width=0,size=0.1)+geom_errorbarh(aes(xmin = exposure-se,xmax = exposure+se),colour="gray60",alpha=1,height=0,size=0.1)+theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x=exposure[4], y=outcome[4])+scale_fill_manual(values=c("#34669a","#9ac43a"))+scale_color_manual(values=c("#34669a","#9ac43a"))
p5=ggplot(tmp4,aes(exposure,outcome,group=type,color=type))+geom_point(size=1,shape=21,fill="gray90")+geom_smooth(method = "lm",size=0.5,linetype=2,se=F)+geom_errorbar(aes(ymin = outcome-se,ymax = outcome+se),colour="gray60",alpha=1,width=0,size=0.1)+geom_errorbarh(aes(xmin = exposure-se,xmax = exposure+se),colour="gray60",alpha=1,height=0,size=0.1)+theme_bw()+theme(axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+labs(x=exposure[4], y=outcome[4])+scale_fill_manual(values=c("#34669a","#9ac43a"))+scale_color_manual(values=c("#34669a","#9ac43a"))
pdf("../figure/figure_MR_micorbiome_meta_crossbars_new.pdf",useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,ncol=4,nrow = 4)
dev.off()
##MR combined new
mr=read.delim("../result/result_mr_microbiome_to_meta_4links_new.txt",header = T,sep = "\t")
result_mr=read.delim("../result/result_mr_microbiome_to_meta_new.txt",header = T,sep = "\t")
result_mr=result_mr[which(result_mr$name%in%mr$name),]
result_mr=result_mr[-which(result_mr$method=="Simple mode"),]
result_mr=result_mr[-which(result_mr$method=="Weighted mode"),]
mr=result_mr
library(meta)
library(ggplot2)
library(gridExtra)
mr$beta_combin_fix=NA
mr$se_combin_fix=NA
mr$P_combin_fix=NA
positions=c("2_replication","1_discovery")
plot=list()
meta=unique(mr$outcome)
for(i in 1:length(meta)){
  tmp=mr[which(mr$outcome==meta[i]),]
  data=data.frame(method=rep(c("Egger","WM","IVW"),2),cohort=c(rep("1_discovery",3),rep("2_replication",3)),beta=as.numeric(tmp[,5]),up=as.numeric(tmp[,5]+tmp[,6]),low=as.numeric(tmp[,5]-tmp[,6]))
  #data=data[which(data$method=="IVW"),]
  #data=data[which(data$method=="WM"),]
  data=data[which(data$method=="Egger"),]
  plot[[i]]=ggplot(data, aes(x=cohort, y=beta, ymin=low, ymax=up, colour=cohort,shape=method))+geom_pointrange(size=0.3, position=position_dodge(width=c(0.5))) + coord_flip()+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),text = element_text(size=5))+scale_color_manual(values=c("#34669a","#9ac43a"))+theme(panel.grid=element_blank())+labs(y="Effect size (beta)", x=meta[i])+geom_hline(yintercept = 0,linetype=3,colour="grey10",alpha=0)+theme(axis.title=element_text(size=6))+scale_x_discrete(limits = positions)
}
pdf("../figure/figure_MR_micorbiome_meta_combined_new_Egger.pdf",useDingbats = F)
do.call(grid.arrange, c(plot,nrow=6,ncol=4))
dev.off()
##MR combined old
library(meta)
library(ggplot2)
library(gridExtra)
mr=read.delim("../result/result_mr_microbial_meta_bidirectional.txt",header = T,sep = "\t")
mr=mr[which(mr$inverse_MR_discovery_pval_IVW>0.05&mr$inverse_MR_rep_pval_IVW>0.05),]
mr$beta_combin_fix=NA
mr$se_combin_fix=NA
mr$P_combin_fix=NA
positions=c("3_combine","2_replication","1_discovery")
plot=list()
for(i in 1:nrow(mr)){
  tmp=metagen(as.numeric(mr[i,c(5,12)]),as.numeric(mr[i,c(6,13)]))
  mr$beta_combin_fix[i]=tmp$TE.fixed
  mr$se_combin_fix[i]=tmp$seTE.fixed
  mr$P_combin_fix[i]=tmp$pval.fixed
  data=data.frame(cohort=c("1_discovery","2_replication","3_combine"),beta=as.numeric(mr[i,c(5,12,22)]),up=as.numeric(mr[i,c(5,12,22)]+mr[i,c(6,13,23)]),low=as.numeric(mr[i,c(5,12,22)]-mr[i,c(6,13,23)]))
  plot[[i]]=ggplot(data, aes(x=cohort, y=beta, ymin=low, ymax=up, colour=cohort))+geom_pointrange(shape=15, size=0.3, position=position_dodge(width=c(0.1))) + coord_flip()+theme_bw()+theme(axis.line = element_line(colour = "black"),panel.border = element_blank(),legend.position="none",axis.text.y=element_blank(),axis.ticks.y=element_blank(),text = element_text(size=5))+scale_color_manual(values=c("blue4","brown3","gray70"))+theme(panel.grid=element_blank())+labs(y="Effect size (beta IVW)", x=mr$outcome[i])+geom_hline(yintercept = 0,linetype=3,colour="grey10",alpha=0)+theme(axis.title=element_text(size=6))+scale_x_discrete(limits = positions)
}
pdf("../figure/figure_MR_micorbiome_meta_combined_new.pdf",useDingbats = F)
do.call(grid.arrange, c(plot,nrow=7,ncol=4))
dev.off()

###medication
annotation=read.delim("./key_lld_1183meta_annotation.txt",header = T,quote = "",sep = "\t",check.names = F)
annotation$name=gsub("\"","",annotation$name)
annotation$super_class=gsub("\"","",annotation$super_class)
annotation$class=gsub("\"","",annotation$class)
annotation$sub_class=gsub("\"","",annotation$sub_class)
annotation$direct_parent=gsub("\"","",annotation$direct_parent)
diet=read.delim("~/Documents/UMCG/research-projects/microbiome_data/meta_data/LLD_207pheno_log_imp_1135samples.txt",header = T,sep = "\t")
diet=diet[,48:125]
vsv=read.delim("../data/data_1054samples_LLD_baseline_1782vsv_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
diet=diet[row.names(vsv),]
dsv=read.delim("~/Documents/UMCG/research-projects/microbiome_data/mgs_SVs/LLD_base_fup_SVs/01.cleanData/20200615_LLD12_deletionStructuralVariation_1471samples.tsv",header = T,sep = "\t",check.names = F)
dsv=dsv[,colSums(is.na(dsv)==F)>nrow(dsv)*0.1]
dsv=dsv[which(row.names(dsv)%in%row.names(vsv)),]
dsv=merge(vsv[,1:2],dsv,by="row.names",all = T)
row.names(dsv)=dsv[,1]
dsv=dsv[,-c(1:3)]
dsv=dsv[row.names(vsv),]
meta=data.frame(t(read.delim("./mQTL/data_1054samples_LLD_baseline_1183plasma_metabolites_adj_age_sex_smk_oral_qtl.txt",header = T,row.names = 1,sep = "\t")),check.names = F)
meta=meta[row.names(vsv),]
sp=read.delim("../data/data_1054samples_LLD_baseline_156sp_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
path=read.delim("../data/data_1054samples_LLD_baseline_343path_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
bgc=read.delim("../data/data_1054samples_LLD_baseline_3075bgc_adj_age_sex_smk_ppi_anti_lax.txt",header = T,sep = "\t",check.names = F)
#data=cbind(diet,meta,vsv,dsv)
asso_vsv=read.delim("../result/table_cor_meta_vSV_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F)
colnames(asso_vsv)[3:9]=paste(colnames(asso_vsv)[3:9],"_vsv",sep = "")
asso_dsv=read.delim("../result/table_cor_meta_dSV_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F)
colnames(asso_dsv)[3:11]=paste(colnames(asso_dsv)[3:11],"_dsv",sep = "")
asso_sp=read.delim("../result/table_cor_meta_sp_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(asso_sp)[3:9]=paste(colnames(asso_sp)[3:9],"_sp",sep = "")
asso_path=read.delim("../result/table_cor_meta_path_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(asso_path)[3:9]=paste(colnames(asso_path)[3:9],"_path",sep = "")
asso_bgc=read.delim("../result/table_cor_meta_bgc_dis_rep_delta_associations.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(asso_bgc)[3:9]=paste(colnames(asso_bgc)[3:9],"_bgc",sep = "")
asso_diet=read.delim("../result/table_cor_meta_diet_associations_FDR.txt",header = T,sep = "\t",check.names = F,quote = "")
colnames(asso_diet)[9:11]=paste(colnames(asso_diet)[9:11],"_diet",sep = "")
asso_vsv=merge(asso_diet[,c(1,8:11)],asso_vsv[,1:9],by="meta",all = F)
asso_dsv=merge(asso_diet[,c(1,8:11)],asso_dsv[,1:11],by="meta",all = F)
asso_sp=merge(asso_diet[,c(1,8:11)],asso_sp[,1:9],by="meta",all = F)
asso_path=merge(asso_diet[,c(1,8:11)],asso_path[,1:9],by="meta",all = F)
asso_bgc=merge(asso_diet[,c(1,8:11)],asso_bgc[,1:9],by="meta",all = F)
library(mediation)
library(ppcor)
#vsv
result=asso_vsv
for(i in 1:nrow(result)){
  cor=cor.test(diet[,as.character(result$diet[i])],vsv[,as.character(result$vSVs[i])],method = "spearman")
  result$Cor_diet_sv[i]=cor$estimate
  result$Pval_diet_sv[i]=cor$p.value
}
result$Qval_diet_sv=p.adjust(result$Pval_diet_sv,method = "BH")
result=result[which(result$Qval_diet_sv<0.05),]
result=result[which(result$Pval_delta_vsv<0.05),]#87
for(i in 1:nrow(result)){
  data=data.frame(cbind(diet[,as.character(result$diet[i])],meta[,as.character(result$meta[i])],vsv[,as.character(result$vSVs[i])]))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate[i]=summary$d.avg.p
  result$Pval_direct[i]=summary$z.avg.p
  #inverse mediate
  colnames(data)=c("X","M","Y")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}
result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
result=result[which(result$Qval_mediate<0.05&result$Pval_mediate_inverse>0.05|result$Pval_mediate>0.05&result$Qval_mediate_inverse<0.05),]
result=merge(annotation[,c(1,12:19)],result,by="meta",all = F)
write.table(result,file = "../result/medication_diet_impact_meta_through_vSV_FDR.txt",quote = F,sep = "\t",row.names = F)
#dsv
result=asso_dsv
for(i in 1:nrow(result)){
  cor=cor.test(diet[,as.character(result$diet[i])],dsv[,as.character(result$dSVs[i])],method = "spearman")
  result$Cor_diet_sv[i]=cor$estimate
  result$Pval_diet_sv[i]=cor$p.value
}
result$Qval_diet_sv=p.adjust(result$Pval_diet_sv,method = "BH")
result=result[which(result$Qval_diet_sv<0.05),]
result=result[which(result$Pval_delta_change_dsv<0.05),]#3
for(i in 1:nrow(result)){
  data=data.frame(cbind(diet[,as.character(result$diet[i])],meta[,as.character(result$meta[i])],dsv[,as.character(result$dSVs[i])]))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M")
  model.m=glm(M~X,family=binomial(link='logit'),data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate[i]=summary$d.avg.p
  result$Pval_direct[i]=summary$z.avg.p
  #inverse mediate
  colnames(data)=c("X","M","Y")
  model.m=lm(M~X,data)
  model.y=glm(Y~X+M,family=binomial(link='logit'),data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}
result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
result=result[which(result$Qval_mediate<0.05&result$Pval_mediate_inverse>0.05|result$Pval_mediate>0.05&result$Qval_mediate_inverse<0.05),]
result=merge(annotation[,c(1,12:19)],result,by="meta",all = F)
write.table(result,file = "../result/medication_diet_impact_meta_through_dSV_FDR.txt",quote = F,sep = "\t",row.names = F)
#sp
result=asso_sp
for(i in 1:nrow(result)){
  cor=cor.test(diet[,as.character(result$diet[i])],sp[,as.character(result$species[i])],method = "spearman")
  result$Cor_diet_sp[i]=cor$estimate
  result$Pval_diet_sp[i]=cor$p.value
}
result$Qval_diet_sp=p.adjust(result$Pval_diet_sp,method = "BH")
result=result[which(result$Qval_diet_sp<0.05),]
result=result[which(result$Pval_delta_sp<0.05),]#123
for(i in 1:nrow(result)){
  data=data.frame(cbind(diet[,as.character(result$diet[i])],meta[,as.character(result$meta[i])],sp[,as.character(result$species[i])]))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate[i]=summary$d.avg.p
  result$Pval_direct[i]=summary$z.avg.p
  #inverse mediate
  colnames(data)=c("X","M","Y")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}
result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
result=result[which(result$Qval_mediate<0.05&result$Pval_mediate_inverse>0.05|result$Pval_mediate>0.05&result$Qval_mediate_inverse<0.05),]
result=merge(annotation[,c(1,12:19)],result,by="meta",all = F)
write.table(result,file = "../result/medication_diet_impact_meta_through_sp_FDR.txt",quote = F,sep = "\t",row.names = F)
#path
result=asso_path
for(i in 1:nrow(result)){
  cor=cor.test(diet[,as.character(result$diet[i])],path[,as.character(result$pathway[i])],method = "spearman")
  result$Cor_diet_path[i]=cor$estimate
  result$Pval_diet_path[i]=cor$p.value
}
result$Qval_diet_path=p.adjust(result$Pval_diet_path,method = "BH")
result=result[which(result$Qval_diet_path<0.05),]
result=result[which(result$Pval_delta_path<0.05),]#16
for(i in 1:nrow(result)){
  data=data.frame(cbind(diet[,as.character(result$diet[i])],meta[,as.character(result$meta[i])],path[,as.character(result$pathway[i])]))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate[i]=summary$d.avg.p
  result$Pval_direct[i]=summary$z.avg.p
  #inverse mediate
  colnames(data)=c("X","M","Y")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}
result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
result=result[which(result$Qval_mediate<0.05&result$Pval_mediate_inverse>0.05|result$Pval_mediate>0.05&result$Qval_mediate_inverse<0.05),]
result=merge(annotation[,c(1,12:19)],result,by="meta",all = F)
write.table(result,file = "../result/medication_diet_impact_meta_through_path_FDR.txt",quote = F,sep = "\t",row.names = F)
#bgc
result=asso_bgc
for(i in 1:nrow(result)){
  cor=cor.test(diet[,as.character(result$diet[i])],bgc[,as.character(result$bgc[i])],method = "spearman")
  result$Cor_diet_bgc[i]=cor$estimate
  result$Pval_diet_bgc[i]=cor$p.value
}
result$Qval_diet_bgc=p.adjust(result$Pval_diet_bgc,method = "BH")
result=result[which(result$Qval_diet_bgc<0.05),]
result=result[which(result$Pval_delta_bgc<0.05),]#3170
for(i in 1:nrow(result)){
  data=data.frame(cbind(diet[,as.character(result$diet[i])],meta[,as.character(result$meta[i])],bgc[,as.character(result$bgc[i])]))
  data=na.omit(data)
  #diet influence meta through microbiome
  colnames(data)=c("X","Y","M")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate[i]=summary$d.avg.p
  result$Pval_direct[i]=summary$z.avg.p
  #inverse mediate
  colnames(data)=c("X","M","Y")
  model.m=lm(M~X,data)
  model.y=lm(Y~X+M,data)
  summary=summary(mediate(model.m ,model.y,treat = "X", mediator = "M",boot = T,sims = 1000))
  result$Pval_mediate_inverse[i]=summary$d.avg.p
  result$Pval_direct_inverse[i]=summary$z.avg.p
}
result$Qval_mediate=p.adjust(result$Pval_mediate,method = "BH")
result$Qval_mediate_inverse=p.adjust(result$Pval_mediate_inverse,method = "BH")
result=result[which(result$Qval_mediate<0.05&result$Pval_mediate_inverse>0.05|result$Pval_mediate>0.05&result$Qval_mediate_inverse<0.05),]
result=merge(annotation[,c(1,12:19)],result,by="meta",all = F)
write.table(result,file = "../result/medication_diet_impact_meta_through_bgc_FDR.txt",quote = F,sep = "\t",row.names = F)


##
library(stringr)
annotation=read.delim("../data/key_lld_1183meta_annotation_type.txt",header = T,sep = "\t",quote = "",check.names = F)
vsv=read.delim("../result/medication_diet_impact_meta_through_vSV_FDR.txt",header = T,sep = "\t")
dsv=read.delim("../result/medication_diet_impact_meta_through_dSV_FDR.txt",header = T,sep = "\t")
dsv=dsv[,-c(20:21)]
colnames(dsv)=colnames(vsv)
sp=read.delim("../result/medication_diet_impact_meta_through_sp_FDR.txt",header = T,sep = "\t")
colnames(sp)=colnames(vsv)
path=read.delim("../result/medication_diet_impact_meta_through_path_FDR.txt",header = T,sep = "\t")
colnames(path)=colnames(vsv)
bgc=read.delim("../result/medication_diet_impact_meta_through_bgc_FDR.txt",header = T,sep = "\t")
bgc$bgc=str_split_fixed(bgc$bgc,"--Entryname=",2)[,2]
colnames(bgc)=colnames(vsv)
media=data.frame(diet=c(as.character(sp$diet),as.character(bgc$diet),as.character(vsv$diet),as.character(dsv$diet)),microbe=c(as.character(sp$vSVs),as.character(bgc$vSVs),as.character(vsv$vSVs),as.character(dsv$vSVs)),meta=c(as.character(sp$meta),as.character(bgc$meta),as.character(vsv$meta),as.character(dsv$meta)))
media=cbind(media,rbind(sp[,c(11:13,15:30)],bgc[,c(11:13,15:30)],vsv[,c(11:13,15:30)],dsv[,c(11:13,15:30)]))
media=merge(media,annotation[,c(1,12:20)],by="meta",all = F)
#coffee
tmp=media[which(media$Qval_mediate>0.05),]
tmp=tmp[grep("coffee",tmp$diet),]
tmp=tmp[which(tmp$Cor_baseline_diet>0.2&tmp$Cor_baseline_vsv>0.1&tmp$Cor_diet_sv>0.1),]
#alcohol beer
tmp=media[which(media$Qval_mediate>0.05),]
tmp=tmp[which(tmp$Cor_baseline_vsv*tmp$Cor_diet_sv>0),]
#Pipecolic acid
tmp=media[which(media$Qval_mediate<0.05),]
tmp=tmp[which(tmp$Cor_baseline_vsv*tmp$Cor_diet_sv>0),]


#######mediation plot
library(stringr)
annotation=read.delim("../data/key_lld_1183meta_annotation_type.txt",header = T,sep = "\t",quote = "",check.names = F)
vsv=read.delim("../result/medication_diet_impact_meta_through_vSV_FDR.txt",header = T,sep = "\t")
dsv=read.delim("../result/medication_diet_impact_meta_through_dSV_FDR.txt",header = T,sep = "\t")
dsv=dsv[,-c(20:21)]
sp=read.delim("../result/medication_diet_impact_meta_through_sp_FDR.txt",header = T,sep = "\t")
path=read.delim("../result/medication_diet_impact_meta_through_path_FDR.txt",header = T,sep = "\t")
bgc=read.delim("../result/medication_diet_impact_meta_through_bgc_FDR.txt",header = T,sep = "\t")
bgc$bgc=str_split_fixed(bgc$bgc,"--Entryname=",2)[,2]
bgc$bgc=str_split_fixed(bgc$bgc,"--OS=",2)[,1]
media=data.frame(diet=c(as.character(sp$diet),as.character(bgc$diet),as.character(vsv$diet),as.character(dsv$diet)),microbe=c(as.character(sp$species),as.character(bgc$bgc),as.character(vsv$vSVs),as.character(dsv$dSVs)),meta=c(as.character(sp$name),as.character(bgc$name),as.character(vsv$name),as.character(dsv$name)))
media=cbind(media,rbind(sp[,25:30],bgc[,25:30],vsv[,25:30],dsv[,25:30]))
#sankey plot
media$microbe_abb=NA
media$microbe_abb[grep("s__",media$microbe)] =str_split_fixed(media$microbe[grep("s__",media$microbe)],"s__",2)[,2]
media$microbe_abb[grep("s__",media$microbe,invert = T)] =as.character(media$microbe)[grep("s__",media$microbe,invert = T)]
library(ggplot2)
library(ggalluvial)
library(gridExtra)
#diet-microbe-meta
net=media[which(media$Qval_mediate<0.05&media$Pval_mediate_inverse>0.05),c(1,3,10)]
net$fre=1
p=ggplot(net,aes(axis1 = net$diet, axis2 = net$microbe_abb, axis3 = net$meta,y= net$fre))+
  scale_x_discrete(limits = c("Diet", "Microbiome", "Meta")) +
  geom_alluvium(aes(fill = net$meta),cex=0.1)+
  geom_stratum() +theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+geom_text(stat = "stratum",cex=0.8,aes(label = after_stat(stratum))) 
pdf("../figure/fig_mediation_10linkages_sankey.pdf",useDingbats = F)
grid.arrange(p,ncol=2,nrow = 11)
dev.off()
#diet-meta-microbe
net=media[which(media$Qval_mediate>0.05&media$Qval_mediate_inverse<0.05),c(1,3,10)]
net$fre=1
p=ggplot(net,aes(axis1 = net$diet, axis2 = net$meta, axis3 = net$microbe_abb,y= net$fre))+
  scale_x_discrete(limits = c("Diet", "Meta","Microbiome")) +
  geom_alluvium(aes(fill = net$meta),cex=0.1)+
  geom_stratum() +theme_bw()+theme(legend.position="none",axis.text = element_text(size = 5),axis.title = element_text(size = 6),panel.grid=element_blank())+geom_text(stat = "stratum",cex=0.8,aes(label = after_stat(stratum))) 
pdf("../figure/fig_mediation_185linkages_sankey.pdf",useDingbats = F)
grid.arrange(p,ncol=2,nrow = 1)
dev.off()
