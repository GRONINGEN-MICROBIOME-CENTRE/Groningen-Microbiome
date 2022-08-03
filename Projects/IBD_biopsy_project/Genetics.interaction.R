#
library(data.table)
library(ggplot2)
library(foreach)
library(lme4)
library(nlme)
library(factoextra)
library(ggsci)
library(gridExtra)
library(vegan)
library(RColorBrewer)
library(ggalluvial)
require(cowplot)
library(ggpubr)
library("dplyr")
library(edgeR)
library(lme4)
library(limma)
library(mixOmics)
library(robCompositions)
library(microbiome)
source("Microbiome.function.R")
library(lmerTest)
library(jcolors)
library(crayon)

# =========================================================================================================
# files prepare
# =========================================================================================================

# annotate eQTLs
eQTLgene=read.table("eQTLs/eQTL_pairs.txt",header = T,sep = "\t",stringsAsFactors = F)
annotation=read.table("annotation.file.txt",sep = "\t",header = T,stringsAsFactors = F)
eQTLgene=merge(eQTLgene,annotation,by.x = "ExpressionGene",by.y = "Gene",all=F)

# import gene and bacteria data
covariate_bac=read.table("Covariate.bacteria.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
remove=rownames(covariate_bac)[covariate_bac$Diagnosis=="UC" & covariate_bac$Location_rough=="ileum" & covariate_bac$Inflammation=="Yes"]
genus=genus[!rownames(genus) %in% remove,]
covariate_bac=covariate_bac[rownames(covariate_bac) %in% rownames(genus),]
genus=genus[rownames(genus) %in% rownames(covariate_bac),]
bacteria=read.table("OutputTable/CLR.bacteria.txt",sep = "\t",row.names = 1,header = T,check.names = F,stringsAsFactors = F)
bacteria=bacteria[rownames(bacteria) %in% rownames(covariate_bac),]

gene=read.table(file = "OutputTable/Genes.basic.Nocorrection.protein.coding.txt",sep = "\t",row.names = 1,header = T,check.names = F)
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T)
samples_count=intersect(rownames(gene),(rownames(genus)))

gene=gene[rownames(gene) %in% samples_count,]
genus=genus[rownames(genus) %in% samples_count,]
covariate_rna=covariate_rna[rownames(covariate_rna) %in% samples_count,]

# RNA data correction (age, gender, BMI, bacth, inflammation, location (no diagnosis, location == diagnosis), and medication)
covariate_rna=covariate_rna[,c("Inflammation","Location_rough","age_at_biopsy","sex","BMI","Batch","Aminosalicylates","Thiopurines","Steroids","Biological_use")]
covariate_rna=covariate_rna[order(rownames(covariate_rna)),]
gene=gene[order(rownames(gene)),]
genes_correct = apply(gene,2,function(x){
  x.subset = x[!is.na(x)]
  covariate.subset = covariate_rna[!is.na(x),,drop = FALSE]
  subset.data=cbind(x.subset,covariate.subset)
  colnames(subset.data)[1]="tmp.gene"
  
  fit=lm(tmp.gene~age_at_biopsy+sex+BMI+Batch+Inflammation+Location_rough+Aminosalicylates+Thiopurines+Steroids,data=subset.data)
  x.resid = resid(fit)
  x[!is.na(x)] = x.resid
  x[is.na(x)] = 0
  return(x)
})
genes_correct=as.data.frame(genes_correct)

# select inflammed genes + NOD2 gene + ABO + MCM6(LCT) gene
inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation_genes=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
inflammation_genes=append(inflammation_genes,c("NOD2","ABO","MCM6","FUT2"))
genes_correct=genes_correct[,colnames(genes_correct) %in% inflammation_genes]

# final target genes + NOD2 gene + ABO + MCM6(LCT)
target_gene=intersect(colnames(genes_correct),eQTLgene$id)
target_gene=eQTLgene[eQTLgene$id %in% target_gene,]
target_gene=target_gene[,c("rsID","id")]
target_gene=rbind(target_gene,c("rs2066844","NOD2"))
target_gene=rbind(target_gene,c("rs2066845","NOD2"))
target_gene=rbind(target_gene,c("rs2066847","NOD2"))
target_gene=rbind(target_gene,c("rs182549","MCM6"))
target_gene=rbind(target_gene,c("rs8176645","ABO"))
# target_gene=rbind(target_gene,c("rs601338","FUT2")) # already included in eQTL

write.table(target_gene,"eQTLs/FinalPairs.txt",sep = "\t",quote = F,row.names = F)

# =========================================================================================================
# start analysis
# =========================================================================================================
genotype.coupling=read.table("eQTLs/genetics.coupling.file.txt",header = T,stringsAsFactors = F,sep = "\t")
genotype=read.table("eQTLs/dosage.dosages.txt",row.names = 1,sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genotype=as.data.frame(t(genotype))
genotype=genotype[rownames(genotype) %in% genotype.coupling$WES_ID,]
genotype.coupling=genotype.coupling[genotype.coupling$WES_ID %in% rownames(genotype),]
genotype=genotype[order(rownames(genotype)),]
genotype.coupling=genotype.coupling[order(genotype.coupling$WES_ID),]
rownames(genotype)=genotype.coupling$UMCGresearchID
target_gene=target_gene[target_gene$rsID %in% colnames(genotype),]
genes_correct=genes_correct[,colnames(genes_correct) %in% target_gene$id]
target_gene=target_gene[target_gene$id %in% colnames(genes_correct),]
covariate_rna=read.table("Covariate.rna.organized.txt",sep = "\t",stringsAsFactors = F,header = T)

# =========================== mbQTL ===========================
basic_factor=covariate_rna[,c("ResearchID","Inflammation","Location_rough","age_at_biopsy","sex","BMI","Biological_use","resec_ileocec","resec_part_small"),drop=F]
basic_factor$biopsy_ID=rownames(basic_factor)

mm=0
modelfit0 = foreach(i=1:ncol(genotype),.combine = rbind) %do%  {
  tmp.snp=colnames(genotype)[i]
  tmp.data=genotype[,tmp.snp,drop=F]
  tmp.data=merge(tmp.data,basic_factor,by.x = "row.names",by.y = "ResearchID",all = F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 5))
  colnames(tmp.result)=c("SNP","Taxa","beta","se","pvalue")
  
  mm=mm+1
  cat(yellow("===",mm,"===",tmp.snp,"\n"))
  for(n in 1:ncol(bacteria)){
    tmp.bacteria=colnames(bacteria)[n]
    tmp.test=merge(tmp.data,bacteria[,tmp.bacteria,drop=F],all=F,by.x ="biopsy_ID",by.y="row.names")
    colnames(tmp.test)=c("BiopsyID","ResearchID","SNP","inflammation","location","age","sex","BMI","biological","ileocec","small","bacteria")
    tmp.test=tmp.test[tmp.test$SNP!=-1,]
    tmp.test=na.omit(tmp.test)
    tmp.model=lm(bacteria ~ inflammation + location + age +sex + BMI + +biological+ileocec+small+SNP,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=tmp.bacteria
    tmp.result$SNP[n]=tmp.snp
    tmp.result$beta[n]=tmp.model$Estimate[rownames(tmp.model)=="SNP"]
    tmp.result$se[n]=tmp.model$`Std. Error`[rownames(tmp.model)=="SNP"]
    tmp.result$pvalue[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="SNP"]

    
  }
  return.string=as.data.frame(tmp.result)
}
modelfit0$FDR=p.adjust(modelfit0$pvalue,method = "BH")
modelfit0=merge(modelfit0,target_gene,by.x = "SNP",by.y="rsID",all=F)
write.table(modelfit0,file = "OutputTable/mbQTL.target.txt",sep = "\t",quote = F,row.names = F)

tmp.snp="rs2066847"
tmp.bacteria="Bacteroides"
tmp.data=genotype[,tmp.snp,drop=F]
tmp.data=merge(tmp.data,basic_factor,by.x = "row.names",by.y = "ResearchID",all = F)
tmp.test=merge(tmp.data,bacteria[,tmp.bacteria,drop=F],all=F,by.x ="biopsy_ID",by.y="row.names")
colnames(tmp.test)=c("BiopsyID","ResearchID","SNP","inflammation","location","age","sex","BMI","biological","ileocec","small","bacteria")
tmp.test=tmp.test[tmp.test$SNP!=-1,]
tmp.test$SNP=as.factor(tmp.test$SNP)
tmp.model=lm(bacteria ~ inflammation + location + age +sex + BMI + +biological+ileocec+small,data = tmp.test)
tmp.test$residual=resid(tmp.model)

ggplot(tmp.test, aes(x=SNP, y=bacteria,fill=SNP)) + 
  geom_boxplot(position=position_dodge(2),outlier.shape = NA)+
  theme_bw()+
  geom_jitter(width = 0.2,size=1,shape=21)+
  scale_fill_manual(values = c("#117733", "#88CCEE", "#BB5566"))+
  guides(color=FALSE)+guides(fill=FALSE)+
  theme(strip.background =element_rect(fill="white"),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+xlab("")
ggsave("OutputPlot/test.test.pdf",width = 5,height = 5)


# =========================== interaction (eQTL) ===========================
mm=0
modelfit = foreach(i=1:ncol(genes_correct),.combine = rbind) %do%  {
  tmp.gene=colnames(genes_correct)[i]
  tmp.snp=target_gene$rsID[target_gene$id==tmp.gene]
  tmp.data=merge(genes_correct[,tmp.gene,drop=F],covariate_rna[,"ResearchID",drop=F],by="row.names",all=F)
  tmp.data=merge(tmp.data,genotype[,tmp.snp,drop=F],by.x="ResearchID",by.y="row.names",all=F)
  
  tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 8))
  colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_snp","Beta_interact","P_taxa","P_snp","P_interact")
  
  mm=mm+1
  cat(yellow("===",mm,"===",tmp.gene,"\n"))
  for(n in 1:ncol(bacteria)){
    tmp.bacterria=colnames(bacteria)[n]
    tmp.test=merge(tmp.data,bacteria[,tmp.bacterria,drop=F],all=F,by.x ="Row.names",by.y="row.names")
    colnames(tmp.test)=c("BiopsyID","ResearchID","Gene","SNP","Bac")
    tmp.test=tmp.test[tmp.test$SNP!=-1,]
    tmp.model=lm(Gene ~ Bac + SNP + Bac * SNP,data = tmp.test)
    tmp.model=as.data.frame(summary(tmp.model)$coefficients)
    
    tmp.result$Taxa[n]=tmp.bacterria
    tmp.result$Gene[n]=tmp.gene
    tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
    tmp.result$Beta_snp[n]=tmp.model$Estimate[rownames(tmp.model)=="SNP"]
    tmp.result$Beta_interact[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac:SNP"]
    tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
    tmp.result$P_snp[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="SNP"]
    tmp.result$P_interact[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac:SNP"]
    
  }
  return.string=as.data.frame(tmp.result)
}
modelfit$FDR=p.adjust(modelfit$P_interact,method = "BH")

# test NOD2
tmp.gene="NOD2"
tmp.genotype=genotype[,colnames(genotype) %in% target_gene$rsID[target_gene$id=="NOD2"],drop=F]
tmp.genotype=tmp.genotype[tmp.genotype$rs2066844!=-1,]
tmp.genotype=tmp.genotype[tmp.genotype$rs2066845!=-1,]
tmp.genotype=tmp.genotype[tmp.genotype$rs2066847!=-1,]
tmp.genotype$defect=rowSums(tmp.genotype)
tmp.genotype$defect[tmp.genotype$defect>0]=1
tmp.snp="defect"
tmp.data=merge(genes_correct[,tmp.gene,drop=F],covariate_rna[,"ResearchID",drop=F],by="row.names",all=F)
tmp.data=merge(tmp.data,tmp.genotype[,tmp.snp,drop=F],by.x="ResearchID",by.y="row.names",all=F)

tmp.result=as.data.frame(matrix(nrow = ncol(bacteria),ncol = 8))
colnames(tmp.result)=c("Gene","Taxa","Beta_taxa","Beta_snp","Beta_interact","P_taxa","P_snp","P_interact")
for(n in 1:ncol(bacteria)){
  tmp.bacterria=colnames(bacteria)[n]
  tmp.test=merge(tmp.data,bacteria[,tmp.bacterria,drop=F],all=F,by.x ="Row.names",by.y="row.names")
  colnames(tmp.test)=c("BiopsyID","ResearchID","Gene","SNP","Bac")
  tmp.test=tmp.test[tmp.test$SNP!=-1,]
  tmp.model=lm(Gene ~ Bac + SNP + Bac * SNP,data = tmp.test)
  tmp.model=as.data.frame(summary(tmp.model)$coefficients)
  
  tmp.result$Taxa[n]=tmp.bacterria
  tmp.result$Gene[n]=tmp.gene
  tmp.result$Beta_taxa[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac"]
  tmp.result$Beta_snp[n]=tmp.model$Estimate[rownames(tmp.model)=="SNP"]
  tmp.result$Beta_interact[n]=tmp.model$Estimate[rownames(tmp.model)=="Bac:SNP"]
  tmp.result$P_taxa[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac"]
  tmp.result$P_snp[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="SNP"]
  tmp.result$P_interact[n]=tmp.model$`Pr(>|t|)`[rownames(tmp.model)=="Bac:SNP"]
  
}














