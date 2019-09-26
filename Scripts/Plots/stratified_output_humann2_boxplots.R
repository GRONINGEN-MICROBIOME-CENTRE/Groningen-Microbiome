#Creator: Arnau Vich 
#Year: 2019
#Example of barplots to illustrate the bacterial contribution to each pathway
#For this plot you need to subset the pathway of interest including species stratified levels (e.g. PW1, PW1|Ecoli, PW1|Fprausnitzii)



lld <- read.delim("./fig2/lld.txt", row.names=1, check.names = F)
mibs<- read.delim("./fig2/mibs.txt", row.names=1, check.names = F)
ibd<- read.delim("./fig2/ibd.txt", row.names=1, check.names = F)

rel_lld=sweep(lld, 2, colSums(lld), '/')
rel_ibd=sweep(ibd, 2, colSums(ibd), '/')
rel_mibs=sweep(mibs, 2, colSums(mibs), '/')

pheno_ibd=read.delim("./IBD_phenos.txt", row.names=1, check.names = F)
pheno_lld=read.delim("./LLD_phenos.txt", row.names=1, check.names = F)
pheno_mibs=read.delim("./MIBS_phenos.txt", row.names=1, check.names = F)

lld2=merge(pheno_lld, as.data.frame(t(rel_lld)), by = "row.names")
ibd2=merge(pheno_ibd, as.data.frame(t(rel_ibd)), by = "row.names")
mibs2=merge(pheno_mibs, as.data.frame(t(rel_mibs)), by = "row.names")

ibd3=as.data.frame(t(ibd2[,c(1,19,38,49:ncol(ibd2))]))
mibs3=as.data.frame(t(mibs2[,c(1,19,38,49:ncol(mibs2))]))
lld3=as.data.frame(t(lld2[,c(1,19,38,49:ncol(lld2))]))


x1=merge(ibd3, mibs3, by = "row.names", all.x = T)
row.names(x1)=x1$Row.names
x1$Row.names=NULL
x2=merge(x1, lld3, by = "row.names", all.x = T)

x3=as.data.frame(t(x2))

colnames(x3)=unlist(x3[1,])
x3=x3[-1,]

x4=melt(x3,id.vars = c("Row.names","cohort","PPI"))

x4$variable2=gsub(".*s__","",x4$variable)
x4[x4$variable=="PWY0-1297:_superpathway_of_purine_deoxyribonucleosides_degradation|unclassified",]$variable2="unclassified"
ggplot(x4,aes(variable2,as.numeric(as.character(value2)),fill=PPI)) + geom_boxplot(outlier.size = 0.2)+facet_grid(cohort2~.) + theme_minimal()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("Relative contribution to Purine deoxyribonucleosides degradation (PWY0-1297)") + xlab("")

