# re-normalize RNAseq using clr

# RNA-seq data
gene1=read.table("RNAseq/NewRelease.GeneCount.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
gene2=read.table("RNAseq/OldRelease.ExpressionTable.txt",sep = "\t",header = T,stringsAsFactors = F,row.names = 1,check.names = F)
rownames(gene1)=gsub("\\..*","",rownames(gene1))
overlap=intersect(rownames(gene1),rownames(gene2))
gene1=gene1[rownames(gene1) %in% overlap,]
gene2=gene2[rownames(gene2) %in% overlap,]
colnames(gene1)=gsub(".count","",colnames(gene1))
coupling=read.table("MetaData/couplingFile.txt",header = T,stringsAsFactors = F,check.names = F,sep = "\t")
coupling=coupling[coupling$`RNA-seqID` %in% colnames(gene1),]
gene1=gene1[,colnames(gene1) %in% coupling$`RNA-seqID`]
coupling=coupling[order(coupling$`RNA-seqID`),]
gene1=gene1[,order(colnames(gene1))]
colnames(gene1)=coupling$ChangedID
gene=cbind(gene2,gene1)

# normalize and transform
count1=as.data.frame(t(gene1))
count1=count1[rowSums(count1)!=0,colSums(count1)!=0]
count1=apply(count1,1,function(x){
  x=x/sum(x)
})
count1=as.data.frame(t(count1))
expression1 <- zCompositions::cmultRepl(count1, method="CZM", label=0)
expression1=compositions::clr(expression1)
expression1=as.data.frame(expression1)

count2=as.data.frame(t(gene2))
count2=count2[,colnames(count2) %in% colnames(expression1)]
count2=count2[rowSums(count2)!=0,colSums(count2)!=0]
count2=apply(count2,1,function(x){
  x=x/sum(x)
})
count2=as.data.frame(t(count2))
expression2 <- zCompositions::cmultRepl(count2, method="CZM", label=0)
expression2=compositions::clr(expression2)
expression2=as.data.frame(expression2)

expression1=expression1[,colnames(expression1) %in% colnames(expression2)]
expression2=expression2[,colnames(expression2) %in% colnames(expression1)]

expression=rbind(expression2,expression1)
write.table(expression,file = "RNAseq/Merged.RNAseq.CLR.txt",quote = F,row.names = T,sep = "\t")

# annotation
annot=read.table("annotation.file.txt",sep = "\t",header = T,stringsAsFactors = F)
annot=annot[annot$Gene %in% colnames(expression),]

expression=expression[,order(colnames(expression))]
annot=annot[order(annot$Gene),]

colnames(expression)=annot$id
write.table(expression,file = "RNAseq/Input.RNAseq.CLR.annot.txt",quote = F,row.names = T,sep = "\t")
















































