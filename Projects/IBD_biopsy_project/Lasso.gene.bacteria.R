stabs_stability <- function( x, y, gene_name){
  
  ## Extract the expression for this gene (response variable)
  y_i <- y[,grep(paste0("^",gene_name,"$"),colnames(y))]
  
  ## Make sure y_i is numeric before model fitting 
  stopifnot(class(y_i) == "numeric")
  
  ## perform stability selection using glmnet lasso
  stab.glmnet <- stabsel(x = x, y = y_i,
                         fitfun = glmnet.lasso, cutoff = 0.6,
                         PFER = 1)
  
  taxa.selected <- names(stab.glmnet$selected)
  if(length(taxa.selected) == 0) taxa.selected <-"None"
  
  taxa.selected.df <- as.data.frame(cbind(gene_name,taxa.selected))
  
  return(taxa.selected.df)
  
}
library("stabs","glmnet","methods","doParallel")
# lasso gene~bacteria selection

genes_CD.correct=read.table("RNAseq.CLR.CD.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genes_UC.correct=read.table("RNAseq.CLR.UC.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
genes_Control.correct=read.table("RNAseq.CLR.Control.correction.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)

inflammation=read.table("OutputTable/RNAseq.inflammation.compare.txt",sep = "\t",header = T)
inflammation=(inflammation[inflammation$FDR<0.05,])
inflammation=intersect(inflammation$Gene[inflammation$group==3],intersect(inflammation$Gene[inflammation$group==1],inflammation$Gene[inflammation$group==2]))
genes_CD.correct=genes_CD.correct[,colnames(genes_CD.correct) %in% inflammation]
genes_UC.correct=genes_UC.correct[,colnames(genes_UC.correct) %in% inflammation]
genes_Control.correct=genes_Control.correct[,colnames(genes_Control.correct) %in% inflammation]

genus=read.table("OutputTable/CLR.genus.txt",sep = "\t",header = T,stringsAsFactors = F,check.names = F)
dim(genus)

input.genes.list="GUCA2A"

# CD
genus_CD=genus[rownames(genus) %in% rownames(genes_CD.correct),]
genes_CD.correct=genes_CD.correct[rownames(genes_CD.correct) %in% rownames(genus_CD),]

genus_CD=genus_CD[order(rownames(genus_CD)),]
genes_CD.correct=genes_CD.correct[order(rownames(genes_CD.correct)),]
stopifnot(all(rownames(genes_CD.correct) == rownames(genus_CD)))

y <- as.matrix(genes_CD.correct) #response
x <- as.matrix(genus_CD) #predictors

stabs_stability(x,y,input.genes.list)







































