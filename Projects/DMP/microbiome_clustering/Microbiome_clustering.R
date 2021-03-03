# ==========================================
# By: Weersma Group, UMCG, 20202
#
# Script used for enterotype analysis / microbiome clustering
# on DMP microbiome data
# 
# ==========================================

# load libraries
library(vegan)
library(ggplot2)
library(cluster)
library(factoextra)
library(ade4)
library(clusterSim)

#============================================================================
# calculator for JSD distance
# ===========================================================================
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
 }

# ===========================================================================
# partition around medoid (PAM) clustering
# ===========================================================================
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
}

# ==============================================
# MAIN
# ==============================================

# ==========================================================================================
# import metaphlan data of dag3, filter & re-scale the data (re-calculate rel. abundance)
# ==========================================================================================
metaphlan=read.table("DMP_metaphlan_merged.txt", 
                     header = T,
                     stringsAsFactors = F,sep = "\t",check.names = F,row.names = 1)
taxa=as.data.frame(t(metaphlan))

# I) DATA PREP
# =================================
genus=taxa[,grep("g__",colnames(taxa))]
genus=genus[,grep("s__",colnames(genus),invert = T)]
colnames(genus)=lapply(colnames(genus),function(x){
  strsplit(x,"g__")[[1]][2]
})
genus = genus[rowSums(genus > 0) > 0,]
genus = genus[,colSums(genus > 0) > 0]
genus=as.data.frame(t(genus))
genus.dist=dist.JSD(genus)
genus.cluster=pam.clustering(genus.dist, k=3)

# II) cluster the data using PAM clustering
# ==========================================================================================
# > perform clustering
nclusters = index.G1(t(genus), genus.cluster, d = genus.dist, centrotypes = "medoids")
nclusters=NULL
# > calculate CH index to identify optimal number of clusters
# ==========================================================================================
for (k in 1:10) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    genus.cluster_temp=pam.clustering(genus.dist, k)
    nclusters[k]=index.G1(t(genus),genus.cluster_temp,  d = genus.dist,
                          centrotypes = "medoids")
  }
}
# III) output
# ==========================================================================================
# > CH index plot
pdf("CHindex.pdf")
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
dev.off()
cluster=data.frame(row.names = colnames(genus),Cluster=genus.cluster)
# > table output 
write.table(cluster,file = "Clusters.txt",sep = "\t",row.names = T,quote = F)
# > PCoA colored by clusters
pdf("PCoA.pdf")
beta_diversity=vegdist(as.data.frame(t(genus)),method = "bray")
pcoa_analysis=as.data.frame(cmdscale(beta_diversity,k=4))
table=merge(cluster,pcoa_analysis,by="row.names",all = F)
ggplot (table, aes(V1,V2)) + geom_point(aes(colour = as.factor(Cluster)),size=1) + theme_bw()
dev.off()

# =================================================
# IV) calculation of enterotype composition
# =================================================
library(RColorBrewer)
library(ggplot2)
library(ggsci)

# Data prep function
# > changes metaphlan result to a composition table, selects top n most abundant features
# =========================================================================================
CompositionTable <- function(x,n){ 
  require(foreach)
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  return(composition_table)
}

# > SPECIES COMPOSITION 
# =================================
cluster=read.table("Clusters.species.txt",header = T,sep = "\t",stringsAsFactors = F,row.names = 1)
metaphlan=read.table("DAG3_metaphlan_merged.txt",header = T,
                     stringsAsFactors = F,sep = "\t",check.names = F,row.names = 1)
taxa=as.data.frame(t(metaphlan))

species=taxa[,grep("s__",colnames(taxa))]
species=species[,grep("t__",colnames(species),invert = T)]
colnames(species)=lapply(colnames(species),function(x){
  strsplit(x,"s__")[[1]][2]
})
species = species[rowSums(species > 0) > 0,]
species = species[,colSums(species > 0) > 0]

cluster1=species[rownames(species) %in% rownames(cluster[cluster$Cluster==1,,drop=F]),]
cluster2=species[rownames(species) %in% rownames(cluster[cluster$Cluster==2,,drop=F]),]

cluster1_table=CompositionTable(cluster1,5)
ggplot (cluster1_table, aes(x=ID, y=Relative,fill=Level)) + 
        geom_bar (stat = "identity") + 
        theme_classic() + xlab("") + ylab("Relative_abundance")+
        theme(axis.text.x = element_text(size=0),
              legend.position="bottom")+
        scale_fill_brewer(palette = "Paired")+
        guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("Species.cluster1.pdf")

cluster2_table=CompositionTable(cluster2,5)
ggplot (cluster2_table, aes(x=ID, y=Relative,fill=Level)) + 
  geom_bar (stat = "identity") + 
  theme_classic() + xlab("") + ylab("Relative_abundance")+
  theme(axis.text.x = element_text(size=0),
        legend.position="bottom")+
  scale_fill_lancet()+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("Species.cluster2.pdf")


# > GENUS COMPOSITION 
# =================================
cluster=read.table("Clusters.genus.txt",header = T,sep = "\t",stringsAsFactors = F,row.names = 1)
metaphlan=read.table("DAG3_metaphlan_merged.txt",header = T,
                     stringsAsFactors = F,sep = "\t",check.names = F,row.names = 1)
taxa=as.data.frame(t(metaphlan))

genus=taxa[,grep("g__",colnames(taxa))]
genus=genus[,grep("s__",colnames(genus),invert = T)]
colnames(genus)=lapply(colnames(genus),function(x){
  strsplit(x,"g__")[[1]][2]
})
genus = genus[rowSums(genus > 0) > 0,]
genus = genus[,colSums(genus > 0) > 0]

cluster1=genus[rownames(genus) %in% rownames(cluster[cluster$Cluster==1,,drop=F]),]
cluster2=genus[rownames(genus) %in% rownames(cluster[cluster$Cluster==2,,drop=F]),]

cluster1_table=CompositionTable(cluster1,5)
ggplot (cluster1_table, aes(x=ID, y=Relative,fill=Level)) + 
  geom_bar (stat = "identity") + 
  theme_classic() + xlab("") + ylab("Relative_abundance")+
  theme(axis.text.x = element_text(size=0),
        legend.position="bottom")+
  scale_fill_brewer(palette = "Paired")+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("Genus.cluster1.pdf")

cluster2_table=CompositionTable(cluster2,5)
ggplot (cluster2_table, aes(x=ID, y=Relative,fill=Level)) + 
  geom_bar (stat = "identity") + 
  theme_classic() + xlab("") + ylab("Relative_abundance")+
  theme(axis.text.x = element_text(size=0),
        legend.position="bottom")+
  scale_fill_lancet()+
  guides(fill=guide_legend(nrow=3, byrow=TRUE))
ggsave("Genus.cluster2.pdf")






