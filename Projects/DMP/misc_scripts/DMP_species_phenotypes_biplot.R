# ==================================
# By: Weersma Group, UMCG (2020)
# DMP project, plotting function for
# biplot of species bray-curtis PCoA
# with phenotypes
# ==================================

# == load libraries ==
library(ggplot2)
library (vegan)

# 
#boxy:/groups/umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/forbiplots/Bray_Curtis.RData

# ==========================================
#    ============ MAIN ================
# ==========================================
# set wd
setwd('umcg-lifelines/tmp03/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/forbiplots')

# load bray-curtis dissimilarity matrix
load("Bray_Curtis.RData")
# load phenotypes data
phe=read.csv("phenotypes_arnau.csv")
# load microbiome data (species)
sp=read.table("microbiome_8208samples_S_filtered_prab.csv")

# tidy up column names
colnames(sp)=sapply(strsplit(colnames(sp), ".", fixed=TRUE), tail, 1)
colnames(sp)=gsub(".*s__","",colnames(sp))

# # perform NMDS clustering (not used in DMP manuscript)
# nmds = metaMDS(sp2, distance = "bray")
# 
# # align data between tables
 phe3=phe[row.names(phe)%in%row.names(sp),]
phe3$SEX=1
phe3$SEX[phe3$ANTHRO.SEX=="F"]=2
phe3$HEALTH2[ phe3$HEALTH_STATUS=="poor" ]=2
phe3$HEALTH2[ phe3$HEALTH_STATUS=="mediocre" ]=3
phe3$HEALTH2[ phe3$HEALTH_STATUS=="good" ]=4
phe3$HEALTH2[ phe3$HEALTH_STATUS=="very good" ]=5
phe3$HEALTH2[ phe3$HEALTH_STATUS=="excellent" ]=6
phe3=phe3[match(rownames(sp2),rownames(phe3)),]
# 
# # fit phenotype vectors on the NMDS data (not used in DMP manuscript)
# en=envfit(nmds,phe3, na.rm=T, permutations=9999)

data.scores$SEX = phe3$SEX

# calculate bray-curtis dissimilarity
dist_sp <- vegdist(sp2,  method = "bray")
# do PCOA
myPCOA <- pcoa(dist_sp)
phe4=phe3
phe4$BRISTOL[is.na(phe4$BRISTOL)]=median(phe4$BRISTOL)
#en=envfit(myPCOA,phe4, na.rm=T, permutations=9999)
myPCOAv2= cmdscale(dist_sp)

# fit phenotype vectors on the PCOA data
en2=envfit(myPCOAv2,phe4, na.rm=T, permutations=9999)
en_coord_cont = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)
en_coord_sp = as.data.frame(scores(en2, display = "species"))
coord2=en_coord_cont/2
row.names(coord2)=c("AGE","BMI","BRISTOL STOOL SCALE","SEX","HEALTH STATUS")

# plot PCOA with phenotype vectors
g <- ggplot(data = en_coord_sp, aes(x = Dim1, y = Dim2)) + 
  geom_point(data = en_coord_sp, size = 2, alpha = 0.1, colour="black")  +
  geom_segment(aes(x = 0, y = 0, xend = Dim1, yend = Dim2), data = coord2, size =1, colour = "firebrick") +
  geom_text(data = coord2, aes(x = Dim1, y = Dim2), colour = "firebrick", fontface = "bold", label = row.names(coord2)) + 
  theme_bw( ) +xlab ("PCoA1") + ylab ("PCoA2")  + xlim(-0.75,0.32) 

sp3=sp2[,"Prevotella_copri", drop=F]

en_coord_sp=merge(en_coord_sp,sp3, by="row.names")

ggsave(plot = g,filename = 'DMP_biplot.png')
