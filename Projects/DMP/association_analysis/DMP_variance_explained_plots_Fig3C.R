# ==================================
# By: Weersma Group, UMCG (2020)
# DMP project, plotting functions
# ==================================

# ============================================================================================
# ============================================================================================
# ===================================== ADONIS PLOTS =======================================
# ============================================================================================
# ============================================================================================
library(ggplot2)

#set wd
#NOTE: current path is placeholder, it should be folder with this code and ./group_adonis sub-folder
setwd('D:/Vbox/shared/dag/git_14_05/DMP/microbiome_variance/') 

# input folder
inFldr <- "./group_adonis"

# phenotype grouping file
inGroups <- read.table(paste0(inFldr,'/groups_v2.csv'),sep=',',header = T)
s2R2 <- c()

# read results folder for multivariate adonis, taxa
inFolder <- paste0(inFldr,'/output_adonis_multivariate/step2/taxa/')
t <- dir(path = inFolder,pattern = '*.RData')
t <- t[order(t)]
for (cnt in c(1:nrow(inGroups))) {
  inR <- readRDS(paste0(inFolder,'/',t[cnt]))
  R2sum <- sum(inR$aov.tab$R2[seq(1,length(inR$aov.tab$R2)-2)])
  s2R2 <- c(s2R2,R2sum)
}

allRes <- cbind.data.frame(inGroups$Group,s2R2)
colnames(allRes) <- c("Group","Adonis.seq")

# big groups
allResBigg <- allRes[-c(3,4,5,6,7),]
# per disease
allResDis <- allRes
allResDis <- allResDis[c(3,4,5,6,7),]
allResDis[[2]] <- allResDis[[2]]/sum(allResDis[[2]])*allResBigg[[2]][10]

# load & process (pwys)
# ====================
s2R2 <- c()
inFolder <-  paste0(inFldr,'/output_adonis_multivariate/step2/pwys/')
t <- dir(path = inFolder,pattern = '*.RData')
t <- t[order(t)]
for (cnt in c(1:nrow(inGroups))) {
  inR <- readRDS(paste0(inFolder,'/',t[cnt]))
  R2sum <- sum(inR$aov.tab$R2[seq(1,length(inR$aov.tab$R2)-2)])
  s2R2 <- c(s2R2,R2sum)
}
allResP <- cbind.data.frame(inGroups$Group,s2R2)
colnames(allResP) <- c("Group","Adonis.seq")

# big groups
allResPBigg <- allResP[-c(3,4,5,6,7),]
# per disease
allResPDis <- allResP
allResPDis <- allResPDis[c(3,4,5,6,7),]
allResPDis[[2]] <- allResPDis[[2]]/sum(allResPDis[[2]])*allResPBigg[[2]][10]

# prep output folder
if (!dir.exists('./group_adonis/plots_adonis/')) {
  dir.create('./group_adonis/plots_adonis/')
}

# ========== do plots (taxa) =============
# (Fig 3/c, part I)
inDFss <- allResBigg[,c("Group","Adonis.seq")]
colnames(inDFss) <- c("Group","R2")
inDFss$Group <- factor(as.character(inDFss$Group),levels = inDFss$Group[order(inDFss$R2,decreasing = F)])
inDFss <- inDFss[order(inDFss$R2, decreasing = T),]
inDFss$csumR2 <- inDFss$R2[1]
for (rn in c(2:nrow(inDFss))) {inDFss$csumR2[rn] <- inDFss$csumR2[rn-1] + inDFss$R2[rn]}
inDFss$ypos <- inDFss$csumR2
for (rn in c(1:nrow(inDFss))) {inDFss$ypos[rn] <- inDFss$ypos[rn] - inDFss$R2[rn]/2}
inDFss$R2lbl <- paste0(format(inDFss$R2*100,digits = 1),'%')
inDFss$Data <- "Taxa"
rownames(inDFss) <- inDFss$Group

#colrs <- c("red","blue","orange")

g <- ggplot(inDFss,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col() + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + ylim(0,0.165) + theme_classic()
g

# ========== prep for plots (pwys) ===================
# (Fig 3/c, part II)
inDFssP = allResPBigg
colnames(inDFssP) <- c("Group","R2")
inDFssP$Group <- factor(as.character(inDFssP$Group),levels = levels(inDFss$Group) )
rownames(inDFssP) <- inDFssP$Group
inDFssP <- inDFssP[rownames(inDFss),]

inDFssP$csumR2 <- inDFssP$R2[1]
for (rn in c(2:nrow(inDFssP))) {inDFssP$csumR2[rn] <- inDFssP$csumR2[rn-1] + inDFssP$R2[rn]}
inDFssP$ypos <- inDFssP$csumR2
for (rn in c(1:nrow(inDFssP))) {inDFssP$ypos[rn] <- inDFssP$ypos[rn] - inDFssP$R2[rn]/2}
inDFssP$R2lbl <- paste0(format(inDFssP$R2*100,digits = 1),'%')
inDFssP$Data <- "PWYs"

g <- ggplot(inDFssP,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col() + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + ylim(0,0.165) + theme_classic()
g

# merged plot (Fig 3/c)
# ======================
inDFssPT <- rbind.data.frame(inDFss,inDFssP)
inDFssPT$Data <- factor(inDFssPT$Data,levels = c("Taxa","PWYs"))
g <- ggplot(inDFssPT,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col(width=0.85) + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + theme_bw()
g <- g + ylab("Beta-diversity variance explained (R2)") + theme(text = element_text(size = 14)) + ggtitle("") + xlab("Data layer")# + ylim(0,0.3)
g
ggsave(plot = g,filename = "./group_adonis/plots_adonis/Fig_3c_adonis_taxa_pwys_seq_pergroup.png",width = 5,height = 8) 
write.table(inDFssPT,file = "./group_adonis/plots_adonis/Fig_3c_adonis_taxa_pwys_seq_pergroup_datatable.csv",sep=',',row.names = F)


# =========== for diseases only =============
# ========================================== 

# ========== do plots (taxa) =============
inDFss <- allResDis[,c("Group","Adonis.seq")]
colnames(inDFss) <- c("Group","R2")
inDFss$Group <- factor(as.character(inDFss$Group),levels = inDFss$Group[order(inDFss$R2,decreasing = F)])
inDFss <- inDFss[order(inDFss$R2, decreasing = T),]
inDFss$csumR2 <- inDFss$R2[1]
for (rn in c(2:nrow(inDFss))) {inDFss$csumR2[rn] <- inDFss$csumR2[rn-1] + inDFss$R2[rn]}
inDFss$ypos <- inDFss$csumR2
for (rn in c(1:nrow(inDFss))) {inDFss$ypos[rn] <- inDFss$ypos[rn] - inDFss$R2[rn]/2}
inDFss$R2lbl <- paste0(format(inDFss$R2*100,digits = 1),'%')
inDFss$Data <- "Taxa"
rownames(inDFss) <- inDFss$Group
g <- ggplot(inDFss,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col() + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + ylim(0,0.03) + theme_classic()
g

# ========== prep for plots (pwys) ===================
inDFssP = allResPDis
colnames(inDFssP) <- c("Group","R2")
inDFssP$Group <- factor(as.character(inDFssP$Group),levels = levels(inDFss$Group) )
rownames(inDFssP) <- inDFssP$Group
inDFssP <- inDFssP[rownames(inDFss),]

inDFssP$csumR2 <- inDFssP$R2[1]
for (rn in c(2:nrow(inDFssP))) {inDFssP$csumR2[rn] <- inDFssP$csumR2[rn-1] + inDFssP$R2[rn]}
inDFssP$ypos <- inDFssP$csumR2
for (rn in c(1:nrow(inDFssP))) {inDFssP$ypos[rn] <- inDFssP$ypos[rn] - inDFssP$R2[rn]/2}
inDFssP$R2lbl <- paste0(format(inDFssP$R2*100,digits = 1),'%')
inDFssP$Data <- "PWYs"

g <- ggplot(inDFssP,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col() + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + ylim(0,0.03) + theme_classic()
g


# both
inDFssPT <- rbind.data.frame(inDFss,inDFssP)
inDFssPT$Data <- factor(inDFssPT$Data,levels = c("Taxa","PWYs"))
g <- ggplot(inDFssPT,aes(x = Data, y=R2,col=Group,fill=Group))  + geom_col(width=0.85) + 
  geom_text(aes(y = ypos, label = R2lbl), color = "black") + theme_bw()
g <- g + ylab("Beta-diversity variance explained (R2)") + theme(text = element_text(size = 14)) + ggtitle("") + xlab("Data layer")# + ylim(0,0.3)
g
ggsave(plot = g,filename = "./group_adonis/plots_adonis/adonis_taxa_pwys_seq_diseases.png",width = 5,height = 8) 
write.table(inDFssPT,file = "./group_adonis/plots_adonis/adonis_taxa_pwys_seq_pergroup_diseases.csv",sep=',',row.names = F)