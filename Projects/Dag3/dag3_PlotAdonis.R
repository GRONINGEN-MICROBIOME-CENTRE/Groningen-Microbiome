# ===============================================
# make pretty plots for adonis analysis for DAG3
# DAG3 phenotypes version v26
# 07/09/2020, R. Gacesa (UMCG)

# ===============================================
# libs
library(ggplot2)

setwd('D:/UMCG/DAG3_statistics_V2/')
inDF <- read.table('D:/UMCG/DAG3_results_v26/adonis/adonis_taxa_v26.csv',sep = ',',header=T)
# inDF$FDR.BH <- p.adjust(inDF$pval)
# inDF$Significant <- "No"
# inDF$Significant[inDF$FDR.BH <= 0.05] <- "Yes"
# inDF <- inDF[order(inDF$FDR.BH),]
# 
# write.table(inDF,'D:/UMCG/DAG3_results_v26/adonis/adonis_taxa_v26_rdy.csv',sep=',',row.names = F)
# 
# # load grouped file
# inDF <- read.table('D:/UMCG/DAG3_results_v26/adonis/adonis_taxa_v26_rdy_grouped.csv',sep=',',header=T,stringsAsFactors = T)
# inDF <- inDF[order(inDF$R2),]
# inDF$Variable <- factor(as.character(inDF$Variable),levels = inDF$Variable[order(inDF$R2,decreasing = F)])

# plot config
# yl <- 0.02
# cntplot <- 0
# colrs <- c("red","blue")
# # do plots per group
# for (grp in unique(inDF$Group)) {
#   print(paste0(' >> PLOTTING ',grp))
#   cntplot <- cntplot + 1
#   dftp <- inDF[inDF$Group==grp,]
#   if ("Yes" %in% dftp$Significant & "No" %in% dftp$Significant) {
#     dftp$Significant <- factor(as.character(dftp$Significant),levels=c("No","Yes"))
#     g <- ggplot(dftp,aes_string(y="R2",x="Variable",col="Significant",fill="Significant")) + 
#       scale_color_manual(values=colrs) + scale_fill_manual(values=colrs)
#   } else if ("Yes" %in% dftp$Significant & !("No" %in% dftp$Significant)) {
#     g <- ggplot(dftp,aes_string(y="R2",x="Variable",col="Significant",fill="Significant")) + 
#       scale_color_manual(values=colrs[2]) + scale_fill_manual(values=colrs[2])
#   } else if ("No" %in% dftp$Significant & !("Yes" %in% dftp$Significant)) {
#     g <- ggplot(dftp,aes_string(y="R2",x="Variable",col="Significant",fill="Significant")) + 
#       scale_color_manual(values=colrs[1]) + scale_fill_manual(values=colrs[1])
#   }
#   g <- g + geom_col() + coord_flip() +
#     theme(legend.position="bottom") + theme(text = element_text(size = 14)) + ggtitle(grp) + xlab('') + ylab('Beta-diversity variance explained (R^2)') + 
#     theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,yl) + theme(axis.text.y = element_text(face="bold"))
#   ggsave(plot = g,filename = paste0('plots_adonis/adonis_plot_',grp,'.png'),height = 2*(1.25+8/50*nrow(dftp)),width = 8,limitsize = F)
# }

# variant 2 (now also with nominal significance)
# load grouped file
inDF <- read.table('D:/UMCG/DAG3_results_v26/adonis/adonis_taxa_v26_rdy_grouped.csv',sep=',',header=T,stringsAsFactors = T)
inDF <- inDF[order(inDF$R2),]
inDF$Significant <- as.character(inDF$Significant)
inDF$Significant[inDF$FDR.BH > 0.05 & inDF$pval < 0.05] <- "Nominal"
#inDF$Variable <- factor(as.character(inDF$Variable),levels = inDF$Variable[order(inDF$R2,decreasing = F)])
inDF$Variable2 <- factor(as.character(inDF$Variable2),levels = inDF$Variable2[order(inDF$R2,decreasing = F)])

# plot config
yl <- 0.013
cntplot <- 0
colrs <- c("red","blue","orange")
# do plots per group
for (grp in unique(inDF$Group.smaller)) {
  print(paste0(' >> PLOTTING ',grp))
  cntplot <- cntplot + 1
  dftp <- inDF[inDF$Group.smaller==grp,]
  dftp$Significant <- factor(as.character(dftp$Significant),levels=c("No","Yes","Nominal"))
  g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) + 
         scale_color_manual(values=colrs) + scale_fill_manual(values=colrs)
  
  if ("Yes" %in% dftp$Significant & "No" %in% dftp$Significant & "Nominal" %in% dftp$Significant) {
    dftp$Significant <- factor(as.character(dftp$Significant),levels=c("No","Yes","Nominal"))
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs) + scale_fill_manual(values=colrs)
  } else if ("Yes" %in% dftp$Significant & "No" %in% dftp$Significant & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1,2)]) + scale_fill_manual(values=colrs[c(1,2)])
  } else if ("Yes" %in% dftp$Significant & !("No" %in% dftp$Significant) & "Nominal" %in% dftp$Significant) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(2,3)]) + scale_fill_manual(values=colrs[c(2,3)])
  } else if (!("Yes" %in% dftp$Significant) & ("No" %in% dftp$Significant) & "Nominal" %in% dftp$Significant) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1,3)]) + scale_fill_manual(values=colrs[c(1,3)])
  } else if ("Yes" %in% dftp$Significant & !("No" %in% dftp$Significant) & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(2)]) + scale_fill_manual(values=colrs[c(2)])
  } else if (!("Yes" %in% dftp$Significant) & ("No" %in% dftp$Significant) & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1)]) + scale_fill_manual(values=colrs[c(1)])
  } else if (!("Yes" %in% dftp$Significant) & !("No" %in% dftp$Significant) & ("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(3)]) + scale_fill_manual(values=colrs[c(3)])
  }
  
  
  g <- g + geom_col() + coord_flip() + ggtitle(paste0("Taxa VS ",grp,"")) + 
    theme(legend.position="bottom") + theme(text = element_text(size = 14)) + xlab('') + ylab('Beta-diversity variance explained (R^2)') + 
    theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,yl) + theme(axis.text.y = element_text(face="bold"))
  ggsave(plot = g,filename = paste0('D:/UMCG/DAG3_results_v26/plots/adonis/adonis_taxa_all_plot_',grp,'.png'),height = 2*(1.25+8/50*nrow(dftp)),width = 8,limitsize = F)
}
# extra: meds, only significant and nominally sig
toPlot <- inDF[inDF$Group.smaller=="Medication.Use" & inDF$Significant != "No",]

g <- ggplot(toPlot,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
  scale_color_manual(values=colrs[c(3,2)]) + scale_fill_manual(values=colrs[c(3,2)]) 
g <- g + geom_col() + coord_flip() +
  theme(legend.position="bottom") + theme(text = element_text(size = 14)) + ggtitle(paste0("Taxa VS Medications (sig. results)")) +
  xlab('') + ylab('Beta-diversity variance explained (R^2)') + 
  theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,yl) + theme(axis.text.y = element_text(face="bold"))
ggsave(plot = g,filename = paste0('D:/UMCG/DAG3_results_v26/plots/adonis/adonis_taxa_signom_plot_meds.png'),height = 2*(1.25+8/50*nrow(toPlot)),width = 8,limitsize = F)

# ====================================================================
# ====================================================================
#     ==================== PATHWAYS ==============================  
# ====================================================================
# ====================================================================
# variant 2 (now also with nominal significance)
# load grouped file
inDF <- read.table('D:/UMCG/DAG3_results_v26/adonis/adonis_pathways_rdy_grouped_v26.csv',sep=',',header=T,stringsAsFactors = T)
inDF <- inDF[order(inDF$R2),]
inDF$Significant <- as.character(inDF$Significant)
inDF$Significant[inDF$FDR.BH > 0.05 & inDF$pval < 0.05] <- "Nominal"
#inDF$Variable <- factor(as.character(inDF$Variable),levels = inDF$Variable[order(inDF$R2,decreasing = F)])
inDF$Variable2 <- factor(as.character(inDF$Variable2),levels = inDF$Variable2[order(inDF$R2,decreasing = F)])

# plot config
yl <- 0.013
cntplot <- 0
colrs <- c("red","blue","orange")
# do plots per group
for (grp in unique(inDF$Group.smaller)) {
  print(paste0(' >> PLOTTING ',grp))
  cntplot <- cntplot + 1
  dftp <- inDF[inDF$Group.smaller==grp,]
  dftp$Significant <- factor(as.character(dftp$Significant),levels=c("No","Yes","Nominal"))
  g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) + 
    scale_color_manual(values=colrs) + scale_fill_manual(values=colrs)
  
  if ("Yes" %in% dftp$Significant & "No" %in% dftp$Significant & "Nominal" %in% dftp$Significant) {
    dftp$Significant <- factor(as.character(dftp$Significant),levels=c("No","Yes","Nominal"))
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs) + scale_fill_manual(values=colrs)
  } else if ("Yes" %in% dftp$Significant & "No" %in% dftp$Significant & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1,2)]) + scale_fill_manual(values=colrs[c(1,2)])
  } else if ("Yes" %in% dftp$Significant & !("No" %in% dftp$Significant) & "Nominal" %in% dftp$Significant) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(2,3)]) + scale_fill_manual(values=colrs[c(2,3)])
  } else if (!("Yes" %in% dftp$Significant) & ("No" %in% dftp$Significant) & "Nominal" %in% dftp$Significant) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1,3)]) + scale_fill_manual(values=colrs[c(1,3)])
  } else if ("Yes" %in% dftp$Significant & !("No" %in% dftp$Significant) & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(2)]) + scale_fill_manual(values=colrs[c(2)])
  } else if (!("Yes" %in% dftp$Significant) & ("No" %in% dftp$Significant) & !("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(1)]) + scale_fill_manual(values=colrs[c(1)])
  } else if (!("Yes" %in% dftp$Significant) & !("No" %in% dftp$Significant) & ("Nominal" %in% dftp$Significant)) {
    g <- ggplot(dftp,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
      scale_color_manual(values=colrs[c(3)]) + scale_fill_manual(values=colrs[c(3)])
  }
  
  
  g <- g + geom_col() + coord_flip() + ggtitle(paste0("PWYs VS ",grp,"")) + 
    theme(legend.position="bottom") + theme(text = element_text(size = 14)) + xlab('') + ylab('Beta-diversity variance explained (R^2)') + 
    theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,yl) + theme(axis.text.y = element_text(face="bold"))
  ggsave(plot = g,filename = paste0('D:/UMCG/DAG3_results_v26/plots/adonis/adonis_pathways_all_plot_',grp,'.png'),height = 2*(1.25+8/50*nrow(dftp)),width = 8,limitsize = F)
}
# extra: meds, only significant and nominally sig
toPlot <- inDF[inDF$Group.smaller=="Medication.Use" & inDF$Significant != "No",]

g <- ggplot(toPlot,aes_string(y="R2",x="Variable2",col="Significant",fill="Significant")) +
  scale_color_manual(values=colrs[c(3,2)]) + scale_fill_manual(values=colrs[c(3,2)]) 
g <- g + geom_col() + coord_flip() +
  theme(legend.position="bottom") + theme(text = element_text(size = 14)) + ggtitle(paste0("PWYs VS Medications (sig. results)")) +
  xlab('') + ylab('Beta-diversity variance explained (R^2)') + 
  theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(0,yl) + theme(axis.text.y = element_text(face="bold"))
ggsave(plot = g,filename = paste0('D:/UMCG/DAG3_results_v26/plots/adonis/adonis_pathways_signom_plot_meds.png'),height = 2*(1.25+8/50*nrow(toPlot)),width = 8,limitsize = F)