sp=data.frame(t(read.delim("~/Documents/UMCG/research-projects/microbiome_data/merged_microbime_data/GMHI_index/4347_final_relative_abundances.txt",header = T,row.names = 1,sep = "\t",check.names = F)),check.names = F)
colnames(sp)=gsub("s__","",colnames(sp))
result=data.frame(sp=as.character(colnames(sp)),presence_rate=colSums(sp>0)/nrow(sp)*100,mean_abu=colMeans(sp))
result$median_abu=NA
for(i in 1:ncol(sp)){
  result$median_abu[i]=median(sp[,i])
}
list=read.delim("~/Documents/UMCG/research-projects/microbiome_data/merged_microbime_data/GMHI_index/50sp.txt",header = T,sep = "\t",check.names = F)
result$GMHI=NA
result$GMHI=ifelse(result$sp%in%as.character(list$Species),"yes","no")
length(which(result$GMHI=="yes"))
result$select_SNV=NA
length(which(result$presence_rate>20&result$mean_abu>0.5|result$presence_rate>20&result$median_abu>0.5))
result$select_SNV[which(result$presence_rate>20&result$mean_abu>0.5|result$presence_rate>20&result$median_abu>0.5)]="yes"
sum(result$mean_abu[which(result$select_SNV=="yes")])#81%
write.table(result,file = "../table/species_list_based_on_GMHI_NatComPaper.txt",quote = F,sep = "\t",row.names = F)



result=result[order(result$mean_abu,decreasing = T),]
result$select_SNV[1:40]="yes"
sum(result$mean_abu[1:40])#76%
result=result[order(result$median_abu,decreasing = T),]
result$select_SNV[1:40]="yes"
sum(result$mean_abu[which(result$select_SNV=="yes")])#82%
length(which(result$select_SNV=="yes"))


##select several sp for real data
result=read.delim("../data/selected_sp_for_real.tsv",header = T,sep = "\t")
result$Sensitivity=result$TP/(result$TP+result$FN)
result$Precision=result$TP/(result$TP+result$FP)
select=data.frame(sp=unique(as.character(result$Taxa)),mean_Sensitivity=NA,mean_Precision=NA)
for(i in 1:nrow(select)){
  select$mean_Sensitivity[i]= mean(result$Sensitivity[which(result$Taxa==as.character(select$sp[i]))],na.rm = T)
  select$mean_Precision[i]= mean(result$Precision[which(result$Taxa==as.character(select$sp[i]))],na.rm = T)
}

plot(select$mean_Sensitivity,select$mean_Precision,xlim = c(0,1),ylim = c(0,1))
library(cluster)
row.names(select)=select$sp
select=select[,2:3]
dist=as.matrix(dist(as.matrix(select),method = "manhattan"))
tmp=data.frame(pam(dist,k=10)$clustering)
select=cbind(select,tmp[row.names(select),])
colnames(select)[3]="cluster"
write.table(select,file = "../table/select_10sp_for_hmp.txt",quote = F,sep = "\t")










