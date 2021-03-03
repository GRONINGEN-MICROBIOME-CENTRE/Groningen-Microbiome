####main codes for core microbiome selection and microbial co-abundance network analysis
#1.core microbiome selection
#2.compute FDR for edges


###100x sampling based core microbiome selection
##core micorbes
sp=read.delim("../data/dag3_8268samples_721species.txt",header = T,sep = "\t",row.names = 1)
sp=read.delim("../data/dag3_8268samples_569pathway_reads.txt",header = T,sep = "\t",row.names = 1,check.names = F,quote = "")
sp=sp[rowSums(sp)>0,colSums(sp>0)>nrow(sp)*0.05]#140 species 347 pathways
library(matrixStats)
result=NULL
sd=NULL
percent=seq(0.01,1,by=0.01)
for (i in percent) {
  tmp_percent=NULL
  for(j in 1:100){
    permutation=sample(row.names(sp),round(nrow(sp)*i,0))
    per_sp=sp[permutation,]
    tmp_percent=rbind(tmp_percent,colSums(per_sp>0)/round(nrow(sp)*i,0))
  }
  result=cbind(result,colMeans(tmp_percent))
  sd=cbind(sd,colSds(tmp_percent))
}
colnames(result)=1:100
colnames(sd)=1:100
library(stringr)
#row.names(result)=str_split_fixed(row.names(result),".s__",2)[,2]
row.names(result)=str_split_fixed(row.names(result),": ",2)[,1]
result=cbind(result,rowMeans(result))
colnames(result)[101]="mean_rate"
pdf("../figure/core_microbiome_path_sampling_threshold.pdf",width = 10,height = 10,useDingbats = F)
par(mfrow=c(5,2),mgp=c(2,1,-0.2))
for(i in 1:nrow(result)){
  plot(as.numeric(colnames(result)[1:100]),result[i,1:100],type = "b",cex=0.4,frame = FALSE, pch = 19, col = "lightpink", xlab = "Sampling percentage (%)", ylab = "Presence rate",main = row.names(result)[i],cex.main=0.8,xlim = c(0,100),ylim = c(min(result[i,1:100])-max(sd[i,1:100]),max(result[i,1:100])+max(sd[i,1:100])))
  segments(as.numeric(colnames(result)[1:100])-0.02,result[i,1:100]-sd[i,1:100],as.numeric(colnames(result)[1:100])+0.02,result[i,1:100]+sd[i,1:100],col = "gray60",lty = 1,cex=0.02)
}
dev.off()
library(ggplot2)
library(ggrepel)
library(reshape2)
data=data.frame(t(result[,1:100]))
data$id=row.names(data)
data=melt(data,id="id")
data$value=data$value*100
data$id=as.numeric(as.character(data$id))
p1=ggplot(data, aes(id, value, group=variable)) +geom_line(color="black", size=0.05,alpha=0.8)+xlab(label = "Sampling percentage (%)") + ylab(label = "Presence rate (%)")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank())+ theme(legend.position = "none")+scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)) 
data=data.frame(result)
data=data[order(data$mean_rate,decreasing = T),]
data$order=1:nrow(data)
data$mean_rate=data$mean_rate*100
data_tmp=data
data_tmp$name=row.names(data_tmp)
data_tmp$name[which(data_tmp$mean_rate<80)]=NA
p2=ggplot(data_tmp,aes(order,mean_rate))+geom_point(shape=20,alpha=1,cex=0.2)+geom_hline(yintercept = 80,linetype=3,colour="mediumpurple1",alpha=0.5)+labs(x="Rank of core microbiome",y="Presence rate (%)")+theme_bw()+theme(panel.grid=element_blank(),axis.title=element_text(size=8))+theme(legend.position="none")+geom_text_repel(aes(order,mean_rate,label=factor(data_tmp$name)),size=1,alpha=0.85,colour="black",direction = "both",segment.size=0.1,segment.alpha = 0.7,segment.color="red")
result=data.frame(result)
result$mean_rate=result$mean_rate*100
data=data.frame(percent=1:100,number=NA)
for(i in 1:100){
  data$number[i]=length(which(result$mean_rate>=i))
}
p3=ggplot(data, aes(percent, number)) +geom_line(color="black", size=0.5,alpha=0.8)+xlab(label = "Presence rate (%)") + ylab(label = "Number of core microbiome")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank())+ theme(legend.position = "none")+scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)) 
write.table(result,file = "../table/result_dag3_core_microbiome_path.txt",quote = F,sep = "\t")
library(gridExtra)
pdf(file="../figure/core_microbiome_path.pdf",useDingbats=F)
grid.arrange(p1,p3,p2,ncol=1,nrow = 3)
dev.off()

###2.sparcc (Python) for network inference
## sparcc
python SparCC.py $INPUT_PATH --cor_file=$OUTPUT_PATH/real_cors_sparcc_sp.txt
python MakeBootstraps.py $INPUT_PATH -n $BOOT_ITER -t permutation_#.txt -p $OUTPUT_PATH/Resamplings/
# compute sparcc on permutated datasets
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99
do
python SparCC.py $OUTPUT_PATH/Resamplings/permutation_$i.txt --cor_file=$OUTPUT_PATH/Bootstraps/sim_cor_$i.txt
done
# compute p-value from bootstraps
python PseudoPvals.py $OUTPUT_PATH/real_cors_sparcc_sp_fg.txt $OUTPUT_PATH/Bootstraps/sim_cor_#.txt $BOOT_ITER -o $OUTPUT_PATH/pvals_two_sided_sp_fg.txt -t two_sided
## spiec-easi
