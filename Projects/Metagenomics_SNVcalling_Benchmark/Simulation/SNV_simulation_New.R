library(tidyverse)
library(viridis)
.libPaths("R_packages")
library(patchwork)
MS = read_tsv("Simulation_Paper_2Refs/Benchmark.tsv")
SS = read_tsv("Simulation_Paper/Benchmark.tsv")
SS$Benchmark_Design = "uniref"


wsjPal <- c('#65C1E8',
            '#D85B63',
            '#D680AD',
            '#5C5C5C',
            '#C0BA80',
            '#FDC47D',
            '#EA3B46')

for (Data in list(MS, SS)){

  Data %>% select(-X9) %>% mutate(Specificity = TN/(TN+FN), Sensitivity = TP/(TP+FN), Precision=TP/(TP+FP) ) -> Data
   
  Data %>% drop_na() %>% group_by(Tool)  %>% summarise(M_Precision = mean(Precision), SD_Precision = sd(Precision),M_sensitivity = mean(Sensitivity), SD_sensitivity = sd(Sensitivity) ) -> Summary_data
  Summary_data %>% ggplot(aes(x=M_Precision, y=M_sensitivity, col=Tool)) + geom_point(size=3,alpha=0.7) + geom_errorbar(aes(ymin=M_sensitivity-SD_sensitivity, ymax=M_sensitivity+SD_sensitivity ), alpha=0.2) +
    geom_errorbarh(aes(xmin=M_Precision-SD_Precision, xmax=M_Precision+SD_Precision), alpha=0.2) + theme_bw() + scale_color_manual(values=wsjPal) + xlab("Mean Precision (TP/TP+FP)") + ylab("Mean sensitivity (TP/TP+FN)") -> Summary_plot
  
  
  ggplot(Data,aes(x=Tool, y=TP/1000)) + geom_boxplot()  + theme_bw() + coord_flip() +  scale_fill_manual(values=wsjPal)+ theme(legend.position = "none") -> TP_plot
  ggplot(Data,aes(x=Tool, y=FN/1000)) + geom_boxplot() +  theme_bw() + coord_flip() +  scale_fill_manual(values=wsjPal)+ theme(legend.position = "none") -> FN_plot
  ggplot(Data,aes(x=Tool, y=FP/1000)) + geom_boxplot() +  theme_bw() + coord_flip() +  scale_fill_manual(values=wsjPal)+ theme(legend.position = "none") -> FP_plot
  
  ggplot(Data,aes(x=Tool, y=Sensitivity)) + geom_boxplot(aes(fill=Tool),alpha=0.7) + geom_point() + theme_bw() + coord_flip() + scale_fill_manual(values=wsjPal)+ theme(legend.position = "none")  -> Sensitivity_plot
  ggplot(Data,aes(x=Tool, y=Precision)) + geom_boxplot(aes(fill=Tool), alpha=0.7) + geom_point() + theme_bw() + coord_flip() + scale_fill_manual(values=wsjPal) + theme(legend.position = "none")  -> Precision_plot
  
  
  Sensitivity_plot / Precision_plot / (TP_plot + FN_plot + FP_plot) | Summary_plot +  plot_layout(guides = "collect") +   plot_annotation(title =unique(Data$Benchmark_Design) )  -> Composition
  print(Composition)

  ggsave(paste(c("Plots/",unique(Data$Benchmark_Design),"General.pdf"),collapse=""),plot = Composition, height=13, width=13, units="in")
  ggsave(paste(c("Plots/",unique(Data$Benchmark_Design),"General.png"), collapse=""),plot = Composition, height=13, width=13, units="in")
  
}

MS_R = read_tsv("Simulation_Paper_2Refs/ROC_curves.tsv") %>% mutate(Benchmark_Design="multiref")
SS_R =  read_tsv("Simulation_Paper/ROC_curves.tsv") %>% mutate(Benchmark_Design="uniref")


for (ROC in list(MS_R, SS_R)){
	
	ROC %>% drop_na() %>% group_by(Tool,Quantile_Phred)  %>% summarise(M_Precision = mean(Precision), SD_Precision = sd(Precision),M_sensitivity = mean(Sensitivity), SD_sensitivity = sd(Sensitivity) ) -> ROC_sum
	ROC_sum %>% mutate(Up_Precision = M_Precision+SD_Precision, Down_Precision = M_Precision-SD_Precision, UP_Sensitivity = M_sensitivity+SD_sensitivity, Down_sensitivity = M_sensitivity-SD_sensitivity ) -> ROC_sum
	ROC_sum %>% mutate(Down_Precision = ifelse(Down_Precision<0,0,Down_Precision), Down_sensitivity=ifelse(Down_sensitivity<0,0,Down_sensitivity)) -> ROC_sum

	#Recall (TP / TP+FN ) , Sensitivity, True Positive Rate
	ROC_sum %>% ggplot(aes(y=M_sensitivity, x=Quantile_Phred)) + geom_line(aes( col=Tool)) + theme_bw() + scale_color_manual(values = wsjPal)  + ylab("Mean sensitivity (TP/TP+FN)") +
	  geom_ribbon(aes(ymin=Down_sensitivity, ymax=UP_Sensitivity, fill=Tool), linetype=2, alpha=0.05)+ scale_fill_manual(values = wsjPal)  -> RSensitivity_plot

	#Precision (TP/ TP+FP) , Positive Prediction value
	ROC_sum %>% ggplot(aes(y=M_Precision, x=Quantile_Phred)) + geom_line(aes(col=Tool)) + theme_bw()+ scale_color_manual(values=wsjPal)    + ylab("Mean Precision (TP/TP+FP)") +
	  geom_ribbon(aes(ymin=Down_Precision, ymax=Up_Precision, fill=Tool), linetype=2, alpha=0.05) + scale_fill_manual(values = wsjPal)  -> RPrecision_plot

	RSensitivity_plot | RPrecision_plot  +  plot_layout(guides = "collect") -> Composition2
	ggsave(paste(c("Plots/",unique(ROC$Benchmark_Design),"ROC_probabilistic.pdf"),collapse=""),plot = Composition2, height=13, width=13, units="in")
	ggsave(paste(c("Plots/",unique(ROC$Benchmark_Design),"ROC_probabilistic.png"), collapse=""),plot = Composition2, height=13, width=13, units="in")
}
