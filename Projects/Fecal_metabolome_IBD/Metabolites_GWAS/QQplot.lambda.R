# plot QQ plot and calculate lamda of genomics inflation index
library(qqman)
args<-commandArgs(T)
meta = args[1]

summary=read.table(meta,sep=" ",header = F,stringsAsFactors = F)
colnames(summary)=c("P","SNP","CHR","BP")
summary=na.omit(summary)

pdf(paste0(meta,".QQplot.pdf"))
qq(summary$P)
dev.off()

z = qnorm(summary$P/2)
lambda = round(median(z^2) / qchisq(0.5, 1), 3)

cat(meta,"Lambda---",lambda)
