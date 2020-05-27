# plot QQ plot and calculate lamda of genomics inflation index
args<-commandArgs(T)
meta = args[1]

require(data.table)
library(qqman)
summary=fread(meta,sep=" ",header = F,stringsAsFactors = F)
colnames(summary)=c("P","SNP","CHR","BP")
summary=na.omit(summary)

png(paste0(meta,".QQplot.png"))
qq(summary$P)
dev.off()

png(paste0(meta,".MahattanPlot.png"))
manhattan(summary, chr="CHR", bp="BP", snp="SNP", p="P" )
dev.off()

z = qnorm(summary$P/2)
lambda = round(median(z^2) / qchisq(0.5, 1), 3)

cat(meta,"Lambda---",lambda)
