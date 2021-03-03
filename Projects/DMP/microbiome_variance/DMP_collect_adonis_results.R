# ======================================================================
# By: Weersma Group, UMCG (2020)
#
# DMP Variance explained analysis,
#
# helper script, collects all results from big adonis run
#
# Note: Script has to be called with input folder (with adonis results)
# ======================================================================

# === MAIN ===

# === collect command line parameters =======
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  stop("requires input folder")
}
# ===========================================

# go through the folder, collect all results
inFld <- args[1]
fls <- list.files(inFld)

# result collector
resDF <- NULL

# iterate over files
for (f in fls) {
   if (grepl('csv$',f)) {
      print(paste0(' collecting ',f))
      inT <- paste0(inFld,'/',f)
      #print(inT)
      inR <- read.table(inT,sep=',',header=T,quote='"')
      inR$file <- f
      resDF <- rbind.data.frame(resDF,inR)
   }
}

# output files to collect
outF <- "adonis_out.csv"
resDF <- resDF[order(resDF$Var),]
write.table(resDF,outF,sep=',',row.names=F)


