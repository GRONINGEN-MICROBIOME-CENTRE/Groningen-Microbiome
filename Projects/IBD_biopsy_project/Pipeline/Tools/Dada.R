
args<-commandArgs(T)

run_dada2_pipeline = function(
  project.folder,output.folder,
  R1.tag="1_paired.fq.gz",
  R2.tag="2_paired.fq.gz",
  unite.species="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/database/silva_species_assignment_v132.fa",
  unite.ref="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/database/silva_nr_v132_train_set.fa",
  visualize = T
) {

library(dada2)

files = list.files(project.folder)

fnFs = paste(sep="/",project.folder,grep(R1.tag,files,value = T))
fnRs = paste(sep="/",project.folder,grep(R2.tag,files,value = T))

if(!all(sub(R1.tag,"",fnFs)==sub(R2.tag,"",fnRs))) {
  stop("Something is wrong with the names of your sequences, please check\n")
} else {
  sample.names = sub(paste0(project.folder,"/"),"",sub(R1.tag,"",fnFs))
}
FWD <- "CCTACGGGNGGCWGCAG"
REV <- "GACTACHVGGGTATCTAATCC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(project.folder, "filtN", basename(fnFs))
fnRs.filtN <- file.path(project.folder, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)

primerHits <- function(primer, fn) {
  require(ShortRead)
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
aa=rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
write.table(aa, file =paste0(output.folder,"Primer.check.before.cutadapt.txt"),sep="\t",quote=F)

cutadapt <- "/home/umcg-hushixian/.local/bin/cutadapt"

path.cut <- file.path(project.folder, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
# Trim FWD and the reverse-complement of REV off of R1 (forward reads)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
# Trim REV and the reverse-complement of FWD off of R2 (reverse reads)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}
bb=rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))
write.table(bb, file =paste0(output.folder,"Primer.check.after.cutadapt.txt"),sep="\t",quote=F)

# Forward and reverse fastq filenames have the format:
cutFs <- sort(list.files(path.cut, pattern = "_1_paired.fq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2_paired.fq.gz", full.names = TRUE))

if(!all(sub(R1.tag,"",cutFs)==sub(R2.tag,"",cutRs))) {
  stop("Something is wrong with the names of your sequences, please check\n")
} else {
  sample.names =  paste("S",sapply(strsplit(basename(cutFs), "_"), `[`, 1),sep = "")
}

#extra filtering within data2
filtFs <- file.path(path.cut, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path.cut, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) = sample.names
names(filtRs) = sample.names

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, minLen = 160,
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE

#training error models
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


#run dada2 clustering
dadaFs <- dada(filtFs, err=errF,pool=TRUE, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR,pool=TRUE, multithread=TRUE)
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

#extra filtering by length 
nchars.seq = nchar(colnames(seqtab))
maxreads = as.integer(names(which.max(table(nchars.seq))))
seqtab.filtered = seqtab[, nchars.seq > (maxreads -5) & nchars.seq < (maxreads + 5)]

#checking for chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#generate filtering summary statistics
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged",
                     "nonchim")
rownames(track) <- sample.names

#table(nchar(getSequences(seqtab.nochim)))
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
species <- assignSpecies(seqtab.nochim, unite.species, tryRC = TRUE)


#length.sequences = nchar(rownames(taxa))
write.table(track, file =paste0(output.folder,"dada2.read.statistics.txt"),sep="\t")
write.table(seqtab.nochim, file = paste0(output.folder,"SV.table.txt"),sep="\t")
write.table(taxa,file = paste0(output.folder,"SV.taxa.txt"),sep="\t")
write.table(species,file = paste0(output.folder,"SV.species.txt"),sep="\t")

#summary taxa
aa=t(seqtab.nochim)
summarize_taxa <- function(counts, taxonomy) {
  require('plyr')
  if(is.matrix(taxonomy)) {
    #message('multiple taxonomies')
    alply(taxonomy, 2, summarize_taxa, counts = counts, .dims = TRUE)
  } else if(is.matrix(counts)) {
    #message('multiple counts')
    apply(counts, 2, summarize_taxa, taxonomy = taxonomy)
  } else {
    #message('summarize')
    tapply(counts, taxonomy, sum)
  }
}
mm=summarize_taxa(aa, taxa)

write.table(as.data.frame(mm$Kingdom), file =paste0(output.folder,"Kingdom.level.txt"),sep="\t",quote=F)
write.table(as.data.frame(mm$Phylum), file =paste0(output.folder,"Phylum.level.txt"),sep="\t",quote=F)
write.table(as.data.frame(mm$Class), file =paste0(output.folder,"Class.level.txt"),sep="\t",quote=F)
write.table(as.data.frame(mm$Order), file =paste0(output.folder,"Order.level.txt"),sep="\t",quote=F)
write.table(as.data.frame(mm$Family), file =paste0(output.folder,"Family.level.txt"),sep="\t",quote=F)
write.table(as.data.frame(mm$Genus), file =paste0(output.folder,"Genus.level.txt"),sep="\t",quote=F)


}


run_dada2_pipeline(
  project.folder = args[1],
  output.folder = args[2],
  visualize = F
)
