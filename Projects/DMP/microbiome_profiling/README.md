## DMP microbiome profiling, example jobs

This folder shows examples of SLURM jobs used for microbiome profiling of DMP samples (example is for 1 sample, equivalent jobs were executed for all samples):

- LL11_A01_871_q1.sh : QC of raw data
- LL11_A01_871_p1.sh : trimming, filtering & removal of human reads
- LL11_A01_871_q2.sh : QC of cleaned data
- LL11_A01_871_p2.sh : taxonomy profiling (Metaphlan2)
- LL11_A01_871_p3.sh : functional profiling (Humann2)
- LL11_A01_871_vs.sh : profiling of virulence factors
- LL11_A01_871_ac.sh : profiling of antibiotic resistance genes
