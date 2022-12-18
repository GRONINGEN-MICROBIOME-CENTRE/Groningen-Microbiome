
# RNAseq raw data processing

*Samples input folder check*
---
- Create a folder contains all samples
- One sample has one folder with two fastq files inside. sample_1.fastq and sample_2.fastq


*Build reference index for STAR mapping*
---
- Check STAR genome index.
- Reference genome is from GTEx V7 (https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq)

```
sbatch Build.index.sh
```


*Create job file for each sample*
---
- the pipeline in *Generate.process.sh* contains the following:

1. 26 million paired-end 150-bp reads were generated per sample. The quality of the raw reads was checked using FastQC with default parameters (v0.11.7). 
2. The adaptors and low-quality reads were clipped using Trimmomatic (v0.36) with settings length <50 nucleotides, quality <25. 
3. Then reads quality was checked agian using FastQC.
4. Reads were aligned to the human genome (Homo_sapiens_assembly19.fasta) using STAR (v2.7.3).
5. Reads sorting and mapping statistics were obtained using SAMtools (v0.1.19), sambamba(v0.7.0) and picard (v2.20.5.). 
6. Gene expression was estimated through HTSeq (0.9.1)based on the annotation from GTEx v7(gencode.v19.annotation.patched_contigs.gtf), resulting in an RNA expression dataset of 57825 genes.

Note, the reference genome and annotation files are derived from GTEx v7 (https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq)

```
bash Generate.process.sh $SAMPLE_PATH/
```


*Submit all jobs*
---
- sbatch all jobs.
