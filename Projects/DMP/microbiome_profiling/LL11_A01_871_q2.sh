#!/bin/bash
#SBATCH --job-name=q2_LL11_A01_871
#SBATCH --error=__q2_LL11_A01_871.err
#SBATCH --output=__q2_LL11_A01_871.out
#SBATCH --mem=4gb
#SBATCH --time=3:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load FastQC/0.11.7-Java-1.8.0_162
# --- RUN QC --- 
echo "Running FastQC on cleaned reads" 
fastqc -t 8 -q -o LL11_A01_871/qc_postclean LL11_A01_871/clean_reads/LL11_A01_871_kneaddata_merged.fastq
