#!/bin/bash
#SBATCH --job-name=q1_LL11_A01_871
#SBATCH --error=__q1_LL11_A01_871.err
#SBATCH --output=__q1_LL11_A01_871.out
#SBATCH --mem=4gb
#SBATCH --time=3:59:00
#SBATCH --cpus-per-task=4
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load FastQC/0.11.7-Java-1.8.0_162
# --- RUN QC --- 
echo "Running FastQC on unclean reads" 
fastqc -t 8 -q -o LL11_A01_871/qc_preclean LL11_A01_871_1.fq.gz
fastqc -t 8 -q -o LL11_A01_871/qc_preclean LL11_A01_871_2.fq.gz
