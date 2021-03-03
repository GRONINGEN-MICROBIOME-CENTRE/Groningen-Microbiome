#!/bin/bash
#SBATCH --job-name=h_LL11_A01_871
#SBATCH --error=__h_LL11_A01_871.err
#SBATCH --output=__h_LL11_A01_871.out
#SBATCH --mem=32gb
#SBATCH --time=23:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load Python/3.6.4-intel-2018a
export PATH="$PATH":/data/umcg-rgacesa/tools/metaphlan2/utils:/data/umcg-rgacesa/tools_install/FastQC:/data/umcg-rgacesa/tools/metaphlan2:/data/umcg-rgacesa/tools_install/bin
echo "Running Humann2" 
# --- RUN HUMANN2 --- 
humann2 --input LL11_A01_871/clean_reads/LL11_A01_871_kneaddata_merged.fastq --output LL11_A01_871/humann2/ --taxonomic-profile LL11_A01_871/metaphlan/LL11_A01_871_metaphlan.txt --threads 16 --o-log LL11_A01_871/LL11_A01_871_humann2.log --remove-temp-output
# --- CLEANUP --- 
echo "Cleaning redundant data"
echo " --> ALL DONE !!! <-- "
