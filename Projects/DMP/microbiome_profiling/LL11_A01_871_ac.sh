#!/bin/bash
#SBATCH --job-name=cs_LL11_A01_871
#SBATCH --error=__cs_LL11_A01_871.err
#SBATCH --output=__cs_LL11_A01_871.out
#SBATCH --mem=8gb
#SBATCH --time=3:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load Python/2.7.11-foss-2016a
echo "Running CARD/SB analysis"
# --- RUN shortBRED vs CARD DB Markers --- 
/data/umcg-rgacesa/tools/shortBRED/shortbred_quantify.py --markers /data/umcg-rgacesa/tools/DB_CARD_2018_11/markers_CARD_2018_11.fa --wgs LL11_A01_871/clean_reads/LL11_A01_871_kneaddata_merged.fastq --threads 16 --results LL11_A01_871/CARD_SB/LL11_A01_871_CARD_sb.txt --tmp LL11_A01_871/CARD_SB/tmp --usearch /data/umcg-rgacesa/tools/usearch/usearch10.0.240_i86linux32
rm -r LL11_A01_871/CARD_SB/tmp
