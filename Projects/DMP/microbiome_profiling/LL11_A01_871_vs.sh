#!/bin/bash
#SBATCH --job-name=vs_LL11_A01_871
#SBATCH --error=__vs_LL11_A01_871.err
#SBATCH --output=__vs_LL11_A01_871.out
#SBATCH --mem=8gb
#SBATCH --time=3:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load Python/2.7.11-foss-2016a
module load Biopython/1.65-foss-2016a-Python-2.7.11
echo "Running VFDB_SB analysis"
# --- RUN shortBRED vs VFDB DB Markers --- 
/data/umcg-rgacesa/tools/shortBRED/shortbred_quantify.py --markers /data/umcg-rgacesa/tools/shortBRED/markers_VFDB/markers_core_VFDB.fa --wgs LL11_A01_871/clean_reads/LL11_A01_871_kneaddata_merged.fastq --threads 16 --results LL11_A01_871/VFDB_SB/LL11_A01_871_vfdb_sb.txt --tmp LL11_A01_871/VFDB_SB/tmp --usearch /data/umcg-rgacesa/tools/usearch/usearch10.0.240_i86linux32
NLINES=$(grep -c '' LL11_A01_871/clean_reads/LL11_A01_871_1_kneaddata_paired_1.fastq)
NREADS=$(( NLINES / 4 ))
