#!/bin/bash
#SBATCH --job-name=m_LL11_A01_871
#SBATCH --error=__m_LL11_A01_871.err
#SBATCH --output=__m_LL11_A01_871.out
#SBATCH --mem=16gb
#SBATCH --time=7:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load Biopython/1.65-foss-2016a-Python-2.7.11
module load Bowtie2/2.2.8-foss-2018a
echo "Running Metaphlan" 
# --- RUN METAHPLAN --- 
/data/umcg-rgacesa/tools/metaphlan2/metaphlan2.py LL11_A01_871/clean_reads/LL11_A01_871_kneaddata_merged.fastq --input_type multifastq --nproc 16 --bowtie2db /data/umcg-rgacesa/tools/metaphlan2/metaphlan_databases/mpa_v20_m200 -o LL11_A01_871/metaphlan/LL11_A01_871_metaphlan.txt --tmp_dir LL11_A01_871/metaphlan_tmp 2>&1 | tee LL11_A01_871/LL11_A01_871_metaphlan.log
# --- SUBMIT HUMANN2 JOB ---- 
echo "Submitting humann2 job" 
sbatch LL11_A01_871_p3.sh 
