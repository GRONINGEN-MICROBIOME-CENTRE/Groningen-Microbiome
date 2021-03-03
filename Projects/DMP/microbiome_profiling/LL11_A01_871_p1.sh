#!/bin/bash
#SBATCH --job-name=kn_LL11_A01_871
#SBATCH --error=__kn_LL11_A01_871.err
#SBATCH --output=__kn_LL11_A01_871.out
#SBATCH --mem=16gb
#SBATCH --time=7:59:00
#SBATCH --cpus-per-task=8
#SBATCH --open-mode=truncate
# --- LOAD MODULES --- 
module load Biopython/1.65-foss-2016a-Python-2.7.11
module load Bowtie2/2.2.8-foss-2018a
# --- MAKE FOLDERS ---- 
mkdir LL11_A01_871
mkdir ./LL11_A01_871/qc_preclean 
mkdir ./LL11_A01_871/qc_postclean 
mkdir ./LL11_A01_871/filtering_data 
mkdir ./LL11_A01_871/clean_reads 
mkdir ./LL11_A01_871/humann2 
mkdir ./LL11_A01_871/metaphlan 
mkdir ./LL11_A01_871/VFDB_SB 
mkdir ./LL11_A01_871/CARD_SB 
# --- SUBMIT FASTQC JOB (PRE-knead) ---- 
echo "Submitting FastQC (pre-kneaddata)" 
sbatch LL11_A01_871_q1.sh 
# --- RUN KNEADDATA ---- 
echo "Running Kneaddata" 
kneaddata --input LL11_A01_871_1.fq.gz --input LL11_A01_871_2.fq.gz --threads 16 --processes 16 --output LL11_A01_871/filtering_data/ --log LL11_A01_871/LL11_A01_871_kneaddata.log -db /data/umcg-rgacesa/tools_install/trimmomatic_DB/
#   -->  clean kneaddata results: 
cat LL11_A01_871/filtering_data/LL11_A01_871_1_kneaddata_paired_1.fastq > LL11_A01_871/filtering_data/LL11_A01_871_kneaddata_merged.fastq
cat LL11_A01_871/filtering_data/LL11_A01_871_1_kneaddata_paired_2.fastq >> LL11_A01_871/filtering_data/LL11_A01_871_kneaddata_merged.fastq
mv LL11_A01_871/filtering_data/*kneaddata_paired_1.fastq LL11_A01_871/clean_reads
mv LL11_A01_871/filtering_data/*kneaddata_paired_2.fastq LL11_A01_871/clean_reads
mv LL11_A01_871/filtering_data/*kneaddata_merged.fastq LL11_A01_871/clean_reads
rm -r LL11_A01_871/filtering_data
# --- SUBMIT FASTQC JOB (POST-knead) ---- 
echo "Submitting FastQC job (post-knead)" 
sbatch LL11_A01_871_q2.sh 
# --- SUBMIT Virulence finder SB job ---- 
echo "Submitting Virulence finder SB job" 
sbatch LL11_A01_871_vs.sh 
# --- SUBMIT CARD resistome SB job ---- 
echo "Submitting CARD SB job" 
sbatch LL11_A01_871_ac.sh 
# --- SUBMIT METAPHLAN JOB ---- 
echo "Submitting Metaphlan job" 
sbatch LL11_A01_871_p2.sh 
