#!/bin/bash
#SBATCH --job-name=preQC.test           
#SBATCH --output=preQC.test.out
#SBATCH --error=preQC.test.err
#SBATCH --time=23:59:00
#SBATCH --mem=10gb

ml FastQC

input="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Input/test"
run_name=$(basename $input)

mkdir -p /groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/preQC/$run_name
out_preQC="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/preQC/"

for file in $input/*

do

name=${file%_001*}

fastqc $file -o $out_preQC/$run_name/

done      
         
ml multiqc
         
multiqc $out_preQC/$run_name/ -o $out_preQC/$run_name/
