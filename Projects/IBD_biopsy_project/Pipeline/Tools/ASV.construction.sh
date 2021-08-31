#!/bin/bash
#SBATCH --job-name=ASV.test           
#SBATCH --output=ASV.test.out
#SBATCH --error=ASV.test.err
#SBATCH --time=23:59:00
#SBATCH --mem=10gb

ml R
input="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Trimmomatic/test/"

mkdir -p "/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Dada2/test"
output="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Dada2/test/"

Rscript /groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/pipeline/Dada2.R $input $output

