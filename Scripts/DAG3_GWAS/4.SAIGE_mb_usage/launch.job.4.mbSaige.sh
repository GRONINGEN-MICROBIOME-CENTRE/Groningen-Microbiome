#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=06:59:59
#SBATCH --mem=32gb
### 02-20-2020 changed memory from 40 to 10gb
### 24-02-2020 changed memory from 10 to 32gb
## added full path to rscript
ml RPlus
hostname

echo "outcome: ${1}"
echo "snpfile: ${2}"
echo "testfolder: ${3}"
echo "outputfolder: ${4}"
echo "features: ${5}"
echo "chromosome: ${6}"


Rscript /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/4.Saige_GWAS_permb_per_chromosome.R -m ${1} \
                                                                                                 -s ${2} \
                                                                                                 -t ${3} \
                                                                                                 -o ${4} \
                                                                                                 -f ${5} \
                                                                                                 -c ${6}
