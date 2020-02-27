#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=8:59:59
#SBATCH --mem=40gb

ml RPlus
hostname

echo "workdir: ${1}"
echo "chromosome: ${2}"
echo "root: ${3}"


Rscript /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/vcftoGDS.R -w ${1} \
                                                                        -c ${2} \
                                                                        -r ${3}
