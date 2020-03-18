#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --ntasks=3
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=32G


if [[ -z "$1" ]]
then
echo "Error: set input"
exit
fi
if [[ -z "$2" ]]
then
echo "Error: set output"
exit
fi

#ugli_samples="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/filter_ugli_ID"
#pairing_samples="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/replace_UGLIID"
module load BCFtools

hostname
echo "input: ${1}"
echo "output: ${2}"

###sample_filter

###filter for imputation quality
bcftools filter -i 'INFO>0.4' -Oz --threads 3 $1 > $2
#done




