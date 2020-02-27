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



ugli_samples="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/filter_ugli_ID"
pairing_samples="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/replace_UGLIID"

module load BCFtools


hostname
echo "input: ${1}"
echo "output: ${2}"

###sample_filter
bcftools view --force-samples --samples-file ${ugli_samples} --thread 3 -Oz $1 > $2
###sample_rename
bcftools reheader --samples ${pairing_samples} --thread 3 -Oz $2 > $3
###

#bcftools reheader --samples replace_UGLIID ${chr}.DAG3_imputed_UGLI_ID.vcf.gz -o ${chr}.DAG3_imputed_DAG_ID.vcf.gz


#done




