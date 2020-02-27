jobpath="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/bcf_filter_info.job.sh"

for chr in {1..22} "X" 
do

input="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf/${chr}.DAG3_imputed_DAG_ID.vcf.gz"
output="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf/${chr}.DAG3_ID_imputed_info_filtered.vcf.gz"
log="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/logs/"

sbatch -J "$chr.info_filter" -o "${log}${chr}.DAG3_ID_imputed_info_filtered.out" -e "${log}${chr}.DAG3_ID_imputed_info_filtered.err" -v ${jobpath} ${input} ${output} 

sleep 0.5

done




chr=22
