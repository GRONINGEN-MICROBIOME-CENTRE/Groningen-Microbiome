jobpath="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/vct2gds.job.sh"
chr=22
for chr in {1..22} "X" 
do

workdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf"
chr=$chr
root=".DAG3_ID_imputed_info_filtered.vcf.gz"
logs="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/logs/"


sbatch -J "${chr}_vcf2gds" -o "${logs}${chr}_vcf2gds.out" -e "${logs}${chr}_vcf2gds.err" -v ${jobpath} ${workdir} ${chr} ${root}

sleep 0.5

done


#for X
chr="X"
workdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf/X_proc"
chr=$chr
root="_DS_dip_merged.vcf.gz"
logs="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/logs/"
