### 07-02-2020 - added chromosome X
jobpath="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/launch.job.4.mbSaige.sh"
list="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/3rd.list"
logs="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/logs"
i=242 ####change for the log file naming and job namning
cat $list | while read line 
do
  let i+=1
  for chr in {1..22} "X"
    do

   mb=${line}
   snpfile="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/genotypes/UGLI_SNPs_for_GRM.gds"
   testgenotypes="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_gds/"
   out="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa"
   features="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/features_and_phenotypes/DAG3_gen_BMI_merged_log_v2.txt"
   chr=${chr}
   mkdir -p $out



   sbatch -J "${i}.saigemb.chr${chr}" -o "${logs}/saigemb.${i}.chr${chr}.out" -e "${logs}/saigemb.${i}.chr${chr}.err" -v ${jobpath} ${mb} ${snpfile} ${testgenotypes} ${out} ${features} ${chr}

   sleep 0.5

   done

done



chr=22
mb="k_Bacteria.p_Proteobacteria.c_Betaproteobacteria.o_Burkholderiales.f_Sutterellaceae.g_Sutterella" ## test v1
mb="k__Archaea.p__Euryarchaeota.c__Methanobacteria" ## test v2

mb="k__Bacteria.p__Proteobacteria.c__Betaproteobacteria.o__Burkholderiales.f__Sutterellaceae.g__Sutterella"

q