#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=35:59:59
#SBATCH --mem=32gb
## 24-02-2020 changed wkdirs and files
##changed required memory to 32gb
list="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa/to.clump.list"
wkdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa"
ml PLINK/1.9-beta6-20190617

cat $list | while read line 
do
  for chr in {1..22} "X"
    do
    mb=${line}
    wkfile="${wkdir}/${mb}_gwas/p-values.chr.${chr}"
  
    ## create the corresponding  SNP colum in the p-values files
    awk -F'\t' -vOFS='\t' '{ $(NF+1)=$2 ":" $3"_"$5"_"$6 ; print}' $wkfile  > ${wkfile}.v2
    sed  -i -e '1s/chr:pos_ref_alt/SNP/g' ${wkfile}.v2
    plink --bfile /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/imputed_vcf/plink/chr_${chr} \
          --clump ${wkfile}.v2 \
          --clump-field  pval  \
          --clump-p1 0.00000001 \
          --clump-r2 0.20 \
          --out ${wkdir}/${mb}_gwas/plink.clumped.${chr}
    
    done
done

