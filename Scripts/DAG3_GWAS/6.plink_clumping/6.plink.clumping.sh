#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=5:59:59
#SBATCH --mem=32gb
#SBATCH --output=/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/scripts/clumping.out
#SBATCH --error=/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/scripts/clumping.err

## changelog
## 24-02-2020 changed wkdirs and files
## changed required memory to 32gb
## 04-03-2020
## added filter to only clump choromosome with significant signals
## added flag at the end of each feature
## changed pvalues clumping to 5x10-8 instead of 1x10-8

#wkdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa"
list="/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2/to.plot.2"
wkdir="/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2/"
refdir="/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/reference_plink_GSA"

#ml PLINK/1.9-beta6-20190617
ml plink


cat $list | while read line 
do
  mb=${line}
  for chr in {1..22} "X"
    do
    wkfile="${wkdir}/${mb}_gwas/p-values.chr.${chr}"
    
    ### try every chromosome for significant signals
    nsig=$(awk -F'\t' -vOFS='\t' '$12<0.00000005'  ${wkfile}|wc -l)
    ##in significant then clump
    if [ ${nsig} -gt 0 ]; then
    ## create the corresponding  SNP colum in the p-values files
    awk -F'\t' -vOFS='\t' '{ $(NF+1)=$2 ":" $3"_"$5"_"$6 ; print}' $wkfile  > ${wkfile}.v2
    sed  -i -e '1s/chr:pos_ref_alt/SNP/g' ${wkfile}.v2
    plink --bfile ${refdir}/chr_${chr} \
          --clump ${wkfile}.v2 \
          --clump-field  pval  \
          --clump-p1 0.00000005 \
          --clump-r2 0.20 \
          --out ${wkdir}/${mb}_gwas/signal.${chr}
    fi
    done
    echo "$mb has finished clumping eval"
    
done

