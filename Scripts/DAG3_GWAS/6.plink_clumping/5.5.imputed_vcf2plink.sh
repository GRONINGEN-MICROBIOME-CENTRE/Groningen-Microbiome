#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=09:00:59
#SBATCH --mem=40gb


wkdir="/groups/umcg-ugli/tmp01/umcg-elopera/testDAG/"

ml PLINK/2.0-alpha2-20191006

for chr in {1..22} "X"
do
###first create the files with the reference SNP positions form the original .vgz file

plink2 --vcf ${wkdir}/imputed_vcf/${chr}.DAG3_imputed_DAG_ID.vcf.gz \
--const-fid  \
--set-all-var-ids  @:#_\$r_\$a  \
  --make-bed \
--out ${wkdir}/imputed_vcf/chr_${chr}
done    
