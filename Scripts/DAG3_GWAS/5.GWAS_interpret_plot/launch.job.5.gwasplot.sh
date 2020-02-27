#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=60:59:59
#SBATCH --mem=48gb

## added full path to rscript
ml RPlus
hostname


Rscript /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/scripts/5.mb_GWAS_merge_plot.R -i /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa \
                                                                     -m /groups/umcg-ugli/tmp01/umcg-elopera/testDAG/output/out_v2_taxa/gwas.list 
