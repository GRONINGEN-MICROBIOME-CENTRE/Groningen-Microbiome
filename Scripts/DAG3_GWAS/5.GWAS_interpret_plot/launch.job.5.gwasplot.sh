#!/bin/bash
#SBATCH --nodes=1
#SBATCH --time=59:59:59
#SBATCH --mem=48gb
#SBATCH --output=/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/scripts/graph_0.5n_filter.out
#SBATCH --error=/groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/scripts/graph_0.5n_filter.err

## added full path to rscript
ml RPlus
hostname


Rscript /groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/scripts/5.mb_GWAS_merge_plot.R -i /groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2 \
                                                                                     -m /groups/umcg-wijmenga/tmp04/umcg-elopera/DAG3/output/out_v2/to.plot.2 
