for batch in {01..87}

do

echo '#!/bin/bash'  > Batch.${batch}.plot.sh
echo "#SBATCH --job-name=QQplot.batch.${batch}" >> Batch.${batch}.plot.sh
echo "#SBATCH --error=QQplot.batch.${batch}.err" >> Batch.${batch}.plot.sh
echo "#SBATCH --output=QQplot.batch.${batch}.out" >> Batch.${batch}.plot.sh
echo "#SBATCH --mem=10gb" >> Batch.${batch}.plot.sh
echo "#SBATCH --time=5:00:00" >> Batch.${batch}.plot.sh
echo "#SBATCH --cpus-per-task=6" >> Batch.${batch}.plot.sh

echo "ml R" >> Batch.${batch}.plot.sh
echo "awk 'NR!=1{print \$1}' /groups/umcg-gastrocol/tmp04/Metabolic_Project/Metabolics/Numeric/CD_numeric.metabolic.000${batch}.txt | sort | uniq > Merged_out/Features.${batch}.txt" >> Batch.${batch}.plot.sh

echo "cat Merged_out/Features.${batch}.txt | while read line" >> Batch.${batch}.plot.sh
echo "do" >> Batch.${batch}.plot.sh
echo "Rscript QQplot.lambda.R Merged_out/\$line.Batch.${batch}.summary.txt" >> Batch.${batch}.plot.sh
echo "done" >> Batch.${batch}.plot.sh

done

