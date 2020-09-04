

for i in Batch*
do
cd /groups/umcg-gastrocol/tmp04/Metabolic_Project/Numeric_meta/$i

for chr in {1..22}
do
sbatch MetaAnalysis.chr.${chr}.sh
done

done



for i in Batch*
do
cd /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/$i

for chr in {1..22}
do
sbatch MetaAnalysis.chr.${chr}.sh
done

done

