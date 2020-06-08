mkdir -p Merged_out

for batch in {01..87}
do

for chr in {1..22}
do

zcat Batch.${batch}/Output_${chr}/eQTL.*.summary.txt.gz | awk 'FNR==1 && NR!=1 { while (/^PValue/) getline; } 1 {print $1,$2,$3,$4,$5}' >> Merged_out/Batch.$batch.summary.txt

echo -e "batch $batch and Chr $chr is merged"
done

awk 'NR!=1{print $1}' /groups/umcg-gastrocol/tmp04/Metabolic_Project/Metabolics/Numeric/CD_numeric.metabolic.000${batch}.txt | sort | uniq > Merged_out/Features.${batch}.txt
cat Merged_out/Features.${batch}.txt | while read line
do
awk -v var="$line" '{if($5==var)print $1,$2,$3,$4}' Merged_out/Batch.$batch.summary.txt > Merged_out/$line.Batch.${batch}.summary.txt
done

echo -e "$line is merged"
rm Merged_out/Batch.$batch.summary.txt
rm Merged_out/Features.${batch}.txt 

done
