##if you have any questions please contact: lianminchen@yeah.net
##the panphlan including two parts, first part for creating genome scanning file (present or absent of genes) from MGS for certain species, second part will merge all genome scanning files  (each samples one marker file) you show interest and do downstream analysis
##run panphlan for .bam files as an example in Boxy cluster

################Part1
##### step_1: extract samples name

for sample in *.bam

 do i=$(echo ${sample%.bam})

##### step_2: basic settings
  echo "#!/bin/bash" >> ./$i.panphlan.sh
  echo "#SBATCH --job-name=$i.panphlan" >> ./$i.panphlan.sh
  echo "#SBATCH --error=$i.panphlan.err" >> ./$i.panphlan.sh
  echo "#SBATCH --output=$i.panphlan.out" >> ./$i.panphlan.sh
  echo "#SBATCH --mem=32gb" >> ./$i.panphlan.sh
  echo "#SBATCH --time=59:00:00" >> ./$i.panphlan.sh
  echo "#SBATCH --cpus-per-task=6" >> ./$i.panphlan.sh

##### step_3: path
  echo "export PATH=\$PATH:/home/umcg-lchen/bin/" >> ./$i.panphlan.sh
  echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/" >> ./$i.panphlan.sh

##### step_4: create files  
  echo "mkdir -p ./$i" >> ./$i.panphlan.sh
  echo "mkdir -p ./$i/filtering_data/" >> ./$i.panphlan.sh
  echo "mkdir -p ./$i/clean_reads/" >> ./$i.panphlan.sh

##### step_5: picard 
  echo "module load picard" >> ./$i.panphlan.sh
  echo "module load Python/2.7.11-foss-2015b" >> ./$i.panphlan.sh
  echo "echo Starting Picard" >> ./$i.panphlan.sh
  echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=/groups/umcg-dag3/tmp03/LLD_followup_microbiome/panphlan_analysis/data/$i.bam F=./$i/filtering_data/$i.fastq1 F2=./$i/filtering_data/$i.fastq2" >> ./$i.panphlan.sh
  echo "echo Picard finished" >> ./$i.panphlan.sh

##### step_6: kneaddata
  echo "module load Python/2.7.11-foss-2015b" >> ./$i.panphlan.sh
  echo "module load kneaddata" >> ./$i.panphlan.sh
  echo "kneaddata --input ./$i/filtering_data/$i.fastq1 -t 6 -p 7 --input ./$i/filtering_data/$i.fastq2 -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$i/filtering_data/ --log ./$i/clean_reads/$i.log" >> ./$i.panphlan.sh
  echo "echo kneaddata finished" >> ./$i.panphlan.sh
  
##### step_7: move files
  echo "echo Moving files" >> ./$i.panphlan.sh
  echo "cat ./$i/filtering_data/$i\_kneaddata_paired_1.fastq > ./$i/filtering_data/$i\_kneaddata_merged.fastq" >> ./$i.panphlan.sh
  echo "cat ./$i/filtering_data/$i\_kneaddata_paired_2.fastq >> ./$i/filtering_data/$i\_kneaddata_merged.fastq" >> ./$i.panphlan.sh
  echo "mv ./$i/filtering_data/*kneaddata_merged.fastq ./" >> ./$i.panphlan.sh
  echo "rm -r ./$i/filtering_data/" >> ./$i.panphlan.sh
  echo "rm -r ./$i/" >> ./$i.panphlan.sh
  echo "rm -r ./$i.bam" >> ./$i.panphlan.sh

##### step_8: fastq to specific strain gene in each metagenomic sample
  echo "module load Python/2.7.11-foss-2015b" >> ./$i.panphlan.sh
  echo "module load SAMtools/1.9-foss-2015b" >> ./$i.panphlan.sh
  echo "module load Bowtie2/2.3.4.1-foss-2015b" >> ./$i.panphlan.sh
  echo "/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/panphlan_map.py -c rtorques16 -i ./$i\_kneaddata_merged.fastq --i_bowtie2_indexes /groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan_DB -o ../map_results/rtorques16/$i.csv --verbose" >> ./$i.panphlan.sh

  echo "/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/panphlan_map.py -c sparasanguinis16 -i ./$i\_kneaddata_merged.fastq --i_bowtie2_indexes /groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan_DB -o ../map_results/sparasanguinis16/$i.csv --verbose" >> ./$i.panphlan.sh

  echo "/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/panphlan_map.py -c svestibularis16 -i ./$i\_kneaddata_merged.fastq --i_bowtie2_indexes /groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan_DB -o ../map_results/svestibularis16/$i.csv --verbose" >> ./$i.panphlan.sh

  echo "/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/panphlan_map.py -c fprausnitzii16 -i ./$i\_kneaddata_merged.fastq --i_bowtie2_indexes /groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan_DB -o ../map_results/fprausnitzii16/$i.csv --verbose" >> ./$i.panphlan.sh

  echo "/groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan/panphlan/panphlan_map.py -c ssalivarius16 -i ./$i\_kneaddata_merged.fastq --i_bowtie2_indexes /groups/umcg-gastrocol/tmp03/metagenomic_tools/panphlan_DB -o ../map_results/ssalivarius16/$i.csv --verbose" >> ./$i.panphlan.sh
  
done



################Part2
module load Python/2.7.11-foss-2015b
/groups/umcg-gastrocol/tmp04/metagenomic_tools/panphlan/panphlan/panphlan_profile.py -c r39BFAA -i ./ --o_dna result_gene_presence_absence_r39FAA.txt --i_bowtie2_indexes /groups/umcg-gastrocol/tmp04/metagenomic_tools/panphlan_DB --add_strains --verbose

