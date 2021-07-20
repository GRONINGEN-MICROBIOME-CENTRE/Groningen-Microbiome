
## Strainphlan 2.0 for version 3.0 please check: 
## For questions please contact: lianminchen@yeah.net/t.sinha@rug.nl
##the strainphlan including two parts, first part for creating marker file from MGS, second part will merge all maker files (each samples one marker file) you show interest and do phylogenetic analysis
##run strainphlan for .bam files as an example in Boxy cluster

###part1
# step_1: extract samples name

for sample in *_.bam

 do i=$(echo ${sample%_.bam})

# step_2: basic settings
  echo "#!/bin/bash" >> ./$i.strainphlan.sh
  echo "#SBATCH --job-name=$i.strainphlan" >> ./$i.strainphlan.sh
  echo "#SBATCH --error=$i.strainphlan.err" >> ./$i.strainphlan.sh
  echo "#SBATCH --output=$i.strainphlan.out" >> ./$i.strainphlan.sh
  echo "#SBATCH --mem=30gb" >> ./$i.strainphlan.sh
  echo "#SBATCH --time=59:00:00" >> ./$i.strainphlan.sh
  echo "#SBATCH --cpus-per-task=6" >> ./$i.strainphlan.sh
echo "##if you have any questions please contact: lianminchen@yeah.net The strainphlan including two parts, first part for creating marker file from MGS, second part will merge all maker files (each samples one marker file) you show interest and do phylogenetic analysis; Run strainphlan for Trishla's MGS" >> ./$i.strainphlan.sh

# step_3: path  
  echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/" >> ./$i.strainphlan.sh
  echo "export PATH=\$PATH:/home/umcg-lchen/bin/" >> ./$i.strainphlan.sh
  echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/strainphlan_src/" >> ./$i.strainphlan.sh
  echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/" >> ./$i.strainphlan.sh

# step_4: create files  
  echo "mkdir -p ./$i" >> ./$i.strainphlan.sh
  echo "mkdir -p ./$i/filtering_data/" >> ./$i.strainphlan.sh
  echo "mkdir -p ./$i/clean_reads/" >> ./$i.strainphlan.sh

# step_4.1: move files

 echo "echo gz file" >> ./$i.strainphlan.sh
 echo "gunzip $i\_1.fq.gz" >> ./$i.strainphlan.sh
 echo "gunzip $i\_2.fq.gz" >> ./$i.strainphlan.sh
 echo "mv ./$i\_1.fq ./$i/filtering_data/" >> ./$i.strainphlan.sh
 echo "mv ./$i\_2.fq ./$i/filtering_data/" >> ./$i.strainphlan.sh

# step_5: picard 
 echo "module load picard" >> ./$i.strainphlan.sh
 echo "module load Python/2.7.11-foss-2015b" >> ./$i.strainphlan.sh
 echo "echo Starting Picard" >> ./$i.strainphlan.sh
 echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=/groups/umcg-dag3/tmp03/LLD_followup_microbiome/lld_fup/$i.bam F=./$i/filtering_data/$i.fastq1 F2=./$i/filtering_data/$i.fastq2" >> ./$i.strainphlan.sh
 echo "echo Picard finished" >> ./$i.strainphlan.sh

# step_6: kneaddata
  echo "module load Python/2.7.11-foss-2015b" >> ./$i.strainphlan.sh
  echo "module load kneaddata" >> ./$i.strainphlan.sh
  echo "kneaddata --input ./$i/filtering_data/$i\_1.fq -t 6 -p 7 --input ./$i/filtering_data/$i\_2.fq -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$i/filtering_data/ --log ./$i/clean_reads/$i.log" >> ./$i.strainphlan.sh
  echo "echo kneaddata finished" >> ./$i.strainphlan.sh

# step_7: move files
  echo "echo Moving files" >> ./$i.strainphlan.sh
  echo "cat ./$i/filtering_data/$i\_1_kneaddata_paired_1.fastq > ./$i/filtering_data/$i\_kneaddata_merged.fastq" >> ./$i.strainphlan.sh
  echo "cat ./$i/filtering_data/$i\_1_kneaddata_paired_2.fastq >> ./$i/filtering_data/$i\_kneaddata_merged.fastq" >> ./$i.strainphlan.sh
  echo "mv ./$i/filtering_data/*kneaddata_merged.fastq ./$i/clean_reads/" >> ./$i.strainphlan.sh
  echo "rm -r ./$i/filtering_data/" >> ./$i.strainphlan.sh

# step_8: metaphlan
  echo "module load Python/2.7.11-foss-2015b" >> ./$i.strainphlan.sh
  echo "module load MetaPhlAn" >> ./$i.strainphlan.sh
  echo "module load Pysam/0.9.0-foss-2015b-Python-2.7.11" >> ./$i.strainphlan.sh
  echo "module load numpy" >> ./$i.strainphlan.sh
  echo "python /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/metaphlan2.py --bowtie2db /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200 ./$i/clean_reads/$i\_kneaddata_merged.fastq --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl --input_type multifastq --nproc 6 -s ./$i/clean_reads/$i.sam.bz2 --bowtie2out ./$i/clean_reads/$i.bowtie2_out.bz2 -o ./$i/clean_reads/$i.profile" >> ./$i.strainphlan.sh
  echo "echo Metaphlan sam.bz2 finished" >> ./$i.strainphlan.sh

# step_9: strainphlan
  echo "python /groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/strainphlan_src/sample2markers.py --ifn_samples ./$i/clean_reads/$i.sam.bz2 --input_type sam --output_dir ./$i/clean_reads --nprocs 6 &> ./$i/clean_reads/log.txt --samtools_exe /home/umcg-lchen/bin/samtools  --bcftools_exe /home/umcg-lchen/bin/bcftools" >> ./$i.strainphlan.sh
  echo "echo strainphlan markers finished" >> ./$i.strainphlan.sh
  echo "rm -r ./$i/clean_reads/*.fastq" >> ./$i.strainphlan.sh
done


###part2
#!/bin/bash
#SBATCH --job-name=strainphlan
#SBATCH --error=strain_claades.err
#SBATCH --output=strain_clades.out
#SBATCH --mem=100gb
#SBATCH --time=59:00:00
#SBATCH --cpus-per-task=10

export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/strainphlan_src/
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/
export PATH=$PATH:/home/umcg-lchen/bin/
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/

module load MetaPhlAn
module load PythonPlus/2.7.11-foss-2015b-v17.06.1

python /groups/umcg-gastrocol/tmp03/metagenomic_tools/strainphlan_new_version/strainphlan.py --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl --ifn_samples /groups/umcg-dag3/tmp03/trishla_strain/step_2/*.markers --output_dir /groups/umcg-dag3/tmp03/trishla_strain/strain_result/ --nprocs_main 10 --print_clades_only > /groups/umcg-dag3/tmp03/trishla_strain/strain_result/all_clades_fullname.txt
echo clades finished
