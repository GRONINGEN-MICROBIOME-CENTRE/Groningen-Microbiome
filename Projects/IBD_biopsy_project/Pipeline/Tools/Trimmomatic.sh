#!/bin/bash
#SBATCH --job-name=Trimmomatic.test           
#SBATCH --output=Trimmomatic.test.out
#SBATCH --error=Trimmomatic.test.err
#SBATCH --time=23:59:00
#SBATCH --mem=10gb

input="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Input/test"
run_name=$(basename $input)

mkdir -p /groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Trimmomatic/$run_name
output="/groups/umcg-weersma/tmp01/Shixian/Biopsy_16S/Trimmomatic/"

for file in $input/*

do

base=$(basename $file)
name=${base%_L001*}

echo $name >> $run_name.sample.list

done

ml Java
cat $run_name.sample.list | uniq | while read line

do

java -jar /home/umcg-hushixian/Trimmomatic-0.36/trimmomatic-0.36.jar PE \
  -phred33 /$input/$line\_L001_R1_001.fastq.gz /$input/$line\_L001_R2_001.fastq.gz \
   $output/$run_name/$line\_1_paired.fq.gz $output/$run_name/$line\_1_unpaired.fq \
   $output/$run_name/$line\_2_paired.fq.gz $output/$run_name/$line\_2_paired.fq \
   ILLUMINACLIP:/home/umcg-hushixian/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 \
   LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 HEADCROP:8 MINLEN:50

done


rm $run_name.sample.list
