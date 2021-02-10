Author: Trishla Sinha, Febuary 2021

## In progress 

#!/bin/bash

for i in *_1.fq.gz; do 
	sid=${i%_1.fq.gz}

	echo "#!/bin/bash" > "$sid".sh
	echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
	echo "#SBATCH --error=$sid.err" >> "$sid".sh
	echo "#SBATCH --output=$sid.out" >> "$sid".sh
	echo "#SBATCH --mem=50gb" >> "$sid".sh
	echo "#SBATCH --time=05:50:00" >> "$sid".sh
	echo "#SBATCH --cpus-per-task=8" >> "$sid".sh

	#Loading required modules 
	echo "ml picard" >> "$sid".sh
	echo "ml kneaddata" >> "$sid".sh
	echo "ml BioPerl"  >> "$sid".sh
	echo "module load Bowtie2"  >> "$sid".sh
	
	# Making the relavent directories to store all outputs 
	echo "mkdir ./$sid/"  >> "$sid".sh
	echo "mkdir ./$sid/clean_reads/"  >> "$sid".sh
	echo "mkdir ./$sid/test_new/"  >> "$sid".sh
	echo "mkdir ./$sid/consensus_markers"  >> "$sid".sh
	
	# Defining path to directory 
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_3/MetaPhlAn-3.0/metaphlan" >> "$sid".sh

	#BAM to FASTQ (if your files are not in fast q format) 
	#echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$sid.bam F=./$sid/clean_reads/"$sid".fastq1 F2=./$sid/clean_reads/"$sid".fastq2" >> "$sid".sh
	
	#Clean 1: Remove human reads and triming
	echo "kneaddata --input "$sid"_1.fq.gz -t 6 -p 7 --input "$sid"_2.fq.gz -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$sid/clean_reads/ --log ./$sid/clean_reads/$sid.log" >> "$sid".sh
	
	#Clean 2: Extra cleaning step: Remove reads with 50 of k-mers matching the human genome or the UniVec database
	echo "kraken2 --db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/decon/ --use-names --confidence 0.5 --threads 5 --report ./$sid/"$sid"_kraken.txt --unclassified-out ./$sid/clean_reads/clean#.fq --paired ./$sid/clean_reads/"$sid"_1_kneaddata_paired_1.fastq ./$sid/clean_reads/"$sid"_1_kneaddata_paired_2.fastq > ./$sid/clean_reads/summary_kraken.txt" >> "$sid".sh


	echo "cat ./"$sid"/clean_reads/clean_1.fq >> ./"$sid"/clean_reads/merged.fastq" >> "$sid".sh

	echo "cat ./"$sid"/clean_reads/clean_2.fq >> ./"$sid"/clean_reads/merged.fastq" >> "$sid".sh

	# Metaphlan 

	echo "metaphlan.py ./"$sid"/clean_reads/merged.fastq --input_type fastq --unknown_estimation --bowtie2out ./$sid/metagenome.bowtie2.bz2 --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_3/MetaPhlAn-3.0/metaphlan/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.pkl --nproc 6 -o ./"$sid"/test_new/"$sid"_default.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

	
	echo "rm -r ./$sid/clean_reads/*.fastq" >> "$sid".sh
	echo "rm -r ./$sid/clean_reads/*.fq" >> "$sid".sh
  
  
done
