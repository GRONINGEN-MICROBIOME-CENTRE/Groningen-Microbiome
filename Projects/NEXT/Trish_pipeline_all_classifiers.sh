#!/bin/bash

for i in *_1.fq.gz; do 
	sid=${i%.fq.gz}

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
	echo "mkdir ./$sid/test_05/"  >> "$sid".sh
	echo "mkdir ./$sid/test_08/"  >> "$sid".sh
	echo "mkdir ./$sid/test_old/"  >> "$sid".sh
	echo "mkdir ./$sid/test_new/"  >> "$sid".sh
	echo "mkdir ./$sid/mOTUS/" >> "$sid".sh
	
	# Defining paths to various directories (Kraken2, Bracken, MetaPhlan 3 and MOTUs) 
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2" >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/Bracken" >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/Bracken/src/" >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_3/MetaPhlAn-3.0/metaphlan" >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/mOTUS/mOTUs_v2.5/mOTUs_v2" >> "$sid".sh

	#BAM to FASTQ (if your files are not in fast q format) 
	#echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$sid.bam F=./$sid/clean_reads/"$sid".fastq1 F2=./$sid/clean_reads/"$sid".fastq2" >> "$sid".sh
	
	#Clean 1: Remove human reads and triming
	echo "kneaddata --input ./$sid/clean_reads/"$sid"_1.fq.gz -t 6 -p 7 --input ./$sid/clean_reads/"$sid"_2.fq.gz -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$sid/clean_reads/ --log ./$sid/clean_reads/$sid.log" >> "$sid".sh
	
	#Clean 2: Extra cleaning step: Remove reads with 50 of k-mers matching the human genome or the UniVec database
	echo "kraken2 --db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/decon/ --use-names --confidence 0.5 --threads 5 --report ./$sid/"$sid"_kraken.txt --unclassified-out ./$sid/clean_reads/clean#.fq --paired ./$sid/clean_reads/"$sid"_kneaddata_paired_1.fastq ./$sid/clean_reads/"$sid"_kneaddata_paired_2.fastq > ./$sid/clean_reads/summary_kraken.txt" >> "$sid".sh
	
	# Classify using Kraken
	echo "kraken2 --db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/standard_DB/ --use-names --confidence 0.5 --threads 8 --report ./$sid/test_05/"$sid"_0.5_kraken.txt --paired ./$sid/clean_reads/clean_1.fq ./$sid/clean_reads/clean_2.fq > ./$sid/test_05/summary_kraken.txt" >> "$sid".sh
	echo "kraken2 --db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/standard_DB/ --use-names --confidence 0.8 --threads 8 --report ./$sid/test_08/"$sid"_0.8_kraken.txt --paired ./$sid/clean_reads/clean_1.fq ./$sid/clean_reads/clean_2.fq > ./$sid/test_08/summary_kraken.txt" >> "$sid".sh


	#Estimate using Bracken

	echo "bracken -d /groups/umcg-gastrocol/tmp03/metagenomic_tools/kraken2/standard_DB/ -i ./$sid/test_08/"$sid"_0.8_kraken.txt -o ./Bracken_reports/"$i"_bracken.txt -w ./Bracken_reports/"$sid"_sp -l S " >> "$sid".sh

	#echo "rm -r ./$sid/clean_reads/" >> "$sid".sh

	# mOTUS

	#echo ""
	echo "motus profile -f ./$sid/clean_reads/clean_1.fq -r ./$sid/clean_reads/clean_2.fq -o ./$sid/mOTUs/"$sid"_motus_profile.txt -t 8 -C parenthesis -not_renormalise_cami -g 3" >> "$sid".sh
	echo "motus profile -f ./$sid/clean_reads/clean_1.fq -r ./$sid/clean_reads/clean_2.fq -o ./$sid/mOTUs/"$sid"_motus_strict_profile.txt -t 8 -C parenthesis -not_renormalise_cami -g 6 -l 90" >> "$sid".sh
	echo "motus profile -f ./$sid/clean_reads/clean_1.fq -r ./$sid/clean_reads/clean_2.fq -o ./$sid/mOTUs/"$sid"_motus_count_profile.txt -p -t 8 -c -g 5 -y insert.scaled_counts" >> "$sid".sh
	#echo "rm -r ./$sid/clean_reads/" >> "$sid".sh


	echo "ml Biopython"  >> "$sid".sh

	echo "cat ./"$sid"/clean_reads/clean_1.fq >> ./"$sid"/clean_reads/merged.fastq" >> "$sid".sh

	echo "cat ./"$sid"/clean_reads/clean_2.fq >> ./"$sid"/clean_reads/merged.fastq" >> "$sid".sh

	# Metaphlan 

	echo "metaphlan.py ./"$sid"/clean_reads/merged.fastq --input_type fastq --unknown_estimation --bowtie2out ./$sid/metagenome.bowtie2.bz2 --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_3/MetaPhlAn-3.0/metaphlan/metaphlan_databases/mpa_v30_CHOCOPhlAn_201901.pkl --nproc 6 -o ./"$sid"/test_new/"$sid"_default.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

	echo "metaphlan.py ./"$sid"/metagenome.bowtie2.bz2 --input_type bowtie2out --legacy-output  --nproc 6 -o ./"$sid"/test_old/"$sid"_old.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

	#echo "metaphlan.py ./"$sid"/metagenome.bowtie2.bz2 --input_type bowtie2out --add_viruses --unknown_estimation --nproc 6 -o ./"$sid"/test_viruses/"$sid"_virus.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

	#echo "metaphlan.py ./"$sid"/metagenome.bowtie2.bz2 --input_type bowtie2out --CAMI_format_output --unknown_estimation --nproc 6 -o ./"$sid"/test_CAMI/"$sid"_CAMI.txt --tmp_dir ./"$sid"/clean_reads/" >> "$sid".sh

	#echo "rm -r ./$sid/clean_reads/" >> "$sid".sh
done

