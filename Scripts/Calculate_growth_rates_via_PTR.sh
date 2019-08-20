#!/bin/bash

#Creator: Arnau Vich
#Year: 2015
#Software: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5087275/

for i in *bam; do 
	sid=${i%.bam}

	echo "#!/bin/bash" > "$sid".sh
	echo "#SBATCH --job-name="$sid"_metagenomes" >> "$sid".sh
	echo "#SBATCH --error=$sid.err" >> "$sid".sh
	echo "#SBATCH --output=$sid.out" >> "$sid".sh
	echo "#SBATCH --mem=60gb" >> "$sid".sh
	echo "#SBATCH --time=23:50:00" >> "$sid".sh
	echo "#SBATCH --cpus-per-task=3" >> "$sid".sh

	echo "ml picard" >> "$sid".sh
	echo "ml kneaddata" >> "$sid".sh

	echo "mkdir ./$sid/"  >> "$sid".sh
	echo "mkdir ./$sid/clean_reads/"  >> "$sid".sh
	echo "mkdir ./$sid/PTR/"  >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/GEM/" >> "$sid".sh
	echo "export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/" >> "$sid".sh

	echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$sid.bam F=./$sid/clean_reads/"$sid".fastq1 F2=./$sid/clean_reads/"$sid".fastq2" >> "$sid".sh
	echo "kneaddata --input ./$sid/clean_reads/"$sid".fastq1 -t 6 -p 7 --input ./$sid/clean_reads/"$sid".fastq2 -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$sid/clean_reads/ --log ./$sid/clean_reads/$sid.log" >> "$sid".sh

	echo "ml Python/2.7.9-foss-2015b" >> "$sid".sh
	echo "ml pandas/0.16.2-foss-2015b-Python-2.7.9" >> "$sid".sh
	echo "ml numpy/1.9.2-foss-2015b-Python-2.7.9" >> "$sid".sh
	echo "ml lmfit" >> "$sid".sh
	echo "ml dill" >> "$sid".sh

	echo "python /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/PTRC.py -db_path_name /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/index -pe -i1 /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/runs/$sid/clean_reads/*kneaddata_paired_1.fastq -i2 /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/runs/$sid/clean_reads/*kneaddata_paired_2.fastq -outfol /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/runs/$sid/PTR/ -m /groups/umcg-gastrocol/tmp03/metagenomic_tools/PTR/PTRC1.1/runs/$sid/PTR/mapping.map CA" >> "$sid".sh

	echo "rm -r ./$sid/clean_reads/" >> "$sid".sh
	echo "rm ./$sid/PTR/mapping.map" >> "$sid".sh
	echo "mv ./$sid/PTR/mapping.ptr ./$sid/PTR/$sid.ptr" >> "$sid".sh

done
