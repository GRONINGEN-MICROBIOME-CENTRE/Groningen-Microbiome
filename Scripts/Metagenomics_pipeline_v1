Metagenomics protocol 
===================

#Creator: Arnau Vich
#Pipeline used for taxa and functional characterization.
#Paths follow the configurations of Boxy cluster of the genetic department of the UMCG
#Year: 2016

1.Load tools and set variables for SLURM environment / Boxy cluster UMCG
-------------------------------------------------------------------------

```
#!/bin/bash
#SBATCH --job-name=$SAMPLE_ID_metagenomes
#SBATCH --error=$SAMPLE_ID.err
#SBATCH --output=$SAMPLE_ID.out
#SBATCH --mem=70gb
#SBATCH --time=5:59:00
#SBATCH --cpus-per-task=6
```

```
module load picard 
module load Python/3.4.1-foss-2015b
module load Bowtie2
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/
export PATH=$PATH:/home/umcg-avich/.local/bin
mkdir ./$SAMPLE_ID
mkdir ./$SAMPLE_ID/filtering_data/
mkdir ./$SAMPLE_ID/clean_reads/
mkdir ./$SAMPLE_ID/metaphlan/
mkdir ./$SAMPLE_ID/humann2/
mkdir ./$SAMPLE_ID/DUDes/
```

2.QC step: 
---------------

Currently our samples are sequenced in the Broad, so the the BAM files are already QCed by removing low quality reads and trimming adapters

Convert BAM to FASTQ with PICARD:
```
java -jar ${EBROOTPICARD}/picard.jar SamToFastq I=$SAMPLE_ID.bam F=./$SAMPLE_ID/filtering_data/$SAMPLE_ID.fastq1 F2=./$SAMPLE_ID/filtering_data/$SAMPLE_ID.fastq2
```

Remove reads mapping the human genome and quality filtering with Kneaddata:
```
kneaddata --input ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1.fq -t 6 -p 7 --input ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_2.fq -db /groups/umcg-gastrocol/tmp03/metagenomic_tools/kneaddata-0.5.4/Homo_sapiens_Bowtie2_v0.1/ --output ./$SAMPLE_ID/filtering_data/ --log ./$SAMPLE_ID/clean_reads/$SAMPLE_ID.log
```

Reorganize data and remove intermediate files:
```
cat ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1_kneaddata_paired_1.fastq > ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_kneaddata_merged.fastq 
cat ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_1_kneaddata_paired_2.fastq >> ./$SAMPLE_ID/filtering_data/$SAMPLE_ID_kneaddata_merged.fastq 
mv ./$SAMPLE_ID/filtering_data/*kneaddata_paired_1.fastq ./$SAMPLE_ID/clean_reads/
mv ./$SAMPLE_ID/filtering_data/*kneaddata_paired_2.fastq ./$SAMPLE_ID/clean_reads/ 
mv ./$SAMPLE_ID/filtering_data/*kneaddata_merged.fastq ./$SAMPLE_ID/clean_reads/
rm -r ./$SAMPLE_ID/filtering_data/
```

3.Taxonomic classification: You can use Metaphlan2 or DUDes or both
----------------------------------------------------------------------

Metaphlan2 
```
echo Starting taxonomy classification using Metaphlan
metaphlan2.py ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_kneaddata_merged.fastq  --input_type multifastq --mpa_pkl /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200.pkl --nproc 6 -o ./$SAMPLE_ID/metaphlan/$SAMPLE_ID_metaphlan.txt --tmp_dir ./$SAMPLE_ID/clean_reads/
```

DUDes (2 steps)
```
echo Starting taxonomy profiling using DUDes
bowtie2 -x /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/custom_db/db_refseq_20052017 --no-unal --fast -p 6 -k 50 -1 ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_1_kneaddata_paired_1.fastq -2 ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_1_kneaddata_paired_2.fastq -S ./$SAMPLE_ID/DUDes/$SAMPLE_ID_output.sam 
python3 /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/DUDes.py -s ./$SAMPLE_ID/DUDes/$SAMPLE_ID_output.sam -d /groups/umcg-gastrocol/tmp03/metagenomic_tools/dudes_v0_07/custom_db/DUDES_refseq_db.npz -t 6 -m 50 -a 0.0005 -l strain -o ./$SAMPLE_ID/DUDes/$SAMPLE_ID 
```

4. Taxonomic classification using mOTUSv2
------------------------------------------

mOTUs2

Load dependencies

```
echo Starting taxonomy profiling using mOTUs2

ml SAMtools/1.5-foss-2015b
ml BWA/0.7.15-foss-2015b
ml Python
mkdir ./"$sid"/mOTUs
export PATH=\$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/mOTUS/mOTUs_v2
```

Profiling 

```
# Output mOTUs counts
motus profile -f ./"$sid"/"$sid"_1_kneaddata_paired_1.fastq -r ./"$sid"/"$sid"_1_kneaddata_paired_2.fastq -o ./"$sid"/"$sid"_counts_profile.txt -q -t 6 -c -g 3 -y insert.scaled_counts

# Output relative abundances

motus profile -f ./"$sid"/"$sid"_1_kneaddata_paired_1.fastq -r ./"$sid"/"$sid"_1_kneaddata_paired_2.fastq -o ./"$sid"/"$sid"_abundance_profile.txt -q -t 6 -g 3 -y insert.scaled_counts
```


5.Functional profiling / pathway identification
------------------------------------------------

You can run this part as a separate job, if so you can use the following parameters in the cluster: 
```
#!/bin/bash
#SBATCH --job-name=$SAMPLE_ID_metagenomes
#SBATCH --error=$SAMPLE_ID.err2
#SBATCH --output=$SAMPLE_ID.out2
#SBATCH --mem=40gb
#SBATCH --time=23:50:00
#SBATCH --cpus-per-task=6
module load picard 
module load Python/3.4.1-foss-2015b
module load Bowtie2
export PATH=$PATH:/groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/
export PATH=$PATH:/home/umcg-avich/.local/bin
```

Run Humann2, WARNING: Check Humann2 configuration!!! (We currently use Uniref90 + Chocophlan db)
```
echo Starting pathways prediction using Humann2
humann2 --input ./$SAMPLE_ID/clean_reads/$SAMPLE_ID_kneaddata_merged.fastq --output ./$SAMPLE_ID/humann2/ --taxonomic-profile ./$SAMPLE_ID/metaphlan/$SAMPLE_ID_metaphlan.txt --threads 6 --o-log ./$SAMPLE_ID/clean_reads/$SAMPLE_ID.full.humann2.log --remove-temp-output
mv ./$SAMPLE_ID/clean_reads/*.log ./$SAMPLE_ID/
```
