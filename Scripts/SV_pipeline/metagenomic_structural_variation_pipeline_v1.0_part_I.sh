#################################################
# Title:  Metagenomic structural variation pipeline v1.0 - part I
# Date:   Nov. 2019
# Author: Daoming Wang, wangdaoming94@outlook.com
# Cite:   Zeevi D, Korem T, Godneva A, Bar N, Kurilshikov A, Lotan-Pompan M, Weinberger A, Fu J, Wijmenga C, Zhernakova A & Segal E (2019) Structural variation in the gut microbiome associates with host health. Nature 568, 43–48.
#################################################
# Usage: 
# sh metagenomic_structural_variation_pipeline_v1.0_part_I.sh
# cd $wd/$job_pre.sbatch
# for i in *.sh;do sbatch $i;done
#################################################
# Note:
# Beware of the changes of software and database paths
#################################################

wd=/groups/umcg-lld/tmp04/umcg-dwang/SVs/lld # absolute path of work space
input_list=$wd/input/LLD_s10.list # a list contains all absolute paths of the raw sequencing data (*.bam), see 'lld_s10.list' as a refernce
job_pre=s01.SVs_perSample.LLD_s10 # a job prefix you prefer
keanddata_db=/groups/umcg-lld/tmp04/umcg-dwang/database/keaddata-0.5.4 # database of keaddata

# rm -r $wd/$job_pre.sbatch
# rm -r $wd/$job_pre.res

mkdir -p  $wd/$job_pre.sbatch
mkdir -p  $wd/$job_pre.res

cat $input_list | while read line
do 
	sample_id=$(basename $line .bam)
	oo="$wd/$job_pre.sbatch/$job_pre.$sample_id.sh"

	echo "#!/bin/bash" > $oo
	echo "#SBATCH --job-name=$job_pre.$sample_id" >> $oo
	echo "#SBATCH --error=$job_pre.$sample_id.err" >> $oo
	echo "#SBATCH --output=$job_pre.$sample_id.out" >> $oo
	echo "#SBATCH --mem=40gb" >> $oo
	echo "#SBATCH --time=24:00:00" >> $oo
	echo "#SBATCH --cpus-per-task=6" >> $oo
	echo "#SBATCH --constraint=\"tmp04\"" >> $oo
	echo "#SBATCH --export=NONE" >> $oo
	echo "#SBATCH --get-user-env=L" >> $oo

#	>> step 1: reads cleaning <<
	echo "echo Step 1 - reads cleaning started on \`date\`" >> $oo

	echo "module load picard" >> $oo
	echo "module load Python/2.7.11-foss-2015b" >> $oo
	echo "module load Bowtie2/2.3.4.1-foss-2015b" >> $oo

	echo "mkdir -p $wd/$job_pre.res/clean" >> $oo
	echo "mkdir -p $wd/$job_pre.res/clean/$sample_id" >> $oo
	echo "mkdir -p $wd/$job_pre.res/clean/$sample_id/filtering_data/" >> $oo
	echo "mkdir -p $wd/$job_pre.res/clean/$sample_id/clean_reads/" >> $oo
  
	echo "echo 1.1 Picard started on \`date\`" >> $oo
	echo "java -jar \${EBROOTPICARD}/picard.jar SamToFastq I=$wd/rawdata/$sample_id.bam F=$wd/$job_pre.res/clean/$sample_id/filtering_data/${sample_id}_1.fastq F2=$wd/$job_pre.res/clean/$sample_id/filtering_data/${sample_id}_2.fastq" >> $oo
	echo "echo 1.1 Picard finished on \`date\`" >> $oo

  	echo "echo 1.2 kneaddata started on \`date\`" >> $oo
	echo "kneaddata -i $wd/$job_pre.res/clean/$sample_id/filtering_data/${sample_id}_1.fastq -i $wd/$job_pre.res/clean/$sample_id/filtering_data/${sample_id}_2.fastq -db $keanddata_db -t 6 -p 7 -o $wd/$job_pre.res/clean/$sample_id/filtering_data/ --log $wd/$job_pre.res/clean/$sample_id/clean_reads/$sample_id.log" >> $oo
	echo "echo 1.2 kneaddata finished finished on \`date\`" >> $oo
	echo "echo Step 1 - reads cleaning finished on \`date\`" >> $oo

#	>> step 2: icra <<
	echo "echo Step 2 - icra started on \`date\`" >> $oo
	echo "module load Cython/0.24.1-foss-2015b-Python-2.7.11" >> $oo
	echo "mkdir -p $wd/$job_pre.res/icra" >> $oo
	echo "python /groups/umcg-lld/tmp04/umcg-dwang/software/SGVFinder/src/ICRA_cmd.py $wd/$job_pre.res/icra $wd/$job_pre.res/clean/$sample_id/filtering_data/${sample_id}_1_kneaddata_paired --pe " >> $oo
	echo "echo Step 2 - icra finished on \`date\`" >> $oo

#	>> step 3: SGVF_PerFile <<
	echo "echo Step 3 - SGVF_PerFile started on \`date\`" >> $oo
	echo "mkdir -p $wd/$job_pre.res/sgvf_perfile/$sample_id" >> $oo
	echo "python /groups/umcg-lld/tmp04/umcg-dwang/software/SGVFinder/src/SGVF_PerFile_cmd.py $wd/$job_pre.res/icra/${sample_id}_1_kneaddata_paired.jsdel $wd/$job_pre.res/sgvf_perfile/$sample_id.jsdel 100 --x_coverage 0.01 --rate_param 10" >> $oo # KEY command! you can change the parameters to make it suitable with your project 
	echo "echo Step 3 - SGVF_PerFile finished on \`date\`" >> $oo

#	>> step 4: remove intermediate files <<
	echo "rm -rf $wd/$job_pre.res/clean/$sample_id" >> $oo
	echo "rm $wd/$job_pre.res/icra/${sample_id}_1_kneaddata_paired.jsdel $wd/$job_pre.res/icra/${sample_id}_1_kneaddata_paired.pmp" >> $oo
done
