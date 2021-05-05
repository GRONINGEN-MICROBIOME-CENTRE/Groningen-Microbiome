#!/bin/sh


###############################
###########PART I ############
###############################


#Run in a directory with Kneaddata cleaned FQ files, after concatenation of Forward and reverse reads

for sample in *_kneaddata_merged.fastq

 do i=$(echo ${sample%_kneaddata_merged.fastq})

  echo "#!/bin/bash" > ./$i.strainphlan.sh
  echo "#SBATCH --job-name=$i.strainphlan" >> ./$i.strainphlan.sh
  echo "#SBATCH --error=$i.strainphlan.err" >> ./$i.strainphlan.sh
  echo "#SBATCH --output=$i.strainphlan.out" >> ./$i.strainphlan.sh
  echo "#SBATCH --mem=40gb" >> ./$i.strainphlan.sh
  echo "#SBATCH --time=15:00:00" >> ./$i.strainphlan.sh
  echo "#SBATCH --cpus-per-task=6" >> ./$i.strainphlan.sh

 #Module loading; this is specific of Calculon server (Groningen)
 #Binaries were downloaded from https://www.dropbox.com/sh/m4na8wefp53j8ej/AAAsshd6cVhuZAlaMKdVLTa0a/bin?dl=0&subfolder_nav_tracking=1
  echo "ml PythonPlus/2.7.11-foss-2015b-v17.06.1" >> ./$i.strainphlan.sh
  echo 'export PATH=$PATH:/groups/umcg-lld/tmp04/umcg-sandreusanchez/Oral_microbiome_PD_MGS/Strain_resolved_analysis/Run_Strainphlan/Binaries'  >> ./$i.strainphlan.sh
  #SAM and BCF are actually not needed with the binaries, but is fine while running MP. Not fine for StrainPhlan. 
  SAM=/apps/software/SAMtools/0.1.19-foss-2015b/bin/samtools
  Meta=/apps/software/MetaPhlAn/2.7.2-foss-2015b/metaphlan2.py
  markers=/apps/software/MetaPhlAn/2.7.2-foss-2015b/strainphlan_src/sample2markers.py
  BCF=/apps/software/BCFtools/1.7-foss-2015b/bin/bcftools
###### Run metaphlan  to generate read alignments to the markers
 echo "python $Meta --bowtie2db /groups/umcg-gastrocol/tmp03/metagenomic_tools/metaphlan_2/db_v20/mpa_v20_m200 ./$i\_kneaddata_merged.fastq --mpa_pkl /apps/software/MetaPhlAn/2.7.2-foss-2015b/databases/mpa_v20_m200.pkl --input_type multifastq --nproc 6 -s ./$i.sam.bz2 --bowtie2out ./$i.bowtie2_out.bz2 -o ./$i.profile" >> ./$i.strainphlan.sh

###### Generate sample-specific markers (by a majority rule)
  echo "python $markers --ifn_samples ./$i.sam.bz2 --input_type sam --output_dir markers_$i/ --nprocs 6" >> ./$i.strainphlan.sh

###Node 014 is messed up in cluster
  sbatch --exclude=umcg-node014 ./$i.strainphlan.sh
done


###############################
###########PART II ############
###############################


#Run this part once the jobs from the first part are completed (!)

ml PythonPlus/2.7.11-foss-2015b-v17.06.1
export PATH=$PATH:/groups/umcg-lld/tmp04/umcg-sandreusanchez/Oral_microbiome_PD_MGS/Strain_resolved_analysis/Run_Strainphlan/Binaries
#Step 1 get Clades available from the markers generated in the first loop
python /apps/software/MetaPhlAn/2.7.2-foss-2015b/strainphlan.py  --mpa_pkl /apps/software/MetaPhlAn/2.7.2-foss-2015b/databases/mpa_v20_m200.pkl --ifn_samples markers*/*.markers --output_dir Clades --nprocs_main 5 --print_clades_only > Clades/all_clades_fullname.txt

#Step 2 Use strainphlan to generate a multiple sequence analysis of the markes extracted in I, on the clades extracted from Step1
while read line; do

First=${line:0:1}
if [ $First == "(" ]; then
arrIN=(${line//,/ })
arrIN2=${arrIN#"('"}
line=${arrIN2%"'"}
fi

echo Clades_$line\.sh
echo "#!/bin/bash" > Clades_$line\.sh
echo "#SBATCH --job-name=$line.strainphlan" >> Clades_$line\.sh
echo "#SBATCH --error=$line.strainphlan.err" >> Clades_$line\.sh
echo "#SBATCH --output=$line.strainphlan.out" >> Clades_$line\.sh
echo "#SBATCH --mem=20gb" >> Clades_$line\.sh
echo "#SBATCH --time=02:00:00" >> Clades_$line\.sh
echo "#SBATCH --cpus-per-task=5" >> Clades_$line\.sh

echo 'ml PythonPlus/2.7.11-foss-2015b-v17.06.1' >> Clades_$line\.sh
echo 'export PATH=$PATH:/groups/umcg-lld/tmp04/umcg-sandreusanchez/Oral_microbiome_PD_MGS/Strain_resolved_analysis/Run_Strainphlan/Binaries' >> Clades_$line\.sh
echo "python /apps/software/MetaPhlAn/2.7.2-foss-2015b/strainphlan.py  --mpa_pkl /apps/software/MetaPhlAn/2.7.2-foss-2015b/databases/mpa_v20_m200.pkl --ifn_samples markers*/*.markers --output_dir Clades/$line --nprocs_main 5 --clades $line --keep_alignment_files"  >> Clades_$line\.sh

sbatch --exclude=umcg-node014 Clades_$line\.sh
rm Clades_$line\.sh
done <Clades/all_clades_fullname.txt





