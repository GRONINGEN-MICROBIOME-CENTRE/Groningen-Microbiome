#!/bin/bash
#SBATCH --job-name=s02.dist
#SBATCH --error=s02.dist.err
#SBATCH --output=s02.dist.out
#SBATCH --mem=5gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

mash=/data/p290416/opt/Anaconda3/bin/mash
wd=/data/p290416/project/SNV_benchmark/mash

mkdir -p $wd/s02.dist.res

for i in /data/p290416/project/SNV_benchmark/test_genome_anno/SEQS/*.fa
do
	$mash dist $wd/s01.sketch.res/merged.msh $i
done


