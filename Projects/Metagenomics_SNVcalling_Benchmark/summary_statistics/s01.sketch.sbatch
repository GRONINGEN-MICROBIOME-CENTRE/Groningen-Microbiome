#!/bin/bash
#SBATCH --job-name=s01.sketch
#SBATCH --error=s01.sketch.err
#SBATCH --output=s01.sketch.out
#SBATCH --mem=5gb
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --export=NONE
#SBATCH --get-user-env=L

mash=/data/p290416/opt/Anaconda3/bin/mash
wd=/data/p290416/project/SNV_benchmark/mash

mkdir -p $wd/s01.sketch.res

$mash sketch  -s 10000 -k 17 -o $wd/s01.sketch.res/merged /data/p290416/project/SNV_benchmark/test_genome_anno/SEQS/*.fa

