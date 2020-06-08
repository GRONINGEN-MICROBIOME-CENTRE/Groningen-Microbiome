#!/bin/bash
#BATCH --job-name=Numeric.batchtmp.chrtmp
#SBATCH --error=Numeric.batchtmp.chrtmp.err
#SBATCH --output=Numeric.batchtmp.chrtmp.out
#SBATCH --mem=12gb
#SBATCH --time=2:00:00
#SBATCH --cpus-per-task=6


module load Java
java -XX:ParallelGCThreads=5 -Xmx10G -jar /groups/umcg-gastrocol/tmp04/Metabolic_Project/Settings/eqtl-mapping-pipeline-1.4nZ/eqtl-mapping-pipeline.jar --mode metaqtl --settings Marc_XML.xml

cd Output_chrtmp
zcat eQTLs.txt.gz | awk '{OFS="\t"}{print $1,$2,$3,$4,$5,$9,$10,$12,$13,$14,$18}' > eQTL.batchtmp.chrtmp.summary.txt

gzip eQTL.batchtmp.chrtmp.summary.txt
rm SNPQCLog.txt.gz
rm eQTLs.txt.gz
                 
