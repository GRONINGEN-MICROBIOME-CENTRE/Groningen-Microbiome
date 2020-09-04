
nohup sh -c 'for chr in {1..22}; do java -jar ~/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i Metabolic.IBD."${chr}" -I PLINK_BED -O TRITYPER -o trityper_"${chr}"; done' &
