lld_genetic_input="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/HC/"
cd_genetic_input="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/CD/"
uc_genetic_input="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Genetics/Combined/UC/"

coupling="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Metabolics/"
trait="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Metabolics/Binary/"


for batch in {01..52}

do

mkdir -p /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch

cp /groups/umcg-gastrocol/tmp04/Metabolic_Project/MetaAnalysis_binary_template/Marc_XML.chr.tmp.xml /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/
cp /groups/umcg-gastrocol/tmp04/Metabolic_Project/MetaAnalysis_binary_template/MetaAnalysis.chr.tmp.sh /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/

for chr in {1..22}

do

mkdir -p /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Output_$chr
final_out="/groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Output_$chr/"

sed "s|final_out|${final_out}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.tmp.xml \
> /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|lld_genetic_input|${lld_genetic_input}\/trityper_${chr}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|cd_genetic_input|${cd_genetic_input}\/trityper_${chr}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|uc_genetic_input|${uc_genetic_input}\/trityper_${chr}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|lld_coupling|${coupling}\/LLD.coupling.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|cd_coupling|${coupling}\/CD.coupling.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|uc_coupling|${coupling}\/UC.coupling.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|lld_trait|${trait}\/CT_binary.metabolic.000${batch}.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|cd_trait|${trait}\/CD_binary.metabolic.000${batch}.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|uc_trait|${trait}\/UC_binary.metabolic.000${batch}.txt|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|lld_annot|${trait}\/CT_binary.metabolic.txt.annot|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|cd_annot|${trait}\/CD_binary.metabolic.txt.annot|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed -i "s|uc_annot|${trait}\/UC_binary.metabolic.txt.annot|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/Marc_XML.chr.${chr}.xml

sed "s|batchtmp|${batch}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/MetaAnalysis.chr.tmp.sh \
> /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/MetaAnalysis.chr.${chr}.sh

sed -i "s|chrtmp|${chr}|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/MetaAnalysis.chr.${chr}.sh
sed -i "s|Marc_XML.xml|Marc_XML.chr.${chr}.xml|g" /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/MetaAnalysis.chr.${chr}.sh

echo -e "batch $batch ---- chr $chr is done"

done

rm /groups/umcg-gastrocol/tmp04/Metabolic_Project/Binary_meta/Batch.$batch/*tmp.*

done
~         
