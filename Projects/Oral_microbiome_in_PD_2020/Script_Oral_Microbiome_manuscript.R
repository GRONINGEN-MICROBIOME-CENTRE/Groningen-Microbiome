setwd("~/../Resilio Sync/Transfer/Collaborations/Nijmegen/Final_version/")
source(file = "Functions_Oral_Microbiome_manuscript.R")
source("~/../Resilio Sync/Transfer/Collaborations/Debby/PPI/Data/ancom_scripts/ancom_v2.1.R")

set.seed(2020)

###############################################
################DATA LOAD######################
###############################################
#####################
#Abundance data MP2##
#####################
matrix_taxonomy = read_tsv("../Merged_batch.tsv")
matrix_taxonomy[-1,] -> matrix_taxonomy
matrix_taxonomy %>% select(-"Pipelines_feature_counts") -> matrix_taxonomy
###########Separate on phylogenetic levels#############
Taxon = vector()
list_split = strsplit(matrix_taxonomy$ID, "\\|")
for (i in seq(1:length(list_split))){
  element =  list_split[[i]]
  Taxon = c(Taxon, element[length(element)])
}
matrix_taxonomy %>% mutate(ID = Taxon) -> matrix_taxonomy
#Species <  Genus < Family < Order < Class < Phylum
matrix_taxonomy %>% filter(grepl("t__",ID)) -> strain_matrix ; matrix_taxonomy %>% filter(grepl("s__",ID)) -> species_matrix ; matrix_taxonomy %>% filter(grepl("g__",ID)) -> genus_matrix ; matrix_taxonomy %>% filter(grepl("f__",ID)) -> family_matrix ; matrix_taxonomy %>% filter(grepl("o__",ID)) -> order_matrix ; matrix_taxonomy %>% filter(grepl("c__",ID)) -> class_matrix ;  matrix_taxonomy %>% filter(grepl("p__",ID)) -> phylum_matrix
#######################
###Metadatada##########
#######################
sample_metadata = read_delim("Metadata2.csv", delim =";")
sample_metadata %>% mutate(StudieID = str_replace(sample_metadata$StudieID,"PI*","D")) %>% mutate(StudieID = str_replace(StudieID,"D0","D")) -> sample_metadata
sample_metadata %>% mutate(ID_patient = StudieID ,StudieID = paste0(sample_metadata$StudieID,"P"), Source="pocket")-> metadata_P  ; sample_metadata %>% mutate(ID_patient = StudieID, StudieID = paste0(sample_metadata$StudieID,"S"), Source="saliva")-> metadata_S  ; sample_metadata = rbind(metadata_P,metadata_S)
sample_metadata %>% mutate(Tissue_n_group = paste(Group,Source)) -> sample_metadata
sample_metadata %>% mutate(Disease= `Periodontitis score (continious)`) -> sample_metadata
#########################
####Humann2 abundance####
#########################
pathway_abundance = read_tsv("../Merged_batch_pathabundance.tsv") %>% filter(! grepl("\\|",`# Pathway`) )
pathway_abundance %>% select(-"Pipelines_merged_Abundance") -> pathway_abundance

#Immuno markers
read_tsv("Immuno_cells.txt") -> immuno_cell

######################
######################
#Exploration of data##
######################
######################

############################
####Taxonomy independent####
############################

###Contamination
contamination_stats(File = "../Reads_contaminated.txt.txt", metadata = sample_metadata)
#No signiciant difference in the number of reads between saliva an pocket
#No significant effect of tissue source on contamination
#Significant effect between thte total number of reads and the percentage of contamination

##Taxonomic distribution
distribution_taxonomy()
#Try to colour by groups

##########################
####Taxonomy dependent####
##########################

################
###FUNCTIONS####
################
Do_ancom = function(Taxonomy, sample_metadata){
  Results_ancom_final = tibble()
  for (GROUP_n in 1:length(unique(sample_metadata$Tissue_n_group))){
    for (GROUP2_n in (GROUP_n+1):length(unique(sample_metadata$Tissue_n_group))){
      GROUP = unique(sample_metadata$Tissue_n_group)[GROUP_n]
      GROUP2 = unique(sample_metadata$Tissue_n_group)[GROUP2_n]
      if (! paste(c(GROUP,GROUP2),collapse="-") %in% c("0 pocket-1 pocket", "0 saliva-1 saliva","1 pocket-0 pocket", "1 saliva-0 saliva" )){ next }
      if (is.na(GROUP2)){ break }
      print(paste(c(GROUP,GROUP2),collapse="-"))
      sample_metadata %>% filter(Tissue_n_group %in% c(GROUP, GROUP2)) -> metadata_groups
      Perform_ancom(filter(Taxonomy, ID %in% metadata_groups$StudieID), metadata_groups) %>% mutate(TEST = paste(c(GROUP,GROUP2),collapse="-")) -> Results_ancom
      Results_ancom_final = rbind(Results_ancom_final,Results_ancom)
    }
  }
  return(Results_ancom_final)
}
Analysis_of_taxonomic_level = function(Taxa_matrix, Directory, sample_metadata){
  print("Preprocessing data")
  #Get data ready
  Preprocess_f(Taxa_matrix, paste(c(Directory,"/Summary.csv"), collapse="")) -> Taxa_matrix2
  print("Checking top species")
  #Check top species in each Group
  Richness_data = tibble()
  for (GROUP in unique(sample_metadata$Tissue_n_group)){
    Sub_data = filter(Taxa_matrix2, ID %in% filter(sample_metadata, Tissue_n_group== GROUP)$StudieID)
    apply(select(Sub_data, -ID) , 2, FUN = function(x){ sum(x) > 0 }) -> columns_present
    names(columns_present[columns_present == T]) -> To_keep
    Sub_data %>% select(To_keep) -> Sub_data
    
    print(c(GROUP, dim(Sub_data)))
    Sub_data  %>% summarise_if(is.numeric, median) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% arrange(desc(V1)) %>% head(n=10) %>% print()
    apply(Sub_data , 1, FUN = function(x){ length(as_vector(x)[as_vector(x) > 0]) }) %>% as.vector() -> Rich
    Richness_data = rbind(Richness_data, select(mutate(Sub_data, Group = GROUP, Richness = Rich), Group, Richness))
  }
  for (GROUP in unique(sample_metadata$Source)){
    print(GROUP)
    Sub_data = filter(Taxa_matrix2, ID %in% filter(sample_metadata, Source== GROUP)$StudieID)
    Sub_data  %>% summarise_if(is.numeric, median) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% arrange(desc(V1)) %>% head(n=10) %>% print()
    Sub_data  %>% summarise_if(is.numeric, median) %>% t() %>% as.data.frame() %>% rownames_to_column() %>% arrange(desc(V1)) %>% head(n=10) %>% 
      ggplot() + geom_bar(aes(x=rowname, y =V1), stat="identity") + ggtitle(GROUP) +theme_bw()+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) -> FIG
    print(FIG)
  }

  print("Taxa richness (number of taxa per group)")
  ggplot(Richness_data, aes(x=Group, y=Richness))+  geom_boxplot() + geom_point() + theme_bw() -> Richness_plot ; print(Richness_plot)
  dunn.test(g = Richness_data$Group, x = Richness_data$Richness)

  
  
  print("Alpha and Beta diversity analysis")
  #Alpha and Beta diversity
  Diversity_analysis(Taxa_matrix2,sample_metadata) -> Figures_diversity
  Figures_diversity[[5]] %>% ggplot(aes(x=Pheno, y=R2, fill=Tissue, col=Pvalue<0.05)) + geom_bar(stat="identity",position = "dodge") + coord_flip() + ylab("Phenotype")+ 
    theme_bw() +  scale_fill_manual(values = wesanderson::wes_palette("Royal1")) + scale_color_manual(values = c("white","black")) -> Var_explained
  SAVE(PLOT = Var_explained, NAME = "manuscript/Plots/PERMANOVA_R2.pdf")
  
  #No differences in diversity betwen saliva/pocket or healthy/sick
  print("Filtering taxa")
  #Differential abundance
  Filter_abundance(Taxa_matrix2, threshold = 0.2, min_abundance = 0.0001) -> Taxa_matrix2
  
  print("Preparing for differential abundance")
  sample_metadata %>% arrange(StudieID) -> sample_metadata
  sample_metadata %>% filter(StudieID %in% Taxa_matrix2$ID) -> sample_metadata
  Taxa_matrix2 %>% arrange(ID) -> Taxa_matrix2 
  
  
  ####################ANCOM
  print("Using ANCOM function")
  

  Do_ancom(Taxa_matrix2, sample_metadata) -> Results_ancom_final
  Results_ancom_final %>% filter(TEST == "0 pocket-1 pocket") %>% filter(detected_0.6 == T) -> DA_pocket_s
  Results_ancom_final %>% filter(TEST == "0 saliva-1 saliva") %>% filter(detected_0.6 == T) -> DA_saliva_s
  return(list(Taxa_matrix2, Results_ancom_final))
  
  
}
study_taxon = function(Meta, abundance, taxa, Stats){
  Meta %>% filter(StudieID %in% abundance$ID) %>% arrange(StudieID) -> Meta
  abundance %>% arrange(ID) %>% select(taxa) %>% as_vector() %>% as.vector() -> abundance_taxa
  tibble(Group = Meta$Group, Source = Meta$Source, Tissue_n_group = Meta$Tissue_n_group, Abundance = abundance_taxa) -> abundance_taxon
  abundance_taxon %>% ggplot(aes(x=Tissue_n_group, y = Abundance)) + geom_boxplot() +geom_point() + theme_bw() + labs(title=taxa) -> FIG
  print(FIG)
  abundance_taxon %>% group_by(Tissue_n_group) %>% summarise(n()) %>% print()
  filter(abundance_taxon,Abundance > 0) %>% group_by(Tissue_n_group) %>% summarise(n()) %>% print()
  Stats %>% filter(taxa_id == taxa) %>% print()
}
Correlate_with_phenotypes = function(Ancom_results, sample_metadata, abundance_matrix, tissue, cells="no"){
  filter(Ancom_results,detected_0.6 == T)$taxa_id -> To_check
  if (length(To_check) == 0){ print("No significant ANCOM results")}
  sample_metadata %>% filter(StudieID %in% abundance_matrix$ID) %>% filter(Source == tissue) %>% arrange(StudieID) -> sample_metadata
  abundance_matrix %>% arrange(ID) %>% filter(ID %in% sample_metadata$StudieID) %>% select(To_check)  -> abundance_matrix

  if (cells == "no"){ 
    Variables = c("IL-1Ra (pg/mL)", "IL-6 (pg/mL)","hsCRP (ug/mL)", "WBC (10^6/mL)", "Periodontium right (SUVmean)")
    sample_metadata = sample_metadata
  } else {
    Variables =  c("NEUT#(10^3/uL)","LYMPH#(10^3/uL)","MONO#(10^3/uL)", "NEUT%(%)", "LYMPH%(%)", "MONO%(%)")  
    immuno_cell %>% filter(StudieID %in% sample_metadata$StudieID) %>% arrange(ID_patient) %>% mutate(`Age (years)` = sample_metadata$`Age (years)`) -> sample_metadata
  }
  
  Final_result = tibble()
  for (Species in colnames(abundance_matrix)){
    abundance_matrix %>% select(Species) %>% as_vector %>% as.vector() -> Regressor
    for (Dependent in Variables){
      Age = sample_metadata$`Age (years)`
      dependent = sample_metadata%>% select(Dependent) %>% as_vector() %>% as.vector()
      Data = tibble(Regressor = Regressor, Dependent = dependent, Age = Age)
      summary(lm(scale(dependent)[,1] ~ scale(Regressor)[,1] + Age)) -> Model
      Result = tibble(Beta = Model$coefficients[2], Pval=Model$coefficients[11], Biomarker=Dependent, Taxa =Species)
      rbind(Result, Final_result) -> Final_result
    }
  }
  return(Final_result)
}
SAVE = function(NAME, PLOT, Width = 6.3, Height = 6.3, Dpi = 300, Scale = 1){
ggsave(filename = NAME, PLOT,
       width = Width, height = Height, dpi = Dpi, units = "in", scale = Scale)
}


######################
######################
#####Analysis########
######################
######################

#############
###Strain####
#############
M = tibble()
for (ID in strain_matrix$ID){
  S = str_split(ID, "\\|")[[1]]
  strain = S[length(S)]
  Species = S[length(S)-1]
  M = rbind(M, tibble(strain, Species))
}
M %>% group_by(Species) %>% summarise(N = n()) %>% arrange(desc(N)) %>% filter(N>1) -> Multiple_Strain
Analysis_of_taxonomic_level(strain_matrix, Directory = "Straub_results", sample_metadata = sample_metadata) -> strain_results
strains_matrix2 <- strain_results[[1]] ; ancom_strains <- strain_results[[2]]
lapply(ancom_strains$taxa_id , FUN = function(x){ str_split(x, "\\|")[[1]][7] } )  %>% as_vector() -> Spc
ancom_strains %>% mutate(S = Spc)  %>% filter(S %in% Multiple_Strain$Species) %>% arrange(TEST)

ancom_strains %>% filter(TEST == "0 pocket-1 pocket") %>% filter(detected_0.6==T) -> pocket_strains
ancom_strains %>% filter(TEST == "0 saliva-1 saliva") %>% filter(detected_0.6==T) -> saliva_strains


###################
#######Species#####
###################

Analysis_of_taxonomic_level(species_matrix, Directory = "Species_results", sample_metadata = sample_metadata) -> species_results
species_matrix2 <- species_results[[1]] ; ancom_species <- species_results[[2]]

ancom_species %>% filter(TEST == "0 pocket-1 pocket") %>% filter(detected_0.6==T) -> pocket_species
write_tsv(path="manuscript/Supplement/Species_ancom.tsv",x = ancom_species)


if ("0 pocket-1 pocket" %in% ancom_species$TEST){ ancom_species %>% mutate(TEST = ifelse(TEST == "0 pocket-1 pocket","1 pocket-0 pocket", ifelse(TEST == "0 saliva-1 saliva","1 saliva-0 saliva" ,TEST))) -> ancom_species}
ancom_species %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)
ancom_species %>% filter(detected_0.6 == T) %>% arrange(desc(W)) 
#####################################################
###########Systemic inflammation markers#############
#####################################################
Do_heatmap = function(DF, OUT, size=0.7){
  DF %>% mutate(log10Pval = -log10(Pval)) %>% select(Biomarker, Taxa, log10Pval) %>% spread(Taxa, log10Pval) %>%filter(! is.na(Biomarker)) %>% as.data.frame() %>% column_to_rownames(var = "Biomarker") -> wide_results
  DF %>% select(Biomarker, Taxa, Beta) %>% spread(Taxa, Beta) %>% filter(! is.na(Biomarker)) %>% as.data.frame() %>% column_to_rownames(var = "Biomarker") -> wide_results2
  
  qvals <- matrix("", nrow=nrow(wide_results), ncol=ncol(wide_results))
  qvals[wide_results > -log10(0.05)]= "*"
  
  my_palette = colorRampPalette(c("steelblue3","white","red"))(n = 101)
  pdf(file=OUT,height = 10,width = 5)
  gplots::heatmap.2(as.matrix(wide_results2), scale="none",dendrogram="both",cellnote = qvals, cexRow = 0.7,cexCol = size,trace="none",key = T,notecol="black", col=my_palette, keysize=2,notecex=2.0, 
                    distfun=function(x) dist(x, method="euclidean"),Colv=T, hclustfun=function(x) hclust(x, method="ward.D2"),margins=c(10,10)) -> Fig_pheno
  dev.off()
}

Correlate_with_phenotypes(dplyr::filter(ancom_species, TEST == "1 pocket-0 pocket"), sample_metadata, mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="pocket") -> Inflammation_pocket
Correlate_with_phenotypes(dplyr::filter(ancom_species, TEST == "1 saliva-0 saliva"), sample_metadata, mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="saliva") -> Inflammation_saliva
mutate(Inflammation_pocket, Tissue="pocket") -> Inflammation_pocket
mutate(Inflammation_pocket, Tissue="saliva") -> Inflammation_saliva
write_tsv(path="manuscript/Supplement/Species_Inflammation.tsv",x = rbind(Inflammation_pocket,Inflammation_saliva))

Do_heatmap(DF = Inflammation_pocket ,OUT= "manuscript/Supplement/Figures/Pocket_species.pdf")
Do_heatmap(DF = Inflammation_saliva ,OUT= "manuscript/Supplement/Figures/Saliva_species.pdf")

Inflammation_pocket_all %>% mutate(FDR = p.adjust(Pval, "fdr")) %>% ggplot(aes(x=Beta, y=-log10(Pval), col=FDR<0.05)) + geom_point() + facet_wrap(~Biomarker, scales = "free") + theme_bw()

make_pheatmap = function(Inflammation){
  Inflammation %>% mutate(log10Pval = -log10(Pval)) %>% dplyr::select(c(Biomarker, Taxa, log10Pval)) %>% spread(Biomarker, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Taxa") -> wide_results
  Inflammation %>% select(Biomarker, Taxa, Beta) %>% spread(Biomarker, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Taxa") -> wide_results2
  wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
  pheatmap::pheatmap(wide_results, display_numbers = wide_results_Annotation, fontsize_number = 10, fontsize_col = 9, fontsize_row = 5 ) -> Fig_0
  print(Fig_0)
}
make_pheatmap(Inflammation_pocket)
make_pheatmap(Inflammation_saliva)


#Check which blood cells are correlated
Inflammation_pocket_all %>% filter( Biomarker == "WBC (10^6/mL)") %>% filter(Pval < 0.05) -> Species_blood_count
Correlate_with_phenotypes(mutate(Species_blood_count, detected_0.6=T, taxa_id=Taxa), sample_metadata, mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="pocket", "cells") -> cells_pocket
Inflammation_saliva_all %>% filter( Biomarker == "WBC (10^6/mL)") %>% filter(Pval < 0.05) -> Species_blood_count_saliva
Correlate_with_phenotypes(Ancom_results = mutate(Species_blood_count_saliva, detected_0.6=T, taxa_id=Taxa), sample_metadata = sample_metadata, abundance_matrix =mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="saliva", cells="cells") -> cells_saliva
######################################
#Species of interest in bibliography##
######################################

study_taxon(taxa= "s__Treponema_denticola", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Tannerella_forsythia", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Anaeroglobus_geminatus", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Porphyromonas_gingivalis", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)

study_taxon(taxa= "s__Anaeroglobus_geminatus", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Capnocytophaga_gingivalis", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Capnocytophaga_granulosa", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Parvimonas_micra", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Parvimonas_unclassified", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Propionibacterium_propionicum", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)


study_taxon(taxa= "s__Bulleidia_extructa", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Olsenella_uli", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Olsenella_unclassified", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)

study_taxon(taxa= "s__Bulleidia_extructa", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Eubacterium_brachy", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Eubacterium_infirmum", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Eubacterium_yurii", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Filifactor_alocis", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Tannerella_forsythia", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)
study_taxon(taxa= "s__Treponema_denticola", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)

colnames(species_matrix2)[grepl("Streptococcus", colnames(species_matrix2))]
Neisseria
Fusobacterium
Lactobacillus

study_taxon(taxa= "s__Streptococcus_parasanguinis", abundance= species_matrix2, Meta= sample_metadata, Stats= ancom_species)

####################################
###Comparison of specific groups####
####################################

#Red cluster
compositional_transformation(select(species_matrix2, -ID)) %>% mutate(ID = species_matrix2$ID) -> transformed_species_matrix 
transformed_species_matrix %>% mutate(Red_complex = s__Porphyromonas_gingivalis + s__Tannerella_forsythia + s__Treponema_denticola) -> transformed_species_matrix
transformed_species_matrix %>% arrange(ID) -> transformed_species_matrix ; sample_metadata %>% filter(StudieID %in% transformed_species_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2
transformed_species_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_species_matrix
wilcox.test(Red_complex ~ Group, filter(transformed_species_matrix, Source=="pocket"))
wilcox.test(Red_complex ~ Group, filter(transformed_species_matrix, Source=="saliva"))
transformed_species_matrix %>% ggplot(aes(x=factor(Source), y= Red_complex, col=factor(Group))) + geom_boxplot(outlier.shape = NA) + geom_sina()  + theme_bw() +
  scale_color_manual(values = wes_palette("Royal1"),labels =c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab("Location") + ylab("CLR-abundance Red complex") -> Fig_Redcomplex

ancom_species %>% mutate(taxa_id = as.character(taxa_id)) -> ancom_species
rbind(ancom_species,c("Red_complex", 200, T, T, T, T, 5, "0 pocket-1 pocket"),c("Red_complex", 200, T, T, T, T, 5, "0 saliva-1 saliva")) -> ancom_species2

Correlate_with_phenotypes(dplyr::filter(ancom_species2, TEST == "0 pocket-1 pocket"), sample_metadata, transformed_species_matrix, tissue ="pocket") %>% filter(Taxa == "Red_complex") %>% arrange(Pval)
Correlate_with_phenotypes(dplyr::filter(ancom_species2, TEST == "0 saliva-1 saliva"), sample_metadata, transformed_species_matrix, tissue ="saliva") %>% filter(Taxa == "Red_complex") %>% arrange(Pval)


#Propionibacterium +Parvimonas -> s__Propionibacterium_propionicum, "s__Parvimonas_micra",  "s__Parvimonas_unclassified"
transformed_species_matrix %>% mutate(ClusterB = s__Parvimonas_micra + s__Parvimonas_unclassified + s__Propionibacterium_propionicum) -> transformed_species_matrix
transformed_species_matrix %>% arrange(ID) -> transformed_species_matrix ; sample_metadata %>% filter(StudieID %in% transformed_species_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2
transformed_species_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_species_matrix
wilcox.test(ClusterB ~ Group, filter(transformed_species_matrix, Source=="pocket"))
wilcox.test(ClusterB ~ Group, filter(transformed_species_matrix, Source=="saliva"))
#transformed_species_matrix %>% ggplot(aes(x=factor(Group), y= ClusterB)) + geom_boxplot() + facet_wrap(~Source) + theme_bw()
transformed_species_matrix %>% ggplot(aes(x=factor(Source), y= ClusterB, col=factor(Group))) + geom_boxplot(outlier.shape = NA) + geom_sina()  + theme_bw() +
  scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab("Location") + ylab("CLR-abundance bacteria atherosclerotic plaque") -> Fig_atherosclerosis

rbind(ancom_species,c("ClusterB", 200, T, T, T, T, 5, "0 pocket-1 pocket"),c("ClusterB", 200, T, T, T, T, 5, "0 saliva-1 saliva")) -> ancom_species2

Correlate_with_phenotypes(dplyr::filter(ancom_species2, TEST == "0 pocket-1 pocket"), sample_metadata, transformed_species_matrix, tissue ="pocket") %>% filter(Taxa == "ClusterB") %>% arrange(Pval)
Correlate_with_phenotypes(dplyr::filter(ancom_species2, TEST == "0 saliva-1 saliva"), sample_metadata, transformed_species_matrix, tissue ="saliva") %>% filter(Taxa == "ClusterB") %>% arrange(Pval)


(Fig_Redcomplex + labs(title = "A")) + (Fig_atherosclerosis+labs(title = "B")+ theme(legend.position = "none")) + plot_layout(guides = "collect",ncol = 1) -> Combined_complexes
SAVE(PLOT = Combined_complexes, NAME = "manuscript/Plots/combined_complexes_abundance.pdf")
SAVE(PLOT = Combined_complexes, NAME = "manuscript/Plots/combined_complexes_abundance.pdf")


#Check if they correlate with systemic inflammation markers
Correlate_with_phenotypes(tibble(taxa_id = c("ClusterB", "Red_complex"), detected_0.6 =c(T,T)), abundance_matrix = transformed_species_matrix, sample_metadata = sample_metadata , tissue="pocket") -> Association_Clusters_inflammation
Association_Clusters_inflammation %>% arrange(Pval) -> Results_groups_pocket_inflammation
Correlate_with_phenotypes(tibble(taxa_id = c("ClusterB", "Red_complex"), detected_0.6 =c(T,T)), abundance_matrix = transformed_species_matrix, sample_metadata = sample_metadata , tissue="saliva") -> Association_Clusters_inflammation_saliva
Association_Clusters_inflammation_saliva %>% arrange(Pval) -> Results_groups_saliva_inflammation



#Cluster A (mild PD):
transformed_species_matrix %>% mutate(ClusterA = s__Campylobacter_concisus + s__Campylobacter_gracilis + s__Campylobacter_rectus+ s__Campylobacter_showae + 
                                        s__Corynebacterium_durum + s__Corynebacterium_matruchotii +
                                        s__Fusobacterium_nucleatum + s__Fusobacterium_periodonticum +
                                        s__Leptotrichia_buccalis +  s__Leptotrichia_hofstadii + 
                                        s__Leptotrichia_shahii+ s__Leptotrichia_unclassified + s__Leptotrichia_wadei +
                                        s__Tannerella_forsythia) -> transformed_species_matrix

transformed_species_matrix %>% arrange(ID) -> transformed_species_matrix ; sample_metadata %>% filter(StudieID %in% transformed_species_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2
transformed_species_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_species_matrix
wilcox.test(ClusterA ~ Group, filter(transformed_species_matrix, Source=="pocket"))
wilcox.test(ClusterA ~ Group, filter(transformed_species_matrix, Source=="saliva"))
transformed_species_matrix %>% ggplot(aes(x=factor(Source), y= ClusterA, col=factor(Group))) + geom_boxplot(outlier.shape = NA) + geom_sina()  + theme_bw() +
  scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab("Location") + ylab("CLR-abundance bacteria Cluster A") -> Fig_clusterA

#Cluster B (severe PD):

transformed_species_matrix %>% mutate(ClusterB_ade = Red_complex + s__Filifactor_alocis + s__Treponema_denticola+
                                        s__Treponema_lecithinolyticum+s__Treponema_maltophilum+s__Treponema_medium+
                                        s__Treponema_socranskii+s__Treponema_vincentii) -> transformed_species_matrix

transformed_species_matrix %>% arrange(ID) -> transformed_species_matrix ; sample_metadata %>% filter(StudieID %in% transformed_species_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2
transformed_species_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_species_matrix
wilcox.test(ClusterB_ade ~ Group, filter(transformed_species_matrix, Source=="pocket"))
wilcox.test(ClusterB_ade ~ Group, filter(transformed_species_matrix, Source=="saliva"))
transformed_species_matrix %>% ggplot(aes(x=factor(Source), y= ClusterB_ade, col=factor(Group))) + geom_boxplot(outlier.shape = NA) + geom_sina()  + theme_bw() +
  scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab("Location") + ylab("CLR-abundance bacteria Cluster B") -> Fig_clusterB


(Fig_clusterA + labs(title = "A")) + (Fig_clusterB+labs(title = "B")+ theme(legend.position = "none")) + plot_layout(guides = "collect",ncol = 1) -> Combined_complexes2
SAVE(PLOT = Combined_complexes2, NAME = "manuscript/Plots/combined_complexes_abundance_AB_S.png")
SAVE(PLOT = Combined_complexes2, NAME = "manuscript/Plots/combined_complexes_abundance_AB_S.pdf")



rbind(ancom_species,c("ClusterB_ade", 200, T, T, T, T, 5, "0 pocket-1 pocket"),c("ClusterB_ade", 200, T, T, T, T, 5, "0 saliva-1 saliva")) -> ancom_species2


Correlate_with_phenotypes(tibble(taxa_id = c("ClusterB_ade"), detected_0.6 =c(T,T)), abundance_matrix = transformed_species_matrix, sample_metadata = sample_metadata , tissue="pocket") -> Association_Clusters_inflammation

##Other clusters
transformed_species_matrix %>% mutate(Koren = s__Veillonella_atypica +s__Veillonella_dispar+s__Veillonella_parvula+
                                        s__Veillonella_parvula + s__Veillonella_unclassified + s__Streptococcus_anginosus+ s__Streptococcus_australis +s__Streptococcus_cristatus + s__Streptococcus_gordonii+               
                                        s__Streptococcus_infantis + s__Streptococcus_intermedius + s__Streptococcus_mitis_oralis_pneumoniae+
                                        s__Streptococcus_mutans + s__Streptococcus_oligofermentans + s__Streptococcus_parasanguinis + s__Streptococcus_salivarius +
                                        s__Streptococcus_sanguinis + s__Streptococcus_tigurinus + s__Streptococcus_vestibularis) -> transformed_species_matrix
transformed_species_matrix %>% arrange(ID) -> transformed_species_matrix ; sample_metadata %>% filter(StudieID %in% transformed_species_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2 ; transformed_species_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_species_matrix
wilcox.test(Koren ~ Group, filter(transformed_species_matrix, Source=="pocket")) ; wilcox.test(Koren ~ Group, filter(transformed_species_matrix, Source=="saliva"))
transformed_species_matrix %>% ggplot(aes(x=factor(Source), y= ClusterB_ade, col=factor(Group))) + geom_boxplot(outlier.shape = NA) + geom_sina()  + theme_bw() +
  scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab("Location") + ylab("CLR-abundance Koren") 



########################
####Other phenotypes####
########################
New_ancom = function(M, Meta, Pheno){
  select(M, -ID) %>% t() %>% as_tibble() %>% mutate(Taxa = colnames(select(M, -ID))) %>%
    `colnames<-`(c(M$ID, "Taxa")) %>% as.data.frame() %>% column_to_rownames("Taxa") -> S_DF
  List_processing <- feature_table_pre_process(S_DF, Meta, "StudieID", group_var = Pheno , out_cut = 0.05, zero_cut = 0.90, lib_cut=0, neg_lb=F)
  
  Meta_corrected = List_processing[[2]]
  Abundance_corrected = List_processing[[1]]
  
  ANCOM(Abundance_corrected, Meta_corrected, List_processing[[3]], Pheno, "fdr", "0.1") -> ANCOM_results
  as_tibble(ANCOM_results[[1]]) %>% arrange(desc(W)) -> Results_1
  as_tibble(ANCOM_results[[2]]) -> Results_2
  colnames(Results_2) = c("taxa_id", "CLR_mean_diff", "W", "Structural_0")
  left_join(Results_1, Results_2, by="taxa_id") -> Results_1
  
  return(Results_1)
}

Filter_abundance(species_matrix2, threshold = 0.2, min_abundance = 0.0001) %>% arrange(ID) -> Taxa_matrix2
sample_metadata %>% arrange(StudieID) %>% filter(StudieID %in% Taxa_matrix2$ID) -> sample_metadata_s
sample_metadata_s %>% mutate(PPDdeepest=`PPDdeepest (mm)`, PPDmean = `PPDmean (mm)`) -> sample_metadata_s

results_phenos_sp = tibble()
ancom_analysis = T
ancom_total = tibble()
#Do for each pheno and tissue
for (Tissue in c("pocket", "saliva")){
  for (Pheno in c("PPDdeepest", "PPDmean")){
    ##Select tissue
    sample_metadata_s %>% filter(Source == Tissue) %>% select(c("StudieID", Pheno)) -> Phenos
    if (ancom_analysis == T){
      #Do ancom
      New_ancom(Taxa_matrix2,Phenos, Pheno) -> Ancom_result
      ancom_total = rbind(ancom_total,mutate(Ancom_result,Pheno = Pheno, Tissue = Tissue))
    } else {
      #Do regression
      compositional_transformation(select(Taxa_matrix2, -ID)) %>% mutate(ID = Taxa_matrix2$ID) %>% filter(ID %in% Phenos$StudieID) %>% select(-ID) -> Taxa_matrix3
      Dependent = Phenos %>% select(Pheno) %>% as_vector() %>% as.vector()
      apply(Taxa_matrix3, 2, FUN= function(x){ summary(lm(Dependent ~ x)) -> Results; return(c(Results$coefficients[2], Results$coefficients[8]))}) -> Associations
      as.data.frame(Associations) %>% t() %>% as_tibble() %>% mutate(Taxa = colnames(as.data.frame(Associations))) %>% arrange(V2) %>% mutate(FDR = p.adjust(V2, "fdr")) %>% mutate(Source = Tissue, Phenotype = Pheno) -> Associations
      Associations %>% ggplot() + geom_point(aes(x= V1, y = -log10(V2), col = FDR<0.05)) + labs(title = paste(c(Tissue,"-",Pheno), collapse="")) + theme_bw() -> Plot_assoc
      print(Plot_assoc)
      PATH = paste(c("Species_results/Results_", Pheno, "_", Tissue, ".csv" ), collapse="")
      
      write_csv(x = Associations , path = PATH)
      results_phenos_sp = rbind(results_phenos_sp, mutate(Associations, Tissue = Tissue, Pheno = Pheno))
    }  
    #Associations %>% filter(FDR < 0.1) %>% mutate(taxa_id = Taxa) -> Assoc_i
    #Correlate_with_phenotypes(mutate(Assoc_i, detected_0.6=T), sample_metadata_s, mutate(Taxa_matrix3, ID=Phenos$StudieID), tissue ="pocket") -> Inflammation
    #Inflammation %>% arrange(Pval) -> Inflammation
    
    
  }
}

#results_phenos_sp %>% group_by(Source, Phenotype, FDR<0.05) %>% summarise(n())
ancom_total %>% filter(detected_0.6 == T) %>% filter( Structural_0 == "No")
write_tsv(path="manuscript/Supplement/Species_otherphenos.tsv",x = ancom_total)



################################################
###Other taxonomic levels...####################
################################################

#Genus
Analysis_of_taxonomic_level(genus_matrix, Directory = "Genus_results", sample_metadata = sample_metadata) -> genus_results
genus_matrix2 <- genus_results[[1]] ; ancom_genus <- genus_results[[2]]
ancom_genus %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)
study_taxon(taxa= "g__Propionibacterium", abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
study_taxon(taxa= "g__Parvimonas", abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
study_taxon(taxa= "g__Streptococcus",abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
study_taxon(taxa= "g__Neisseria",abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
study_taxon(taxa= "g__Fusobacterium",abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
study_taxon(taxa= "g__Lactobacillus",abundance= genus_matrix2, Meta= sample_metadata, Stats= ancom_genus)
#Family
Analysis_of_taxonomic_level(family_matrix, Directory = "Family_results", sample_metadata = sample_metadata) -> family_results
family_matrix2 <- family_results[[1]] ; ancom_family <- family_results[[2]]
ancom_family %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)
compositional_transformation(select(family_matrix2, -ID)) %>% mutate(ID = species_matrix2$ID) -> transformed_family_matrix 

transformed_family_matrix %>% mutate(Mirtra = f__Porphyromonadaceae + f__Micrococcaceae + f__Streptococcaceae) -> transformed_family_matrix
transformed_family_matrix %>% arrange(ID) -> transformed_family_matrix ; sample_metadata %>% filter(StudieID %in% transformed_family_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2 ; transformed_family_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_family_matrix
wilcox.test(Mirtra ~ Group, filter(transformed_family_matrix, Source=="pocket")) ; wilcox.test(Mirtra ~ Group, filter(transformed_family_matrix, Source=="saliva"))

#Order
Analysis_of_taxonomic_level(order_matrix, Directory = "Order_results", sample_metadata = sample_metadata) -> order_results
order_matrix2 <- order_results[[1]] ; ancom_order <- order_results[[2]]
ancom_order %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)
#Phylum
Analysis_of_taxonomic_level(phylum_matrix, Directory = "Order_results", sample_metadata = sample_metadata) -> phylum_results
phylum_matrix2 <- phylum_results[[1]] ; ancom_phylum <- phylum_results[[2]] ; mixed_phylum <- phylum_results[[3]]
ancom_phylum %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)
compositional_transformation(select(phylum_matrix2, -ID)) %>% mutate(ID = phylum_matrix2$ID) -> transformed_phylum_matrix 

transformed_phylum_matrix %>% mutate(Jonsson = p__Proteobacteria + p__Actinobacteria) -> transformed_phylum_matrix
transformed_phylum_matrix %>% arrange(ID) -> transformed_phylum_matrix ; sample_metadata %>% filter(StudieID %in% transformed_phylum_matrix$ID) %>% arrange(StudieID) -> sample_metadata_2 ; transformed_phylum_matrix %>% mutate(Group=sample_metadata_2$Group,Source=sample_metadata_2$Source,Medication_usage=sample_metadata_2$Medication_usage) ->transformed_phylum_matrix
wilcox.test(Jonsson ~ Group, filter(transformed_phylum_matrix, Source=="pocket")) ; wilcox.test(Jonsson ~ Group, filter(transformed_phylum_matrix, Source=="saliva"))
######


##########################
###Pathways###############
##########################
####ANCOM on clinical group#######

Prepare_pathways = function(pathway_abundance){
  new_i = vector()
  for (i in colnames(pathway_abundance)){
    if (i=="# Pathway"){ 
      i = "Pathway"
    }else{
      i = str_replace(string = i,pattern = "_merged_Abundance",replacement = "")
    }
    new_i = c(new_i,i)
  }
  colnames(pathway_abundance) = new_i
  pathway_abundance %>% mutate_if(is_numeric,make_relative) -> pathway_abundance
  new_i = colnames(pathway_abundance)
  Pathway = pathway_abundance$Pathway ; pathway_abundance %>% select(-Pathway) %>% t() -> pathway_abundance2 
  colnames(pathway_abundance2) = Pathway ; pathway_abundance2 %>% as_tibble() %>% mutate(ID = new_i[2:length(new_i)]) -> pathway_abundance2
  pathway_abundance2 %>% select(-c(UNMAPPED,UNINTEGRATED)) -> pathway_abundance2
  pathway_abundance2 %>% select(c("ID",colnames(select(pathway_abundance2,-ID)))) -> pathway_abundance2
  arrange(pathway_abundance2,ID) -> pathway_abundance2
  return(pathway_abundance2)
}
Prepare_pathways(pathway_abundance) -> pathway_abundance2
Diversity_analysis(pathway_abundance2,sample_metadata) -> Figures_diversity
Figures_diversity[[5]]
Filter_abundance(pathway_abundance2, threshold = 0.2, min_abundance = 0.001 ) -> pathway_abundance2
Do_ancom(pathway_abundance2, sample_metadata) -> Results_ancom_pathway
#write_tsv(Results_ancom_pathway, "Pathway_results.csv")
write_tsv(path="manuscript/Supplement/Pathway_ancom.tsv",x = Results_ancom_pathway)


Results_ancom_pathway %>% ggplot() + geom_point(aes(x=CLR_mean_diff, y=W, col=detected_0.6==T)) +
  theme_bw() + facet_wrap(~ TEST)

mutate(compositional_transformation(select(pathway_abundance2, -ID)), ID = pathway_abundance2$ID) -> pathway_abundance3

#######Inflammation markers############

Correlate_with_phenotypes(dplyr::filter(Results_ancom_pathway, TEST == "1 pocket-0 pocket"), sample_metadata,pathway_abundance3 , tissue ="pocket") -> Inflammation_pocket
Correlate_with_phenotypes(Ancom_results = dplyr::filter(Results_ancom_pathway, TEST == "1 saliva-0 saliva"), sample_metadata = sample_metadata,abundance_matrix = pathway_abundance3, tissue ="saliva") -> Inflammation_saliva

write_csv(path = "Pathways_targeted_pocket.csv",x=Inflammation_pocket)
write_csv(path = "Pathways_targeted_saliva.csv",x=Inflammation_saliva)

Do_heatmap(DF = Inflammation_pocket ,OUT= "manuscript/Supplement/Figures/Pocket_pathway.pdf", size = 0.3)


make_pheatmap(Inflammation_pocket)
make_pheatmap(Inflammation_saliva)

#########Other phenotypes#############
sample_metadata %>% arrange(StudieID) %>% filter(StudieID %in% pathway_abundance2$ID) -> sample_metadata_s
sample_metadata_s %>% mutate(PPDdeepest=`PPDdeepest (mm)`, PPDmean = `PPDmean (mm)`) -> sample_metadata_s
ancom_total_path = tibble()
#Do for each pheno and tissue
for (Tissue in c("pocket", "saliva")){
  for (Pheno in c("PPDdeepest", "PPDmean")){
    ##Select tissue
    sample_metadata_s %>% filter(Source == Tissue) %>% select(c("StudieID", Pheno)) -> Phenos
      #Do ancom
      New_ancom(pathway_abundance2,Phenos, Pheno) -> Ancom_result
      ancom_total_path = rbind(ancom_total_path,mutate(Ancom_result,Pheno = Pheno, Tissue = Tissue))
}
}
#results_phenos_sp %>% group_by(Source, Phenotype, FDR<0.05) %>% summarise(n())
ancom_total_path %>% filter(detected_0.6 == T) %>% filter( Structural_0 == "No")
write_tsv(path="manuscript/Supplement/Pathway_otherphenos.tsv",x = ancom_total_path)


##Other Phenotypes
"""
Filter_abundance(pathway_abundance2, threshold = 0.2, min_abundance = 0.0001) %>% arrange(ID) -> Taxa_matrix2
sample_metadata %>% arrange(StudieID) %>% filter(StudieID %in% Taxa_matrix2$ID) -> sample_metadata_s
sample_metadata_s %>% mutate(PPDdeepest=`PPDdeepest (mm)`, PPDmean = `PPDmean (mm)`) -> sample_metadata_s

results_phenos = tibble()
for (Tissue in c("pocket", "saliva", "all")){
  for (Pheno in c("PPDdeepest", "PPDmean")){
    if (Tissue == "all"){
      sample_metadata_s %>% select(c("StudieID", Pheno, "Source","ID_patient")) -> Phenos
    } else { sample_metadata_s %>% filter(Source == Tissue) %>% select(c("StudieID", Pheno)) -> Phenos }
    pathway_abundance3 %>% filter(ID %in% Phenos$StudieID) %>% select(-ID) -> Taxa_matrix3
    Dependent = Phenos %>% select(Pheno) %>% as_vector() %>% as.vector()
    
    if (Tissue == "all"){
      Source = Phenos %>% select("Source") %>% as_vector() %>% as.vector()
      apply(Taxa_matrix3, 2, FUN= function(x){INPUT = tibble(Dependent=Dependent, Source = Phenos$Source, Taxa=x, ID_patients = Phenos$ID_patient) ; tryCatch({summary(lme(Dependent ~ Taxa + Source, random=~ 1 | ID_patients, INPUT)) -> Results}, error = function(e){ return(c(NA,NA))}); Results = as.data.frame(Results$tTable)["Taxa",] ; return(c(Results$Value, Results$`p-value`))}) -> Associations
    } else {
      apply(Taxa_matrix3, 2, FUN= function(x){ summary(lm(Dependent ~ x)) -> Results; return(c(Results$coefficients[2], Results$coefficients[8]))}) -> Associations
    }
    as.data.frame(Associations) %>% t() %>% as_tibble() %>% mutate(Taxa = colnames(as.data.frame(Associations))) %>% arrange(V2) %>% mutate(FDR = p.adjust(V2, "fdr")) %>% mutate(Source = Tissue, Phenotype = Pheno) -> Associations
    Associations %>% ggplot() + geom_point(aes(x= V1, y = -log10(V2), col = FDR<0.05)) + labs(title = paste(c(Tissue,"-",Pheno), collapse="")) + theme_bw() -> Plot_assoc
    print(Plot_assoc)
    
    PATH = paste(c("Results_", Pheno, "_", Tissue, "_pathway.csv" ), collapse="")
    
    write_csv(x = Associations , path = PATH)
    results_phenos = rbind(results_phenos, mutate(Associations, Tissue = Tissue, Pheno = Pheno))
  }
}


results_phenos %>% group_by(Tissue, Pheno, FDR< 0.05) %>% summarise(n())
results_phenos %>% filter(Pheno == "PPDmean" & FDR < 0.05) %>% arrange(V2)
write_csv()
""" 
####



Inflammation_pocket %>% filter( Biomarker == "WBC (10^6/mL)") %>% filter(Pval < 0.05) -> Species_blood_count
Correlate_with_phenotypes(mutate(Species_blood_count, detected_0.6=T, taxa_id=Taxa), sample_metadata, mutate(compositional_transformation(select(pathway_abundance2, -ID)), ID = pathway_abundance2$ID), tissue ="pocket", "cells") -> cells_pocket
Inflammation_saliva %>% filter( Biomarker == "WBC (10^6/mL)") %>% filter(Pval < 0.05) -> Species_blood_count_saliva
Correlate_with_phenotypes(mutate(Species_blood_count_saliva, detected_0.6=T, taxa_id=Taxa), sample_metadata, mutate(compositional_transformation(select(pathway_abundance2, -ID)), ID = pathway_abundance2$ID), tissue ="saliva", "cells") -> cells_saliva


#########################
########################
##Figures.################
#######################
#########################
library(wesanderson)
library(patchwork)


#1. Phylum level composition
Preprocess_f(phylum_matrix, paste(c("Phylum_results","/Summary.csv"), collapse="")) -> Taxa_matrix2

metadata_phylum = arrange(filter(sample_metadata, StudieID %in% Taxa_matrix2$ID),StudieID)
#Check top species in each Group
Taxa_matrix2 %>% arrange(ID) %>% mutate(Location = metadata_phylum$Source, Group = as.factor(metadata_phylum$Group)) %>% gather(Taxa, Abundance, 2:12, factor_key=TRUE) %>%
  group_by(Group, Location, Taxa) %>% summarise(Median_abundance = median(Abundance), Mean_abundance=mean(Abundance), Sd_abundance=sd(Abundance)) %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) %>%
  arrange(Mean_abundance) %>% 
  ggplot(aes(x= reorder(Taxa, -Mean_abundance), y = Mean_abundance, fill=Clinical_Group)) + geom_bar(stat="identity",position = "dodge",col="black")+
  geom_errorbar(aes(ymin=Mean_abundance, ymax=Median_abundance+Sd_abundance), width=.2,position=position_dodge(.9)) +facet_wrap(~Location) + ylab("Mean realtive abundance (%)") + xlab("Phylum") +
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size= 5, face="italic")) + scale_fill_manual(values = wes_palette("Royal1")) -> Fig_phylum_composition

Preprocess_f(species_matrix, paste(c("Species_results","/Summary.csv"), collapse="")) -> Taxa_matrix2
Taxa_matrix2 %>% arrange(ID) %>% mutate(Location = metadata_phylum$Source, Group = as.factor(metadata_phylum$Group)) %>% gather(Taxa, Abundance, 2:236, factor_key=TRUE) -> Taxa_matrix3
Taxa_matrix3 %>% group_by(Taxa) %>% summarise(MA = median(Abundance)) %>% arrange(desc(MA)) %>% head(10) -> Top_species
Taxa_matrix3 %>% filter(Taxa %in% Top_species$Taxa) %>% group_by(Group, Location, Taxa) %>% summarise(Median_abundance = median(Abundance),Sd_abundance=sd(Abundance), Mean_abundance=mean(Abundance)) %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) %>%
  arrange(Mean_abundance) %>% 
  ggplot(aes(x= reorder(Taxa, -Mean_abundance), y = Mean_abundance, fill=Clinical_Group)) + geom_bar(stat="identity",position = "dodge", color="black")+ facet_wrap(~Location) +
  geom_errorbar(aes(ymin=Mean_abundance, ymax=Mean_abundance+Sd_abundance), width=.2,position=position_dodge(.9)) + ylab("Mean relative abundance (%)") + xlab("Species")+
  theme_bw() + theme(axis.text.x = element_text(angle = 90, hjust = 1, size= 5,face="italic")) + scale_fill_manual(values = wes_palette("Royal1")) -> Fig_sp_composition


(Fig_phylum_composition + labs(title = "A")) + (Fig_sp_composition+labs(title = "B")) + plot_layout(guides = "collect",ncol = 1) -> Combined_composition

SAVE(PLOT=Combined_composition, NAME = "manuscript/Plots/combined_composition.pdf")


#1. Reviewer's comments, do boxplots:
Preprocess_f(phylum_matrix, paste(c("Phylum_results","/Summary.csv"), collapse="")) -> Taxa_matrix2
Taxa_matrix2 %>% arrange(ID) %>% mutate(Location = metadata_phylum$Source, Group = as.factor(metadata_phylum$Group)) %>% gather(Taxa, Abundance, 2:12, factor_key=TRUE) %>%
  group_by(Group, Location, Taxa) %>% summarise(Median_abundance = median(Abundance), Mean_abundance=mean(Abundance), Sd_abundance=sd(Abundance)) %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) %>%
  arrange(Mean_abundance) -> Info_median_phylum
Taxa_matrix2 %>% arrange(ID) %>% mutate(Location = metadata_phylum$Source, Group = as.factor(metadata_phylum$Group)) %>% gather(Taxa, Abundance, 2:(dim(Taxa_matrix2)[2]), factor_key=TRUE) -> Abundances_table
left_join(Abundances_table, select(Info_median_phylum, c(Group, Location,Taxa, Median_abundance ))  ,by = c("Location", "Group", "Taxa") ) %>% drop_na() -> Abundances_table

Abundances_table %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) %>%
  arrange(Median_abundance) %>% filter(!Taxa %in% c("p__Ascomycota", "p__Tenericutes")) %>%
  ggplot(aes(x= reorder(Taxa, -Median_abundance), y = log10(Abundance), fill=Clinical_Group)) + geom_boxplot(col="black",outlier.shape = NA )+
  facet_wrap(~Location) + ylab("Log10(Relative abundance (%))") + xlab("Phylum") +
  theme_bw() + theme(axis.text.y = element_text(size= 10, face="italic")) + scale_fill_manual(values = wes_palette("Royal1")) +
  coord_flip() + geom_sina() -> Fig_phylum_composition_v2


#Species
Preprocess_f(species_matrix, paste(c("Species_results","/Summary.csv"), collapse="")) -> Taxa_matrix2
Taxa_matrix3 %>% filter(Taxa %in% Top_species$Taxa) %>% group_by(Group, Location, Taxa) %>% summarise(Median_abundance = median(Abundance),Sd_abundance=sd(Abundance), Mean_abundance=mean(Abundance)) %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) -> Info_median
Taxa_matrix2 %>% arrange(ID) %>% mutate(Location = metadata_phylum$Source, Group = as.factor(metadata_phylum$Group)) %>% gather(Taxa, Abundance, 2:(dim(Taxa_matrix2)[2]), factor_key=TRUE) -> Abundances_table
left_join(Abundances_table, select(Info_median, c(Group, Location,Taxa, Median_abundance ))  ,by = c("Location", "Group", "Taxa") ) %>% drop_na() -> Abundances_table

Abundances_table %>% mutate(Clinical_Group = ifelse(Group == 0, "PD mild", "PD severe")) %>%
  arrange(Median_abundance) %>% 
  ggplot(aes(x= reorder(Taxa, -Median_abundance), y = log10(Abundance), fill=Clinical_Group)) + geom_boxplot(col="black",outlier.shape = NA)+
  facet_wrap(~Location) + ylab("log10(Relative abundance (%))") + xlab("Species") +
  theme_bw() + theme(axis.text.y = element_text(size= 10, face="italic")) + scale_fill_manual(values = wes_palette("Royal1")) +
  coord_flip() + geom_sina() -> Fig_sp_composition_v2

(Fig_phylum_composition_v2 + labs(title = "A")) + (Fig_sp_composition_v2+labs(title = "B")) + plot_layout(guides = "collect",ncol = 1) -> Combined_composition_v2

SAVE(PLOT=Combined_composition_v2, NAME = "manuscript/Plots/combined_composition_v2.pdf", Height = 10, Width=10)





#ggsave(Fig_phylum_composition, filename = "manuscript/Plots/Phylum_composition.png")
#ggsave(Fig_sp_composition, filename = "manuscript/Plots/Sp_composition.png")
#ggsave(Combined_composition, filename = "manuscript/Plots/combined_composition.png",scale = 2)



#2 Alpha and Beta diversity + Heatmap

#Richness diversity
Taxa_matrix3 %>% mutate(Richness = ifelse(Abundance==0, 0, 1)) %>% group_by(ID,Location,Group) %>% summarise(Richness=sum(Richness)) %>%
ggplot(aes(x=Location, y=Richness, col=Group))+  geom_boxplot(outlier.shape = NA)+ geom_sina() + theme_bw() +
   scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") -> Richness_plot 


Prepare_phylobject = function(M, metadata, Tree=NA){
  
  #Prepare abundance object
  M %>% select(-ID) %>% t() %>% as_tibble() %>% mutate(label=colnames(select(M,-ID))) %>% `colnames<-`(c(M$ID, "label")) %>% select(c("label",M$ID)) -> trans_M
  trans_M %>% as.data.frame() %>% column_to_rownames("label") -> ASV_table
  otu_table(ASV_table, taxa_are_rows = TRUE) -> ASV_table
  #Prepare sample name object
  taxmat = matrix(rownames(ASV_table))
  rownames(taxmat) = rownames(ASV_table)
  SampleData = metadata %>% as.data.frame() %>% column_to_rownames("StudieID")
  SampleData = sample_data(SampleData)
  
  tax_table(taxmat) -> taxmat
  phyloseq::phyloseq(ASV_table, taxmat, SampleData) -> phylobjetc
  
  if (! is.na(Tree)){ phyloseq::merge_phyloseq(phylobjetc,Tree) -> phylobjetc }
  
  
  return(phylobjetc)
  
}
Prepare_phylobject(Taxa_matrix2,mutate(sample_metadata, Group = as.factor(Group))) -> phyloobject

ordinate(phyloobject,"PCoA", "bray") -> ord
plot_ordination(phyloobject, ord, color="Group", shape="Source") + theme_bw() + geom_point(size = 4) +
  stat_ellipse(type = "t", linetype = 2) + scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group", shape="Location") + theme_bw() + xlab("PCoA1(34.5%)") + ylab("PCoA2(16.2%") +
  theme(legend.title = element_text(size = 14),legend.text = element_text(size = 13)) -> PCoA_plot
#plot_ordination(phyloobject, ord, shape="Group", type="split")

Diversity_analysis(Taxa_matrix2,sample_metadata) -> Figures_diversity
Figures_diversity[[1]] + scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group", shape="Location")+ 
  theme(legend.title = element_text(size = 14),legend.text = element_text(size = 13))-> clustering_plot
#Figures_diversity[[4]] + scale_color_manual(values = wes_palette("Royal1"),labels = c("Mid PD", "Severe PD")) + labs(color = "Clinical_Group", shape="Location")  -> PCoA_plot
Figures_diversity[[2]] + scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + labs(color = "Clinical_Group") + xlab(label = "Location") + ylab(label = "Diversity") +
  scale_shape_discrete(labels = c("Mild PD", "Severe PD")) -> Shannon_plot
(Richness_plot + labs(title = "A")) + (Shannon_plot+labs(title = "B")) + plot_layout(guides = "collect",ncol = 1) -> Combined_richness
#ggsave(Combined_richness, filename = "manuscript/Plots/combined_richness.png",scale = 2)
SAVE(PLOT=Combined_richness, NAME="manuscript/Plots/combined_richness.pdf")

(PCoA_plot + labs(title = "A")) + (clustering_plot+labs(title = "B")+ theme(legend.position = "none")) + plot_layout(guides = "collect",ncol = 1) -> Combined_beta
#ggsave(Combined_beta, filename = "manuscript/Plots/combined_beta.png",scale = 4)
SAVE(PLOT=Combined_beta, NAME="manuscript/Plots/combined_beta.pdf", Width = 12.6 , Height = 12.6)

#3. Heatmap association with inflammation

make_pheatmap = function(Inflammation){
  Inflammation %>% mutate(log10Pval = -log10(Pval)) %>% dplyr::select(c(Biomarker, Taxa, log10Pval)) %>% spread(Biomarker, log10Pval) %>% as.data.frame() %>% column_to_rownames(var = "Taxa") -> wide_results
  Inflammation %>% select(Biomarker, Taxa, Beta) %>% spread(Biomarker, Beta) %>% as.data.frame() %>% column_to_rownames(var = "Taxa") -> wide_results2
  
  
  wide_results_Annotation = round(wide_results,3)
  wide_results_Annotation[wide_results_Annotation > -log10(0.05)] <- paste(wide_results_Annotation[wide_results_Annotation > -log10(0.05)], "*")
  
  breaksList = seq(-1, 1, by = 0.2)
  
  #wide_results_Annotation = wide_results2 ; wide_results_Annotation[wide_results_Annotation < "0"] <- "-" ;wide_results_Annotation[wide_results_Annotation > "0"] <- "+"
  pheatmap::pheatmap(wide_results2, display_numbers = wide_results_Annotation, fontsize_number = 13, fontsize_col = 13, fontsize_row = 13, color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(length(breaksList)), breaks = breaksList, tsize = 10) -> Fig_0
  return(Fig_0)
}
library(RColorBrewer)
sample_metadata %>% mutate(`hsCRP (ug/mL)` = as.vector(scale(`hsCRP (ug/mL)`)), `WBC (10^6/mL)` = as.vector(scale(`WBC (10^6/mL)`)), `Periodontium right (SUVmean)`=as.vector(scale(`Periodontium right (SUVmean)`)), `Periodontium (SUVmax)`= as.vector(scale(`Periodontium (SUVmax)`)), `IL-1Ra (pg/mL)` = as.vector(scale(`IL-1Ra (pg/mL)`)), `IL-6 (pg/mL)`=as.vector(scale(`IL-6 (pg/mL)`))) -> sample_metadata2
Correlate_with_phenotypes(dplyr::filter(ancom_species, TEST == "1 pocket-0 pocket"), sample_metadata2, mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="pocket") -> Inflammation_pocket
Correlate_with_phenotypes(dplyr::filter(ancom_species, TEST == "1 saliva-0 saliva"), sample_metadata2, mutate(compositional_transformation(select(species_matrix2, -ID)), ID = species_matrix2$ID), tissue ="saliva") -> Inflammation_saliva
Correlate_with_phenotypes(dplyr::filter(Results_ancom_pathway, TEST == "1 pocket-0 pocket"), sample_metadata2,pathway_abundance3 , tissue ="pocket") -> Inflammation_pocket_path

make_pheatmap(Inflammation_pocket) -> heatmap_sp_pocket
make_pheatmap(Inflammation_saliva) -> heatmap_sp_saliva
make_pheatmap(Inflammation_pocket_path) -> heatmap_pathway_pocket

ggsave(heatmap_sp_pocket, filename = "manuscript/Plots/Heatmap_Inflammation_sp_pocket.png",scale = 4)
ggsave(heatmap_sp_saliva, filename = "manuscript/Plots/Heatmap_Inflammation_sp_saliva.png", scale= 4)
ggsave(heatmap_pathway_pocket, filename = "manuscript/Plots/Heatmap_Inflammation_pathway_pocket.png", scale=4)

SAVE(PLOT = heatmap_sp_pocket, NAME= "manuscript/Plots/Heatmap_Inflammation_sp_pocket.png" )
SAVE(PLOT = heatmap_sp_saliva, NAME= "manuscript/Plots/Heatmap_Inflammation_sp_saliva.png" )
SAVE(PLOT = heatmap_pathway_pocket, NAME= "manuscript/Plots/Heatmap_Inflammation_pathway_pocket.png" )


#Pathway plots
Results_ancom_pathway
compositional_transformation(select(pathway_abundance2, -ID)) %>% mutate(ID = pathway_abundance2$ID) -> transformed_ptw_matrix 

select(transformed_ptw_matrix, c(ID, colnames(pathway_abundance2)[grepl("PWY-7111", colnames(pathway_abundance2))])) -> To_plot
left_join(To_plot, select(sample_metadata, c(ID,Source,Group ))) -> To_plot

To_plot %>% ggplot(aes(y=`PWY-7111: pyruvate fermentation to isobutanol (engineered)`, x=factor(Group))) + geom_boxplot() + facet_wrap(~Source) + theme_bw()

#Indv taxa plot
transformed_species_matrix %>% select(ID,s__Anaeroglobus_geminatus) -> To_plot
species_matrix2 %>% select(ID,s__Anaeroglobus_geminatus) -> To_plot

left_join(To_plot, select(sample_metadata, c(ID,Source,Group ))) -> To_plot
To_plot %>% ggplot(aes(y=s__Anaeroglobus_geminatus, x=factor(Group))) + geom_boxplot()+ geom_point() + facet_wrap(~Source) + theme_bw()


#Summary stats
Prev = function(x){
  x[x != 0] = 1
  sum(x)/length(x)
  
}
Summary_t = tibble()
for (Taxa in seq(1,nrow(matrix_taxonomy))){
  Taxa = matrix_taxonomy[Taxa,]
  ID = Taxa$ID
  Taxa %>% select(-ID) %>% t() %>% as.numeric() -> Abundance
  sapply(colnames(Taxa)[2:ncol(Taxa)], FUN = function(x){ str_split(x,"_")[[1]][1]}) %>% as.vector() -> IDs
  DF = tibble(Ab = Abundance, ID = IDs)
  Meta = sample_metadata %>% mutate(ID=StudieID) %>% select(ID, Group, Source) 
  left_join(DF, Meta, by ="ID") -> DF
  DF %>% group_by(Group, Source) %>% summarise(Mean=mean(Ab), SD=sd(Ab), Prevalence = Prev(Ab)) -> SU
  SU %>% mutate(Taxa= ID) -> SU
  Summary_t = bind_rows(Summary_t, SU)
}
write_tsv(x = Summary_t, path = "manuscript/Supplement/Summary_taxa.tsv")

Summary_p = tibble()
for (Taxa in seq(1,nrow(pathway_abundance))){
  Taxa = pathway_abundance[Taxa,]
  colnames(Taxa)[1] = "ID"
  ID = Taxa$ID
  Taxa %>% select(-ID) %>% t() %>% as.numeric() -> Abundance
  sapply(colnames(Taxa)[2:ncol(Taxa)], FUN = function(x){ str_split(x,"_")[[1]][1]}) %>% as.vector() -> IDs
  DF = tibble(Ab = Abundance, ID = IDs)
  Meta = sample_metadata %>% mutate(ID=StudieID) %>% select(ID, Group, Source) 
  left_join(DF, Meta, by ="ID") -> DF
  DF %>% group_by(Group, Source) %>% summarise(Mean=mean(Ab), SD=sd(Ab), Prevalence = Prev(Ab)) -> SU
  SU %>% mutate(Pathway= ID) -> SU
  Summary_p = bind_rows(Summary_p, SU)
}
write_tsv(x = Summary_p, path = "manuscript/Supplement/Summary_ptw.tsv")

################################################
#Strainphlan strain analysis####################
################################################
filename <- system("dir Strainphlan_trees/",intern=TRUE)
All_result = tibble()
tip_name_new = function(x){
  str_split(x, "_")[[1]][1]
  
  
}


List_plots = list()
I = 0
Get_plots = function(Tree_file){
  Tree = read.tree(Tree_file)
  Tree$tip.label = as_vector(lapply(Tree$tip.label,  tip_name_new))
  sample_metadata %>% filter(StudieID %in% Tree$tip.label) -> meta_tree
  meta_tree[match(Tree$tip.label,meta_tree$StudieID),] %>% mutate(Group = as.factor(Group)) -> meta_tree
  cophenetic(Tree) -> Dist_matrix
  pcoa(Dist_matrix) -> PCOA_strains
  Variance = PCOA_strains$values$Relative_eig
  PCOA_strains$vectors %>% as.data.frame()  %>% rownames_to_column("StudieID") %>% as_tibble()  -> PCOA_data
  left_join(PCOA_data, meta_tree, by= "StudieID") -> PCOA_data
  ggplot(PCOA_data, aes(label=StudieID,x=Axis.1,y=Axis.2,color= Group)) + geom_point(aes(shape=Source)) + theme_bw() +
    scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) + ggrepel::geom_text_repel() + 
    xlab(paste(c("PCoA1 (",as.character(round(Variance[1],2)*100),"%)"), collapse="")) + ylab(paste(c("PCoA2 (",as.character(round(Variance[1],2)*100),"%)"), collapse="")) -> PCoA_tree
  meta_tree %>% mutate(Col = ifelse(Group == 0, wes_palette("Royal1")[1], wes_palette("Royal1")[2] ) ) -> meta_tree
  plot(Tree,tip.color = meta_tree$Col ) -> Phylo_tree
  return(list(PCoA_tree, Phylo_tree))
}

for ( Tree_file in list.files(path="Strainphlan_trees/", pattern='*.tree', full.names = TRUE) ){

  Tree = read.tree(Tree_file)
  Tree$tip.label = as_vector(lapply(Tree$tip.label,  tip_name_new))
  sample_metadata %>% filter(StudieID %in% Tree$tip.label) -> meta_tree
  meta_tree[match(Tree$tip.label,meta_tree$StudieID),] %>% mutate(Group = as.factor(Group)) -> meta_tree
  cophenetic(Tree) -> Dist_matrix
  #pcoa(Dist_matrix) -> PCOA_strains
  #PCOA_strains$vectors %>% as.data.frame()  %>% rownames_to_column("StudieID") %>% as_tibble()  -> PCOA_data
  #left_join(PCOA_data, meta_tree, by= "StudieID") -> PCOA_data
  #ggplot(PCOA_data) + geom_point(aes(x=Axis.1,y=Axis.2, shape=Source, color= as.factor(Group))) + theme_bw() +
  #  scale_color_manual(values = wes_palette("Royal1"),labels = c("Mild PD", "Severe PD")) -> PCoA_tree
  #meta_tree %>% mutate(Col = ifelse(Group == 0, wes_palette("Royal1")[1], wes_palette("Royal1")[2] ) ) -> meta_tree
  #plot(Tree,tip.color = meta_tree$Col ) -> Phylo_tree
  L1 = length(unique(meta_tree$Group))
  L2 = length(unique(meta_tree$Source))
  if (L1 > 1 && L2 > 1){
      adonis2(Dist_matrix ~ Source + Group ,meta_tree, strata=StudieID) -> Result
  } else if (L1 >1){
    adonis2(Dist_matrix ~ Group, meta_tree) -> Result
  } else {
    adonis2(Dist_matrix ~ Source, meta_tree) -> Result
  }
  Result %>% as.data.frame() %>% rownames_to_column("Feature") %>% mutate(Tree=Tree_file) %>%
    as_tibble() %>% mutate(Index= I) %>% mutate(N = dim(meta_tree)[1],
  N_pocket =dim(filter(meta_tree, Source=="pocket"))[1], N_saliva=N-N_pocket, 
                N_disease=dim(filter(meta_tree, Group==1))[1], N_healthy=N-N_disease ) -> Result
  All_result = rbind(All_result, Result)
  #Plots = list(PCoA_tree, Phylo_tree)
  I = I+1
}

All_result %>% arrange(`Pr(>F)`) -> All_result
lapply(All_result$Tree,FUN=function(x){str_split(x,"\\.")[[1]][2]}) %>% as_vector() -> S
All_result %>% mutate(Taxon=S) %>%  select(Taxon, Feature, `Pr(>F)`, N, N_saliva, N_pocket, N_disease) %>% gt::gt()
All_result %>% mutate(Taxon=S) %>% select(-c(Tree, Index)) -> To_save
write_tsv(To_save,path = "manuscript/Supplement/Strain_permanova.tsv")

Get_plots(All_result$Tree[1]) -> Plots1
Plots1[[1]]
Get_plots(All_result$Tree[2])
Get_plots(All_result$Tree[3])
