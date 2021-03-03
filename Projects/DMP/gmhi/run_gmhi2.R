

# Parameters identified in original GMHI paper
MH_species <- c("s__Alistipes_senegalensis","s__Bacteroidales_bacterium_ph8","s__Bifidobacterium_adolescentis","s__Bifidobacterium_angulatum","s__Bifidobacterium_catenulatum","s__Lachnospiraceae_bacterium_8_1_57FAA","s__Sutterella_wadsworthensis")

MN_species <- c("s__Anaerotruncus_colihominis","s__Atopobium_parvulum","s__Bifidobacterium_dentium","s__Blautia_producta","s__candidate_division_TM7_single_cell_isolate_TM7c","s__Clostridiales_bacterium_1_7_47FAA","s__Clostridium_asparagiforme","s__Clostridium_bolteae","s__Clostridium_citroniae","s__Clostridium_clostridioforme","s__Clostridium_hathewayi","s__Clostridium_nexile","s__Clostridium_ramosum","s__Clostridium_symbiosum","s__Eggerthella_lenta","s__Erysipelotrichaceae_bacterium_2_2_44A","s__Flavonifractor_plautii","s__Fusobacterium_nucleatum","s__Gemella_morbillorum","s__Gemella_sanguinis","s__Granulicatella_adiacens","s__Holdemania_filiformis","s__Klebsiella_pneumoniae","s__Lachnospiraceae_bacterium_1_4_56FAA","s__Lachnospiraceae_bacterium_2_1_58FAA","s__Lachnospiraceae_bacterium_3_1_57FAA_CT1","s__Lachnospiraceae_bacterium_5_1_57FAA","s__Lachnospiraceae_bacterium_9_1_43BFAA","s__Lactobacillus_salivarius", "s__Peptostreptococcus_stomatis","s__Ruminococcaceae_bacterium_D16","s__Ruminococcus_gnavus","s__Solobacterium_moorei","s__Streptococcus_anginosus","s__Streptococcus_australis","s__Streptococcus_gordonii","s__Streptococcus_infantis","s__Streptococcus_mitis_oralis_pneumoniae","s__Streptococcus_sanguinis","s__Streptococcus_vestibularis","s__Subdoligranulum_sp_4_3_54A2FAA","s__Subdoligranulum_variabile","s__Veillonella_atypica")


# Pre-processing data matrix of species' relative abundances:
# Remove eukaryota and virus species
# Re-normalization of species' relative abundances after removing unclassified and virus species

abun <- read.csv("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/DAG3_data_ready/microbiome_v27/microbiome_8208samples_rdy.csv", row.names=1939, header=T)
abun <- abun[,!grepl(x=colnames(abun), pattern=c("k__Viruses|k__Eukaryota"))]
abun <- sweep(abun, 1, rowSums(abun),"/")
abun[abun < 0.00001] <- 0

abun <- t(abun[,grep(x=colnames(abun), pattern="s__")])
rownames(abun) <- stringr::str_split_fixed(rownames(abun), "\\.", 7)[,7] 

# Extracting Health-prevalent species present in metagenome
# Extracting Health-scarce species present in metagenome
MH_species_metagenome <- abun[row.names(abun) %in% MH_species, ]
MN_species_metagenome <- abun[row.names(abun) %in% MN_species, ]

# Diversity among Health-prevalent species
# Diversity among Health-scarce species
alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
MH_shannon <- apply((MH_species_metagenome), 2, alpha)
MN_shannon <- apply((MN_species_metagenome), 2, alpha) 

# Richness of Health-prevalent species
# Richness of Health-scarce species
R_MH <- apply(MH_species_metagenome, 2, function(i) (sum(i > 0)))
R_MN <- apply(MN_species_metagenome, 2, function(i) (sum(i > 0)))

# Median RMH from 1% of the top-ranked samples (see Methods)
# Median RMN from 1% of the bottom-ranked samples (see Methods)
MH_prime <- 7
MN_prime <- 31

# Collective abundance of Health-prevalent species
# Collective abundance of Health-scarce species
psi_MH <- ((R_MH/MH_prime)*MH_shannon)
psi_MN <- ((R_MN/MN_prime)*MN_shannon)

GMHI <- data.frame(log10((psi_MH+0.00001)/(psi_MN+0.00001))) # 0.00001 added to avoid having the denominator as 0
colnames(GMHI) <- c("GMHI")

# metadata
mdat <- read.csv("/groups/umcg-lifelines/tmp01/projects/dag3_fecal_mgs/DAG3_data_ready/phenotypes/DAG3_metadata_merged_ready_v27.csv", header=T)
mdat$DAG3_sampleID <- as.character(mdat$DAG3_sampleID)

GMHI <- GMHI %>%
     	rownames_to_column("DAG3_sampleID") %>%
        left_join(select(mdat, DAG3_sampleID, MED.DISEASES.None.No.Diseases), by="DAG3_sampleID")
	