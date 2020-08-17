#Functions script manuscript Oral microbiome in PD

library(tidyverse) #Used in all functions
library(vegan) # Used in Diversity_analysis
library(ggforce)
library(ape) #Used in Diversity_analysis
library(ggtree) #Used in Diversity_analysis
#library(cluster) #Used in Diversity_analysis
#library(gamlss) #Used in fit_bezi
library(dunn.test)
select = dplyr::select




#"do_Metrics" generates a summary on Path N of a given matrix. It counts for each feature (taxa) the number of 0 we observe and some statistics about them. 
do_Metrics = function(MData,N){
  #Get column where the ID is, the last column, and remove it
  ID = MData[dim(MData)[2]] 
  MData %>% dplyr::select(-one_of(c("ID"))) -> MData
  #Function taht counts the number of 0s in a specific Taxa (V)
  number_0 = function(V){
    N_0 = 0
    for (i in V){
      if(i == 0){ N_0 = N_0 + 1}
    }
    return(N_0)
  }
  
  Total_summary = tibble()
  #Iterate though taxa and get SUMMARY stats, minimum abundance, maximum, median, standard deviation and number of 0
  for (C in seq(1:dim(MData)[2])){
    MData[,C] %>% mutate_if(is.character, as.numeric) %>% summarise_all(list(min, max,median,mean,sd,number_0)) %>% mutate(Name = colnames(MData[,C])) -> SUMMARY
    colnames(SUMMARY) = c("Minimum_abundance", "Maximum_abundance", "Median_abundance", "Mean_abundance","Standard_deviation","Number_of_0_counts","Organism")
    rbind(Total_summary, SUMMARY) -> Total_summary
  }
  write_csv(x = Total_summary, path = N)
}
#"contamination_stats" creates some plots regarding human contamination. It takes the metadata information and the contamination statistics. It checks on: Contamination vs nReads, Tissue vs Contamination, and the distribution of contamination.
#Host contamination is a common problem on metagenomic experiments, since even though the sequencing input might be high, majority of reads might be comming from the host.
contamination_stats = function(File = "Reads_contaminated.txt.txt", metadata=sample_metadata){
  #Read contamination stats and process
  Contaminated_metrics = read_delim(File, delim=" ")
  Contaminated_metrics %>% mutate(Total= Human_reads+Paired_Reads, Cont_percentage=Human_reads/Total) -> Contaminated_metrics
  #Put contamination metrics with the rest of covaraites
  metadata %>% filter(StudieID %in% Contaminated_metrics$Sample) %>% arrange(StudieID) -> metadata
  Contaminated_metrics %>% filter(Sample %in% metadata$StudieID) %>% arrange(Sample) -> Contaminated_metrics
  Contaminated_metrics %>% mutate(Source = metadata$Source) -> Contaminated_metrics
  #Look at the distribution of contamination
  Contaminated_metrics %>% ggplot() + geom_density(aes(x=Cont_percentage, col=Source)) + theme_bw() -> contamination_distribution
  ggplot(Contaminated_metrics,aes(x=Total, y=Cont_percentage)) + geom_point(aes(col=Source)) + theme_bw() + geom_smooth(method='lm', formula= y~x) ->contamination_vs_reads
  
  #Check the effect of the source in the amount of human reads
  Contaminated_metrics %>%
    ggplot(aes(x=Source, y=Cont_percentage)) + geom_boxplot()+ geom_sina()  + theme_bw() -> plot_source_contamination
  
  print("Doing statistical test on the % of contamination in Saliva vs Pocket")
  print(plot_source_contamination)
  wilcox.test(filter(Contaminated_metrics, Source=="pocket")$Cont_percentage, filter(Contaminated_metrics, Source=="saliva")$Cont_percentage) -> Saliva_vs_pocket
  print(Saliva_vs_pocket)
  print("Distribution of proportion of Human reads")
  print(contamination_distribution)
  print("Effect of the number of reads in the abundance of Human reads")
  print(contamination_vs_reads)
  wilcox.test(Contaminated_metrics$Total, Contaminated_metrics$Cont_percentage) -> Stats_plot
  print(Stats_plot)
  
}
#"distribution_taxonomy", it takes the overall taxonomy matrix and counts the different taxonomic classifications identified, overall comparing the two sources (saliva and pocket), by individual, and whether the number of reads plays a role in the taxons identified.
distribution_taxonomy = function(taxa_matrix = matrix_taxonomy, Contamination_file = "../Reads_contaminated.txt.txt",metadata=sample_metadata){
  t_class = taxa_matrix %>% filter(grepl("t__",ID)) %>% mutate(Level="strain")
  s_class = taxa_matrix %>% filter(grepl("s__",ID)) %>% filter(!grepl("t__",ID)) %>% mutate(Level="species")
  g_class = taxa_matrix %>% filter(grepl("g__",ID) & !grepl("s__",ID)) %>% mutate(Level="genus")
  f_class = taxa_matrix %>% filter(grepl("f__",ID) & !grepl("g__",ID)) %>% mutate(Level="family")
  o_class = taxa_matrix %>% filter(grepl("o__",ID) & !grepl("f__",ID)) %>% mutate(Level="order")
  c_class = taxa_matrix %>% filter(grepl("c__",ID) & !grepl("o__",ID)) %>% mutate(Level="clade")
  p_class = taxa_matrix %>% filter(grepl("p__",ID) & !grepl("c__",ID)) %>% mutate(Level="phylum")
  k_class = taxa_matrix %>% filter(grepl("k__",ID) & !grepl("p__",ID)) %>% mutate(Level="domain")
  
  print("Counting number the species identified per taxonomic level")
  
  raw_numbers = tibble(strain = dim(t_class)[1], species = dim(s_class)[1], 
                       genus = dim(g_class)[1], family = dim(g_class)[1], order = dim(o_class)[1],
                       clade = dim(c_class)[1], phylum = dim(p_class)[1], domain = dim(k_class)[1])
  raw_numbers %>% gather(Taxonomic_level, Number, strain:domain, factor_key=TRUE) -> tax_numbers
  ggplot(tax_numbers) + geom_bar(aes(x=Taxonomic_level,y=Number), stat="identity") + theme_bw() -> taxonomy_numbers
  print(taxonomy_numbers)
  
  
  print("Counting how many items are there per taxonomic level in each individual")
  classified = tibble()
  for (t_s in list(t_class,s_class,g_class,f_class,o_class,c_class,p_class,k_class)){
    taxa = t_s$ID
    transposed = as_tibble(t(select(t_s,-c(ID,Level))))
    colnames(transposed) = taxa
    transposed  %>% mutate_if(is.character, as.numeric) -> transposed
    (! as.matrix(transposed) == 0)*1 -> transposed
    
    classified_i = tibble()
    for (Indv in seq(1:dim(transposed)[1])){
      presences = sum(as_vector(transposed[Indv,]))
      presences = tibble(N=presences, Level=t_s$Level[1], ID =colnames(t_s)[Indv+1])
      classified_i = rbind(classified_i,presences)
    }
    
    classified = rbind(classified,classified_i)
  }
  
  
  classified %>% mutate(ID= str_replace(ID,"_feature_counts","")) %>% filter(!ID=="Pipelines") ->classified
  
  classified %>% arrange(desc(N)) -> ordered
  classified %>% mutate(Level = factor(Level,levels= c("strain","species","genus","family","order","clade","phylum","domain"))) -> classified
  classified %>% mutate(ID = factor(ID,levels= unique(ordered$ID))) -> classified
  
  
  
  classified %>% ggplot(aes(x=ID, y=N)) + geom_bar(stat="identity") + facet_wrap(~Level) + theme_bw() + theme(axis.text.x = element_text(angle = 90)) -> number_classified_by_ID
  print(number_classified_by_ID)
  
  
  print("Loading reads per sample, assesing whether the number of reads affect the taxons identified in each sample")
  Reads = read_delim(Contamination_file," ")
  Reads %>% filter(!Sample== "Pipelines") -> Reads
  Read_number = vector()
  for (entry in arrange(Reads,Sample)$Paired_Reads){
    for (i in c("strain","species","genus","family","order","clade","phylum","domain")){ 
      Read_number = c(Read_number,entry)
    }
  }
  arrange(classified,ID) %>% mutate(Reads= Read_number) -> Classified_reads
  ggplot(Classified_reads, aes(x=N,y=Reads,col=ID)) + geom_point() + facet_wrap(~Level, scales="free") + theme_bw()  + theme(legend.position = "none") -> reads_vs_class
  print(reads_vs_class)
}
#
do_bars = function(taxa_matrix, metadata){
  #MAtch metadata and counts
  metadata %>% filter(StudieID %in% taxa_matrix$ID) %>% arrange(StudieID) -> metadata
  taxa_matrix %>% filter(ID %in% metadata$StudieID) %>% arrange(ID) -> taxa_matrix
  taxa_matrix = select(taxa_matrix, -ID)
  
  metadata %>% mutate(N_ID=as.numeric(gsub("[^0-9.-]", "", StudieID))) %>% arrange(N_ID) -> test
  metadata$StudieID = factor(metadata$StudieID,levels(factor(metadata$StudieID))[match(test$StudieID,metadata$StudieID)])
  
  cbind(metadata,taxa_matrix) -> combined_matrix
  
  #tissue
  print("Study taxa composition by tissue")
  combined_matrix %>% gather(Taxa, Abundance, (dim(metadata)[2]+1):dim(combined_matrix)[2], factor_key=TRUE) %>%
    ggplot() + geom_bar(aes(x=ID_patient,y=Abundance, fill=Taxa), stat = "identity",position="stack") + facet_wrap(~Source, dir="v")+
    theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",legend.title = element_text(size = 3), legend.text = element_text(size = 6))+
    theme(legend.key.size = unit(0.1, "cm")) -> BAR1
  print(BAR1)
  
  print("Study taxa composition in Patients vs Controls")
  combined_matrix %>% gather(Taxa, Abundance, (dim(metadata)[2]+1):dim(combined_matrix)[2], factor_key=TRUE) %>%
    ggplot() + geom_bar(aes(x=StudieID,y=Abundance, fill=Taxa), stat = "identity",position="stack") + facet_wrap(~factor(Group), dir="v",scales="free")+
    theme_bw() + theme(axis.text.x = element_text(angle = 90), legend.position = "bottom",legend.title = element_text(size = 3), legend.text = element_text(size = 6))+
    theme(legend.key.size = unit(0.1, "cm")) -> BAR2
  print(BAR2)
  
  
  
  combined_matrix %>% filter(Source=="saliva") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) -> top
  
  combined_matrix %>% filter(Source=="pocket") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) -> top2
  
  
  
  combined_matrix %>% filter(as.character(Group)=="0") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance)) -> top3
  
  combined_matrix %>% filter(as.character(Group)=="1") %>% select(colnames(taxa_matrix)) %>% summarise_all(median) %>%
    gather(Taxa, Median_relative_abundance, 1:length(colnames(taxa_matrix)), factor_key=TRUE) %>%
    arrange(desc(Median_relative_abundance))  -> top4
  
  top %>% mutate(Sample_group="Saliva") -> top
  top2 %>% mutate(Sample_group="Pocket") -> top2
  top3 %>% mutate(Sample_group="Healthy") -> top3
  top4 %>% mutate(Sample_group="Patient") -> top4
  TOP = rbind(top,top2,top3,top4)
  return(TOP)
  
}
#First, we need to change a bit the taxa matrix. It has each bacteria as a row and each sample as column. We want it the other way around. We also call do_metrics on each of the samples that are preprocessed
Preprocess_f = function(Mtx,Summary_name){
  
  #Remove "_feature_counts" from the Sample name and the column ID (first column)
  Sample = gsub("_feature_counts","",colnames(Mtx))
  Sample = Sample[2:length(Sample)]
  
  #Remove higher orders from the Order name (a Species has info about Genera, Family... Separeated by | and we only one the last one)
  #Taxon = vector()
  #list_split = strsplit(Mtx$ID, "\\|")
  #for (i in seq(1:length(list_split))){
  #  element =  list_split[[i]]
  #  Taxon = c(Taxon, element[length(element)])
  #}
  
  #Transposition. Remove First row, which are the Taxon names
  Mtx2 = as_tibble(t(Mtx)[-1,])
  #Make Abundacnes as numeric
  Mtx2%>% mutate_if(is.character, as.numeric) -> Mtx2
  #Now we can add the taxon name that we produced by taking the last field after |
  colnames(Mtx2) = Mtx$ID
  
  # We count 0s, and other abundance stats and save the file in Summary
  do_Metrics(mutate(Mtx2, ID = Sample), N= Summary_name)
  
  #taxa2 = Mtx2[2:length(Mtx2)]
  taxa_filt= Mtx2
  
  #Add ID again
  mutate(taxa_filt, ID = Sample) -> taxa_filt
  #Move ID to fist position
  taxa_filt[c(dim(taxa_filt)[2], seq(1:(dim(taxa_filt)[2]-1)))] -> taxa_filt
  
  return(taxa_filt)
}
#Next, we want to check Alpha-diversity (within samples microbial composition) differences and Beta-diversity (differences between samples) differences.
#For that, we need to match metadata and abundance data. Then, we explore the data similarity by means of a hierarchical clustering on the abundance matrix. Alpha diversity is calculated (as Shanon entrophy) and compared between Pocket samples (Cases and Controls) and Saliva samples (Cases and Controls) separately by means of a Wilcoxon rank test (non-parametric). Beta diversity (Brais-Curtis distance) is decomposed as a PCoA and then we use a MANOVA test to look for differences between Cases and Controls in three different models. One including both sample sources, and the other two in pocket and saliva respectively. We control for the influence of the number of reads in each sample.
Diversity_analysis = function(taxa_matrix,metadata){
  print("Matching data and metadata IDs")
  metadata %>% filter(StudieID %in% taxa_matrix$ID ) %>% arrange(StudieID) -> metadata  
  taxa_matrix %>% filter(ID %in% metadata$StudieID ) %>% arrange(ID) -> taxa_matrix 
  
  ###Clustering of the data
  ## Hierarchical clustering
  print("Clustering based on disimilarity matrix")
  dist <- vegdist(taxa_matrix[,2:dim(taxa_matrix)[2]],  method = "bray")
  fit = hclust(dist)
  tree = as.phylo(fit)
  tree$tip.label = metadata$StudieID
  ggtree(tree)  %<+%  metadata + geom_tiplab(aes(col=factor(Group)),hjust=-0.3) + geom_tippoint(aes(shape=Source),size=3) + #aes(col=factor(Group)) + 
    xlim(NA,0.7) +  theme(legend.position=c(.9, .6)) -> Cluster_plot
  print(Cluster_plot)
  
  ### Alpha diversity
  print("Computing in sample alpha diversity and plotting")
  shannon_diversity = vegan::diversity(taxa_matrix[,2:dim(taxa_matrix)[2]], index = "shannon")
  
  alpha_diversities =  mutate(metadata, Shannon=shannon_diversity)
  gather(alpha_diversities, Diversity_index, value, c(Shannon), factor_key=TRUE) -> alpha_diversities
  
  #alpha_diversities %>% ggplot(aes(x=factor(Group), y=value,col=Source)) + geom_boxplot() +geom_sina() +
  #  facet_wrap(~Diversity_index, scales= "free") + theme_bw() -> Fig_alpha
  alpha_diversities %>% ggplot(aes(x=factor(Source), y=value,col=factor(Group))) + geom_boxplot(outlier.shape = NA) +geom_sina() +
  facet_wrap(~Diversity_index, scales= "free") + theme_bw() -> Fig_alpha
  print(Fig_alpha)
  
  ### Beta Diversity
  print("Calculating Beta diversity and computing PCoA")
  dist <- vegdist(taxa_matrix[,2:dim(taxa_matrix)[2]],  method = "bray")
  PCoA = pcoa(dist)
  
  ggplot(data=as_tibble(PCoA$values)) + geom_bar(aes(x=seq(1:dim(PCoA$values)[1]), y=Relative_eig), col="black", stat = "identity") + 
    theme_bw() -> Fig_Scree_plot
  print(Fig_Scree_plot)
  Eig_vector = PCoA$vectors[,1:3]
  metadata %>% mutate(PCoA1 = Eig_vector[,1], PCoA2 = Eig_vector[,2], PCoA3=Eig_vector[,3]) -> ordination_data
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA1, y=PCoA2, col=factor(Group), shape=Source)) +
    theme_bw() + labs(x=paste(c("PCoA1(",as.character(round(PCoA$values[1,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse="")) -> Fig_PCoA1
  ggplot(data= ordination_data) + geom_point(aes(x = PCoA2, y=PCoA3, col=Source, shape=factor(Group))) +
    theme_bw() + labs(x=paste(c("PCoA2(",as.character(round(PCoA$values[2,2],3)*100),"%)"), collapse=""), y=paste(c("PCoA3(",as.character(round(PCoA$values[3,2],3)*100),"%)"), collapse="")) -> Fig_PCoA2
  print("First two Components of the PCoA")
  print(Fig_PCoA1)
  #print("Components 2 and 3 of the PCoA")
  #print(Fig_PCoA2)
  
  ##Getting contamination stats 
  Reads = read_delim("../Reads_contaminated.txt.txt",delim=" ")
  metadata %>% mutate(n_reads = arrange(Reads,Sample)$Paired_Reads) -> metadata
  
  ###Alpha diversity statistical testing
  print("Using Dunn test to compare groups")
  
  dunn.test(alpha_diversities$value, alpha_diversities$Tissue_n_group)
  
  
  
  #######
  
  ###Beta diversity Statistical testing
  print("Doing Permanova test on Beta diversity")
  permanova_1 = adonis2(dist~ n_reads+ Source*factor(Group),data=metadata, permutations = 2000,strata=Source, method="bray")
  
  #Divide by tissue
  metadata %>% filter(Source=="pocket") -> meta_pocket
  metadata %>% filter(Source=="saliva") -> meta_saliva
  taxa_matrix %>% filter(ID %in% meta_saliva$StudieID) -> saliva_matrix
  taxa_matrix %>% filter(ID %in% meta_pocket$StudieID) -> pocket_matrix

  Results_total = tibble()
  for (Explanatory in c("Group","Periodontium (SUVmax)", "Periodontium (SUVmean)","IL-1Ra (pg/mL)" ,"IL-6 (pg/mL)" ,"hsCRP (ug/mL)" ,"WBC (10^6/mL)")){
    meta_saliva %>% select(Explanatory) %>% as_vector() -> Expl
    meta_saliva %>% mutate(Regressor = Expl) -> meta_saliva
    meta_pocket %>% select(Explanatory) %>% as_vector() -> Expl
    meta_pocket %>% mutate(Regressor = Expl) -> meta_pocket
    
    meta_pocket %>% drop_na(Explanatory) -> exp_pocket
    meta_saliva %>% drop_na(Explanatory) -> exp_saliva
    saliva_matrix %>% filter(ID %in% exp_saliva$StudieID) %>% select(-ID) -> SM
    pocket_matrix %>% filter(ID %in% exp_pocket$StudieID) %>% select(-ID) -> SP
    saliva_dist <- vegdist(SM,  method = "bray")
    pocket_dist <- vegdist(SP,  method = "bray")
    
    permanova_2_saliva = adonis2(saliva_dist~ n_reads + Regressor ,data=exp_saliva, permutations = 2000, method="bray")
    permanova_2_pocket = adonis2(pocket_dist~ n_reads + Regressor ,data=exp_pocket, permutations = 2000, method="bray")

    Results_pheno = tibble(R2 = permanova_2_saliva$R2[2], Pvalue = permanova_2_saliva$`Pr(>F)`[2], Tissue = "saliva", Pheno=Explanatory)
    rbind(Results_pheno, tibble(R2 = permanova_2_pocket$R2[2], Pvalue = permanova_2_pocket$`Pr(>F)`[2], Tissue = "pocket", Pheno=Explanatory)) -> Results_pheno
    rbind(Results_total,Results_pheno) -> Results_total
  }
  print("Total permanova")
  print(permanova_1)
  print("In saliva permanova")
  print(permanova_2_saliva)
  print("In pocket permanova")
  print(permanova_2_pocket)
  ###########
  
  return(list(Cluster_plot, Fig_alpha, Fig_Scree_plot, Fig_PCoA1,Results_total))
  
}

#Differential abundant bacteria.

#ANCOM
Perform_ancom = function(M, Meta){
  
  #Transmute abundance matrix
  select(M, -ID) %>% t() %>% as_tibble() %>% mutate(Taxa = colnames(select(M, -ID))) %>%
    `colnames<-`(c(M$ID, "Taxa")) %>% as.data.frame() %>% column_to_rownames("Taxa") -> S_DF
  
  
  #The idea is that 0s can have three sources: 
  #They can be outlayers due to some technical factor, in which case they are removed (whole taxon?)
  #out_cut is a value to detect outlayers. It is assumed a mixture of different distributions. If one of the clusters, which has fewer samples than the others is detected really far than the others, is removed. out_Cut is the proportion used to asses that a cluster is smaller.
  #They can be structural. So that they are biologically meaningful. It has two criteria in their paper, one is that is never observed and the other uses a lower bound (?). neg_lb = F means that it does not use the lower bound criteria
  #The other 0s are believed to be not detected due to sampling randomness and are pseudocounted.
  #Zero cut just removes by prevalance of at least 10%
  List_processing <- feature_table_pre_process(S_DF, Meta, "StudieID", group_var = "Tissue_n_group" , out_cut = 0.05, zero_cut = 0.90, lib_cut=0, neg_lb=F)
  
  Structural_0_per_group = as.data.frame(List_processing[[3]]) %>% rownames_to_column("Taxa") %>% as_tibble()
  colnames(Structural_0_per_group) = c("Taxa", "Group_0", "Group_1")
  
  Meta_corrected = List_processing[[2]]
  Abundance_corrected = List_processing[[1]]
  #Difference are the removed by prevalence or outlayers.
  
  
  #This fits the actual ancom, the structural 0s are assumed to have a W=Inf.
  #W is a value that represents how many tests can reject null hypothesis (fdr < 0.05)
  #Each taxon performes n_Taxon - 1 tests. Each of them comparing a ratio of that taxon/another taxon. The idea of clr is similar, but in clr there is a "representative" taxon that is the geometric mean.
  #The test depends: If no covariates it fits a wilcox test, if covariates it fits anova, if mixed it fits a linear miexed model
  #To visualize the effect size, it also uses this "representative ratio" clr and makes some volcano plots.
  #The effect is calculated fitting a linear model with the CLR transformed data (abunce  ~ Group). The beta (difference between group 1 and 0) is reported.
  #WATCH OUT! CLR transformation may give incorrect effects for low abundant bacteria: I have added a pseudocount of half the minimum that seems to 
  ANCOM(Abundance_corrected, Meta_corrected, List_processing[[3]], "Tissue_n_group", "fdr", "0.1") -> ANCOM_results
  as_tibble(ANCOM_results[[1]]) %>% arrange(desc(W)) -> Results_1
  #Changed the script so instead of returning plot returns data
  as_tibble(ANCOM_results[[2]]) -> Results_2
  colnames(Results_2) = c("taxa_id", "CLR_mean_diff", "W", "Structural_0")
  
  Figures = names(ANCOM_results[[3]])
  #ANCOM_results[[3]][[which(Figures == "s__Prevotella_veroralis")]]
  
  arrange(Results_2, taxa_id) -> Results_2
  arrange(Results_1, taxa_id) %>% mutate(CLR_mean_diff = Results_2$CLR_mean_diff) %>% arrange(desc(W)) -> ANCOM_table
  return(ANCOM_table)
  
}
#Bezi
gm_mean = function(x, na.rm=TRUE){ exp(sum(log(x), na.rm=na.rm) / length(x)) }
add_pseudocount = function(x, na.rm=TRUE){
  as_vector(x) -> y
  min(y[!y==0]) -> y
  Res = x+y
  return(Res)
}
fit_bezi = function(abundance_matrix,metadata, Report="Group", All = T){
  #Matrix of abundance --> Add pseudocount
  #Metadata 
  #Covariates
  #Calculate geometric mean per ID, include as covariate
  add_pseudocount(abundance_matrix) -> abundance_matrix
  Geometric_mean = apply(abundance_matrix, 1,FUN=gm_mean)
  metadata %>% mutate(Geom_mean = Geometric_mean) -> metadata
  result_bugg = tibble()
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    metadata %>% mutate(Dependent = dependent/100, Tissue_n_group = factor(Tissue_n_group)) -> regression_input
    distribution_bug = ggplot(regression_input) + geom_density(aes(Dependent, col=Tissue_n_group)) + theme_bw() + facet_wrap(~Source)
    if (length(unique(regression_input$Source)) > 1){
      Null = as.formula("Dependent ~ Source + log(Geom_mean)")
      H1 = as.formula("Dependent ~ Source + log(Geom_mean) + Tissue_n_group")
      return("NOT AVAILABLE")
    }else{
      Null = as.formula("Dependent ~  log(Geom_mean)")
      H1 = as.formula("Dependent ~  log(Geom_mean) + Tissue_n_group")
      #H1 = as.formula("Dependent ~  Tissue_n_group")
    }
    
    if (All == T){
      full.model = gamlss(H1, sigma.formula = H1, nu.formula = H1, data = regression_input,family = "BEZI", method=RS(1000))
    } else{
      full.model = gamlss(H1, data = regression_input,family = "BEZI", method=RS(1000))
    }
    Summary_model = summary(full.model)
    regression_input %>% group_by(Tissue_n_group) %>% summarise(N=mean(Dependent)) -> Mean_values
    paste(Mean_values$N, collapse=" ") -> Mean_0_vs_1
    if (All == T){
      tibble(Taxa=colnames(abundance_matrix)[N],Mu_Effet_Group_1 = Summary_model[3], Pvalue_mu=  Summary_model[30] ,Sigma_effect = Summary_model[5], Pvalue_Sigma=Summary_model[33], Nu_effect = Summary_model[9], Pvalue_Nu=Summary_model[36], Mean_0_vs_1 = Mean_0_vs_1, DIFF=diff(Mean_values$N)) -> output
    } else {
      tibble(Taxa=colnames(abundance_matrix)[N],Mu_Effet_Group_1 = Summary_model[3], Pvalue_mu=  Summary_model[18], Mean_0_vs_1 = Mean_0_vs_1, DIFF=diff(Mean_values$N)) -> output
      #tibble(Taxa=colnames(abundance_matrix)[N],Mu_Effet_Group_1 = Summary_model[2], Pvalue_mu= Summary_model[14], Mean_0_vs_1 = Mean_0_vs_1, DIFF=diff(Mean_values$N)) -> output
    }
    result_bugg = rbind(result_bugg, output)
  }
  
  return(result_bugg) 
}
#Logistic
fit_logistic = function(abundance_matrix, metadata){
  result_bugg = tibble()
  Predictors = compositional_transformation(abundance_matrix)
  for (N in seq(1:dim(Predictors)[2])){
    predictor = as_vector(Predictors[,N])
    metadata %>% mutate(predictor = as.vector(predictor)) -> Model_data
    H1 = as.formula("Group ~  predictor")
    Model = glm(H1, data=Model_data)
    summary(Model)$coefficients -> COEF
    
    Model_data %>% group_by(Tissue_n_group) %>% summarise(N=mean(predictor)) -> Mean_values
    paste(Mean_values$N, collapse=" ") -> Mean_0_vs_1
    
    tibble(Taxa=colnames(abundance_matrix)[N],Logs_odds = COEF[2], Pvalue=COEF[8] , Mean_0_vs_1 = Mean_0_vs_1, DIFF=diff(Mean_values$N)) -> output
    result_bugg = rbind(result_bugg, output)
  }
  return(result_bugg) 
  
}


#Parametric linear model
fit_linear_model = function(abundance_matrix,metadata, covariates,   Report="Group"){
  result_bugg = tibble()
  for (N in seq(1:dim(abundance_matrix)[2])){
    dependent = as_vector(abundance_matrix[,N])
    if (sum(dependent)==0){next
    }else if (sum((!dependent==0)*1)/sum(dependent) < 0.2){ next }
    
    metadata %>% mutate(Dependent = dependent) -> regression_input
    distribution_bug = ggplot(regression_input) + geom_density(aes(Dependent)) + theme_bw() + facet_wrap(~Source)
    
    if ("source" %in% covariates & length(covariates) == 1){ 
      lm(Dependent ~ Source + factor(Group) + n_reads, data= regression_input ) -> regression_output
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      if (Report == "Group"){ pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[3] }
      else if (Report == "Source"){pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[2] }
    }else if(length(covariates)==0){
      lm(Dependent ~ factor(Group) + n_reads, data= regression_input ) -> regression_output
      S = summary(regression_output)
      fold_change = mean(filter(regression_input, Group == 1)$Dependent) /  mean(filter(regression_input, Group == 0)$Dependent)
      pvalue = as_tibble(S$coefficients)$`Pr(>|t|)`[2]
    }
    
    shapiro.test(S$residuals) -> normality
    Normality_pval = normality$p.value
    
    f_t = as_tibble(t(c(colnames(abundance_matrix)[N], pvalue, fold_change, Normality_pval)))
    colnames(f_t) = c("Taxa","P-value","Fold_change(Case/Control)","Normality_pvalue")
    result_bugg = rbind(result_bugg, f_t)
  }
  result_bugg %>% mutate(FDR=p.adjust(`P-value`,"fdr"), FDR_normality=p.adjust(Normality_pvalue,"fdr")) -> result_bugg
  return(result_bugg)
}


#Others
make_relative = function(Column){ return(Column/sum(Column)*100) }
Filter_abundance = function(M, threshold=0.2, min_abundance = 0){
  #Filter species not seen in at least treshold*100% of samples
  apply(select(M, -ID), MARGIN=2, FUN = function(x){ (sum(x != 0)/length(x)) >= 0.2 }) -> Filter_vector
  F_names = colnames(select(M, -ID))[!Filter_vector]
  M %>% select(! F_names) -> M
  #Filter species with mean relative abundance lower than min_abundance
  apply(select(M, -ID), MARGIN=1, FUN = function(x){ sum(x) }) -> Total_reads_sample
  as_tibble(select(M, -ID)/Total_reads_sample) %>% summarise_all(mean) -> All_means
  colnames(select(M, -ID))[as_vector(All_means) > min_abundance] -> Taxa_not_to_filter
  M %>% select(c(ID, Taxa_not_to_filter)) -> M
  return(M)
}
compositional_transformation = function(Count_table){
  #Clr transformation
  library(microbiome)
  
  Count_table -> Counts
  #Transform counts
  colnames(Counts) -> N
  Counts %>% t() -> Check
  Counts_transformed = as_tibble(t(abundances(x=as.data.frame(Check), transform="clr")))
  ###
  Counts_transformed %>% select(colnames(Counts_transformed)[1:(length(colnames(Counts_transformed)))]) -> Counts_transformed
  return(Counts_transformed)
  
}

