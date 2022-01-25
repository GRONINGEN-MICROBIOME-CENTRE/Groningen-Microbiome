library(tidyverse) #Standard packages for data load (readr), visualization (ggplot) and processing (dplyr)
library(patchwork)


############################################################################
#########################METAANALYSIS#######################################
############################################################################


get_level = function(x){
  x = str_split(x, "\\.")[[1]][1]
  return(x)
}


Meta_analyze = function(Gender, Results, BMI=F, Bonferroni= 134*6){
    Results %>% filter(! Metabolite == "Butyrobetain.Carnitine" ) -> Results
    Results %>% filter(! Metabolite == "Deoxycarnitine.Carnitine" ) -> Results
    if (Gender == "Male"){
      Results %>% filter(grepl("\\.0\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_Jul21_2021_MALES.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_BMI_Jul21_2021_MALES.csv") }
    
    } else if (Gender == "Female"){
      Results %>% filter(grepl("\\.1\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_Jul21_2021_FEMALES.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_BMI_Jul21_2021_FEMALES.csv") }
      
    } else{
      Results %>% filter(grepl("\\.3\\.",Source)) %>% select(-c(FDR, Source)) -> Results_LLD
      if (BMI == F){ Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_All_Jul21_2021.csv")
      } else {Results_Rotter = read_csv("Metanalysis/Rotterdam_2021_3/TMAO_RS16s_All_BMI_Jul21_2021.csv") }
    }
  
  #Results_Rotter %>% mutate(Metabolite = ifelse(Metabolite == "Deoxycarnitine",  "y-butyrobetaine", Metabolite))
  Results_LLD %>% mutate(Metabolite = ifelse(Metabolite == "y-butyrobetaine",  "Deoxycarnitine", Metabolite)) %>%
                  mutate(Metabolite = ifelse(Metabolite == "L-Carnitine",  "Carnitine", Metabolite) ) %>% 
    mutate(Metabolite = ifelse(Metabolite == "Butyrobetain.Carnitine", "Deoxycarnitine.Carnitine", Metabolite) )  %>%
    mutate(Metabolite = ifelse(Metabolite == "TMAO.Butyrobetaine", "TMAO.Deoxycarnitine", Metabolite) ) -> Results_LLD
  Results_Rotter %>% mutate(Metabolite = ifelse(Metabolite == "Butyrobetaine.Carnitine", "Deoxycarnitine.Carnitine", Metabolite) ) %>%
    mutate(Metabolite = ifelse(Metabolite == "TMAO.Butyrobetaine", "TMAO.Deoxycarnitine", Metabolite) ) -> Results_Rotter
  if ("N_samples" %in% colnames(Results_LLD)){ Results_LLD %>% select(-N_samples) ->  Results_LLD}
  
  Results_LLD$Level =  sapply(Results_LLD$Bug, get_level) ; Results_Rotter$Level =  sapply(Results_Rotter$Bug, get_level)
  Results_Rotter %>% select(-X1) %>% filter(!Level == "SampleID") -> Results_Rotter
  
  print("Number of species test")
  length(unique(Results_LLD$Bug)) ; length(unique(Results_Rotter$Bug))
  
  print("Plot summary stats")
  Results_LLD %>% ggplot() + geom_point(aes(x=Beta, y=-log10(Pvalue), col = Pvalue<0.05 )) +
    facet_grid(Metabolite~Level, scales="free") + theme_bw() + ggtitle('LLD') -> Summary_LLD
  Results_Rotter %>% ggplot() + geom_point(aes(x=Beta, y=-log10(Pvalue),  col = Pvalue<0.05)) +
    facet_grid(Metabolite~Level, scales="free") + theme_bw() + ggtitle('Rotterdam') -> Summary_rot
  Summary_LLD + Summary_rot + plot_layout(guides = "collect") -> Summary_cohorts
  
  
    ###########################################
    ##############Meta analysis################
    ###########################################
    print("Iteration through taxa")
    #size1 = 1032
    #size2 = 1427 
    
    Meta_stats = tibble()
    for (BUG in unique(Results_LLD$Bug)){
      if (! BUG %in% Results_Rotter$Bug){ next }
      Results_LLD %>% filter(Bug == BUG) -> Bug_LL
      Results_Rotter %>% filter(Bug == BUG)  -> Bug_R
      
      for (Metabolite_i in Bug_LL$Metabolite){
        Bug_LL %>% filter(Metabolite == Metabolite_i) -> Entry_L ; Bug_R %>% filter(Metabolite == Metabolite_i) -> Entry_R
        
        #p1.z	=	abs(qnorm(Entry_L$Pvalue/2)) ; p2.z	=	abs(qnorm(Entry_R$Pvalue/2))
        #if (Entry_L$Beta<0) p1.z = -1*p1.z
        #if (Entry_R$Beta<0) p2.z	=	-1*p2.z
        #size		= size1+size2
        #w.meta.z		=	(p1.z*size1+p2.z*size2)/sqrt(size1^2+size2^2)
        #w.meta.p		=	(1-pnorm(abs(w.meta.z)))*2
        
        ##Meta function
        #meta::metacont(n.e = size1, mean.e = Entry_L$Beta, sd.e = Entry_L$SE, n.c=size2, mean.c=Entry_R$Beta, sd.c=Entry_R$Pvalue) -> meta_value
        rbind(Entry_L, Entry_R) -> Input_data
        meta::metagen(TE=Beta, seTE=SE, data=Input_data, comb.fixed = T, comb.random= T ) -> meta_value
        
        Agreement = ifelse(sign(Entry_L$Beta) == sign(Entry_R$Beta), "Agree", "Desagree")
        Sub_result = tibble(Bug = BUG, Metabolite = Metabolite_i, MetaP=meta_value$pval.fixed, Treatment_effect=meta_value$TE.fixed,SE_meta=meta_value$seTE.fixed,Heterogeneity_Q=meta_value$Q,Heterogeneity_Pvalue=meta_value$pval.Q,  Concordance = Agreement , Meta_random_P=meta_value$pval.random, Treatment_effect_random=meta_value$TE.random, SE_meta_random=meta_value$seTE.random)
        
        left_join(x=Entry_L, y=Entry_R, by=c("Metabolite", "Bug"), suffix = c("_LLD", "_Rott")) -> Daa
        Sub_result =  left_join(x=Sub_result, y=Daa, by= c("Metabolite", "Bug") )
        
        Meta_stats = rbind(Meta_stats, Sub_result)
      }  
    }
    
    Meta_stats %>% mutate(meta_FDR = p.adjust(MetaP, "fdr")) %>% mutate(Bonferroni_indp = MetaP*Bonferroni) -> Meta_stats
    Meta_stats %>% mutate(Association_ID = paste(c(Metabolite, Bug),collapse= "-")) -> Meta_stats
    
    
    #####PLOTS
    print("Uncover taxa associated with more than one metabolite")
    
    Meta_stats %>% filter(meta_FDR < 0.07) %>% filter(Metabolite %in% c("Betaine","Carnitine","Choline", "Deoxycarnitine", "TMAO") ) -> Significant
    Final_figure = NULL
    for (B in unique(Significant$Bug)){
      Significant %>% filter(Bug == B) -> Taxa
      if (dim(Taxa)[1] > 1){
        print(B)
       # ggplot(Taxa) + geom_bar(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite), stat="identity") + theme_bw() +ggtitle(B) +theme(axis.title.x=element_blank(), axis.ticks.x=element_blank()) +ylab("Effect") + coord_flip() -> Fig
        ggplot(Taxa) + geom_point(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite), stat="identity") +  
          geom_pointrange(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite,ymin=Taxa$Treatment_effect-2*as.numeric(Taxa$SE_meta) , ymax=Taxa$Treatment_effect + 2*as.numeric(Taxa$SE_meta)),linetype = "dotted" ) +
          theme_bw() +ggtitle(B) +theme(axis.text.x = element_text(size = 5),axis.title.x=element_blank(),axis.title.y=element_blank() , plot.title = element_text(size = 5, face = "bold")) +ylab("Effect") + coord_flip() +  geom_hline(yintercept = 0) -> Fig
        if (is.null(Final_figure)){ Final_figure = Fig
        }else{ Final_figure = Final_figure  + Fig}
      }
    }
    ##Same as above but with Statistically significant bugs
    Meta_stats  %>% filter(Bonferroni_indp < 0.05) %>% filter(Metabolite %in% c("Betaine","Carnitine","Choline", "Deoxycarnitine", "TMAO") ) -> Significant
    Final_figure2 = NULL
    rng = range(0, -log10(min(Meta_stats$MetaP)))
    for (B in unique(Significant$Bug)){
      Meta_stats %>% filter(Bug == B) -> Taxa
      ggplot(Taxa) + geom_point(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite, col= -log10(Taxa$MetaP)), stat="identity") + 
        geom_pointrange(aes_(y=Taxa$Treatment_effect,x=Taxa$Metabolite,ymin=Taxa$Treatment_effect-2*as.numeric(Taxa$SE_meta) , ymax=Taxa$Treatment_effect + 2*as.numeric(Taxa$SE_meta), col= -log10(Taxa$MetaP)),linetype = "dotted" ) +
        theme_bw() + ggtitle(B) + theme(axis.text.y = element_text(size = 5),axis.text.x = element_text(size = 5),axis.title.x=element_blank(),axis.title.y=element_blank() , plot.title = element_text(size = 5, face = "bold"), legend.key.height= unit(0.5, 'cm'), legend.key.width= unit(0.25, 'cm'),legend.title =element_blank()) +ylab("Effect") + coord_flip() +  geom_hline(yintercept = 0) +
        scale_color_gradient2(low="blue", mid="cyan", high="purple", midpoint=mean(rng), breaks=seq(-100,100,2), limits=c(floor(rng[1]), ceiling(rng[2])))  -> Fig
      if (is.null(Final_figure2)){ Final_figure2 = Fig
      }else{ Final_figure2 = Final_figure2  + Fig}
    }
    if (! is.null(Final_figure2)){ Final_figure2 +  plot_layout(guides = 'collect') -> Final_figure2 }
    asso_id = function(x,y){
      paste(x,y, sep = "-")
    }
    Assoc_ID = asso_id(Meta_stats$Bug, Meta_stats$Metabolite)
    Meta_stats %>% mutate(Association_ID =  Assoc_ID) -> Meta_stats
    Fig_conc = Meta_stats %>% ggplot(aes(x=Concordance, y=-log10(MetaP) )) + geom_boxplot(outlier.shape = NA) + geom_jitter(aes(col = (Heterogeneity_Q<0.05), shape= (meta_FDR<0.05) )) + theme_bw() + geom_hline(yintercept = -log10(0.05/Bonferroni), col="grey", linetype="dashed",size=2)
    
    if (BMI == F){ 
        Output_name = paste(c("Manuscript/Additional_material/Summary_stats_", Gender,"_BMI-",as.character(F) ,".tsv"), collapse="")
        Output_LLD = paste(c("Manuscript/Additional_material/Summary_stats_LLD", Gender,"_BMI-",as.character(F) ,".tsv"), collapse="")
        Output_RS = paste(c("Manuscript/Additional_material/Summary_stats_RS", Gender,"_BMI-",as.character(F) ,".tsv"), collapse="")
    } else{
      Output_name = paste(c("Manuscript/Additional_material/Summary_stats_", Gender,".tsv"), collapse="")
      Output_LLD = paste(c("Manuscript/Additional_material/Summary_stats_LLD", Gender, ".tsv"), collapse="")
      Output_RS = paste(c("Manuscript/Additional_material/Summary_stats_RS", Gender, ".tsv"), collapse="")
      }
    print(paste(c("Saving results in ", Output_name), collapse="" ))
    write_tsv(x= Meta_stats, path = Output_name) ;  write_tsv(x= Results_LLD, path = Output_LLD) ;  write_tsv(x= Results_Rotter, path = Output_RS)
    
    return(list(arrange(Meta_stats, MetaP), Summary_cohorts, Final_figure, Fig_conc,Final_figure2) )
   
}


Gender_specific = function(All, Females, Males, Threshold=0.0012 ){
  Significant_all = filter(All, MetaP < Threshold)
  #Search for male associations not in All
  Males %>% filter(! Association_ID %in% Significant_all$Association_ID) %>% filter(MetaP < Threshold) -> Bugs_different_male
  #Search for Female associations not in All
  Females %>% filter(! Association_ID %in% Significant_all$Association_ID) %>% filter(MetaP < Threshold) -> Bugs_different_female
  
  #Put together ; left join with Males/Females to get them all. Left join with Bugs_different_male/Bugs_different_female to get the "significant"
  All %>% filter(Association_ID %in% c(Bugs_different_male$Association_ID, Bugs_different_female$Association_ID)) -> Stats_all_in_diff
  #left_join(x=Stats_all_in_diff, y=Bugs_different_male, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Males")) -> Stats_all_in_diff
  #left_join(x=Stats_all_in_diff, y=Bugs_different_female, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Females")) -> Stats_all_in_diff
  left_join(x=Stats_all_in_diff, y=Males, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Males")) -> Stats_all_in_diff
  left_join(x=Stats_all_in_diff, y=Females, by=c("Bug", "Metabolite", "Association_ID"), suffix=c("","_Females")) -> Stats_all_in_diff
  
  #Format to get effect sizes and SE
  #1. Mean
  Stats_all_in_diff %>% select(Bug, Metabolite, MetaP, Treatment_effect, MetaP_Males, Treatment_effect_Males, MetaP_Females, Treatment_effect_Females ) -> Stats_diff
  Stats_diff %>% select(-c(MetaP,MetaP_Females, MetaP_Males)) %>% gather(Group, Effect, c(Treatment_effect,Treatment_effect_Males, Treatment_effect_Females), factor_key=TRUE) -> Stats_diff_long
  #2.Standard Error
  Stats_all_in_diff %>% select(Bug, Metabolite, SE_meta, SE_meta_Females, SE_meta_Males ) %>% gather(Group, StdE, 3:5, factor_key=TRUE) -> Stats_se_long
  #3. Format "Group" uniformely
  Stats_diff_long %>% mutate(Group = ifelse(grepl("Males",Group), "Males", ifelse(grepl("Females", Group), "Females", "All"))) -> Stats_diff_long
  Stats_se_long  %>% mutate(Group = ifelse(grepl("Males",Group), "Males", ifelse(grepl("Females", Group), "Females", "All"))) -> Stats_se_long
  #4. Add together
  left_join(x=Stats_diff_long, y=Stats_se_long, by=c("Bug", "Metabolite", "Group")) -> Stats_diff_long
  #6. Plot
  Stats_diff_long %>% ggplot(aes(x=Bug, y = Effect, col=Group)) + geom_point(position=position_dodge(width=0.3)) + facet_wrap(~Metabolite, scales = "free") +
    geom_errorbar(position=position_dodge(width=0.3), aes(ymin=Effect-(2*StdE), ymax=Effect+(2*StdE)), width=.2, linetype = "dotted")  + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
    coord_flip() + theme_bw() + theme(axis.text.y = element_text(size = 5),axis.text.x = element_text(size = 5),
      axis.title.x=element_blank(),axis.title.y=element_blank() ) + geom_hline(yintercept = 0) + scale_color_brewer(palette = "Accent")   -> FPlot
  
  return(list(FPlot,Stats_diff_long))
  
}
draw_colnames_2 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
Plot_heatmap = function(Stats){
  
  Stats %>% mutate(Interpretation = ifelse(Bonferroni_indp<0.05, "  #  ", ifelse(meta_FDR<0.05, "·", " " )) ) -> Stats
  
  ##PLOT HEATMAP
  library(grid)
  assignInNamespace(x="draw_colnames", value="draw_colnames_2",
                    ns=asNamespace("pheatmap"))
  ####################3
  unique(Stats[Stats$MetaP<6.218905e-05*6,1:dim(Stats)[2]]$Bug) -> Bug_keep
  #unique(Stats[Stats$meta_FDR<0.05,1:dim(Stats)[2]]$Bug) -> Bug_keep
  filter(Stats, Bug %in% Bug_keep) -> Check
  
  Check %>% select(Bug, Metabolite, Treatment_effect) %>%  spread(Metabolite, Treatment_effect) %>% 
    as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_results
  Check %>% select(Bug, Metabolite, Interpretation) %>%  spread(Metabolite, Interpretation) %>% 
    as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_annotation
  Check %>% select(Bug, Metabolite, Interpretation, Level_LLD) %>%  spread(Metabolite, Interpretation) %>% 
    as.data.frame() %>%mutate(`Taxonomical level` = Level_LLD) %>% column_to_rownames(var = "Bug")  %>% select(`Taxonomical level`) -> Annotation
  Annotation_column = data.frame( `Feature class` = c("Precursor","Precursor","Precursor","Precursor","TMAO","Ratio","Ratio","Ratio","Ratio") )
  rownames(Annotation_column) = colnames(wide_results)
  
  New_names = c()
  for (i in rownames(wide_results)){
    str_replace(i, "genus.", "g.") -> i
    str_replace(i, "order.", "o.") -> i
    str_replace(i, "family.", "f.") -> i
    str_replace(i, "class.", "c.") -> i
    str_replace(i, "phylum.", "p.") -> i
    if ( grepl("unknown", i) ) { New_names = c(New_names, i) ; next }
    strsplit(i, "\\.id\\.")[[1]][1] -> i
    New_names = c(New_names, i)
  }
  rownames(wide_results)= New_names
  paletteLength <- 50
  #myColor <- colorRampPalette(c("Blue", "white", "Red"))(paletteLength)
  myColor = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(paletteLength)
  myBreaks <- c(seq(min(wide_results), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(wide_results)/paletteLength, max(wide_results), length.out=floor(paletteLength/2)))
  Plot = pheatmap::pheatmap(wide_results, display_numbers = wide_annotation,fontsize_number = 8,color=myColor, breaks=myBreaks, annotation_col = Annotation_column, border_color = NA, number_color = "black") #annotation_row = Annotation
  print(Plot)
  return(Plot)
  
  
}

Cohort_forest = function(Stats){
  
  Stats %>% filter( MetaP<6.218905e-05*6) -> Check
  Check %>% gather(Cohort, Beta ,c(Treatment_effect,Beta_LLD, Beta_Rott)) %>% select(Bug, Metabolite, Cohort, Beta,Heterogeneity_Pvalue) %>% mutate(Cohort = ifelse(Cohort == "Treatment_effect", "Meta-analysis", ifelse( Cohort == "Beta_LLD", "LLD", "RS"))) -> Beta_check
  Check %>% gather(Cohort, SE ,c(SE_meta,SE_LLD, SE_Rott )) %>% select(Bug, Metabolite, Cohort, SE ,Heterogeneity_Pvalue) %>% mutate(Cohort = ifelse(Cohort == "SE_meta", "Meta-analysis", ifelse( Cohort == "SE_LLD", "LLD", "RS"))) -> se_check
  left_join(Beta_check, se_check) -> Check2
  Check2 %>% mutate(Cohort = factor(Cohort, levels = c("Meta-analysis", "RS", "LLD") )) -> Check2
  All_plots = list() ; n = 1
  for (Met in unique(Check2$Metabolite)){
  Check2 %>% mutate(`Significant Heterogeneity` =  as.factor(Heterogeneity_Pvalue<0.05) ) %>% filter(Metabolite == Met) %>% mutate(`Significant Heterogeneity` = factor(`Significant Heterogeneity`, levels= c(F,T) ))  %>% ggplot(aes(y=Beta, x= Bug, col=Cohort )) + facet_grid(rows = vars(Metabolite) ,scales = "free") + theme_bw() + coord_flip() +
    geom_point(aes(shape=`Significant Heterogeneity`),position=position_dodge(width=0.5)) + geom_hline(yintercept = 0) + 
    geom_errorbar(aes(ymin=Beta-2*as.numeric(SE) , ymax=Beta + 2*as.numeric(SE)),width=.2,position=position_dodge(width=0.5),linetype = "dotted" ) + scale_shape(drop = FALSE) +
    scale_color_brewer(palette = "Dark2")  -> Plot
    if(! n == 1) { Plot + theme(legend.position = "none") -> Plot }
    All_plots[[n]]  = Plot
    n = n + 1
  }
  return(All_plots)
  
}

cm_to_inch = function(x){
  return(x/2.54)
}
save_pheatmap_pdf <- function(x, filename, width=10.25, height=8.24 ) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

Save_in_one_pdf = function(List_plots, Name, width=10.25, height=8.24){
  pdf(file=Name,width=10.25, height=8.24)
  print(List_plots)
  dev.off()
  
}


####################################
###############NO BMI ##############
####################################

Data_LLD = read_tsv("Results_16S_noBMI.tsv")
Outcome = Meta_analyze(Gender="All", Results = Data_LLD, BMI=F)
Stats =Outcome[[1]] ;  Fig_cohort = Outcome[[2]] ;  Fig_direct = Outcome[[3]] ;Fig_conc = Outcome[[4]] ; Fig_direct2 = Outcome[[5]]
Stats %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(6.2e-5))
Stats %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Outcome_m =  Meta_analyze("Male", Data_LLD, BMI=F)
Stats_Male =Outcome_m[[1]] ;  Fig_cohort_Male = Outcome_m[[2]] ;  Fig_direct_Male = Outcome_m[[3]] ; Fig_conc_Male = Outcome_m[[4]]
Stats_Male %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Male %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Male %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Outcome_f = Meta_analyze("Female", Data_LLD, BMI=F)
Stats_Female =Outcome_f[[1]] ;  Fig_cohort_Female = Outcome_f[[2]] ;  Fig_direct_Female = Outcome_f[[3]] ; Fig_conc_Female = Outcome_f[[4]]
Stats_Female %>% filter(meta_FDR < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Female %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(7.462687e-05))
Stats_Female %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()


Gender_specific(Stats, Stats_Male, Stats_Female, Threshold = 0.001)

####################################
############### BMI ############## THRESHOLD IS 134 (indp bug) * 6 (indp metabolites, including ratios) --> 0.05/n_test = 6.218905e-05
####################################

Data_LLD = read_tsv("Results_16S.tsv")
Outcome = Meta_analyze("All", Data_LLD, BMI=T)
Stats =Outcome[[1]] ;  Fig_cohort = Outcome[[2]] ;  Fig_direct = Outcome[[3]] ;Fig_conc = Outcome[[4]] ; Fig_direct_sig = Outcome[[5]]
ggsave(filename = "Manuscript/Additional_material/DirectionPlot_TMAOvsPrecursors_FDR07.pdf" , plot = Fig_direct, units = "mm",width = 290, height = 120)
ggsave(filename = "Manuscript/Additional_material/DirectionPlot_Significant.pdf" , plot = Fig_direct_sig, units = "mm",width = 290, height = 120)

Stats %>% filter(Bonferroni_indp < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(6.218905e-05))
Stats %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()

Plot_heatmap(Stats) -> Heatmap1
dev.off()
save_pheatmap_pdf(Heatmap1, "Manuscript/Additional_material/Association_Heatmap.pdf",width=12, height=10) #DONT RUN. For saving with the proper ratio, save in Rstudio window (Export PDF,  width=10.25, height=8.24  inches)
Cohort_forest(Stats) -> Forest1 
Save_in_one_pdf(Forest1,  "Manuscript/Additional_material/Association_bycohort_sign.pdf")

Sign_all = filter(Stats, Bonferroni_indp<0.05)



Outcome_m =  Meta_analyze("Male", Data_LLD, BMI=T)
Stats_Male =Outcome_m[[1]] ;  Fig_cohort_Male = Outcome_m[[2]] ;  Fig_direct_Male = Outcome_m[[3]] ; Fig_conc_Male = Outcome_m[[4]] ; Fig_direct_sig_Male = Outcome_m[[5]]
Stats_Male %>% filter(Bonferroni_indp < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Male %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(6.218905e-05))
Stats_Male %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()
#Are Male associations not seen in the overall? No.
Stats_Male %>% filter(Bonferroni_indp< 0.05) %>% filter(! Association_ID %in% Sign_all$Association_ID)


Cohort_forest(Stats_Male) -> Forest2
Save_in_one_pdf(Forest2,  "Manuscript/Additional_material/Association_bycohort_sign_Male.pdf")



Outcome_f = Meta_analyze("Female", Data_LLD, BMI=T)
Stats_Female =Outcome_f[[1]] ;  Fig_cohort_Female = Outcome_f[[2]] ;  Fig_direct_Female = Outcome_f[[3]]; Fig_conc_Female = Outcome_m[[4]]
Stats_Female %>% filter(Bonferroni_indp < 0.05) %>% group_by(Metabolite) %>% summarise(n())
Stats_Female %>% ggplot(aes(x=Treatment_effect, y=-log10(MetaP), col= meta_FDR < 0.05)) + geom_point(aes(shape=Heterogeneity_Pvalue<0.05),size=2.5, alpha=0.5) + facet_wrap(~Metabolite)  + theme_bw() +
  geom_hline(yintercept = -log10(6.218905e-05))
Stats_Female %>% select(Bug, Metabolite, MetaP, Treatment_effect) %>% gt::gt()
#Are Female associations not seen in the overall? No.
Stats_Female %>% filter(Bonferroni_indp< 0.05) %>% filter(! Association_ID %in% Sign_all$Association_ID)

Cohort_forest(Stats_Female) -> Forest3
Save_in_one_pdf(Forest3,  "Manuscript/Additional_material/Association_bycohort_sign_Female.pdf")

##Make Heatmaps for sex-stratified
Plot_heatmap_sex = function(Stats, Stats2){
  Stats %>% arrange(Bug, Metabolite) -> Stats ; Stats2 %>% arrange(Bug, Metabolite) -> Stats2
  Stats %>% mutate(Interpretation = ifelse(Bonferroni_indp<0.05, "  #  ", ifelse(meta_FDR<0.05, "·", " " )) ) -> Stats
  Stats2 %>% mutate(Interpretation = ifelse(Bonferroni_indp<0.05, "  #  ", ifelse(meta_FDR<0.05, "·", " " )) ) -> Stats2
  ##PLOT HEATMAP
  library(grid) #Rotation of x-axis labels
  assignInNamespace(x="draw_colnames", value="draw_colnames_2", ns=asNamespace("pheatmap"))
  ####################
  unique(Stats[Stats$MetaP<6.218905e-05*6,1:dim(Stats)[2]]$Bug) -> Bug_keep
  unique(Stats2[Stats2$MetaP<6.218905e-05*6,1:dim(Stats2)[2]]$Bug)-> Bug_keep2
  Bug_keep = unique(c(Bug_keep, Bug_keep2))
  
  filter(Stats, Bug %in% Bug_keep) -> Check
  filter(Stats2, Bug %in% Bug_keep) -> Check2
    
  Check %>% select(Bug, Metabolite, Treatment_effect) %>%  spread(Metabolite, Treatment_effect) %>% as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_results
  Check %>% select(Bug, Metabolite, Interpretation) %>%  spread(Metabolite, Interpretation) %>% as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_annotation
  Annotation_column = data.frame( `Feature class` = c("Precursor","Precursor","Precursor","Precursor","TMAO","Ratio","Ratio","Ratio","Ratio") )
  rownames(Annotation_column) = colnames(wide_results)
  
  Check2 %>% select(Bug, Metabolite, Treatment_effect) %>%  spread(Metabolite, Treatment_effect) %>% as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_results2
  Check2 %>% select(Bug, Metabolite, Interpretation) %>%  spread(Metabolite, Interpretation) %>% as.data.frame() %>% column_to_rownames(var = "Bug")  -> wide_annotation2
  Annotation_column2 = data.frame( `Feature class` = c("Precursor","Precursor","Precursor","Precursor","TMAO","Ratio","Ratio","Ratio","Ratio") )
  rownames(Annotation_column2) = colnames(wide_results2)
  
  
  New_names = c()
  for (i in rownames(wide_results)){
    str_replace(i, "genus.", "g.") -> i
    str_replace(i, "order.", "o.") -> i
    str_replace(i, "family.", "f.") -> i
    str_replace(i, "class.", "c.") -> i
    str_replace(i, "phylum.", "p.") -> i
    if ( grepl(".unknown.id", i) ) { next }
    strsplit(i, "\\.id\\.")[[1]][1] -> i
    New_names = c(New_names, i)
  }
  
  
  rownames(wide_results)= New_names
  rownames(wide_results2)= New_names
  
  
  paletteLength <- 50
  myColor = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(paletteLength)
  myBreaks <- c(seq(min(wide_results), 0, length.out=ceiling(paletteLength/2) + 1), 
                seq(max(wide_results)/paletteLength, max(wide_results), length.out=floor(paletteLength/2)))
  Plot_male = pheatmap::pheatmap(wide_results, display_numbers = wide_annotation,fontsize_number = 8,color=myColor, breaks=myBreaks, annotation_col = Annotation_column, border_color = NA, number_color = "black",  treeheight_row = 0,  treeheight_col=0)
  
  row.order = Plot_male$tree_row$order
  col.order = Plot_male$tree_col$order
  
  Plot_female = pheatmap::pheatmap(wide_results2[row.order, col.order], 
                     display_numbers = wide_annotation2[row.order, col.order],
                     fontsize_number = 8,color=myColor, breaks=myBreaks, annotation_col = Annotation_column, border_color = NA, number_color = "black", cluster_rows = F, cluster_cols = F)
  
  dev.off()
  save_pheatmap_pdf(Plot_male, "Manuscript/Additional_material/Association_Heatmap_Male.pdf",width=12, height=10) #DONT RUN. For saving with the proper ratio, save in Rstudio window (Export PDF,  width=10.25, height=8.24  inches)
  save_pheatmap_pdf(Plot_female, "Manuscript/Additional_material/Association_Heatmap_Female.pdf",width=12, height=10) #DONT RUN. For saving with the proper ratio, save in Rstudio window (Export PDF,  width=10.25, height=8.24  inches)
  
  
  
}
Plot_heatmap_sex(Stats_Male, Stats_Female)


###Study gender specific associations
Gender_specific(All = Stats, Males= Stats_Male, Females= Stats_Female, Threshold = 0.0012) -> Gender_asso
######Gender hetergoneity test####
rbind(Stats_Male, Stats_Female) -> Input_data
paste(Input_data$Bug, Input_data$Metabolite, sep = "-") ->ID_a
Input_data %>% mutate(ID = ID_a) -> Input_data
Hetero_results = tibble()
for (id in unique(ID_a)){
  Input_data %>% filter(ID == id) -> Input
  meta::metagen(TE=Treatment_effect, seTE=SE_meta, data=Input, comb.fixed = T, comb.random= F ) -> meta_value
  N = tibble(ID=id,Q_sex=meta_value$Q, P_Qsex= meta_value$pval.Q)
  rbind(Hetero_results, N) -> Hetero_results
}
unique(Hetero_results) %>% arrange(P_Qsex) %>% mutate(FDR= p.adjust(P_Qsex, "fdr")) -> Hetero_results
write_tsv(Hetero_results,"Manuscript/Additional_material/Gender_heterogeneity_stats.tsv")

paste(Gender_asso[[2]]$Bug, Gender_asso[[2]]$Metabolite, sep = "-") ->ID_a
Gender_asso[[2]] %>% mutate(ID = ID_a) -> Gender_assoc
left_join(Gender_assoc, Hetero_results) -> Gender_assoc
arrange(Gender_assoc, P_Qsex) %>% 
filter(P_Qsex < 0.05) -> Sig_gender ; unique(Sig_gender$ID) 

Gender_assoc %>% ggplot(aes(x=Bug, y = Effect, col=Group, shape=P_Qsex<0.05)) + geom_point(position=position_dodge(width=0.3)) + facet_wrap(~Metabolite, scales = "free") +
  geom_errorbar(position=position_dodge(width=0.3), aes(ymin=Effect-(2*StdE), ymax=Effect+(2*StdE)), width=.2, linetype = "dotted")  + theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) +
  coord_flip() + theme_bw() + theme(axis.text.y = element_text(size = 5),axis.text.x = element_text(size = 5),
                                    axis.title.x=element_blank(),axis.title.y=element_blank() ) + geom_vline(xintercept = 0) + scale_shape(name = "Heterogeneity (P<0.05)") +
  geom_hline(yintercept = 0) + scale_color_brewer(palette = "Accent") -> FPlot
ggsave(filename = "Manuscript/Additional_material/Gender_heterogeneity.pdf" , plot = FPlot, units = "mm",width = 290, height = 120)



