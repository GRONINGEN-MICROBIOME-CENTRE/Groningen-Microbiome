
# GMHI function to run either "inference" or prediction
run_gmhi <- function(training_set, testing_set=NULL, prediction=FALSE) {
  
  training <- t(training_set) # taxa need to be on rows
  Healthy <- training[,colnames(training) %in% mdat[mdat$MED.DISEASES.None.No.Diseases %in% "Y",]$DAG3_sampleID]
  Nonhealthy <- training[,colnames(training) %in% mdat[mdat$MED.DISEASES.None.No.Diseases %in% "N",]$DAG3_sampleID]
  
  n_samples_H <- ncol(Healthy)
  n_samples_NH <- ncol(Nonhealthy)
  n_samples_tot <- ncol(training)
  n_samples_H+n_samples_NH==n_samples_tot
  
  # PH: species prevalence among healthy
  # PNH: species prevalence among non-healthy
  PH <- apply(Healthy, 1, function(i) (sum(i > 0))*100/n_samples_H) 
  PNH <- apply(Nonhealthy, 1, function(i) (sum(i > 0))*100/n_samples_NH) 
  
  PH_diff <- (PH-PNH)
  PH_fold <- (PH/PNH)
  PNH_fold <- (PNH/PH)
  all_matrix<-data.frame(cbind(training,PH_diff,PH_fold,PNH_fold))
  
  range(PH_diff)
  range(PH_fold)
  #range(PH_fold[!is.infinite(PH_fold)], na.rm=T)
  #range(PNH_fold[!is.infinite(PNH_fold)])

  theta_f <- seq(1.2, 2, by=0.1) # fold change 
  theta_d <- c(5,10,15,20) # difference
  
  ## From paper
  # As it is common in microbiome data to have discrepancies between species’ relative abundances to span several orders of magnitude, the geometric mean, rather than the arithmetic mean, is more appropriate to represent the mean relative abundance of MH species. 
  # More specifically, the Shannon’s diversity index, which is a weighted geometric mean (by definition) and commonly applied in ecological contexts, is used.

  alpha <- function(x){sum((log(x[x>0]))*(x[x>0]))*(-1)}
  
  report <- list()
  
  for(i in 1:length(theta_f)) {
    
    H_signature_sublist <- list()
    NH_signature_sublist <- list()
    report_sublist <- list()
    GMHI_list <- list()
    
    for(j in 1:length(theta_d)) {
      
      H_signature_sublist[[j]] <- data.frame(subset(as.matrix(all_matrix), all_matrix$PH_fold >= theta_f[i] & all_matrix$PH_diff >= theta_d[j]))
      NH_signature_sublist[[j]] <- data.frame(subset(as.matrix(all_matrix), all_matrix$PNH_fold >= theta_f[i] & all_matrix$PH_diff <= -theta_d[j]))
      
      # # remove Inf/-Inf which is created when the denominator/numerator is 0
      # H_signature_sublist[[j]] <- H_signature_sublist[[j]][!is.infinite(H_signature_sublist[[j]]$PH_diff) & !is.infinite(H_signature_sublist[[j]]$PH_fold),]
      # NH_signature_sublist[[j]] <- NH_signature_sublist[[j]][!is.infinite(NH_signature_sublist[[j]]$PH_diff) & !is.infinite(NH_signature_sublist[[j]]$PNH_fold),]
      
      H_shannon <- apply((H_signature_sublist[[j]][,-c(n_samples_tot:(n_samples_tot+3))]/100), 2, alpha)
      NH_shannon <- apply((NH_signature_sublist[[j]][,-c(n_samples_tot:(n_samples_tot+3))]/100), 2, alpha)
      
      H_sig_count <- apply(H_signature_sublist[[j]][,-c(n_samples_tot:(n_samples_tot+3))], 2, function(i) (sum(i > 0)))
      NH_sig_count <- apply(NH_signature_sublist[[j]][,-c(n_samples_tot:(n_samples_tot+3))], 2, function(i) (sum(i > 0)))
      
      constant <- data.frame(cbind(H_sig_count, NH_sig_count))
      
      HC1 <- constant[with(constant, order(-H_sig_count, NH_sig_count)), ]
      H_constant <- median(HC1$H_sig_count[1:floor(n_samples_H*0.01)])
      
      NHC1 <- constant[with(constant, order(H_sig_count, -NH_sig_count)), ]
      NH_constant <- median(NHC1$NH_sig_count[1:floor(n_samples_NH*0.01)])
      
      H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
      NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
      
      GMHI <- data.frame(log10((H_GMHI + 0.00001) / (NH_GMHI + 0.00001)))
      
      # Join our metadata
      GMHI <- GMHI %>% 
        rownames_to_column("DAG3_sampleID") %>%
        left_join(select(mdat, DAG3_sampleID, MED.DISEASES.None.No.Diseases), by="DAG3_sampleID")
      Healthy_GMHI <- GMHI %>% 
        filter(MED.DISEASES.None.No.Diseases=="Y") %>% 
        column_to_rownames("DAG3_sampleID")
      Nonhealthy_GMHI <- GMHI %>% 
        filter(MED.DISEASES.None.No.Diseases=="N") %>% 
        column_to_rownames("DAG3_sampleID")
      
      colnames(Healthy_GMHI)[1] <- "GMHI"
      colnames(Nonhealthy_GMHI)[1] <- "GMHI"
      
      Healthy_accuracy <- sum(Healthy_GMHI$GMHI > 0) * 100 / n_samples_H
      Nonhealthy_accuracy <- sum(Nonhealthy_GMHI$GMHI < 0) * 100 / n_samples_NH
      
      total_accuracy <- (Healthy_accuracy + Nonhealthy_accuracy)
      total_average_accuracy <- (Healthy_accuracy + Nonhealthy_accuracy) / 2
      
      report_sublist[[j]] <- cbind(theta_f[i], theta_d[j], 
                                   nrow(H_signature_sublist[[j]]), nrow(NH_signature_sublist[[j]]), 
                                   H_constant, NH_constant, 
                                   Healthy_accuracy, Nonhealthy_accuracy,
                                   total_accuracy, total_average_accuracy)
      report[[i]] <- report_sublist
    }
  }
  
  accuracy_table <- matrix(unlist(report), ncol=10, byrow=TRUE)
  colnames(accuracy_table) <- c("Fold change", "Difference",
                                "H_count","NH_count",
                                "H_constant", "NH_constant", 
                                "Healthy accuracy", "Non-healthy accuracy",
                                "Total accuracy","Balanced accuracy")
  accuracy_table <- data.frame(accuracy_table) %>% 
    arrange(desc(Balanced.accuracy))
  
  head(accuracy_table)
  
  write.csv(accuracy_table, paste0("Balanced_accuracy_", deparse(substitute(training_set)), ".csv"), row.names=F)
  
  # Generally the first row in the accuracy_table
  H_signature <- data.frame(subset(all_matrix, all_matrix$PH_fold >= accuracy_table[1,]$Fold.change & all_matrix$PH_diff >= accuracy_table[1,]$Difference))
  NH_signature <- data.frame(subset(all_matrix, all_matrix$PNH_fold >= accuracy_table[1,]$Fold.change & all_matrix$PH_diff <= -accuracy_table[1,]$Difference))
  
  H_species <- row.names(H_signature)
  NH_species <- row.names(NH_signature)
  
  H_constant <- accuracy_table[1,]$H_constant
  NH_constant <- accuracy_table[1,]$NH_constant
  
  if(prediction==TRUE) {
    
    testing <- t(testing_set) # taxa need to be on rows
    
    sp_H <- testing[row.names(testing) %in% H_species, ]
    sp_NH <- testing[row.names(testing) %in% NH_species, ]
    
    H_shannon <- apply((sp_H/100), 2, alpha)
    NH_shannon <- apply((sp_NH/100), 2, alpha)
    
    H_sig_count <- apply(sp_H, 2, function(i) (sum(i > 0)))
    NH_sig_count <- apply(sp_NH, 2, function(i) (sum(i > 0)))
    
    H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
    NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
    
    GMHI <- data.frame(log10((H_GMHI + 0.00001) / (NH_GMHI + 0.00001)))
    colnames(GMHI) <- c("GMHI")
    
    # GMHI classification
    result <- data.frame(cbind(GMHI, H_sig_count, NH_sig_count, H_shannon, NH_shannon, H_GMHI, NH_GMHI))
    
    result <- result %>%
      rownames_to_column("DAG3_sampleID") %>% 
      left_join(select(mdat, DAG3_sampleID, MED.DISEASES.None.No.Diseases), by="DAG3_sampleID")
    
    write.csv(result, paste0("GMHI_result_", "train_", deparse(substitute(training_set)), "_test_", deparse(substitute(testing_set)),".csv"), row.names=F)
    
    # Tabling outputting microbial species in the Response-prevalent and NR-prevalent groups.
    prevalence_input <- data.frame(cbind(PH, PNH, PH_diff, PH_fold, PNH_fold))
    
    H_plus_species <- prevalence_input[row.names(prevalence_input) %in% H_species,]
    H_plus_species$PNH_fold <- NULL
    colnames(H_plus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"
    H_minus_species <- prevalence_input[row.names(prevalence_input) %in% NH_species,]
    H_minus_species$PH_fold <- NULL
    colnames(H_minus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"
    
    prevalence_table <- data.frame(rbind(H_plus_species, H_minus_species)) %>% 
      rownames_to_column("Species") %>% 
      mutate(group=c(rep("H", nrow(H_plus_species)), rep("NH", nrow(H_minus_species)))) %>% 
      select(group, Species, PH, PNH, PH_diff, Fold.change..PH.PNH.or.PNH.PH.) %>% 
      mutate(group=relevel(factor(group, ordered=FALSE), ref="H")) %>% 
      group_by(group, PH) %>% 
      arrange(group, desc(PH), PNH)
    
    colnames(prevalence_table) <- c("Group","Species","Prevalence in Healthy samples (%)","Prevalence in Nohealthy samples (%)",
                                    "Difference (PH-PNH)","Fold-change (PH/PNH or PNH/PH)")
    
    write.csv(prevalence_table, paste0("prevalence_table_result_", "train_", deparse(substitute(training_set)), "_test_", deparse(substitute(testing_set)),".csv"), row.names=F)
    
    out <- list(H_species=H_species, NH_species=NH_species, result=result, accuracy_table=accuracy_table)
    return(out)
    
  }else{
    
    testing <- t(training_set) # taxa need to be on the rows
    
    sp_H <- testing[row.names(testing) %in% H_species, ]
    sp_NH <- testing[row.names(testing) %in% NH_species, ]
    
    H_shannon <- apply((sp_H/100), 2, alpha)
    NH_shannon <- apply((sp_NH/100), 2, alpha)
    
    H_sig_count <- apply(sp_H, 2, function(i) (sum(i > 0)))
    NH_sig_count <- apply(sp_NH, 2, function(i) (sum(i > 0)))
    
    H_GMHI <- ((H_sig_count/H_constant)*H_shannon)
    NH_GMHI <- ((NH_sig_count/NH_constant)*NH_shannon)
    
    GMHI <- data.frame(log10((H_GMHI + 0.00001) / (NH_GMHI + 0.00001)))
    colnames(GMHI) <- c("GMHI")
    
    # GMHI classification
    result <- data.frame(cbind(GMHI, H_sig_count, NH_sig_count, H_shannon, NH_shannon, H_GMHI, NH_GMHI))
    
    result <- result %>%
      rownames_to_column("DAG3_sampleID") %>% 
      left_join(select(mdat, DAG3_sampleID, MED.DISEASES.None.No.Diseases), by="DAG3_sampleID")
    
    write.csv(result, paste0("GMHI_result_", "train_", deparse(substitute(training_set)), "_test_", deparse(substitute(testing_set)),".csv"), row.names=F)
    
    # Tabling outputting microbial species in the Response-prevalent and NR-prevalent groups.
    prevalence_input <- data.frame(cbind(PH, PNH, PH_diff, PH_fold, PNH_fold))
    
    H_plus_species <- prevalence_input[row.names(prevalence_input) %in% H_species,]
    H_plus_species$PNH_fold <- NULL
    colnames(H_plus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"
    H_minus_species <- prevalence_input[row.names(prevalence_input) %in% NH_species,]
    H_minus_species$PH_fold <- NULL
    colnames(H_minus_species)[4] <- "Fold-change (PH/PNH or PNH/PH)"
    
    prevalence_table <- data.frame(rbind(H_plus_species, H_minus_species)) %>% 
      rownames_to_column("Species") %>%  
      mutate(group=c(rep("H", nrow(H_plus_species)),rep("NH", nrow(H_minus_species)))) %>% 
      select(group, Species, PH, PNH, PH_diff, Fold.change..PH.PNH.or.PNH.PH.) %>% 
      mutate(group=relevel(factor(group, ordered=FALSE),ref="H")) %>% 
      group_by(group, PH) %>% 
      arrange(group, desc(PH), PNH)
    
    colnames(prevalence_table) <- c("Group","Species","Prevalence in Healthy samples (%)","Prevalence in Nohealthy samples (%)",
                                    "Difference (PH-PNH)","Fold-change (PH/PNH or PNH/PH)")

    write.csv(prevalence_table, paste0("prevalence_table_result_", "train_", deparse(substitute(training_set)), "_test_", deparse(substitute(testing_set)),".csv"), row.names=F)

    
    out <- list(H_species=H_species, NH_species=NH_species, result=result, accuracy_table=accuracy_table)
    return(out)
  }
}