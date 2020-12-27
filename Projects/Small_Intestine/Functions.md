Small Intestine Microbiota: Function Scripts
============================================


### 1.Summary Statistics 


###### To calculate summary statistics of participant metadata, link to original code: https://github.com/pausura/Project_Intestinal_Microbiome/blob/master/Function%20Scripts/summary_statistics_metadata_function.R

```{r}
summary_statistics_metadata <- function (metadata_input, category_table) {
  
  # Packages needed  
  install.packages("psych", repos = 'https://cran.r-project.org/web/packages/psych/index.html')
  library (psych)  #describe r function
  
  # Create other functions to calculate the different parameters
  
  ## Categorical values - create function to calculate the counts and the percentage for categorical variables
  tblFun <- function(x) {
    # Create a table
    tbl <- table(x)
    # Combine columnes/rows to get the counts and percentage (creates new table -> res)
    res <- cbind(tbl,round(prop.table(tbl)*100,2))
    # Give names to the columns
    colnames(res) <- c('Count','Percentage')
    res
  }
  
  ## NA sum function - counts the number of NA
  nzsum <- function(x) {
    sum (is.na(x))
  }
  
  if (missing(category_table)) {
    
    ## Calculate table1 with the whole data:
    
    my_results = matrix(ncol = 6, nrow = ncol(metadata_input))
    
    for (k in 1:ncol(metadata_input)){
      
      if (is.numeric(metadata_input[,k])) {
        # Keep in "x" the result from describe function (done in the columns) - for each factor  
        x = describe(metadata_input[,k])
        z = nzsum(metadata_input[,k])
        # In the new table ("x"): keep different values in the different columns
        my_results[k,1] = "numerical"
        my_results[k,2] = x$median
        my_results[k,3] = x$mean
        my_results[k,4] = x$sd
        my_results[k,5] = x$n
        my_results[k,6] = z
      }
      
      # Condition: if the column values are categorical  
      else {
        # Keep in "x" the result from tblFun function (done in the columns) - for each factor
        x = tblFun(metadata_input[,k])
        z = nzsum(metadata_input[,k])
        # In the new table ("x"): keep different values in the different columns 
        my_results[k,1]="categorical"
        # toString to keep the possible different values/categories in the same vector/column
        my_results[k,2]=toString(rownames(x))
        # First column table x = 'Count'
        my_results[k,3]=toString(x[,1]) 
        # Second column table x = 'Percentage'
        my_results[k,4]=toString(x[,2])
        # Sum of the values on column1 ("x")
        my_results[k,5]=sum(x[,1])
        my_results[k,6]= z
      }
    }
    
    
    # The column names from the original table = row names from the new table 
    rownames(my_results) = colnames(metadata_input)
    # Give names to the columns of the new table 
    colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "%/SD", "Number_non-zeros(n)", "Number_NA") 
    
    # Export the new table
    write.table (my_results, file = "./total_metadata_table.txt" , quote = F, sep = "\t")  
  }
  
  
  else {
    
    # Merge metadata_table with category_table
    metadata_input <- merge(category_table, metadata_input, by ="row.names")  
    # First column as rownames
    rownames(metadata_input) <- metadata_input[,1]
    metadata_input <- metadata_input[,-1]
    # Create a new column to assign a number to each category
    metadata_input$category <- as.integer(as.factor(metadata_input[,1]))      
    category_number <- nlevels(metadata_input[,1])
    
    # Save the different categories files in the global environment  
    matrix_list <- list()
    j = 1
    
    # Loop for to subset the different variables by category
    for (i in 1:category_number) {
      
      category_matrix <- subset(metadata_input, metadata_input$category== i )
      nam <- paste("Category", j , sep = "")
      matrix_list[[j]] <-  assign(nam, category_matrix)
      
      j <- j + 1  
    }
    
    # Save in different matrix the variables for each category 
    for (ii in 1:category_number) {
      
      new_matrix <- as.data.frame(matrix_list[[ii]])
      my_results = matrix(ncol = 6, nrow = ncol(new_matrix))
      
      # For each loop goes to the next column (numerical way)
      for (iii in 1:ncol(new_matrix)) { 
        
        # Condition: if the column values are numerical/continuous
        if (is.numeric(new_matrix[,iii])) {
          # Keep in "x" the result from describe function (done in the columns) - for each factor  
          x = describe(new_matrix[,iii])
          z = nzsum(new_matrix[,iii])
          # In the new table ("x"): keep different values in the different columns
          my_results[iii,1] = "numerical"
          my_results[iii,2] = x$median
          my_results[iii,3] = x$mean
          my_results[iii,4] = x$sd
          my_results[iii,5] = x$n
          my_results[iii,6] = z
        }
        
        # Condition: if the column values are categorical  
        else {
          # Keep in "x" the result from tblFun function (done in the columns) - for each factor
          x = tblFun(new_matrix[,iii])
          z = nzsum(new_matrix[,iii])
          # In the new table ("x"): keep different values in the different columns 
          my_results[iii,1]="categorical"
          # toString to keep the possible different values/categories in the same vector/column
          my_results[iii,2]=toString(rownames(x))
          # First column table x = 'Count'
          my_results[iii,3]=toString(x[,1]) 
          # Second column table x = 'Percentage'
          my_results[iii,4]=toString(x[,2])
          # Sum of the values on column1 ("x")
          my_results[iii,5]=sum(x[,1])
          my_results[iii,6]= z
        }
      }
      
      # The column names from the original table = row names from the new table 
      rownames(my_results) = colnames(new_matrix)
      # Give names to the columns of the new table 
      colnames(my_results) = c("Type", "Categories/Median", "Counts/Mean", "%/SD", "Number_non-zeros(n)", "Number_NA") 
      
      # Save the name of the variable to title the data.frame (table) 
      name_category <- new_matrix[1,1] 
      name_matrix <- paste(name_category, "_metadata_table1.txt", sep = "")
      final_name_matrix <- paste("./", name_matrix, sep = "")
      
      # Export the new table
      write.table (my_results, file = final_name_matrix , quote = F, sep = "\t") 
    }
  }   
}

```

### 2.Remove 'factor' variables 

###### To remove factor variables with 'n' number of levels (default: n = 1)
```{r}
# x = your dataframe
# nlevels = the number of levels the factor variable should have to be removed 

Remove_FactCols <- function(x, nlevels = 1){
  Remove_cols <- vector()
  for (i in 1:ncol(x)) {
    Name <- colnames(x)[i]
    if(nlevels(x[,i]) == 1){
      Remove_cols <- c(Remove_cols, Name)
    }
  }
  x_new <- x %>% select(-Remove_cols)
  print(paste(c("Function removed a total of ",length(Remove_cols), "variables"), collapse= ""))
  return(x_new)
}
```

### 3. Kruskal-Wallis Test 

###### To test a null of no difference in the distribution of multiple variables between two groups 

```{r}
# NB. Replace "CurrentStomaOrPouchType" with the variable name that contains the groups being compared (this variable should be in the first column)

kruskal.test_function <- function(x) {
#create matrix
  my_DF <- matrix(nrow = ncol(x)-1, ncol = 2)
  colnames(my_DF) <- c("Phenotype", "Pval")
  
#Kruskal-Wallis test, all variables in the dataframe tested
  a=1
  for (i in 2:ncol(x)) {
    my_f <- as.formula(paste(colnames(x)[i], "CurrentStomaOrPouchType", sep = " ~ ")) 
    my_test <- kruskal.test(my_f, x)
    for (k in 1:length(my_test$p.value)) {
      my_DF[a,1] <- colnames(x)[i]
      my_DF[a,2] <- my_test$p.value[k]
      a=a+1
    }
  }
  return(my_DF)
}
```

### 4. TopNTaxa 

###### To select the top 'n' (most abundant) taxa 
```{r}
TopNTaxa <- function(x,n) {
  names(x[1:n])
} 
```

### 5. PCoA 

###### (with plot)... link to original code: https://github.com/pausura/Project_Intestinal_Microbiome/blob/master/Function%20Scripts/pcoa_function.R

```{r}
plot_pcoa <- function(tax_level_table, pcoa_elements, variable_table, top_value) {
   ##Required packages
   library(vegan)
   library(ggplot2)
   library(psych)
   library(RColorBrewer)
   
   ## Generate a distance matrix - Bray method
   beta <- vegdist(tax_level_table, method="bray")
   ## cmdscale -> multidimensional scaling of a data matrix (distance matrix) 
   my_pcoa <- as.data.frame(cmdscale(beta, k = pcoa_elements))
   
   # Variable that contains the ".pdf" string to save plots in pdf files 
   a <- ".pdf"  
   
   ## If the variable_table is numeric:
   if (is.numeric(variable_table[,1])) {
     # Variable that contains colors codes to color the plot 
     my_col=c("#0000FF","#62A4D1","#5BE55B","#FFF000", "#FF0000")
     
     # Merge the variable_table with my_pcoa table --> the tax_level_table and the variable_table need to have the same number of rows/samples
     numeric_new_table <- merge(variable_table, my_pcoa, by="row.names") 
     # First column as Rownames
     rownames(numeric_new_table) <- numeric_new_table[,1]
     numeric_new_table <- numeric_new_table[,-1]
     
     # Create new tables with new row number
     new_variable_table <- subset(numeric_new_table, select = 1:ncol(variable_table))
     new_pcoa <- numeric_new_table
     new_pcoa[1:ncol(variable_table)] <- NULL
     
     # If the category table has more than one column:
     if (ncol(variable_table)>1) {
       # Calculate the mean for each column of the variable_table
       summary_table <- describe(new_variable_table)
       # Find the rows with the highest mean values (more abundant variables): top_value
       top_tax <- summary_table[order(summary_table$mean, decreasing = T)[1:top_value],]
       # Save the names of the highest mean values
       total_tax_names <- cbind(rownames(top_tax))
       # First column as Rownames
       rownames(total_tax_names) <- total_tax_names[,1]
       
       # Transpose the category table to merge  
       t_variable_table <- as.data.frame(t(new_variable_table))    
       # Merge the t_variable_table with the total_tax_names (table with the highest mean values) 
       numeric_plot_table <- merge(total_tax_names, t_variable_table, by = "row.names")    
       # First column as Rownames
       rownames(numeric_plot_table) <- numeric_plot_table[,1]
       numeric_plot_table <- numeric_plot_table[,-1]
       # Remove repeated column (rownames)
       numeric_plot_table <- numeric_plot_table[,-1]
       
       # Transpose the new table to have the variables as columns  
       t_plot_table <- as.data.frame(t(numeric_plot_table))
       
       write.table(t_plot_table, file = "./pcoa_table.txt", sep = "\t", quote = F)  
       
       # For loop to get plots as number of variables in t_plot_table 
       for (jj in 1:ncol(t_plot_table)) {
         # Save the name of the variable to title the plot
         name_category <- colnames(t_plot_table)[jj]
         # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file
         name_pdf <- paste(name_category,a, sep = "")
         # Create the plot
         bla = ggplot(new_pcoa, aes(x=V1, y=V2, geom = "blank", colour = t_plot_table[,jj])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(t_plot_table[,jj])) + theme_classic() + labs(x = "PCoA1", y = "PCoA2") + ggtitle(name_category)
         #x_axis: V1 from my_pcoa table (first PCoA element)
         #y_axis: V2 from my_pcoa table (second PCoA element)
         #colour: a different plot is generated for each column in t_plot_table 
         
         # Create the pdf file
         pdf(name_pdf)
         # Print the plot
         print(bla)
         # Empty the current device to create the next plot in the next loop
         dev.off()
       }
     }
     
     # The variable_table has only one column:
     else { 
       # Save the name of the variable to title the plot
       name_category <- colnames(new_variable_table)[1]
       # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file   
       name_pdf <- paste(name_category, a, sep = "")
       # Create the plot
       bla = ggplot(new_pcoa, aes(x=V1, y= V2, geom = "blank", colour = new_variable_table[,1])) + geom_point() + scale_color_gradientn(colours = my_col, colnames(variable_table[,1])) + theme_classic() + labs(x="PCoA1", y = "PCoA2") + ggtitle(name_category)
       #x_axis: V1 from my_pcoa table (first PCoA element)
       #y_axis: V2 from my_pcoa table (second PCoA element)
       #colour: colored by the column in the variable_table
       
       # Create the pdf file
       pdf(name_pdf)
       # Print the plot   
       print(bla)
       # Empty the current device  
       dev.off()
     }
   }
   
   ## The variable_table is categoric:          
   else { 
     # Merge the variable_table with my_pcoa table --> needed the same num of rows in both tables
     categoric_plot_table <- merge(variable_table, my_pcoa, by="row.names") 
     # First column as Rownames
     rownames(categoric_plot_table) <- categoric_plot_table[,1]
     categoric_plot_table <- categoric_plot_table[,-1]
     
     write.table(categoric_plot_table, file = "./pcoa_table.txt", sep = "\t", quote = F)  
     
     # Calculate the number of categories
     category_number <- nlevels(categoric_plot_table[,1])
     # Create a new column to assign a number to each category
     categoric_plot_table$category <- as.integer(categoric_plot_table[,1])
     # Create a new column to colour the plot depending on the category (level)        
     categoric_plot_table$color = "none" 
     # Create a palette of colors depending on the number of categories
     my_palette <- matrix(brewer.pal(category_number,"Set1"))
     
     # For loop to assign one different color to each different category in a new column
     for (i in 1:category_number){
       categoric_plot_table[categoric_plot_table$category == i,]$color = my_palette[i,1]
     }
     
     # Save the name of the variable to title the plot
     name_category <- colnames(variable_table)[1]
     # Add the "a" variable with ".pdf" string to the name_category to save the plot as a pdf file  
     name_pdf <- paste(name_category, a, sep = "")
     
     # Create the plot
     bla = ggplot (categoric_plot_table, aes(x=V1, y=V2, geom="blank", colour=color)) + geom_point () + scale_color_identity("All_categories", breaks=categoric_plot_table$color, labels=categoric_plot_table$category, guide = "legend") + theme_classic() + labs(x="PCoA1", y="PCoA2")
     #x_axis: V1 from my_pcoa table (first PCoA element)
     #y_axis: V2 from my_pcoa table (second PCoA element)
     #colour: colored by the new color column 
     
     # Create the pdf file
     pdf(name_pdf)
     # Print the plot
     print(bla)
     # Empty the current device
     dev.off()
   }
 }
 
```

### 6. Transform & Filter

###### link to original code: 

```{r}
#Asin Transformation function

transform_and_filter_taxa=function(x, samples_row=T, method="asin", missing_filter=0){
  x[x=="NA"]=0
  x[x=="NaN"]=0
  #if samples are in columns transpose
  if (!samples_row){
    x=as.data.frame(t(x))
  }
  #Exclude/keep columns that pass the missigness threshold
  if (missing_filter>100){
    stop("\n Hey! \n Values should be a proportion of missing values allowed per column: a value from 0 to 100")
  }
  x_filt=x[,((colSums(x !=0) / nrow(x)) *100 )>missing_filter]
  my_num_removed=ncol(x)-ncol(x_filt)
  print (paste(my_num_removed, "species removed due to many missing values"))
  if (method=="asin"){
    x_filt=x_filt/100
    x_filt=asin(sqrt(x_filt))
  } else if (method=="log"){
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0]/2)
    x_filt=x_filt+my_min
    x_filt=log10(x_filt)
  }else if (method=="clr"){
    x_filt=x_filt/100
    #replace 0 by the half of the smallest value observed
    my_min=min(x_filt[x_filt>0]/2)
    x_filt=x_filt+my_min
    #clr transformation (adapted from microbiome package function)
    d <- apply(x_filt, 2, function(b) {
      log(b) - mean(log(b))
    })
    x_filt=d
  }
  return(as.data.frame(x_filt))
}

```

### 7. Univariate Analysis

```{r}
#Function to carry out univariate correlation analysis between all categorical and numerical phenotypes and taxonomy in your dataframe

#__Arguments (inputs):
  # x = your dataframe containing phenotype variables followed by taxonomy (sample ID as rownames)
  # Num.Pheno = number of columns in your dataframe containing phenotypes 
  # Num.Taxa = column number of the first taxonomy (i.e. Num.pheno + 1)
  # TotTaxa = Total taxanomy number

#__Output:
  #A table containing the p-values, correlation coefficients and median or mean relative taxa       abundances. Also includes relevant descriptive statistics.

# VERSION 2. (with group means):
univariate.analysis2 <- function(x, Num.Pheno, Num.Taxa, TotTaxa) {
  
  #__Create a matrix where nrows = 'Ncombinations of total category variable groups * N taxa columns'
  #Calculate number of categorical-group combinations and number of non-catergorical (i.e numeric) variables
  Combinations <- vector()
  NonFactors <- vector()
  for (i in 1:Num.Pheno) {
    if(is.factor(x[,i])) {
      Combinations <- c(Combinations, ncol(combn(nlevels(x[,i]),2)))
    }
    else {
      NonFactors <- c(NonFactors, colnames(x)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Combinations)) * ncol(x[c(Num.Taxa:ncol(x))]), ncol = 18)
    colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Group1', 'Group2', 'CategoricalNumerical','NGroups', 'NSamplesG1/Nnonzero', 'NSampleG2/Nzero', 'NA', 'P-value', 'BonferroniAdjust', 'FDRadjust','BelowAdjustPval', 'CorrelationCoeff', 'Mean1', 'Mean2', 'PrevalenceGroup1', 'PrevalenceGroup2')
  }
  
  
  #__Add correlation analysis data to the matrix 
  a=1 
  #Iterate through taxonomy columns
  for (j in Num.Taxa:ncol(x)) {
    #iterate through phenotype columns
    for(k in 1:Num.Pheno){
      #following conditions applied to all factor variables
      if(is.factor(x[,k])){
        #generate matrix with all the group combinations for each factor variable
        combos <- combn(levels(x[,k]),2)
        #pairwise.wilcox.test 
        wilcox <- (pairwise.wilcox.test(x[,j], x[,k], p.adjust.method = 'none'))$p.value
        #condition applied to all factor variables with more than 2 groups
        if(ncol(combos) >2){
          #iterate through the columns of the group combination matrices 
          for (z in 1:ncol(combos)) {
            #for each column in the combiantion matrices fill a row in my_matrix with the taxa name 
            my_matrix[a,1] <- colnames(x)[j]
            #for each column in the combination matrices fill a row in my_matrix with the relevant phenotype  
            my_matrix[a,2] <- colnames(x)[k]
            #adding the group combinations to my_matrix
            my_matrix[a,3] <- combos[1,z]
            my_matrix[a,4] <- combos[2,z]
            my_matrix[a,5] <- "Category"
            my_matrix[a,6] <- nlevels(x[,k])
            my_matrix[a,7] <- sum(x[,k]==combos[1,z], na.rm = T) #Number of samples in group1
            n_group1=sum(x[,k]==combos[1,z], na.rm = T) 
            my_matrix[a,8] <- sum(x[,k]==combos[2,z], na.rm = T) #Number of samples in group2
            n_group2=sum(x[,k]==combos[2,z], na.rm = T)
            my_matrix[a,9] <- colSums(is.na(x[k]))
            my_matrix[a,10] <- wilcox[c(combos[2,z]), c(combos[1,z])] #pairwise.wilcox.test
            my_matrix[a,15] <- mean(subset(x, x[c(k)] == combos[1,z], select = j)[,1]) #Median taxonomy relative abundances group 3
            my_matrix[a,16] <- mean(subset(x, x[c(k)] == combos[2,z], select = j)[,1]) #Median taxonomy relative abundances group 2
            my_matrix[a,17] <- colSums(subset(x, x[c(k)] == combos[1,z], select = j) != 0, na.rm = T)/ n_group1 *100 #Taxonomy prevalences in group1
            my_matrix[a,18] <- colSums(subset(x, x[c(k)] == combos[2,z], select = j) != 0, na.rm = T)/n_group2 *100 #Taxonomy prevalences in group2
            a=a+1
          }
        } else { # i.e. following conditions applied to all factor variables with 2 groups 
          my_matrix[a,1] <- colnames(x)[j] #same as above
          my_matrix[a,2] <- colnames(x)[k]
          my_matrix[a,3] <- combos[1,1] # the combination matrices only have one column hence always insert column 1 
          my_matrix[a,4] <- combos[2,1]
          my_matrix[a,5] <- "Category"
          my_matrix[a,6] <- nlevels(x[,k])
          my_matrix[a,7] = sum(x[,k]==combos[1,1], na.rm = T)
          n_group1=sum(x[,k]==combos[1,1], na.rm = T)
          my_matrix[a,8] = sum(x[,k]==combos[2,1], na.rm = T)
          n_group2=sum(x[,k]==combos[2,1], na.rm = T)
          my_matrix[a,9] <- colSums(is.na(x[k]))
          my_matrix[a,10] <- wilcox[c(combos[2,1]), c(combos[1,1])]
          my_matrix[a,15] <- mean(subset(x, x[c(k)] == combos[1,1], select = j)[,1])
          my_matrix[a,16] <- mean(subset(x, x[c(k)] == combos[2,1], select = j)[,1])
          my_matrix[a,17] <- colSums(subset(x, x[c(k)] == combos[1,1], select = j) != 0, na.rm = T)/n_group1 *100 # Species prevalence per group
          my_matrix[a,18] <- colSums(subset(x, x[c(k)] == combos[2,1], select = j) != 0, na.rm = T)/n_group2 *100
          a=a+1
        }
      }
      else { # i.e. following conditions applied to all non-factor variables 
        my_matrix[a,1] <- colnames(x)[j]
        my_matrix[a,2] <- colnames(x)[k]
        my_matrix[a,3] <- colnames(x)[k] # these variables are numerical i.e. no groups :hence insert phenotype name again
        my_matrix[a,4] <- colnames(x)[k] # ^
        my_matrix[a,5] <- "Numeric"
        my_matrix[a,6] <- 0
        my_matrix[a,7] <- colSums(x[k] != 0, na.rm = T) #Number of nonZeros 
        my_matrix[a,8] <- colSums(x[k] == 0, na.rm = T) #Number of zeros 
        my_matrix[a,9] <- colSums(is.na(x[k])) #Number of NAs 
        my_matrix[a,10] <- (cor.test(x[,j], x[,k], method = 'spearman', exact = F))$p.value #Spearman correlation p-values
        my_matrix[a,14] <- (cor.test(x[,j], x[,k], method = 'spearman', exact = F))$estimate #Spearman correlation coefficient values
        my_matrix[a,15] <- mean(x[,j], na.rm = T) 
        my_matrix[a,16] <- mean(x[,j], na.rm = T)
        my_matrix[a,17] <- colSums(x[j] != 0, na.rm = T) / nrow(x) *100 #Taxonomy prevalences
        my_matrix[a,18] <- colSums(x[j] != 0, na.rm = T) / nrow(x) *100
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,10] <- as.numeric(as.character(my_matrix[,10]))
  my_matrix[,11] <- p.adjust(c(my_matrix[,10]), method = 'bonferroni')
  my_matrix[,12] <- p.adjust(c(my_matrix[,10]), method = 'fdr')
  my_matrix[,13] <- ifelse(my_matrix[,10] < (0.05/TotTaxa), "yes", "no") #0.05/number of bacteria - is the p-value below of higher than this number (yes/no) 
  return(my_matrix)
}

```

### 8. Remove rows with <n samples

```{r}
# This function removes rows (i.e. univariate comparisons) from the univariate.analysis output dataframe. The rows removed correspond to categorical pairwise comparisons with less than 'n' samples in each group.

#__Arguments (Inputs):
  # x = your output dataframe from univariate.analysis function
  # MinSamples = minumum number of samples you want the groups in each categorical pairwise         comparison to have.

remove.rows <- function(x, MinSamples = 20){  
  remove_rows <- vector()
  for (i in 1:nrow(x)) {
    names <- row.names(x)[i] 
    if(x[i,5] == 'Category' & x[i,7] < MinSamples | x[i,5] == 'Category' & x[i,8] < MinSamples){
      remove_rows <- c(remove_rows, names)
    }
  }
  x_new <- subset(x, !(rownames(x) %in% remove_rows))
  print(paste(c("Function removed a total of ",length(remove_rows), "variables"), collapse= ""))
  return(x_new)
}

```

### 9. Prevalence Filter 

```{r}
# x = dataframe
# RowSpecies = column number of the first species in the dataframe
# threshold = prevalence percentage below which taxa will be filtered out

Prevalence.filter <- function(x, RowSpecies, threshold = 15){
  remove_cols <- vector()
  for (i in RowSpecies:ncol(x)) {
  cname <- colnames(x)[i]
  if(colSums(x[i] != 0) / nrow(x) *100 < threshold){
    remove_cols <- c(remove_cols, cname) 
  }
}
x_new <- x %>% select(-remove_cols)
print(paste(c("Function removed a total of ",length(remove_cols), "variables"), collapse= ""))
return(x_new)
}

```

### 10. Remove Columns with >X% NAs

```{r}
#Function to remove N columns with a >= N% of NAs

RemoveColumnsNA = function(DF, Threshold=0.9){
  #packages
  library(dplyr)
  #Create vector where we will store the columns that have more NA than the threshold
  Remove_columns = vector()
  #Iterate through columns
  for (i in 1:ncol(DF)){
    #i is the column number, so we get the Name of the column and the values
    Name = colnames(DF)[i]
    Column = DF[,i]
    #Calculate the % of NA 
    Percentage_NAs = length(Column[is.na(Column)])/length(Column)
    if (Percentage_NAs >= Threshold){
      #Add to our vector the columns that cross the threshold 
      Remove_columns = c(Remove_columns, Name)
      }
  }
  #Finally, using the vector we remove all the columns which name are in the vector
  DF_new <- DF %>% select(-Remove_columns) 
  print(paste(c("Function removed a total of ",length(Remove_columns), " variables"), collapse= ""))
  return(list(DF_new, Remove_columns))
}
```

### 11. Heatmap conversion

```{r}
#Function to convert coeffiecient estimate and Q-value into a factor variable for a heatmap 

# -3 = 0.0001,   -2 <= 0.001,    -1 <= 0.05,    0 > 0.05,    1 <= 0.05,    2 <= 0.001


FactorVariableHeatmapMaaslin <- function(x){
  for (i in 1:nrow(x)) {
    if(x$Coefficient[i] < 0){
      if(x$Q.value[i] <= 0.0001){
        x$Maaslin_[i] <- -3
      }
      else{
        if(x$Q.value[i] <= 0.001){
          x$Maaslin_[i] <- -2
        }
        else{
          if(x$Q.value[i] <= 0.05){
            x$Maaslin_[i] <- -1
          }
          else{
            x$Maaslin_[i] <- 0 
          }
        }
      }
    }
    else{
      if(x$Q.value[i] <= 0.0001){
        x$Maaslin_[i] <- 3
      }
      else{
        if(x$Q.value[i] <= 0.001){
          x$Maaslin_[i] <- 2
        }
        else{
          if(x$Q.value[i] <= 0.05){
            x$Maaslin_[i] <- 1
          }
          else{
            x$Maaslin_[i] <- 0 
          }
        }
      }
    }
  }
  return(x)
}

```

### 12. Logistic Regression

```{r}

LogisticRegression.function <- function(DF, nPheno, colTaxa, my_preds){
  Nlevels <- vector()
  NonFactors <- vector()
  for (i in 1:nPheno) {
    if(is.factor(DF[,i])) {
      Nlevels <- c(Nlevels, (nlevels(DF[,i])-1))
    }
    else {
      NonFactors <- c(NonFactors, colnames(DF)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Nlevels)) * ncol(DF[c(colTaxa:ncol(DF))]), ncol = 9)
    colnames(my_matrix) <- c('Taxonomy', 'Phenotype', 'Estimate', 'Std.Error', 'z_value','Pr(>|z|)','BonferroniAdjust','FDRAdjust', 'Prevalence')
  }
  
  a=1
  for (j in colTaxa:ncol(DF)) {
    my_f <- as.formula(paste(colnames(DF)[j], paste(my_preds, collapse = " + "), sep = " ~ "))
    my_lm <- glm(my_f, family = binomial(link="logit"), data = DF)
    for(k in 2:nrow(summary(my_lm)$coefficients)){
      my_matrix[a,1] <- colnames(DF)[j]
      my_matrix[a,2] <- rownames(summary(my_lm)$coefficients)[k]
      my_matrix[a,3] <- summary(my_lm)$coefficients[k,1]
      my_matrix[a,4] <- summary(my_lm)$coefficients[k,2]
      my_matrix[a,5] <- summary(my_lm)$coefficients[k,3]
      my_matrix[a,6] <- summary(my_lm)$coefficients[k,4]
      my_matrix[a,9] <- sum(DF[,j])/nrow(DF) *100
      a=a+1
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,6] <- as.numeric(as.character(my_matrix[,6]))
  my_matrix[,7] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
  my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'fdr')
  return(my_matrix)
}

```

### 13. Logistic regression descriptive statistics

```{r}

LgRgr.DescriptiveStats <- function(DF, nPheno, colTaxa){
  Nlevels <- vector()
  NonFactors <- vector()
  for (i in 1:nPheno) {
    if(is.factor(DF[,i])) {
      Nlevels <- c(Nlevels, (nlevels(DF[,i])-1))
    }
    else {
      NonFactors <- c(NonFactors, colnames(DF)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Nlevels)) * ncol(DF[c(colTaxa:ncol(DF))]), ncol = 8)
    colnames(my_matrix) <- c('Taxonomy', 'Phenotype','Group', 'nGroup', 'nGroupOther', 'PrevalenceWithinGroup', 'PrevalenceWithinGroupOther', 'TotalPrevalence')
  }
  
  a=1
  for (j in colTaxa:ncol(DF)) {
    for (i in 1:nPheno) {
      if(is.factor(DF[,i])){
        Levels <- levels(DF[,i])[2:nlevels(DF[,i])]
          for (m in 1:length(Levels)) {
          my_matrix[a,1] <- colnames(DF)[j]
          my_matrix[a,2] <- colnames(DF)[i]
          my_matrix[a,3] <- Levels[m]
          my_matrix[a,4] <- sum(DF[,i] == Levels[m])
          n_samples <- sum(DF[,i] == Levels[m])
          my_matrix[a,5] <- sum(DF[,i] != Levels[m], na.rm = T)
          n_samplesother <- sum(DF[,i] != Levels[m], na.rm = T)
          my_matrix[a,6] <- colSums(subset(DF, DF[c(i)] == Levels[m], select = j) == 1)/n_samples *100 #bacteria prevalence within the group/level 'm'
          my_matrix[a,7] <- colSums(subset(DF, DF[c(i)] != Levels[m], select = j) == 1)/n_samplesother *100 #bacteria prevalence within the groups/levels other than m
          my_matrix[a,8] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria within whole cohort
          a=a+1
          }
      }
      else{
        my_matrix[a,1] <- colnames(DF)[j]
        my_matrix[a,2] <- colnames(DF)[i]
        my_matrix[a,3] <- colnames(DF)[i]
        my_matrix[a,4] <- nrow(DF)
        my_matrix[a,5] <- nrow(DF)
        my_matrix[a,6] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria among the non NA individuals per numerical phenotype 
        my_matrix[a,7] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria among the non NA individuals per numerical phenotype 
        my_matrix[a,8] <- sum(DF[,j])/nrow(DF) *100 #prevalence of bacteria within whole cohort
        a=a+1
      }
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  return(my_matrix)
}

```

```{r}
========================================================================================================================================================
REVISED Manucript functions - 17.12.20
========================================================================================================================================================
```

### 14. Linear_model_2
```{r}

Linear_model_2 <- function(x){
Nlevels <- vector()
NonFactors <- vector()
  for (i in 2:ncol(x)) {
    if(is.factor(x[,i])) {
      Nlevels <- c(Nlevels, (nlevels(x[,i])-1))
    }
    else {
      NonFactors <- c(NonFactors, colnames(x)[i])
    }
    #create matrix 
    my_matrix<-  matrix(nrow = (length(NonFactors) + sum(Nlevels)), ncol = 9)
    colnames(my_matrix) <- c('Intercept', 'Phenotype', 'Estimate', 'Std.Error', 't_value','P.value','Adj.r.squared','BonferroniAdjust','FDRAdjust')
  }


  #Add linear model results to the matrix 
  a=1 
  #Iterate through taxonomy columns
  for (j in 1) {
    my_model <- lm(x[,j] ~ ., data = x[,2:ncol(x)])
    for(k in 2:nrow(summary(my_model)$coefficients)){
      my_matrix[a,1] <- colnames(x)[j]
      my_matrix[a,2] <- rownames(summary(my_model)$coefficients)[k]
      my_matrix[a,3] <- summary(my_model)$coefficients[k,1]
      my_matrix[a,4] <- summary(my_model)$coefficients[k,2]
      my_matrix[a,5] <- summary(my_model)$coefficients[k,3]
      my_matrix[a,6] <- summary(my_model)$coefficients[k,4]
      my_matrix[a,7] <- summary(my_model)$adj.r.squared
      a=a+1
    }
  }
  my_matrix <- as.data.frame(my_matrix)
  my_matrix[,6] <- as.numeric(as.character(my_matrix[,6]))
  my_matrix[,8] <- p.adjust(c(my_matrix[,6]), method = 'bonferroni')
  my_matrix[,9] <- p.adjust(c(my_matrix[,6]), method = 'fdr')
  return(my_matrix)
}

```

### 15. Prevalence Function
```{r}

Prev <- function(x) {sum(x)/length(x)*100
}

```

