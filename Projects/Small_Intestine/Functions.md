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
