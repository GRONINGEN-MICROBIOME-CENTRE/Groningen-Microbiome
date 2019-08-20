### Calculate summary statistics of metadata file ###

# Creator: Paula Sureda | Arnau Vich
# Year: 2017

## USAGE ##

# Copy or import function to R #

#Input files
#1.File contaning categorical or numerical phenotypes
#2.(optional)File contaning in the first column sample_ID and in the second column the category name

#Output files
#If category table is provided the script creates a summary table per each category
#If not the script creates a global summary table (for whole table input)


#Example metadata_input
#
#SID        factor1   factor2    factor3
#Sample1     23.5  	   0.9   	   yes
#Sample2     10.9 	   0.01  	   no
#Sample3     50    	   0.3    	 no

#Example category_table
#
#SID 		category
#Sample1	 cat1
#Sample2	 cat1
#Sample3	 cat2

#Example output
#
#factors      Type      Categories/Median   Counts/Mean     %/SD    Number_non-zeros(n)  Number_NA
#factor1    numerical         4.5               4.32        1.22            34              6
#factor2    categorical     UC,CD,IBDU         4,20,5    13.8,69,17.2       29             11
#factor3    numerical         10                 9.34       0.38            38              2


summary_statistics_metadata <- function (metadata_input, category_table) {
  
	# Packages needed       
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
