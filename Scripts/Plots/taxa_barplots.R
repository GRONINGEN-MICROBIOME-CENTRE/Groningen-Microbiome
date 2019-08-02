### Taxonomy composition barplot per category ###

#Creator: Paula Sureda / Arnau Vich

#Input files
#
#1.Taxonomy table at specific level (species, order, family...)
#2.Category table: file contaning in the first column sample_ID and in the second column the category name
#3.Top_tax_value: number of taxonomy levels to be represented in the plot

# Example tax_level_table
#
#SID		SP1		SP2	   SP3
#Sample1	23.5    3.5    5.5
#Sample2    10.9   43.3   34.6
#Sample3    50     10     30

# Example category_table
#
#SID    category
#Sample1  cat1
#Sample2  cat1
#Sample3  cat2

  
barplot_tax_composition <- function(tax_level_table, category_table, top_tax_value) {
 
	# Packages needed
	library(psych)
    library(reshape2) 
    library(ggplot2)
    
    # Merge tax_level_table with category_table file  
	filum_groups <- merge(category_table, tax_level_table, by = "row.names")
	filum_groups <- filum_groups[,-1]
    
    # Split filum_groups table by categories      
    split_categories <- split(filum_groups, filum_groups[,1])
    
    # Save the different categories files in the global environment
    j = 1
    others_position = top_tax_value + 1
    results_table <- matrix(nrow = others_position, ncol= 2 )

    matrix_list <- list()   
         
    # Loop for to calculate the filum mean by category
    for (i in split_categories) {
    	# Calculate mean
    	summary_table <- describe(i)
        # Find the rows with the highest values by the given value
        top_tax <- summary_table[order(summary_table$mean, decreasing = T)[1:top_tax_value],]
        # Sum values and create another category: others
        sum_top_tax <- sum(top_tax$mean)
        others_tax <- as.numeric((100-sum_top_tax))
              
		total_mean <- cbind(top_tax$mean)
		total_mean <- rbind(total_mean, others_tax)

		results_table[,1] <- total_mean
		results_table[,2] <- c(rownames(top_tax), "others")
              
		nam <- paste("Category", j , sep = "")
		matrix_list[[j]] <-  assign(nam, results_table)

		j <- j + 1        
	}
      
    # Change colnames from all matrix in the list
    for (k in seq_along(matrix_list)) {
    	colnames(matrix_list[[k]]) <- c("mean","bacteria")
    }
      
    # Merge function 
    MyMerge <- function(x, y){
    	df <- merge(x, y, by="bacteria" , all.x= TRUE, all.y= TRUE)
    	return(df)
    }
      
    # Merge the matrix saved in the matrix_list in a same table
    composition_table <- Reduce(MyMerge, matrix_list)
    rownames(composition_table) <- composition_table[,1]
    composition_table <- composition_table[,-1]
    
    # For loop to name the columns depending on the groyup/category in the composition_table  
    for (bb in 1:ncol(composition_table)) {
    	nc <- paste("Category", bb , sep = "")
    	colnames(composition_table)[bb] = nc
    	assign(nc,composition_table)
    }
      
    # Export the table with the top4 abundances per each category          
    write.table(composition_table, file = "~/bacteria_composition.txt", sep = "\t", quote = F)

    ## Stacked barplot ##
    composition_table$bacteria = row.names(composition_table) 
    row.names(composition_table) <- NULL

    my_plot_table <- melt(composition_table, id.vars = "bacteria")
      
    filum_plot <- ggplot (my_plot_table, aes(x=variable, y=as.numeric(value))) + geom_bar (aes(fill = bacteria), stat = "identity") + theme_classic() + xlab("Group") + ylab("relative_abundance")
    
    # Save the plots as pdf.file
    
    pdf("taxonomy_composition_plot.pdf")
    print(filum_plot)
    
    dev.off()
}
