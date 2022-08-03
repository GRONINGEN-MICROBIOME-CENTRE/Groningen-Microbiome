
#===========================================================
# Function commonly used in microbial data analysis
#===========================================================

CompositionTable <- function(x,n){ # change the metaphlan result to a composition table, select top n most abundant features
  
  require(foreach)
  x[is.na(x)]=0
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  
  return(composition_table)
}

AverageTable <- function(compositiontable){ # this function is based on ~CompostionTable, to calculate average abundance
  
  averageTable=matrix(nrow = 1,ncol = length(unique(compositiontable$Level)))
  rownames(averageTable)=deparse(substitute(compositiontable))
  colnames(averageTable)=c(as.character(unique(compositiontable$Level)))
  
  for(i in unique(compositiontable$Level)){
    subset=compositiontable$Relative[compositiontable$Level==i]
    avr=mean(subset[!is.na(subset)])
    averageTable[1,i]=avr
  }
  averageTable=as.data.frame(averageTable)
  averageTable$Others=1-sum(averageTable[!is.na(averageTable)])
  averageTable=as.data.frame(t(averageTable))
  averageTable$Taxa=rownames(averageTable)
  
  return(averageTable)
}

# this function is generate barplot connected with microbial composition change
connectedBarplot <- function(dat, color=c(pal_npg("nrc", alpha = 1)(10),"orange"), space=1, alpha=0.2, ...) {  
  b <- barplot(dat, las=2,col=color, space = space, ...)                     
  for (i in seq_len(ncol(dat) - 1)) {     
    lines(c(b[i]+0.5, b[i+1]-0.5), c(0, 0)) 
    
    for (j in seq_len(nrow(dat))) {     
      if (j == 1) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(dat[j,i], dat[j,i+1]))                       
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(0, dat[j,i], dat[j,i+1], 0),               
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j == 2) {                   
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(dat[1,i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], dat[1,i+1]),                       
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
      if (j > 2) {                    
        lines(c(b[i]+0.5, b[i+1]-0.5), c(colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1]))                      
        polygon(c(b[i]+0.5, b[i]+0.5, b[i+1]-0.5, b[i+1]-0.5),                        
                c(colSums(dat[1:(j-1),])[i], colSums(dat[1:j,])[i], colSums(dat[1:j,])[i+1], colSums(dat[1:(j-1),])[i+1]),              
                col=adjustcolor(color[j], alpha.f=alpha))    
      }      
    }          
  }              
}      


# this function is to generate table for plotting composition abundance 
CompositionPlot <- function(metaphlan){
  metaphlan=foreach(i=1:ncol(metaphlan),.combine = rbind) %do% {
    name=colnames(metaphlan)[i]
    tmp.sub=metaphlan[,i,drop=F]
    colnames(tmp.sub)="Abundance"
    tmp.sub$Sample=rownames(tmp.sub)
    tmp.sub$Taxa=name
    
    return.string=tmp.sub
  }
  return(metaphlan)
}


# this function is perform PCoA analysis and return the first two PCo with variance explained
do_PCoA <- function(DistanceMatrix) {
  require(dplyr)
  # Return error if distance matrix is not a distance matrix
  if (class(DistanceMatrix) != "dist") {
    stop("Distance matrix must be an object of type dist. See ?read_dist for a
         way to read a distance matrix from a file into a dist object.")
  }
  
  # Run the PCoA
  PCoA <- cmdscale(DistanceMatrix, k = 2, eig = TRUE)
  
  # Extract PCoA coordinates
  Coordinates <- data.frame(
    Sample = row.names(as.data.frame(PCoA$points)),
    PCoA1 = PCoA$points[,1],
    PCoA2 = PCoA$points[,2],
    row.names = NULL
  ) %>%
    as.tbl()
  
  # Calculate variance explained
  Eigenvalues <- eigenvals(PCoA) 
  Variance <- Eigenvalues / sum(Eigenvalues) 
  Variance1 <- 100 * signif(Variance[1], 2)
  Variance2 <- 100 * signif(Variance[2], 2)
  
  # Return
  Result <- list(
    Coordinates = Coordinates,
    Variance1 = paste(Variance1,"%",sep = ""),
    Variance2 = paste(Variance2,"%",sep = "")
  )
  return(Result)
}

# this function is filter bacterial data based on different levels (phylum, order, class, etc)
SelectLevel <- function(taxa,level){
  require(crayon)
  
  if(level=="phylum"){
    taxa=taxa[,grep("p__",colnames(taxa))]
    taxa=taxa[,grep("c__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"p__")[[1]][2]
    })
    return(taxa)
  }else if(level=="class"){
    taxa=taxa[,grep("c__",colnames(taxa))]
    taxa=taxa[,grep("o__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"c__")[[1]][2]
    })
    return(taxa)
  }else if(level=="order"){
    taxa=taxa[,grep("o__",colnames(taxa))]
    taxa=taxa[,grep("f__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"o__")[[1]][2]
    })
    return(taxa)
  }else if(level=="family"){
    taxa=taxa[,grep("f__",colnames(taxa))]
    taxa=taxa[,grep("g__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"f__")[[1]][2]
    })
    return(taxa)
  }else if(level=="genus"){
    taxa=taxa[,grep("g__",colnames(taxa))]
    taxa=taxa[,grep("s__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"g__")[[1]][2]
    })
    return(taxa)
  }else if(level=="species"){
    taxa=taxa[,grep("s__",colnames(taxa))]
    taxa=taxa[,grep("t__",colnames(taxa),invert = T)]
    colnames(taxa)=lapply(colnames(taxa),function(x){
      strsplit(x,"s__")[[1]][2]
    })
    return(taxa)
  }else{
    cat(yellow("Please choose correct level"))
  }
  
}

SelectLevel.keepName <- function(taxa,level){
  require(crayon)
  
  if(level=="phylum"){
    taxa=taxa[,grep("p__",colnames(taxa))]
    taxa=taxa[,grep("c__",colnames(taxa),invert = T)]
    return(taxa)
  }else if(level=="class"){
    taxa=taxa[,grep("c__",colnames(taxa))]
    taxa=taxa[,grep("o__",colnames(taxa),invert = T)]
    return(taxa)
  }else if(level=="order"){
    taxa=taxa[,grep("o__",colnames(taxa))]
    taxa=taxa[,grep("f__",colnames(taxa),invert = T)]
    return(taxa)
  }else if(level=="family"){
    taxa=taxa[,grep("f__",colnames(taxa))]
    taxa=taxa[,grep("g__",colnames(taxa),invert = T)]
    return(taxa)
  }else if(level=="genus"){
    taxa=taxa[,grep("g__",colnames(taxa))]
    taxa=taxa[,grep("s__",colnames(taxa),invert = T)]
    return(taxa)
  }else if(level=="species"){
    taxa=taxa[,grep("s__",colnames(taxa))]
    taxa=taxa[,grep("t__",colnames(taxa),invert = T)]
    return(taxa)
  }else{
    cat(yellow("Please choose correct level"))
  }
  
}


# this function is to filter taxa based on present rate and relative abundance
Filter <- function(taxa,relative,present){
  taxa = taxa[,colSums(taxa > 0) >= present*nrow(taxa)]
  taxa = taxa[,colMeans(taxa) >= relative]
}

draw_colnames_45 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 45, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_45",
  ns = asNamespace("pheatmap")
)

colors=read.table("taxonomy_color.txt",sep = "\t",stringsAsFactors = F,header = T,comment.char = "")
TaxaColor <- function(taxa,level){
  tmp=colors[colors$name %in% taxa,]
  tmp=tmp[tmp$taxa==level,]
  tmp=tmp[order(match(tmp$name,taxa)),]
  return(setNames(tmp$color,taxa))
}


CountLevel=function(taxa){
  species=SelectLevel(taxa,"species")
  genus=SelectLevel(taxa,"genus")
  family=SelectLevel(taxa,"family")
  order=SelectLevel(taxa,"order")
  class=SelectLevel(taxa,"class")
  phylum=SelectLevel(taxa,"phylum")
  cat(yellow("s__","number is",ncol(species),"\n"))
  cat(yellow("g__","number is",ncol(genus),"\n"))
  cat(yellow("f__","number is",ncol(family),"\n"))
  cat(yellow("o__","number is",ncol(order),"\n"))
  cat(yellow("c__","number is",ncol(class),"\n"))
  cat(yellow("p__","number is",ncol(phylum),"\n"))
}


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







