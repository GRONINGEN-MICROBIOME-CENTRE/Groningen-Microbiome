
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

# this function is to filter taxa based on present rate and relative abundance
Filter <- function(taxa,relative,present){
  taxa = taxa[,colSums(taxa > 0) >= present*nrow(taxa)]
  taxa = taxa[,colMeans(taxa) >= relative]
  taxa = as.data.frame(apply(taxa, 1, function(x){x=x/sum(x)}))
  taxa=as.data.frame(t(taxa))
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
