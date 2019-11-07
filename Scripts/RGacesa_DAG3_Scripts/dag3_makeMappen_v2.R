# library(rgdal)

library(geojsonio) # for geocoding
library(ggmap) # basic stuff for maps
library(leaflet) # for fancy maps
library(htmltools) # required for leaflet pretty maps
library(sf) # for reverse geocoding
library(dplyr) # for 'group by'
library(vegan)
library(Rtsne)

# cute cute function for PCA plotting regression plotting
# takes pcas table and responsevar (has to be in pcas dataset)
# ===========================================
plotPCAs <- function(pcas,nrPCs=3,responseVar,doCentroids=T) {
  # PC 1-2
  g <- ggplot(pcas ,aes_string(x = "PC1",y="PC2",col=responseVar )) + geom_point(size=2.25,alpha=0.8) 
  pcas$rV <- pcas[[responseVar]]
  if (doCentroids) {
    centroids <- aggregate(cbind(PC1,PC2)~rV,pcas,mean)
    colnames(centroids)[[1]] <- responseVar
    g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)
  }
  print(g)
  #ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_12.png',sep=''),width = sX,height = sY)
  if (nrPCs >= 3) {
    # PC 2-3
    g <- ggplot(pcas ,aes_string(x = "PC2",y="PC3",col=responseVar )) + geom_point(size=2.25,alpha=0.8) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC2,PC3)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC2,y=centroids$PC3))
    }
    print(g)
    #ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_23.png',sep=''),width = sX,height = sY)
  }
  if (nrPCs >= 4) {
    # PC 3-4
    g <- ggplot(pcas ,aes_string(x = "PC3",y="PC4",col=responseVar )) + geom_point(size=2.25,alpha=0.4) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC3,PC4)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC3,y=centroids$PC4))
    }
    print(g)
  }
  if (nrPCs >= 5) {
    # PC 4-5
    g <- ggplot(pcas ,aes_string(x = "PC4",y="PC5",col=responseVar )) + geom_point(size=2.25,alpha=0.4) 
    if (doCentroids) {
      centroids <- aggregate(cbind(PC4,PC5)~rV,pcas,mean)
      colnames(centroids)[[1]] <- responseVar    
      g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC4,y=centroids$PC5))
    }
    print(g)
  }
}

# cute cute function for PCA plotting regression plotting
# takes pcas table and responsevar (has to be in pcas dataset)
# ===========================================
plotTSNE <- function(tSNE,responseVar,doCentroids=T) {
  # PC 1-2
  g <- ggplot(pcas ,aes_string(x = "PC1",y="PC2",col=responseVar )) + geom_point(size=2.25,alpha=0.4) 
  pcas$rV <- pcas[[responseVar]]
  if (doCentroids) {
    centroids <- aggregate(cbind(PC1,PC2)~rV,pcas,mean)
    colnames(centroids)[[1]] <- responseVar
    g <- g + geom_point(data=centroids,shape=4,stroke=3,size=5,aes(x=centroids$PC1,y=centroids$PC2),alpha=1)
  }
  print(g)
  #ggsave(plot=g,paste(outFolderName,'/exploratory/diag5_bc_12.png',sep=''),width = sX,height = sY)
}

# cute cute function for linear regression plotting
# takes linear model fit as argument
# ===========================================
ggplotRegression <- function (fit) {
  require(ggplot2)
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point() +
    stat_smooth(method = "lm", col = "blue") +
    labs(title = paste("Adj R2 = ",signif(summary(fit)$adj.r.squared, 2),
                       #"Intercept =",signif(fit$coef[[1]],2 ),
                       " Coef =",signif(fit$coef[[2]], 2),
                       " P =",signif(summary(fit)$coef[2,4], 2)))
}

# =================================================================
# =================================================================
# MAIN 
# =================================================================
# =================================================================
setwd('C:/Users/ranko/Documents/UMCG/DAG3_stats/')
source('leafletPatch.R')
source('../myLibs/R_Microbiome_scripts.R')
source('../myLibs/R_ML_scripts_v3.R')

# # # ============ PREP METAGENOMIC DATA ========================================================================================================
# #  go over our samples, sorted by ???, give them location
# taxRaw <- as.data.frame(t(read.table('data/DAG3_metaphlan_merged.txt',sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"')))
# taxRaw$DAG3_sampleID <- rownames(taxRaw)
# taxRaw$DAG3_sampleID <- gsub("_metaphlan","",taxRaw$DAG3_sampleID)
# # # filter Taxonomy Data (0.1% prevalence, 0.0001 abundance - pretty liberal filter)
# taxRawF <- filterMetaGenomeDF(inDF = taxRaw,keepLevels = c("K","P","C","O","F","G","S"),presPerc = 0.001,minMRelAb = 0.0001,rescaleTaxa = F)
# # # remove "bad samples"(bacteria < 0.95,)
# taxRawF <- taxRawF[taxRawF$k__Bacteria >= 0.95,]
# taxRawF$unclassified <- NULL
# taxRdy <- filterMetaGenomeDF(inDF = taxRawF,keepLevels = c("P","C","F","O","G","S"),presPerc = -1,minMRelAb = -1,rescaleTaxa = T)
# write.table(taxRdy,sep='\t',row.names = F,file = 'data_processed/DAG3_metaphlan_filtered_cleaned_rescaled_nontransformed.tsv')
# # # unfiltered Taxonomy Data
# taxRaw <- as.data.frame(t(read.table('data/DAG3_metaphlan_merged.txt',sep='\t',header = T,row.names = 1,stringsAsFactors = F,quote = '"')))
# taxRaw$DAG3_sampleID <- rownames(taxRaw)
# taxRaw$DAG3_sampleID <- gsub("_metaphlan","",taxRaw$DAG3_sampleID)
# taxRawF <- filterMetaGenomeDF(inDF = taxRaw,keepLevels = c("K","P","C","O","F","G","S"),presPerc = -1,minMRelAb = -1,rescaleTaxa = F)
# # # remove "bad samples"(bacteria < 0.95,)
# taxRawF <- taxRawF[taxRawF$k__Bacteria >= 0.95,]
# taxRawF$unclassified <- NULL
# taxRdy <- filterMetaGenomeDF(inDF = taxRawF,keepLevels = c("P","C","O","F","G","S"),presPerc = -1,minMRelAb = -1,rescaleTaxa = T)
# write.table(taxRdy,sep='\t',row.names = F,file = 'data_processed/DAG3_metaphlan_unfiltered_cleaned_rescaled_nontransformed.tsv')
# # # ============ END OF PREP METAGENOMIC DATA ===================================================================================================

# > load geolocated data
inDF <- read.table(file = 'data_processed/DAG3_geolocated.tsv', header = T,sep='\t',stringsAsFactors = F)
inDF$codeForGoogle <- NULL
#   > remove bad data
rdyDFall <- inDF[!is.na(inDF$lon) & !is.na(inDF$lat),]

# > LOAD METAGENOMIC DATA
taxRdy <- read.table('data_processed/DAG3_metaphlan_filtered_cleaned_rescaled_nontransformed.tsv',header=T,sep='\t',stringsAsFactors = F)
#   > keep only samples for which we have "high quality" MGS
dag3IDsToKeep <- taxRdy$DAG3_sampleID
dag3overlap <- rdyDFall$DAG3_sampleID[rdyDFall$DAG3_sampleID %in% dag3IDsToKeep]
dag3nonmatch <- rdyDFall$DAG3_sampleID[!rdyDFall$DAG3_sampleID %in% dag3IDsToKeep]
dag3extra <- dag3IDsToKeep[!dag3IDsToKeep %in% rdyDFall$DAG3_sampleID]
#   > note: we removed pilot here!
rdyDFall <- rdyDFall[rdyDFall$DAG3_sampleID %in% dag3IDsToKeep,]

# PREP MUNICIPALY DATA
# =====================
dutchTowns <- geojsonio::geojson_read("data_mappen/dutchTownships.geojson",what = "sp")
# EXTRACT CENTERS
lats <- c()
longs <- c()
for (r in c(1:nrow(dutchTowns))) {
  longs <- c(longs,dutchTowns[r,1]@polygons[[1]]@labpt[[1]])
  lats <- c(lats,dutchTowns[r,1]@polygons[[1]]@labpt[[2]])
}
dutchTowns$lats <- lats
dutchTowns$longs <- longs

# get metadata
metaTowns <- read.table('data_mappen/dutch_stats_2015.csv',sep=',',header=T,quote = '"')
tt <- as.character(dutchTowns$name)
rownames(metaTowns) <- metaTowns$NAME
metaTownsJ <- metaTowns[tt,]
dutchTowns2 <- dutchTowns
dutchTowns2$NAME <- metaTownsJ$NAME
dutchTowns2$POP.TOTAL <- metaTownsJ$POP.TOTAL
dutchTowns2$POP.KM2 <- metaTownsJ$POP.KM2
dutchTowns2$URBAN.INDEX <- metaTownsJ$URBAN.INDEX
# define Urban-Index
urban.index.def = c("Very High Density Urban","High Density Urban","Moderate Density Urban","Low Density Urban","Rural","No Data")
dutchTowns2$URBAN.INDEX.F <-  urban.index.def[dutchTowns2$URBAN.INDEX]
dutchTowns2$URBAN.INDEX.F[is.na(dutchTowns2$URBAN.INDEX.F)] <- "No Data"
dutchTowns2$URBAN.INDEX.F <- as.factor(dutchTowns2$URBAN.INDEX.F)
dutchTowns2$PROVINCE <- metaTownsJ$PROVINCE
write.table(metaTownsJ,'data_mappen/dutch_stats_2015_cleaned.tsv',sep='\t',row.names = F)

# do reverse-geocoding
# prep coordinatesrdyDFall.xy
rdyDFall.xy <- rdyDFall
coordinates(rdyDFall.xy) <- ~lon+lat
dutchTownsSF <- st_as_sf(dutchTowns2)
# collector for results
print ("REVERSE GEOCODING !")
smplsAnnot <- data.frame()
#  > iterate over samples
options(warn=-1)
rdyDFall$municipaly <- NA
for (r in c(1:nrow(rdyDFall))) {
  print (paste(' > reverse geocoding sample',r))
  # pick sample
  smp <- rdyDFall[r,]
  # make geometry point from it
  tPt = st_as_sf(smp, coords = c("lon", "lat"), crs = 4326, agr = "constant")
  # intersect with dutch geometry
  isect <- st_intersection(tPt, dutchTownsSF)
  if (nrow(isect) == 1) {
  # extract municipaly and add to sample
    smp$municipaly <- isect$name
    smplsAnnot <- rbind.data.frame(smplsAnnot,data.frame(lon=smp$lon,lat=smp$lat,municipaly=smp$municipaly))
    rdyDFall$municipaly[[r]] <- as.character(smp$municipaly)
  } else {
    print ("WARNING: can't find municipaly - weird stuff or not in NL")
    print (" >> trying quick and dirty fix")
    # try to move it a bit and repeat
    success = F
    for (i in c(1:50)) {
      smp$lon = smp$lon+runif(1,min=0,max=0.01)
      smp$lat = smp$lat+runif(1,min=0,max=0.01)
      tPt = st_as_sf(smp, coords = c("lon", "lat"), crs = 4326, agr = "constant")
      isect <- st_intersection(tPt, dutchTownsSF)
      if (nrow(isect) == 1) {
        print(paste('   >>> SUCCESS <small move - up to 1 km>, assigned to',isect$name))
        smp$municipaly <- isect$name
        smplsAnnot <- rbind.data.frame(smplsAnnot,data.frame(lon=smp$lon,lat=smp$lat,municipaly=smp$municipaly))
        success = T
        rdyDFall$municipaly[[r]] <- as.character(smp$municipaly)
        break
      }
    }
    if (nrow(isect) == 1 & !success) {
      for (i in c(1:50)) {
        smp$lon = smp$lon+runif(1,min=0,max=0.1)
        smp$lat = smp$lat+runif(1,min=0,max=0.1)
        tPt = st_as_sf(smp, coords = c("lon", "lat"), crs = 4326, agr = "constant")
        isect <- st_intersection(tPt, dutchTownsSF)
        if (nrow(isect) == 1) {
          print(paste('   >>> SUCCESS <bigger move - up to 10km>, assigned to',isect$name))
          smp$municipaly <- isect$name
          smplsAnnot <- rbind.data.frame(smplsAnnot,data.frame(lon=smp$lon,lat=smp$lat,municipaly=smp$municipaly))
          success = T
          rdyDFall$municipaly[[r]] <- as.character(smp$municipaly)
          break
        }
      }
    }
    if (!success) {
      print ('  >>> FAILED! assigning NA as municipaly')
      print(smp)
      smplsAnnot <- rbind.data.frame(smplsAnnot,data.frame(lon=smp$lon,lat=smp$lat,municipaly=NA))
    }
  }
}
print (paste(' >> successfully reverse-geocoded',sum(!is.na(smplsAnnot$municipaly)),'samples;',sum(is.na(smplsAnnot$municipaly)),'samples failed!'))
options(warn=0)
write.table(smplsAnnot,'data_processed/DAG3_reverse_geocoded_summary.tsv',sep='\t',row.names = F)
write.table(rdyDFall,'data_processed/DAG3_reverse_geocoded.tsv',sep='\t',row.names = F)

# exploratory plot and merge with map
# =========================================
smplsAnnotS <- smplsAnnot
smplsAnnotS$municipaly <- as.factor(as.character(smplsAnnotS$municipaly))
# summarize per municipaly
smplsAnnotT <- smplsAnnotS %>%
  group_by(municipaly) %>%
  summarise(no_rows = length(municipaly))
smplsAnnotT <- as.data.frame(smplsAnnotT[order(smplsAnnotT$no_rows,decreasing = T),])
smplsAnnotT$municipaly <- factor(smplsAnnotT$municipaly, levels = smplsAnnotT$municipaly[order(smplsAnnotT$no_rows,decreasing = T)])
ggplot(smplsAnnotT[order(smplsAnnotT$no_rows,decreasing = T),][1:20,],
       aes(x=municipaly,y=no_rows,col=municipaly,fill=municipaly)) + geom_col()

# add NR of samples per province:
# >> handle 0 samples case
smplsAnnotT$municipaly <- as.character(smplsAnnotT$municipaly)
uNames <- as.character(unique(dutchTowns2$name))
c = 0
for (u in uNames) {
  c <- c+1
  print(c)
  if (!u %in% smplsAnnotT$municipaly) {
    print (paste(' no ',u,'; adding 0',sep=''))
    smplsAnnotT <- rbind.data.frame(smplsAnnotT,data.frame(municipaly=u,no_rows=0))
  }
}
# reorder stuff to match dutchtowns2
tt <- as.character(dutchTowns2$name)
rownames(smplsAnnotT) <- smplsAnnotT$municipaly
smplsAnnotTJ <- smplsAnnotT[tt,]
# merge it
dutchTowns2$samplesNR <- smplsAnnotTJ$no_rows
# save it
saveRDS(dutchTowns2,file="data_mappen/dutchTownsMapRdy.Rdata")

# ===========================================
# ===========================================
# link real metagenomic data to map !!!
# ===========================================
# ===========================================

# ==================== LOAD DATA ======================
taxRdy <- read.table('data_processed/DAG3_metaphlan_filtered_cleaned_rescaled_nontransformed.tsv',header=T,sep='\t',stringsAsFactors = F)
dutchTowns2 <- readRDS(file="data_mappen/dutchTownsMapRdy.Rdata")
smplsAnnot <- read.table('data_processed/DAG3_reverse_geocoded.tsv',sep='\t',header = T)
smplsAnnot$PSEUDOIDEXT <- NULL
smplsAnnot$POSTCODE <- NULL

source('dag3helper.R')
#  >> load SPECIES
dataS <- taxRdy[,grep('s__',colnames(taxRdy))]
dataS$DAG3_sampleID <- taxRdy$DAG3_sampleID
dataSann <- merge(smplsAnnot,dataS,by="DAG3_sampleID")
#  >> laod GENERA
dataG <- taxRdy[,-grep('s__',colnames(taxRdy))]
dataG <- taxRdy[,grep('g__',colnames(taxRdy))]
dataG$DAG3_sampleID <- taxRdy$DAG3_sampleID
dataGann <- merge(smplsAnnot,dataG,by="DAG3_sampleID")
# > calculate diversity
dataDiv <- getDiversityForGementee(dataSann)
# add to dutchTowns2
dutchTowns2$diversity <- 0.0
dutchTowns2$diversity.sd <- 0.0
c <- 0
for (n in dutchTowns2$NAME) {
    c <- c + 1
    #print(n)
    #print(dataDiv$Diversity[dataDiv$Gementee==n])
    if (n %in% dataDiv$Gementee) {
      dutchTowns2$diversity[c] <- dataDiv$Diversity[dataDiv$Gementee==n]
      dutchTowns2$diversity.sd[c] <- dataDiv$Diversity.SD[dataDiv$Gementee==n]
    }
}
dutchTowns2$diversity.sd[is.na(dutchTowns2$diversity.sd)] <- 0.0

# > calculate most abundant SPECIES
if (!file.exists('data_mappen/abundance_topSpecies.csv')) {
  ab <- getAbundancesForGementee(dataSann,getTopX = 1)
  ab$Taxon <- shortenNames2(ab$Taxon)
  write.table(ab,'data_mappen/abundance_topSpecies.csv',sep=',',row.names = F)
} else {
  ab <- read.table('data_mappen/abundance_topSpecies.csv',sep=',',header = T)
}
ab$Gementee <- as.character(ab$Gementee)
ab$Taxon <- as.character(ab$Taxon)
#  >> add to dutchTowns2
c = 0
dutchTowns2$dominantSpec <- NA
for (n in as.character(dutchTowns2$NAME)) {
  c <- c + 1
  if (n %in% dataDiv$Gementee) {
    dutchTowns2$dominantSpec[c] <- as.character(ab$Taxon[ab$Gementee==n])
  }
}

# > calculate most abundant GENUS
if (!file.exists('data_mappen/abundance_topGenus.csv')) {
  abG <- getAbundancesForGementee(dataGann,getTopX = 1)
  abG$Taxon <- shortenNames2(abG$Taxon)
  write.table(abG,'data_mappen/abundance_topGenus.csv',sep=',',row.names = F)
} else {
  abG <- read.table('data_mappen/abundance_topGenus.csv',sep=',',header = T)
}
abG$Gementee <- as.character(abG$Gementee)
abG$Taxon <- as.character(abG$Taxon)
#  >> add to dutchTowns2
c = 0
dutchTowns2$dominantGen <- NA
for (n in as.character(dutchTowns2$NAME)) {
  c <- c + 1
  if (n %in% dataDiv$Gementee) {
    dutchTowns2$dominantGen[c] <- as.character(abG$Taxon[abG$Gementee==n])
  }
}

# > calculate  ...


# add township HTML link
# dutchTowns2$townshipHtmlLink <- paste('Details/',dutchTowns2$name,'_details',sep='')
# dutchTowns2$townshipHtmlLink <- gsub("[^0-9A-Za-z///' ]","_", dutchTowns2$townshipHtmlLink ,ignore.case = TRUE)
# dutchTowns2$townshipHtmlLink <- gsub(" ","_", dutchTowns2$townshipHtmlLink ,ignore.case = TRUE)
# dutchTowns2$townshipHtmlLink <- gsub("'","_", dutchTowns2$townshipHtmlLink ,ignore.case = TRUE)
# dutchTowns2$townshipHtmlLink <- gsub("__","_", dutchTowns2$townshipHtmlLink ,ignore.case = TRUE)
# dutchTowns2$townshipHtmlLink <- gsub("__","_", dutchTowns2$townshipHtmlLink ,ignore.case = TRUE)
# dutchTowns2$townshipHtmlLink <- tolower(paste(dutchTowns2$townshipHtmlLink,'.html',sep=''))
# # add province HTML links
# dutchTowns2$provinceHtmlLink <- paste('Details/',dutchTowns2$PROVINCE,'_details',sep='')
# dutchTowns2$provinceHtmlLink <- gsub("[^0-9A-Za-z///' ]","_", dutchTowns2$provinceHtmlLink ,ignore.case = TRUE)
# dutchTowns2$provinceHtmlLink <- gsub(" ","_", dutchTowns2$provinceHtmlLink ,ignore.case = TRUE)
# dutchTowns2$provinceHtmlLink <- gsub("'","_", dutchTowns2$provinceHtmlLink ,ignore.case = TRUE)
# dutchTowns2$provinceHtmlLink <- gsub("__","_", dutchTowns2$provinceHtmlLink ,ignore.case = TRUE)
# dutchTowns2$provinceHtmlLink <- gsub("__","_", dutchTowns2$provinceHtmlLink ,ignore.case = TRUE)
# dutchTowns2$provinceHtmlLink <- tolower(paste(dutchTowns2$provinceHtmlLink,'.html',sep=''))

saveRDS(dutchTowns2,file="data_mappen/dutchTownsMapRdy.Rdata")

# make cute cute pretty map of dutch municipalies
# ==================================================
# > colors
palN <- colorNumeric("viridis", domain = NULL,reverse = T,alpha = T)
palNR <- colorNumeric("viridis", domain = NULL,reverse = F,alpha = T)
palF <- colorFactor("viridis",NULL)
palFU <- colorFactor("YlOrRd",NULL)

# > make pretty pretty labels
labsPoly <- lapply(seq(nrow(dutchTowns2)), function(i) {
  divStr <- "No Data"
  if ( (dutchTowns2$diversity[i] > 0.0) & (dutchTowns2$diversity.sd[i] > 0.0) ) {
    divStr <- paste0( format2D(dutchTowns2$diversity[i]),' +/- ',format2D(dutchTowns2$diversity.sd[i]) )
  } else if (dutchTowns2$diversity[i] > 0.0) {
    divStr <- paste0(format2D(dutchTowns2$diversity[i]))
  }
  paste0( '<b>',dutchTowns2$name[i],', ',dutchTowns2$PROVINCE[i],'</b><br>', 
          'Population: ',dutchTowns2$POP.TOTAL[i], '<br>', 
          'Density: ',dutchTowns2$POP.KM2[i],'/km2<br>', 
          'Urban-Index: ',dutchTowns2$URBAN.INDEX[i],' (',dutchTowns2$URBAN.INDEX.F[i],')<br>',
          'Microbiome samples: ',dutchTowns2$samplesNR[i],'<br>',
          'Shannon-diversity: ',divStr,'<br>',
          'Most abundant Species: ',dutchTowns2$dominantSpec[i],'<br>',
          'Most abundant Genus: ',dutchTowns2$dominantGen[i],''
          ) 
})
labsCircle <- lapply(seq(nrow(dutchTowns2)), function(i) {
  paste0( '<b>',dutchTowns2$name[i], '</b><br>', 
          'Microbiome samples:',dutchTowns2$samplesNR[i], '<br>')
          })

# subsample markers by <sqrt>
markerSubS <- data.frame()
for (u in unique(smplsAnnot$municipaly)) {
  nrPick <- ceiling(sqrt(sum(smplsAnnot$municipaly==u)))
  print(paste(' >> grabbing',nrPick,'samples from',u))
  sub <- smplsAnnot[smplsAnnot$municipaly==u,]
  markerSubS <- rbind.data.frame(markerSubS,sample_n(sub, nrPick))
}

# make map
# ==========================
# prep for coloring
dutchTowns2$diversity[dutchTowns2$diversity==0] <- NA
dutchTowns2$samplesNR[dutchTowns2$samplesNR==0] <- NA

l <- leaflet(dutchTowns2) %>% 
  # tiles
  addProviderTiles(providers$Esri.NatGeoWorldMap,options = providerTileOptions(minZoom = 7, maxZoom = 12)) %>%
  addTiles() %>%
  # polygons (urban index)
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.4,
              fillColor = ~palF(URBAN.INDEX),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.75,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$townshipHtmlLink,'>Details</a>',sep=''),
              group = "Urban-Index"
              ) %>%
  # polygons (provinces)
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.4,
              fillColor = ~palF((PROVINCE)),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.75,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$provinceHtmlLink,'>Details</a>',sep=''),
              #popup= '<a href = "https://rstudio.github.io/leaflet/">R</a>',
              group = "Provinces"
  ) %>%
  # polygons (NR samples)
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.6,
              fillColor = ~palN(sqrt(samplesNR)),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.85,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$townshipHtmlLink,'>Details</a>',sep=''),
              group = "Samples"
  ) %>%
  #diversity
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.6,
              fillColor = ~palN(sqrt(diversity)),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.85,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$townshipHtmlLink,'>Details</a>',sep=''),
              group = "Diversity"
  ) %>%
  #most abundant species
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.6,
              fillColor = ~palFU(dominantSpec),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.85,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$townshipHtmlLink,'>Details</a>',sep=''),
              group = "Species"
  ) %>%
  #most abundant Genus
  addPolygons(stroke = T,weight = 0.5,
              smoothFactor = 0.3,
              fillOpacity=0.6,
              fillColor = ~palFU(dominantGen),
              label = lapply(labsPoly, HTML),
              labelOptions = labelOptions(textsize = "15px"),
              highlightOptions = highlightOptions(fillOpacity=0.85,weight=2.0),
              #popup= paste('<a href=',dutchTowns2$townshipHtmlLink,'>Details</a>',sep=''),
              group = "Genus"
  ) %>%
  
  # add legend (urban index)
  addLegend(position=c("topleft"),title="Urban Index:",labels=unique(dutchTowns2$URBAN.INDEX.F),
            colors = unique(palF(dutchTowns2$URBAN.INDEX)),
            group="Urban-Index") %>%
  # addLegend(position=c("topright"),title="Govan:",labels=unique(dutchTowns2$URBAN.INDEX.F),
  #           colors = unique(palN(dutchTowns2$URBAN.INDEX)),
  #           group="Provinces") %>%
  #hideGroup("Provinces") %>% hideGroup("Urban-Index") %>%
  # markers for each sample
  addMarkers(lng = markerSubS$lon,lat = markerSubS$lat,
             group = c("Markers")) %>%
  # poop markers
  addMarkers(lng = smplsAnnot$lon,lat = smplsAnnot$lat,clusterOptions = markerClusterOptions(),
             icon = makeIcon("data_mappen/poo-solid.svg","data_mappen/poo-solid.svg",16,16),
             group = c("Samples","Provinces","Urban-Index")) %>%
  # markers (no poop)
  #addMarkers(lng = smplsAnnot$lon,lat = smplsAnnot$lat,clusterOptions = markerClusterOptions(),
  #           icon = makeIcon("data_mappen/_ionicons_svg_ios-beaker.svg","data_mappen/_ionicons_svg_ios-beaker.svg",16,16),
  #           group = c("Samples","Provinces","Urban-Index")) %>%
  # layer control
  addLayersControl(overlayGroups = c("Urban-Index", "Provinces","Samples","Markers","Diversity","Species","Genus"),
     options = layersControlOptions(collapsed = FALSE)) %>% 
  # hide ticks
  hideGroup("Provinces") %>% hideGroup("Samples") %>% hideGroup("Markers") %>% hideGroup("Diversity") %>% 
  hideGroup("Species") %>% hideGroup("Genus")

setwd('./mappen/')
saveWidget(l, file="mappen_real_v2.html")
setwd('..')

# # awesome icon
# icons <- awesomeIcons(
#   icon = 'beaker',
#   iconColor = 'white',
#   library = 'ion'
#   #markerColor = getColor(df.20),
#   #spin=T
# )

# =================================================
# ========== do some statistics =============
# =================================================
# > load data and prep it
metaTownsJ <- read.table('data_mappen/dutch_stats_2015_cleaned.tsv',sep='\t',header=T)
taxRdy <- read.table('data_processed/DAG3_metaphlan_unfiltered_cleaned_rescaled_nontransformed.tsv',header=T,sep='\t',stringsAsFactors = F)
dutchTowns2 <- readRDS(file="data_mappen/dutchTownsMapRdy.Rdata")
smplsAnnot <- read.table('data_processed/DAG3_reverse_geocoded.tsv',sep='\t',header = T)
smplsAnnot$PSEUDOIDEXT <- NULL
#smplsAnnot$POSTCODE <- NULL
#smplsAnnot$GEBPLTS <- NULL
#smplsAnnot$lon <- NULL
#smplsAnnot$lat <- NULL
metaTownsJf <- metaTownsJ[,c(1,5,8,9,10,11)]
colnames(metaTownsJf) <- c("municipaly","POP.KM2","HOUSES.KM2","URBAN.INDEX","ADDRESS.DENSITY","PROVINCE")
# > prep microbiome
taxRdyF <- filterMetaGenomeDF(taxRdy,keepLevels = c("P","S","G"), rescaleTaxa = T,presPerc = -1,minMRelAb = -1)
taxRdyFann <- merge(taxRdyF,smplsAnnot,by="DAG3_sampleID")
taxRdyFann <- merge(taxRdyFann,metaTownsJf,by="municipaly")
taxRdyFann$URBAN.INDEX <- as.factor(taxRdyFann$URBAN.INDEX)
taxRdyFann[,grep('__',colnames(taxRdyFann))] <- asin(sqrt(taxRdyFann[,grep('__',colnames(taxRdyFann))]))
taxRdyFann <- purgeMGNames(taxRdyFann)
taxRdyFann$LOG.POP.KM2 <- log(taxRdyFann$POP.KM2)
taxRdyFann$LOG.HOUSES.KM2 <- log(taxRdyFann$HOUSES.KM2)
taxRdyFann$LOG.ADDRESS.DENSITY <- log(taxRdyFann$ADDRESS.DENSITY)
# ==========================================================

# > do tests
#    > taxa vs Urban index
responseVar <- "URBAN.INDEX"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  krus <- kruskal.test(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,Krus.CSq=krus$statistic[[1]],Krus.df=krus$parameter[[1]],Krus.p=krus$p.value))
}
resUrbanIndex <- res[order(res$Krus.p),]
resUrbanIndex$FDR <- p.adjust(resUrbanIndex$Krus.p,method="fdr")
resUrbanIndex[c(1:10),]
#tax = "k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Ruminococcaceae.g__Faecalibacterium.s__Faecalibacterium_prausnitzii"
for (tax in resUrbanIndex$TAXON[resUrbanIndex$FDR <= 0.005]) {
  print(ggplot(taxRdyFann,aes_string(x="URBAN.INDEX",y=tax,col="URBAN.INDEX")) + geom_violin() + geom_boxplot(width=0.2))
}
# =================================================

#    > taxa vs Urban index (numeric, linear)
taxRdyFann$URBAN.INDEX.N <- as.numeric(as.character(taxRdyFann$URBAN.INDEX))
responseVar <- "URBAN.INDEX.N"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  #krus <- kruskal.test(frm,taxRdyFann)
  mdl <- lm(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,adj.r.sq=summary(mdl)$adj.r.squared,
                                         resp.coef=summary(mdl)$coefficients[2,1],
                                         resp.pvalue=summary(mdl)$coefficients[2,4]))
}
resUrbanIndexN <- res[order(res$resp.pvalue),]
resUrbanIndexN$FDR <- p.adjust(resUrbanIndexN$resp.pvalue,method="fdr")
resUrbanIndexN[1:10,]

for (tax in resUrbanIndexN$TAXON[resUrbanIndexN$FDR <= 0.01]) {
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  print(ggplotRegression(mdl))
}
# =================================================

#    > taxa vs population density
responseVar <- "LOG.POP.KM2"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  #krus <- kruskal.test(frm,taxRdyFann)
  mdl <- lm(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,adj.r.sq=summary(mdl)$adj.r.squared,
                                         resp.coef=summary(mdl)$coefficients[2,1],
                                         resp.pvalue=summary(mdl)$coefficients[2,4]))
}
resPopDensity <- res[order(res$resp.pvalue),]
resPopDensity$FDR <- p.adjust(resPopDensity$resp.pvalue,method="fdr")
resPopDensity[1:10,]

for (tax in resPopDensity$TAXON[resPopDensity$FDR <= 0.005]) {
  #print(ggplot(taxRdyFann,aes_string(x=responseVar,y=tax)) + geom_point() + geom_smooth(method="lm"))
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  print(ggplotRegression(mdl))
}
# =======================================================

#    > taxa vs house density
responseVar <- "LOG.HOUSES.KM2"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,adj.r.sq=summary(mdl)$adj.r.squared,
                                         resp.coef=summary(mdl)$coefficients[2,1],
                                         resp.pvalue=summary(mdl)$coefficients[2,4]))
}
resHouseDensity <- res[order(res$resp.pvalue),]
resHouseDensity$FDR <- p.adjust(resHouseDensity$resp.pvalue,method="fdr")
resHouseDensity[1:10,]

for (tax in resHouseDensity$TAXON[resHouseDensity$FDR <= 0.005]) {
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  print(ggplotRegression(mdl))
}
# ===================================================

#    > taxa vs population density
responseVar <- "LOG.ADDRESS.DENSITY"
res <- data.frame()
for (c in grep('__',colnames(taxRdyFann))) {
  tax <- colnames(taxRdyFann)[[c]]
  print (paste(' > testing',tax,'VS',responseVar))
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  res <- rbind.data.frame(res,data.frame(TAXON=tax,adj.r.sq=summary(mdl)$adj.r.squared,
                                         resp.coef=summary(mdl)$coefficients[2,1],
                                         resp.pvalue=summary(mdl)$coefficients[2,4]))
}
resAddressDensity <- res[order(res$resp.pvalue),]
resAddressDensity$FDR <- p.adjust(resAddressDensity$resp.pvalue,method="fdr")
resAddressDensity[1:10,]

for (tax in resAddressDensity$TAXON[resAddressDensity$FDR <= 0.05]) {
  frm <- reformulate(responseVar,tax)
  mdl <- lm(frm,taxRdyFann)
  print(ggplotRegression(mdl))
}

# =================================================================
# do PCA species based) # rda, capscale (vegan)
# =================================================================
# ===================================================================
if (file.exists('data_processed/DAG3_uf_PCA_1to5.tsv')) {
  pcas <- read.table('data_processed/DAG3_uf_PCA_1to5.tsv',sep='\t',header=T)
} else {
  taxRdyForPCA <- taxRdyFann[,grep('s__',colnames(taxRdyFann))]
  bcurtis <- vegdist(taxRdyForPCA,method = "bray")
  r.pca <- prcomp(bcurtis, center = T,scale. = F)
  saveRDS(r.pca,file = 'data_processed/DAG3_uf_BC_PCA.RDS')
  # load 
  r.pca <- readRDS(file= 'data_processed/DAG3_uf_BC_PCA.RDS')
  pcas <- as.data.frame(r.pca$rotation[,1:5])
  pcas$sample <- taxRdyFann$DAG3_sampleID
  pcas$sample <- taxRdyFann$DAG3_sampleID
  qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
  qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
  pcas <- merge(pcas,qc1,by="sample")
  pcas$URBAN.INDEX <- taxRdyFann$URBAN.INDEX
  pcas$URBAN.INDEX.N <- as.numeric(as.character(taxRdyFann$URBAN.INDEX))
  pcas$PROVINCE <- taxRdyFann$PROVINCE
  pcas$municipaly <- taxRdyFann$municipaly
  pcas$lon <- taxRdyFann$lon
  pcas$lat <- taxRdyFann$lat
  pcas$log.lon <- log(pcas$lon)
  pcas$log.lat <- log(pcas$lat)
  pcas$LOG.POP.KM2 <- taxRdyFann$LOG.POP.KM2
  pcas$LAND.OF.BIRTH <- taxRdyFann$GEBLND
  pcas$PLACE.OF.BIRTH <- taxRdyFann$GEBPLTS
  qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned_v2.csv',sep=',',header=T)
  qc2$id.external <- NULL
  qc2$sample.nr <- NULL
  pcas <- merge(pcas,qc2,by="sample")
  write.table(pcas,file = 'data_processed/DAG3_uf_PCA_1to5.tsv',sep='\t',row.names = F)
}

plotPCAs(pcas,nrPCs=3,responseVar = "URBAN.INDEX")
plotPCAs(pcas,nrPCs=3,responseVar = "postclean.reads",doCentroids = F)
plotPCAs(pcas,nrPCs=3,responseVar = "age")
plotPCAs(pcas,nrPCs=3,responseVar = "gender")

plotPCAs(pcas[pcas$PROVINCE %in% c("Drenthe","Groningen","Friesland") ,],nrPCs=3, responseVar = "PROVINCE")
plotPCAs(pcas,nrPCs=3,responseVar = "lon",doCentroids=F)

# read NR has effect on PC 1 !!
# =======================================
lmFit <- lm(PC1 ~ postclean.reads, pcas)
ggplotRegression(lmFit)
#    > not on PC 2
lmFit <- lm(PC2 ~ postclean.reads, pcas)
ggplotRegression(lmFit)
# is it outliers (no)
pcas.mid <- pcas[pcas$postclean.reads < 45000000 & pcas$postclean.reads > 15000000,]
lmFit <- lm(PC1 ~ postclean.reads, pcas.mid)
g <- ggplotRegression(lmFit)
ggsave(plot = g,filename = "dag3_uncorr_QC_reads_vs_PC1.png")

lmFit <- lm(PC3 ~ postclean.reads, pcas.mid)
g <- ggplotRegression(lmFit)
ggsave(plot = g,filename = "dag3_uncorr_QC_reads_vs_PC3.png")

# log? (same)
# pcas$log.reads <- log(pcas$postclean.reads)
# lmFit <- lm(PC1 ~ log.reads, pcas)
# ggplotRegression(lmFit)
# pcas.mid <- pcas[pcas$postclean.reads < 45000000 & pcas$postclean.reads > 15000000,]
# lmFit <- lm(PC1 ~ log.reads, pcas.mid)
# ggplotRegression(lmFit)
# lmFit <- lm(PC2 ~ log.reads, pcas.mid)
# ggplotRegression(lmFit)

# ========================================
# GC !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# >> shows strong effect on PCs, and association with NR reads
plotPCAs(pcas,nrPCs=3,responseVar = "postclean.gc",doCentroids=F)
lmFit <- lm(PC1 ~ postclean.gc, pcas)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename="dag3_uncorr_QC_gc_vs_PC1.png")
lmFit <- lm(PC2 ~ postclean.gc, pcas)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename="dag3_uncorr_QC_gc_vs_PC2.png")

lmFit <- lm(PC3 ~ postclean.gc, pcas)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename="dag3_uncorr_QC_gc_vs_PC3.png")

# lmFit <- lm(postclean.gc ~ postclean.reads, pcas)
# ggplotRegression(lmFit)

# lmFit <- lm(postclean.gc ~ log.reads, pcas2.mid)
# ggplotRegression(lmFit)

# check other QC stuff
# qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned.csv',sep=',',header=T)
# qc2$id.external <- NULL
# qc2$sample.nr <- NULL
# pcas2 <- merge(pcas,qc2,by="sample")
# pcas2.mid <- pcas2[pcas2$postclean.reads < 45000000 & pcas2$postclean.reads > 15000000,]

pcas2 <- pcas
pcas2.mid <- pcas2[pcas2$postclean.reads < 45000000 & pcas2$postclean.reads > 15000000,]

# GENDER has some effect
plotPCAs(pcas2,nrPCs=3,responseVar = "gender",doCentroids=T)
lmFit <- lm(PC1 ~ gender, pcas2)
ggplotRegression(lmFit)
ggplot(pcas2,aes(y=PC1,x=as.factor(gender),col=gender)) + geom_boxplot()
wilcox.test(pcas2$PC1[pcas2$gender=="M"],pcas2$PC1[pcas2$gender=="F"])
# AGE has effect
lmFit <- lm(PC1 ~ age, pcas2)
ggplotRegression(lmFit)
lmFit <- lm(PC2 ~ age, pcas2)
ggplotRegression(lmFit)

# isolation date DOES NOT have effect
pcas2$isolation.date.N <- as.numeric(pcas2$isolation.date)
lmFit <- lm(PC2 ~ isolation.date, pcas2)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename = "dag3_uncorr_QC_isolation_vs_PC1.png")

# LOG (DNA concentration has effect on PC 2,3,4,5 (but not really on 1) )
# =================================
pcas2$log.conc.ng.ul <- log(pcas2$conc.ng.ul)
lmFit <- lm(PC1 ~ log.conc.ng.ul, pcas2)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename = "dag3_uncorr_DNAconc_vs_PC1.png")
lmFit <- lm(PC2 ~ log.conc.ng.ul, pcas2)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename = "dag3_uncorr_DNAconc_vs_PC2.png")
lmFit <- lm(PC3 ~ log.conc.ng.ul, pcas2)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename = "dag3_uncorr_DNAconc_vs_PC3.png")

# is conc correlated with reads? NOPE!
lmFit <- lm(postclean.reads ~ log.conc.ng.ul, pcas2)
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename="dag3_uncorr_DNAconc_vs_ReadNR.png")

# volume for isolation ?? doesn't seem relevant
# actually has a tiny effect
plotPCAs(pcas2,nrPCs=3,responseVar = "vol.ul",doCentroids=T)
g <- ggplot(pcas2,aes(y=PC1,x=as.factor(vol.ul),col=vol.ul)) + geom_boxplot()
wc <- wilcox.test(pcas2$PC1[pcas2$vol.ul==100],pcas2$PC1[pcas2$vol.ul==200])
g <- g + ggtitle(paste0("100 vs 200 Wilcox p-value = ",wc$p.value))
#lmFit <- lm(PC2 ~ vol.ul, pcas2)
#ggplotRegression(lmFit)
ggsave(plot=g,filename="dag3_uncorr_DNAisolation_vs_PC1.png")

# nanodrop
#  > has some effect on PC1, PC 2
# pcas2$log.nano.260.280 <- log(pcas2$nano.260.280)
# pcas2.mid <- pcas2[pcas2$nano.260.280 < 2.5 & pcas2$nano.260.280 > 1,]
# lmFit <- lm(PC1 ~ log.nano.260.280, pcas2.mid)
# ggplotRegression(lmFit)
# lmFit <- lm(PC2 ~ log.nano.260.280, pcas2.mid)
# ggplotRegression(lmFit)

# > vs quality (nothing there)
# lmFit <- lm(conc.ng.ul ~ nano.260.280, pcas2)
# ggplotRegression(lmFit)
# > vs NR reads (nothing there)
# lmFit <- lm(log.reads ~ nano.260.280, pcas2.mid)
# ggplotRegression(lmFit)
# ==============

lmFit <- lm(PC1 ~ URBAN.INDEX.N,pcas )
ggplotRegression(lmFit)
lmFit <- lm(PC2 ~ URBAN.INDEX.N,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC1 ~ lat,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC1 ~ lon,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC1 ~ log.lat,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC2 ~ log.lon,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC2 ~ log.lat,pcas )
ggplotRegression(lmFit)

lmFit <- lm(PC2 ~ LOG.POP.KM2,pcas )
ggplotRegression(lmFit)

# =============================================================================
# DO TSNE
# =============================================================================
# if (file.exists('data_processed/DAG3_uf_tSNEtsv')) {
#   tSNE <- read.table('data_processed/DAG3_uf_tSNEtsv',sep='\t',header=T)
# } else {
#   taxRdyForTSNE <- taxRdyFann[,grep('s__',colnames(taxRdyFann))]
#   bcurtis <- vegdist(taxRdyForTSNE,method = "bray")
#   tSNE <- Rtsne(as.matrix(taxRdyForTSNE), check_duplicates=FALSE, pca=TRUE, 
#                 perplexity=30, theta=0.5, dims=2, pca_center = T,pca_scale = F)
#   tSNE.c = as.data.frame(tSNE$Y) 
#   colnames(tSNE.c) <- c("Coord1","Coord2")
#   tSNE.c$sample <- taxRdyFann$DAG3_sampleID
#   qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
#   qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
#   tSNE.c <- merge(tSNE.c,qc1,by="sample")
#   tSNE.c$URBAN.INDEX <- taxRdyFann$URBAN.INDEX
#   tSNE.c$URBAN.INDEX.N <- as.numeric(as.character(taxRdyFann$URBAN.INDEX))
#   tSNE.c$PROVINCE <- taxRdyFann$PROVINCE
#   tSNE.c$municipaly <- taxRdyFann$municipaly
#   tSNE.c$lon <- taxRdyFann$lon
#   tSNE.c$lat <- taxRdyFann$lat
#   tSNE.c$log.lon <- log(tSNE.c$lon)
#   tSNE.c$log.lat <- log(tSNE.c$lat)
#   tSNE.c$LOG.POP.KM2 <- taxRdyFann$LOG.POP.KM2
#   tSNE.c$LAND.OF.BIRTH <- taxRdyFann$GEBLND
#   tSNE.c$PLACE.OF.BIRTH <- taxRdyFann$GEBPLTS
#   qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned.csv',sep=',',header=T)
#   qc2$id.external <- NULL
#   qc2$sample.nr <- NULL
#   tSNE.c <- merge(tSNE.c,qc2,by="sample")
#   write.table(tSNE.c,file = 'data_processed/tSNE_uf.tsv',sep='\t',row.names = F)
# }
# ggplot(data=tSNE.c,aes(x=Coord1,y=Coord2,col=age)) + geom_point()# + theme(legend.position = "none") 

# ===========================================================================
# ======================= CALCULATE SHANNON ================================
# ===========================================================================
shannonRdy <- taxRdyFann[,grep('s__',colnames(taxRdyFann))]
shannonRdy$div.shannon <- diversity(shannonRdy,MARGIN = 1, index="shannon")
taxRdyFann$div.shannon <- shannonRdy$div.shannon
# > correlate stuff with shannon
shanAndMeta <- taxRdyFann[,-grep('__',colnames(taxRdyFann))]
shanAndMeta <- taxRdyFann[,-grep('__',colnames(taxRdyFann))]
colnames(shanAndMeta)[2] <- "sample"
qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned.csv',sep=',',header=T)
qc2$id.external <- NULL
qc2$sample.nr <- NULL
# merging issues, test: 
mismatches <- as.character(shanAndMeta$sample[!(shanAndMeta$sample %in% qc2$sample)])
write.table(mismatches,"qc_mismatches.csv",row.names = F)
qc2$sample[!(qc2$sample %in% shanAndMeta$sample)]

shanAndMeta <- merge(shanAndMeta,qc2,by="sample")

qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
shanAndMeta <- merge(shanAndMeta,qc1,by="sample")

shanAndMeta$sample <- NULL
shanAndMeta$isolation.person <- NULL
inDFthreeProv <- shanAndMeta[shanAndMeta$PROVINCE %in% c("Drenthe","Friesland","Groningen"),]
inDFthreeProv$PROVINCE <- as.factor(as.character(inDFthreeProv$PROVINCE))
inDFthreeProv$municipaly <- as.factor(as.character(inDFthreeProv$municipaly))

for (rVar in colnames(inDFthreeProv)) {
  if (rVar != "div.shannon") {
    print(rVar)
    form <- reformulate(response = "div.shannon",termlabels = rVar)
    lmFit <- lm(form,inDFthreeProv )
    print(summary(lmFit))
    print(ggplotRegression(lmFit))
  }
}

lmFit <- lm(div.shannon ~ URBAN.INDEX.N,shanAndMeta )
ggplotRegression(lmFit)

# lat matters
# pop.km2 a lil bit (same for houses.km); urban index does seem to matter (for 2,3,5)
# province matters ??

# province matters
ggplot(shanAndMeta[shanAndMeta$PROVINCE %in% c("Drenthe","Friesland","Groningen"),],
       aes(y=div.shannon,x=PROVINCE,col=PROVINCE)) + geom_boxplot()

t <- testOneFeature(inDFthreeProv,responseVar = "PROVINCE",feature = "div.shannon",
               doSave = F,retPlot = T,xLab = "Shannon diversity",title = "Shannon diversity across three provinces")
ggsave(plot=t[[2]],filename = 'DAG3_shannon_provinces.png')

# latitude matters (but longitude does not)
shanAndMeta$LOG.LAT <- log(shanAndMeta$lat)
lmFit <- lm(div.shannon ~ LOG.LAT,shanAndMeta[shanAndMeta$LOG.LAT > 3.96,])
g <- ggplotRegression(lmFit)
ggsave(plot=g,filename = 'DAG3_shannon_latitude.png')

shanAndMeta$LOG.LON<- log(shanAndMeta$lon)
lmFit <- lm(div.shannon ~ LOG.LON,shanAndMeta[shanAndMeta$LOG.LON > 1.65,])
ggplotRegression(lmFit)

# muni
t <- testOneFeature(inDFthreeProv,responseVar = "municipaly",feature = "div.shannon",
               doSave = F,retPlot = T,xLab = "Shannon diversity",title = "Shannon diversity across municipalies")
t[[2]] + theme(legend.position = "none") + geom_hline(yintercept = mean(inDFthreeProv$div.shannon))
# do vs average
res <- data.frame()
for (muni in unique(inDFthreeProv$municipaly)) {
  print(paste(' > testing',muni))
  wct <- wilcox.test(inDFthreeProv$div.shannon[inDFthreeProv$municipaly==muni],inDFthreeProv$div.shannon)
  res <- rbind.data.frame(res,data.frame(MUNI=muni,
                                         NR=sum(inDFthreeProv$municipaly==muni),
                                         AVG=mean(inDFthreeProv$div.shannon[inDFthreeProv$municipaly==muni]),
                                         SD=sd(inDFthreeProv$div.shannon[inDFthreeProv$municipaly==muni]),
                                         pV=wct$p.value))
}
res <- res[order(res$pV,decreasing = F),]
res$FDR <- p.adjust(res$pV)
meanDiv <- mean(inDFthreeProv$div.shannon)
toPlot <- inDFthreeProv[inDFthreeProv$municipaly %in% res$MUNI[1:10],]
g <- ggplot(toPlot,aes(x=municipaly,y=div.shannon,col=municipaly)) + geom_boxplot() + geom_hline(yintercept = meanDiv)
ggsave(plot=g,filename = "DAG3_shannon_muni.png")

# gender has no effect
t <- testOneFeature(inDFthreeProv,responseVar = "gender",feature = "div.shannon",
                    doSave = F,retPlot = T,xLab = "Shannon diversity",title = "Shannon diversity vs gender")
# age
lmfit <- lm(div.shannon ~ age,shanAndMeta[shanAndMeta$age > 20,])
ggplotRegression(lmfit)
lmfit <- lm(div.shannon ~ age,shanAndMeta[shanAndMeta$age < 20,])
ggplotRegression(lmfit)
# DNA conc
shanAndMeta$log.conc.ng.ul <- log(shanAndMeta$conc.ng.ul)
lmfit <- lm(div.shannon ~ log.conc.ng.ul,shanAndMeta[shanAndMeta$log.conc.ng.ul > 1 & shanAndMeta$log.conc.ng.ul < 6,] )
g <- ggplotRegression(lmfit)
ggsave(plot=g,filename = "DAG3_shannon_DNAconc.png")

# DNA extraction volume matters !
#lmfit <- lm(div.shannon ~ vol.ul,shanAndMeta[shanAndMeta$conc.ng.ul < 200,])
#ggplotRegression(lmfit)
shanAndMeta$vol.ul <- as.factor(shanAndMeta$vol.ul)
t <- testOneFeature(shanAndMeta,responseVar = "vol.ul",feature = "div.shannon",
                    doSave = F,retPlot = T,yLab = "Shannon diversity",title = "DNA extraction volume")
ggsave(plot=t[[2]],filename = "DAG3_shannon_DNAextractionVolume.png")

# read number DOES MATTER
lmfit <- lm(div.shannon ~ postclean.reads,shanAndMeta[shanAndMeta$postclean.reads < 40000000 & shanAndMeta$postclean.reads > 15000000,])
g <- ggplotRegression(lmfit)
ggsave(plot=g,filename = "DAG3_shannon_readNR.png")

lmfit <- lm(div.shannon ~ postclean.gc,shanAndMeta)
ggplotRegression(lmfit)

# =================================================================
# do adonis
# =================================================================
# ===================================================================
adonisRdy <- taxRdyFann
colnames(adonisRdy)[2] <- "sample"
qc2 <- read.table('data/LLDAG3_Feces_DNA_concentrations_cleaned.csv',sep=',',header=T)
qc2$id.external <- NULL
qc2$sample.nr <- NULL
adonisRdy <- merge(adonisRdy,qc2,by="sample")
qc1 <- read.table('data/DAG3_QC_summary.csv',sep=',',header=T)
qc1$postclean.reads <- as.numeric(gsub(',','',as.character(qc1$postclean.reads)))
adonisRdy <- merge(adonisRdy,qc1,by="sample")
adonisRdy$sample <- NULL
adonisRdy$status <- NULL
adonisRdy$isolation.person <- NULL
adonisRdy$isolation.date <- NULL
adonisRdy$POSTCODE <- NULL

adonisSpec <- adonisRdy[,grep('s__',colnames(adonisRdy))]
adonisVars <- adonisRdy[,-grep('__',colnames(adonisRdy))]
bcurtis <- vegdist(adonisSpec,method = "bray") 

#adonisVarsTouse <- adonisVars[,c(4,5,8,10,11,12,13,15,16,17,20,21)]
adonisVarsTouse <- adonisVars[,c(4,5,8,10,11,12,13,14,15,16,19,21)]#,10,11,12,13,15,16,17,20,21)]

#adonis <- adonis(bcurtis ~ .,data=adonisVarsTouse,permutations=100,parallel=F)

# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  ad<-adonis(bcurtis ~ adonisVars[,i],permutations=2500)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
adonis_meta <- adonis_meta[order(adonis_meta$`FDR(BH)`),]
rownames(adonis_meta) = colnames(adonisVarsTouse)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
ggsave(g,"uncorr_adonis_results.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.table(adonis_meta,"adonis_meta.csv",sep=",",row.names=F)

# do linear correction for confounders
taxRdyFannCorr <- linearCorrectMGPwy(adonisRdy,corrNames=c("vol.ul","postclean.reads","conc.ng.ul"))

# repeat adonis for good measure
# =========================================
taxRdyFannCorr$nano.260.230 <- NULL
adonisSpec <- taxRdyFannCorr[,grep('s__',colnames(taxRdyFannCorr))]
adonisSpec[adonisSpec < 0] <- 0.0
taxRdyFannCorr[,grep('s__',colnames(adonisRdy))] <- adonisSpec

taxRdyFannCorr3P <- taxRdyFannCorr[taxRdyFannCorr$PROVINCE %in% c("Drenthe","Friesland","Groningen"),]
adonisSpec <- taxRdyFannCorr3P[,grep('s__',colnames(taxRdyFannCorr3P))]
adonisVars <- taxRdyFannCorr3P[,-grep('__',colnames(taxRdyFannCorr3P))]
bcurtis <- vegdist(adonisSpec,method = "bray") 

adonisVarsTouse <- adonisVars[,c(1,4,5,8,10,11,12,13,14,15)]
adonisVarsTouse$municipaly <- as.factor(as.character(adonisVarsTouse$municipaly))
adonisVarsTouse$URBAN.INDEX.N <- as.numeric(as.character(adonisVarsTouse$URBAN.INDEX))

# do univariate adonis
adonis_meta <- matrix(ncol = 7, nrow=ncol(adonisVarsTouse))
for (i in 1:ncol(adonisVarsTouse)){
  print (paste(' >> grinding ',colnames(adonisVarsTouse)[[i]]))
  ad<-adonis(bcurtis ~ adonisVars[,i],permutations=2500)
  aov_table <- ad$aov.tab
  #Df
  adonis_meta[i,1]=aov_table[1,1]
  #SumsOfSqs
  adonis_meta[i,2]=aov_table[1,2]
  #MeanSqs
  adonis_meta[i,3]=aov_table[1,3]
  #FModel
  adonis_meta[i,4]=aov_table[1,4]
  #R2
  adonis_meta[i,5]=aov_table[1,5]
  #Pval
  adonis_meta[i,6]=aov_table[1,6]
}
adonis_meta= as.data.frame(adonis_meta)
adonis_meta <- adonis_meta[order(adonis_meta$`FDR(BH)`),]
rownames(adonis_meta) = colnames(adonisVarsTouse)
adonis_meta$V7=p.adjust(adonis_meta$V6, method = "BH")
adonis_meta$V8="No"
adonis_meta$V8[adonis_meta$V7<0.05]="Yes"
# View(adonis_meta)
colnames(adonis_meta)=c("DF", "SumsOfSqs", "MeanSqs", "FModel","R2","pval","FDR(BH)","Significant")
g <- ggplot(adonis_meta, aes(reorder(row.names(adonis_meta), R2), R2, fill=Significant)) + 
  geom_bar(stat = "identity") + coord_flip() + theme_bw() + 
  ylab ("Explained variance (R^2) ") + xlab ("Factor")  + theme(text = element_text(size=11))
print(g)
ggsave(plot=g,filename = "corr_adonis_meta.png")
print(adonis_meta[order(adonis_meta$`FDR(BH)`),])
write.table(adonis_meta,"corr_adonis_meta.csv",sep=",",row.names=F)
write.table(taxRdyFannCorr,"corr_DAG3.tsv",sep="\t",row.names=F)

# ======================================================================================
#      ================= DO STATISTICS ON CORRECTED METAPHLAN =====================
# ======================================================================================
