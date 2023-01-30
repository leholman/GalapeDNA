#############################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====29.02.2020====####
#############################################


### Script 3 - Distance Data Generation ###


#calculate some distances 'as-the-fish-swims' 
####====0.0 Packages====####
library('gdistance') 
library('sf')
library('geosphere')

####====1.0 As-the-fish-swims Distance & Points  ====####


#Pull in galapagos
gebco.crop <- readRDS("mapBuilding/GalapBathy.rds")
gebco.crop2 <- gebco.crop

#set land to zero, sea to 1
gebco.crop2@data@values[gebco.crop2@data@values>-1] <-0
gebco.crop2@data@values[gebco.crop2@data@values<0] <-Inf

#Places
metadatSites <- read.csv("metadata.site.out.csv",row.names=1)
places <- cbind(metadatSites$lon2,metadatSites$lat2)
sitePairwiseDist <- expand.grid(Start=metadatSites$SiteID,
                                End=metadatSites$SiteID)

#
tr1 <- transition(gebco.crop2, transitionFunction=mean, directions=16)
tr2 <- geoCorrection(tr1, type="r", multpl=TRUE)

sitePairwiseDist$Flat <- metadatSites$lat2[match(sitePairwiseDist$Start,metadatSites$SiteID)]
sitePairwiseDist$Flon <- metadatSites$lon2[match(sitePairwiseDist$Start,metadatSites$SiteID)]
sitePairwiseDist$Tlat <- metadatSites$lat2[match(sitePairwiseDist$End,metadatSites$SiteID)]
sitePairwiseDist$Tlon <- metadatSites$lon2[match(sitePairwiseDist$End,metadatSites$SiteID)]
sitePairwiseDist$Calcdistance <- rep(NA,length(sitePairwiseDist$End))



#empty points table 

pathPointsTable <- c()


#Loop over site comparisons and output needed data

for (row in 1:length(sitePairwiseDist$End)){
  
  
  #set dist to zero and skip loop for comparing a site to itself
  if(sitePairwiseDist$Start[row]==sitePairwiseDist$End[row]){sitePairwiseDist$Calcdistance[row] <- 0
  next()}
  #get the shortest path 
  loopPath <- shortestPath(tr2, c(sitePairwiseDist$Flon[row],
                                  sitePairwiseDist$Flat[row]),
                           c(sitePairwiseDist$Tlon[row],
                             sitePairwiseDist$Tlat[row]), output="SpatialLines")
  #calculate the distance of the shortest path
  looplen <- lengthLine(loopPath)
  #sample points along the path every 1000 metres
  loopPathpoints <- spsample(loopPath,looplen/1000,type="regular")
  
  #output lengths to dataframe
  sitePairwiseDist$Calcdistance[row] <- looplen
  
  #create point output table 
  loopPathPointsize <- length(loopPathpoints)
  loopPointsTable <- data.frame("Order"=1:loopPathPointsize,
                                "Start"=rep(sitePairwiseDist$Start[row],loopPathPointsize),
                                "End"=rep(sitePairwiseDist$End[row],loopPathPointsize),
                                "lon"=loopPathpoints@coords[,1],
                                "lat"=loopPathpoints@coords[,2])
  
  
  pathPointsTable <- rbind(pathPointsTable,loopPointsTable)
  
  print(row)
}



write.csv(sitePairwiseDist,"SiteDistance.csv")
write.csv(pathPointsTable,"pathPointsTable.csv")

## Now these lat long points are used to pull the model data in 


####====2.0 Oceanographic Distance  ====####

modeldat <-read.csv("pathPoints/pathPointsTable_Sep.csv")
modeldatLAND <- modeldat[modeldat$VVEL==0,]


points(modeldatLAND$lon,modeldatLAND$lat,pch=16,col="red",cex=0.4)

####====2.1 Calculate sum of vectors ====####

VectorSum <- function(Northing,Easting){
  if (Northing==0|Easting==0){stop("Northing or Easting is zero - land")} 
  outputMagnitude <- sqrt(Northing^2 + Easting^2)
  print(outputMagnitude)
  
  θ <- tan(Northing/Easting)*(180/pi)

  
  
  
  if (Northing>0 & Easting>0){
    print(90-θ)
    print("N+ E+")} else 
      if (Northing<0 & Easting>0){
        print(90+θ)
        print("N- E+")} else 
          if (Northing<0 & Easting<0){
            print(270-θ)
            print("N- E-")} else 
              if (Northing>0 & Easting<0){
                print(270+θ)
                print("N+ E-")} 
  
  
  print("Working")
}


####====2.2 Calculate angle between two lat lon points ====####
##Maybe this https://www.omnicalculator.com/other/azimuth

####====2.3 Calculate resistance vector  ====####







