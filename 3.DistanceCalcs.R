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


plot(modeldatLAND$lon,modeldatLAND$lat,pch=16,col="red",cex=0.4)

####====2.1 Calculate sum of vectors ====####


#First the angle
vectorAngle <- function(Northing,Easting){
  if(Northing == 0 & Easting == 0){return(0)}
  Northing=Northing+0.0000001
  Easting=Easting+0.0000001
resultant_angle <- atan2(sqrt(Northing^2),sqrt(Easting^2))
resultant_angle_degrees <- (180/pi) * resultant_angle
if(Northing>0 & Easting>0){return(90-resultant_angle_degrees)} else 
  if(Northing<0 & Easting>0){return(90+resultant_angle_degrees)} else
    if(Northing<0 & Easting<0){return(270-resultant_angle_degrees)} else
      if(Northing>0 & Easting<0){return(270+resultant_angle_degrees)}
}

#resultant_angle_degrees <- ifelse(resultant_angle_degrees < 0, 360 + resultant_angle_degrees, resultant_angle_degrees)
#return(resultant_angle_degrees)}

vectorAngle(5,5)
vectorAngle(-5,5)
vectorAngle(-5,-5)
vectorAngle(5,-5)


modeldat$resultantAngle <- unlist(mapply(vectorAngle,modeldat$UVEL,modeldat$VVEL))
hist(unlist(Angles),breaks=100)

#Now the magnitude

vectorSum <- function(Northing,Easting){
  if(Northing == 0 & Easting == 0){return(NA)}
  resultantMagnitude <- sqrt(Northing^2+Easting^2)
  return(resultantMagnitude)}

modeldat$magnitudes <- unlist(mapply(vectorSum,modeldat$UVEL,modeldat$VVEL))



####====2.2 Calculate angle between two lat lon points ====####

azimuth <- function(lat1, lon1, lat2, lon2) {
  
  # Convert decimal degrees to radians
  lat1 <- lat1 * pi / 180
  lon1 <- lon1 * pi / 180
  lat2 <- lat2 * pi / 180
  lon2 <- lon2 * pi / 180
  
  # Calculate the difference in longitude
  delta_lon <- lon2 - lon1
  
  # Calculate the azimuth (bearing) using the Haversine formula
  y <- sin(delta_lon) * cos(lat2)
  x <- cos(lat1) * sin(lat2) - sin(lat1) * cos(lat2) * cos(delta_lon)
  azimuth_rad <- atan2(y, x)
  
  # Convert radians to degrees
  azimuth_deg <- azimuth_rad * 180 / pi
  
  # Make sure the result is in the range [0, 360)
  if (azimuth_deg < 0) {
    azimuth_deg <- azimuth_deg + 360
  }
  
  return(azimuth_deg)
}

azimuth(55.693528318045544, 12.61700639326756,55.690141862839695, 12.571430301020584)



####====2.3 Calculate resistance vector  ====####

exampleAngles <- seq(110,200,10)

resistanceAngle <- modeldat$resultantAngle[1:10]
resistanceMagnitude <- modeldat$magnitudes[1:10]

((exampleAngles/360)-(resistanceAngle/360))*360


#1 calculate difference in angle such that direction of travel is 0 

angleCalc <- function(A,B){
  input <- B - A
  if (input > 0){
    return(input)} else 
      if (input < 0) {
        return(input+360)} else
          if (input == 0){return(0)}
}

#2 use cosin to turn this value such that 1 = same direction, -1 = opposite

input <- 1:360

output <- cos((input/360)*(2*pi))

#3 create a resultant scaler based on this value, maybe direction * mangnitude?



#4 sum all scalers to give oceanographic resistance?


