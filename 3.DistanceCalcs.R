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
library("maditr")

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


sitePairwiseDist.d <- dcast(sitePairwiseDist,Start~End,value.var = "Calcdistance")
sitePairwiseDist.d.out <- as.matrix(sitePairwiseDist.d[,2:24])
rownames(sitePairwiseDist.d.out) <- sitePairwiseDist.d$Start

write.csv(sitePairwiseDist.d.out,"distanceData/SiteDistanceMatrix.csv")
write.csv(sitePairwiseDist,"distanceData/SiteDistance.csv")
write.csv(pathPointsTable,"distanceData/pathPointsTable.csv")

## Now these lat long points are used to pull the model data in 


####====2.0 Oceanographic Distance  ====####

modeldat <-read.csv("distanceData/pathPoints/pathPointsTable_Oct.csv")
modeldatLAND <- modeldat[modeldat$VVEL==0,]


plot(modeldatLAND$lon,modeldatLAND$lat,pch=16,col="red",cex=0.4)




####====2.1 Functions for calculating distance ====####


##Function 1 - calculate the resultant angle from a Northing and an Easting 

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

#testing
vectorAngle(5,5)
vectorAngle(-5,5)
vectorAngle(-5,-5)
vectorAngle(5,-5)


##Function 2 Calculate a resultant magnitude of the vector from a Northign and an Easting

vectorSum <- function(Northing,Easting){
  if(Northing == 0 & Easting == 0){return(NA)}
  resultantMagnitude <- sqrt(Northing^2+Easting^2)
  return(resultantMagnitude)}

modeldat$magnitudes <- unlist(mapply(vectorSum,modeldat$UVEL,modeldat$VVEL))


##Function 3 Calculate a the angle (azimuth) between two geographic points with lat lon 


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

#testing
azimuth(55.68517608483797, 12.57629649327822,55.66184699523609, 12.57955805920173)


##Function 4 Calculate the azimuth difference from angle A to angle B

angleCalc <- function(A,B){
  input <- B - A
  if (input > 0){
    return(input)} else 
      if (input < 0) {
        return(input+360)} else
          if (input == 0){return(0)}
}



##Function 5 scale an azimuth angle into resistance with 1 being 0 degrees and -1 being 180 degrees 

cosTrans <- function(input){
  return(cos((input/360)*(2*pi)))
}



####====2.2 Calculate distance on real data ====####


## 1 Calculate angle and magnitude of each lat lon point along the path
#angle
modeldat$resultantAngle <- unlist(mapply(vectorAngle,modeldat$UVEL,modeldat$VVEL))
#magnitude 
modeldat$magnitudes <- unlist(mapply(vectorSum,modeldat$UVEL,modeldat$VVEL))

## 2 Create an output dataframe for each journey

modeldat$journeyID <- paste(modeldat$Start,modeldat$End,sep="_")

journeyOutput <- data.frame("journeyID"=unique(paste(modeldat$Start,modeldat$End,sep="_")),"OceanographicResistance"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))),"OceanographicResistanceSD"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))))

## 3 loop over each journey 

for (journeyIndex in 1:length(journeyOutput$journeyID)){
  
 # journeyIndex <- 8
  
  loopJourney <- journeyOutput$journeyID[journeyIndex]
  
  loopData <- modeldat[modeldat$journeyID==loopJourney,]
  
  loopData$comparisonAngle <- NA
  loopData$comparisonAngleDiff <- NA
  loopData$comparisonAngleMetric <- NA
  loopData$comparisonMetric <- NA
  
  for (loopLocation in 1:length(loopData$Order)){
    if(loopData$Order[loopLocation]==max(loopData$Order)){next()}
    index <- match(loopData$Order[loopLocation],loopData$Order)
    index2 <- match(loopData$Order[loopLocation+1],loopData$Order)
  
  loopData$comparisonAngle[index] <- azimuth(loopData$lat[index],loopData$lon[index],loopData$lat[index2],loopData$lon[index2])

  loopData$comparisonAngleDiff[index] <- angleCalc(loopData$comparisonAngle[index],loopData$resultantAngle[index])
  
  loopData$comparisonAngleMetric[index] <- cosTrans(loopData$comparisonAngleDiff[index])
  
  loopData$comparisonMetric[index] <- loopData$comparisonAngleMetric[index] * loopData$magnitudes[index]
  
  }
  journeyOutput$OceanographicResistance[journeyIndex] <-mean(loopData$comparisonMetric,na.rm = T)
  journeyOutput$OceanographicResistanceSD[journeyIndex] <-sd(loopData$comparisonMetric,na.rm = T)
}


#Now lets transform and output the data 

journeyOutput$start <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 1)

journeyOutput$end <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 2)

JourneyMatrix <- dcast(journeyOutput,start~end,value.var = "OceanographicResistance")
JourneyMatrix2 <- as.matrix(JourneyMatrix[,2:24])
rownames(JourneyMatrix2) <- JourneyMatrix$start


write.csv(JourneyMatrix2,"distanceData/OceanogrphicResistanceMatrix.csv")
write.csv(journeyOutput,"distanceData/OceanogrphicResistancePairwise.csv")

####====2.3 Generalise across many months of data ====####

for (month in list.files(pattern=".*.csv","distanceData/pathPoints")){
  #month <- "pathPointsTable_Apr.csv"
  monthText <- gsub(".csv","",gsub("pathPointsTable_","",month))
  print(monthText)
  
  modeldat <-read.csv(paste0("distanceData/pathPoints/pathPointsTable_",monthText,".csv"))
  
  ## 1 Calculate angle and magnitude of each lat lon point along the path
  #angle
  modeldat$resultantAngle <- unlist(mapply(vectorAngle,modeldat$UVEL,modeldat$VVEL))
  #magnitude 
  modeldat$magnitudes <- unlist(mapply(vectorSum,modeldat$UVEL,modeldat$VVEL))
  
  ## 2 Create an output dataframe for each journey
  
  modeldat$journeyID <- paste(modeldat$Start,modeldat$End,sep="_")
  
  journeyOutput <- data.frame("journeyID"=unique(paste(modeldat$Start,modeldat$End,sep="_")),"OceanographicResistance"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))),"OceanographicResistanceSD"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))))
  
  ## 3 loop over each journey 
  
  for (journeyIndex in 1:length(journeyOutput$journeyID)){
    
    # journeyIndex <- 8
    
    loopJourney <- journeyOutput$journeyID[journeyIndex]
    
    loopData <- modeldat[modeldat$journeyID==loopJourney,]
    
    loopData$comparisonAngle <- NA
    loopData$comparisonAngleDiff <- NA
    loopData$comparisonAngleMetric <- NA
    loopData$comparisonMetric <- NA
    
    for (loopLocation in 1:length(loopData$Order)){
      if(loopData$Order[loopLocation]==max(loopData$Order)){next()}
      index <- match(loopData$Order[loopLocation],loopData$Order)
      index2 <- match(loopData$Order[loopLocation+1],loopData$Order)
      
      loopData$comparisonAngle[index] <- azimuth(loopData$lat[index],loopData$lon[index],loopData$lat[index2],loopData$lon[index2])
      
      loopData$comparisonAngleDiff[index] <- angleCalc(loopData$comparisonAngle[index],loopData$resultantAngle[index])
      
      loopData$comparisonAngleMetric[index] <- cosTrans(loopData$comparisonAngleDiff[index])
      
      loopData$comparisonMetric[index] <- loopData$comparisonAngleMetric[index] * loopData$magnitudes[index]
      
    }
    journeyOutput$OceanographicResistance[journeyIndex] <-mean(loopData$comparisonMetric,na.rm = T)
    journeyOutput$OceanographicResistanceSD[journeyIndex] <-sd(loopData$comparisonMetric,na.rm = T)
  }
  
  
  #Now lets transform and output the data 
  
  journeyOutput$start <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 1)
  
  journeyOutput$end <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 2)
  
  JourneyMatrix <- dcast(journeyOutput,start~end,value.var = "OceanographicResistance")
  JourneyMatrix2 <- as.matrix(JourneyMatrix[,2:24])
  rownames(JourneyMatrix2) <- JourneyMatrix$start
  
  
  write.csv(JourneyMatrix2,paste0("distanceData/monthlyResistance/OceanogrphicResistanceMatrix_",monthText,".csv"))
  write.csv(journeyOutput,paste0("distanceData/monthlyResistance/OceanogrphicResistancePairwise_",monthText,".csv"))
  
}


## A little plot to visualise 

par(mfrow=c(4,3))

ResistanceMonths <- as.list(c())

for (month in list.files(pattern="OceanogrphicResistanceMatrix.*","distanceData/monthlyResistance/")){


oceanResistance <- as.dist(as.matrix(read.csv(paste0("distanceData/monthlyResistance/",month),row.names = 1)))
#nMDS.resist <- metaMDS(oceanResistance,trymax=500)
#plot(nMDS.resist,type = "t",main="Resistance")

ResistanceMonths <- append(ResistanceMonths,list(as.matrix(oceanResistance)))

}

####====2.4 What happens in an ElNino Year??? ====####

for (month in list.files(pattern=".*.csv","distanceData/pathPointsElNino/")){
  #month <- "pathPointsTable_Apr.csv"
  monthText <- gsub(".csv","",gsub("pathPointsTable_","",month))
  print(monthText)
  
  modeldat <-read.csv(paste0("distanceData/pathPointsElNino//pathPointsTable_",monthText,".csv"))
  
  ## 1 Calculate angle and magnitude of each lat lon point along the path
  #angle
  modeldat$resultantAngle <- unlist(mapply(vectorAngle,modeldat$UVEL,modeldat$VVEL))
  #magnitude 
  modeldat$magnitudes <- unlist(mapply(vectorSum,modeldat$UVEL,modeldat$VVEL))
  
  ## 2 Create an output dataframe for each journey
  
  modeldat$journeyID <- paste(modeldat$Start,modeldat$End,sep="_")
  
  journeyOutput <- data.frame("journeyID"=unique(paste(modeldat$Start,modeldat$End,sep="_")),"OceanographicResistance"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))),"OceanographicResistanceSD"=rep(0,length(unique(paste(modeldat$Start,modeldat$End,sep="_")))))
  
  ## 3 loop over each journey 
  
  for (journeyIndex in 1:length(journeyOutput$journeyID)){
    
    # journeyIndex <- 8
    
    loopJourney <- journeyOutput$journeyID[journeyIndex]
    
    loopData <- modeldat[modeldat$journeyID==loopJourney,]
    
    loopData$comparisonAngle <- NA
    loopData$comparisonAngleDiff <- NA
    loopData$comparisonAngleMetric <- NA
    loopData$comparisonMetric <- NA
    
    for (loopLocation in 1:length(loopData$Order)){
      if(loopData$Order[loopLocation]==max(loopData$Order)){next()}
      index <- match(loopData$Order[loopLocation],loopData$Order)
      index2 <- match(loopData$Order[loopLocation+1],loopData$Order)
      
      loopData$comparisonAngle[index] <- azimuth(loopData$lat[index],loopData$lon[index],loopData$lat[index2],loopData$lon[index2])
      
      loopData$comparisonAngleDiff[index] <- angleCalc(loopData$comparisonAngle[index],loopData$resultantAngle[index])
      
      loopData$comparisonAngleMetric[index] <- cosTrans(loopData$comparisonAngleDiff[index])
      
      loopData$comparisonMetric[index] <- loopData$comparisonAngleMetric[index] * loopData$magnitudes[index]
      
    }
    journeyOutput$OceanographicResistance[journeyIndex] <-mean(loopData$comparisonMetric,na.rm = T)
    journeyOutput$OceanographicResistanceSD[journeyIndex] <-sd(loopData$comparisonMetric,na.rm = T)
  }
  
  
  #Now lets transform and output the data 
  
  journeyOutput$start <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 1)
  
  journeyOutput$end <- sapply(strsplit(journeyOutput$journeyID,split = "_"),'[', 2)
  
  JourneyMatrix <- dcast(journeyOutput,start~end,value.var = "OceanographicResistance")
  JourneyMatrix2 <- as.matrix(JourneyMatrix[,2:24])
  rownames(JourneyMatrix2) <- JourneyMatrix$start
  
  
  write.csv(JourneyMatrix2,paste0("distanceData/monthlyResistanceElNino/OceanogrphicResistanceMatrix_",monthText,".csv"))
  write.csv(journeyOutput,paste0("distanceData/monthlyResistanceElNino/OceanogrphicResistancePairwise_",monthText,".csv"))
  
}

# A litlte plot to visualise

par(mfrow=c(4,3))

ResistanceMonths <- as.list(c())

for (month in list.files(pattern="OceanogrphicResistanceMatrix.*","distanceData/monthlyResistance/")){
  
  
  oceanResistance <- as.dist(as.matrix(read.csv(paste0("distanceData/monthlyResistance/",month),row.names = 1)))
  nMDS.resist <- metaMDS(oceanResistance,trymax=500)
  plot(nMDS.resist,type = "t",main="Resistance")
  
  ResistanceMonths <- append(ResistanceMonths,list(as.matrix(oceanResistance)))
  
}



## A little plot to visualise el nino

par(mfrow=c(4,3))

ResistanceMonths <- as.list(c())

for (month in list.files(pattern="OceanogrphicResistanceMatrix.*","distanceData/monthlyResistanceElNino/")){
  
  
  oceanResistance <- as.dist(as.matrix(read.csv(paste0("distanceData/monthlyResistanceElNino/",month),row.names = 1)))
  nMDS.resist <- metaMDS(oceanResistance,trymax=500)
  plot(nMDS.resist,type = "t",main="Resistance")
  
  ResistanceMonths <- append(ResistanceMonths,list(as.matrix(oceanResistance)))
  
}


####====3.0 Temperature Distance ====####

temp <- read.csv("distanceData/Temp_030623.csv",row.names = 1)

sep.temp <- temp$Sep
names(sep.temp) <- temp$SiteID

#temp.dist <- dist(sep.temp,diag = T,upper = T)

temp.dist <- outer(sep.temp,sep.temp,"-")

write.csv(temp.dist,"distanceData/tempDistance.csv")


####====4.0 Code basement ====####



shuffled <- as.matrix(read.csv("distanceData/monthlyResistance/OceanogrphicResistanceMatrix_Jun.csv",row.names = 1))

shuffled <- names(shuffled)[sample(1:23)]

mantel.rtest(as.dist(as.matrix(read.csv("distanceData/monthlyResistance/OceanogrphicResistanceMatrix_Jan.csv",row.names = 1))),
             as.dist(test),nrepet = 9999)



# CORA to SUAR journey 

loopJourney <- "CHAM_PMOR"
checkdat <- modeldat[modeldat$journeyID==loopJourney,]
hist(checkdat$resultantAngle,breaks=180)
hist(checkdat$magnitudes,breaks=180)
hist(loopData$comparisonAngle)


points(loopData$lon,loopData$lat,cex=0.5,lty=16,col="red")


