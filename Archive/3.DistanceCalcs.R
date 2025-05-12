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
library(dplyr)
library(tidyr)
library(corrplot)

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

vectorAngle <- function(Northing, Easting) {
  if (Northing == 0 & Easting == 0) {
    return(0)  # Define a return value (e.g., 0° for no movement)
  }
  
  # Add a tiny offset only when necessary to avoid NaN
  if (Northing == 0) Northing <- .Machine$double.eps
  if (Easting == 0) Easting <- .Machine$double.eps
  
  # Compute angle in radians (flipping arguments to align 0° with North)
  angle_radians <- atan2(Easting, Northing)
  
  # Convert to degrees
  angle_degrees <- (180 / pi) * angle_radians
  
  # Ensure angle is in [0, 360] range
  return((angle_degrees + 360) %% 360)
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




####====4.0 NewData ====####

# Read in the new datasets 
pathPoints2 <- read.csv("distanceData/workingFolderJan25/pathPointsTable.out.csv")
pathPoints2[is.na(pathPoints2)] <- 0
pathPoints2old <- read.csv("distanceData/workingFolderJan25/pathPointsTableOld.out.csv")
pathPoints2old[is.na(pathPoints2old)] <- 0

# List of input datasets
data_list <- list("pathPoints2" = pathPoints2,
                  "pathPoints2old" = pathPoints2old)

# Define the six U/V column pairs
uv_columns <- list(
  "oct_surf" = list(u = "U_oct_surf", v = "V_oct_surf"),
  "oct_20"   = list(u = "U_oct_20",   v = "V_oct_20"),
  "oct_50"   = list(u = "U_oct_50",   v = "V_oct_50"),
  "year_surf"= list(u = "U_year_surf",v = "V_year_surf"),
  "year_20"  = list(u = "U_year_20",  v = "V_year_20"),
  "year_50"  = list(u = "U_year_50",  v = "V_year_50")
)

# Initialize a list to store results for downstream applications
oceanResistResults <- list()

# Loop over each input dataset
for(dataName in names(data_list)){
  
  dat <- data_list[[dataName]]
  
  # Ensure the journeyID exists (combine Start and End)
  dat$journeyID <- paste(dat$Start, dat$End, sep="_")
  
  # Loop over each U/V dataset combination
  for(uvName in names(uv_columns)){
    
    # Get the column names for U and V
    ucol <- uv_columns[[uvName]]$u
    vcol <- uv_columns[[uvName]]$v
    
    # Create a working copy of the data
    modeldat <- dat
    
    # Calculate resultant angle and magnitude using the chosen U/V pair
    modeldat$resultantAngle <- unlist(mapply(vectorAngle, modeldat[[ucol]], modeldat[[vcol]]))
    modeldat$magnitudes     <- unlist(mapply(vectorSum,   modeldat[[ucol]], modeldat[[vcol]]))
    
    # Identify unique journeys and set up an output dataframe
    journey_ids <- unique(modeldat$journeyID)
    journeyOutput <- data.frame(journeyID = journey_ids,
                                OceanographicResistance = rep(0, length(journey_ids)),
                                OceanographicResistanceSD = rep(0, length(journey_ids)),
                                stringsAsFactors = FALSE)
    
    # Loop over each journey to calculate resistance
    for(journeyIndex in 1:length(journey_ids)){
      
      loopJourney <- journey_ids[journeyIndex]
      loopData <- modeldat[modeldat$journeyID == loopJourney, ]
      
      # Initialize additional columns needed for the calculation
      loopData$comparisonAngle <- NA
      loopData$comparisonAngleDiff <- NA
      loopData$comparisonAngleMetric <- NA
      loopData$comparisonMetric <- NA
      
      # Loop over each pair of sequential points along the journey
      # (skip the last point, as it has no next point)
      for(i in 1:(nrow(loopData)-1)){
        # Calculate the azimuth between sequential points
        loopData$comparisonAngle[i] <- azimuth(loopData$lat[i], loopData$lon[i],
                                               loopData$lat[i+1], loopData$lon[i+1])
        # Calculate the difference between the azimuth and the resultant angle
        loopData$comparisonAngleDiff[i] <- angleCalc(loopData$comparisonAngle[i], loopData$resultantAngle[i])
        # Transform the angle difference into a metric (scale from 1 to -1)
        loopData$comparisonAngleMetric[i] <- cosTrans(loopData$comparisonAngleDiff[i])
        # Multiply by the magnitude to get the comparison metric
        loopData$comparisonMetric[i] <- loopData$comparisonAngleMetric[i] * loopData$magnitudes[i]
      }
      
      # Compute the mean and standard deviation of the comparison metric for this journey
      journeyOutput$OceanographicResistance[journeyIndex] <- mean(loopData$comparisonMetric, na.rm = TRUE)
      journeyOutput$OceanographicResistanceSD[journeyIndex] <- sd(loopData$comparisonMetric, na.rm = TRUE)
    } # end journey loop
    
    # Save the journey output dataframe in the results list with a clear key name
    resultKey <- paste(dataName, uvName, sep = "_")
    oceanResistResults[[resultKey]] <- journeyOutput
    
    # Optional: write the output to CSV files (uncomment if desired)
    # outFile <- paste0("distanceData/OceanogrphicResistance_", dataName, "_", uvName, ".csv")
    # write.csv(journeyOutput, outFile, row.names = FALSE)
    
  } # end uv dataset loop
  
} # end input dataset loop

# oceanResistResults now contains 6 dataframes for each input file.
# For example, to access the oceanographic resistance from pathPoints2 using the oct_surf data:
#    oceanResistResults[["pathPoints2_oct_surf"]]
#
# These results can be used for downstream applications.


# Create a list of dataframes with renamed OceanographicResistance columns:
resist_dfs <- lapply(names(oceanResistResults), function(key) {
  df <- oceanResistResults[[key]][, c("journeyID", "OceanographicResistance")]
  colnames(df)[2] <- key  # rename resistance column to the key name
  return(df)
})

# Merge all dataframes by 'journeyID'; using inner join to keep only journeys common to all datasets
combined_df <- Reduce(function(x, y) merge(x, y, by = "journeyID", all = FALSE), resist_dfs)

# Optionally, check the merged data (first few rows)
head(combined_df)

# Remove rows with missing data (if any) to ensure valid correlations
combined_df_complete <- combined_df[complete.cases(combined_df), ]

# Compute the correlation matrix on the OceanographicResistance columns (excluding journeyID)
corr_matrix <- cor(combined_df_complete[ , -1])

# Create a heatmap of the correlation matrix using corrplot
corrplot(corr_matrix, method = "color", 
         addCoef.col = "black",      # add correlation coefficients in black
         number.cex = 0.7,           # coefficient text size
         tl.cex = 0.8,               # label text size
         tl.col = "black",           # label color
         title = "Correlation Across Oceanographic Resistance Datasets",
         mar = c(0,0,1,0))           # adjust margins for title


pairs(combined_df_complete[ , -1], main = "Pairs Plot of Oceanographic Resistance Datasets")


write.csv(combined_df,"distanceData/workingFolderJan25/OceanogrphicResistNew.csv")


####====5.0 Code basement ====####



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


