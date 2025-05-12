################################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====01.05.2025====####
#############################################

### Script $ - Calculating Oceanographic Resistance
### Here we for each journey we process the ocean model extracted data to generate oceanogrphic resistance


## functions for calculations

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



##Function 5 scale an azimuth angle into resistance with -1 being 0 degrees and 1 being 180 degrees 

cosTrans <- function(input){
  return(-cos((input/360)*(2*pi)))
}




## now lets input the data


# Read in the new datasets 
pathPoints2 <- read.csv("OceanModelDataExtraction/pathPointsTable.out.csv")
pathPoints2[is.na(pathPoints2)] <- 0


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
  
  dat <- pathPoints2
  
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
    resultKey <- paste("pathPoints2" , uvName, sep = "_")
    oceanResistResults[[resultKey]] <- journeyOutput
    
    # Optional: write the output to CSV files (uncomment if desired)
    # outFile <- paste0("distanceData/OceanogrphicResistance_", dataName, "_", uvName, ".csv")
    # write.csv(journeyOutput, outFile, row.names = FALSE)
    
  } # end uv dataset loop


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


write.csv(combined_df,"distanceData/OceanographicResistance.csv")

