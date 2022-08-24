#############################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====29.02.2020====####
#############################################


### Script 2 - Maps ###

# Load packages
library(ncdf4)
library(raster)

##Lots of lines hashed out here, original file is 10GB+ so only subset used ;) 
#load GEBCO_2022 netcdf downloaded on 110822 - https://www.gebco.net/data_and_products/gridded_bathymetry_data/
#gebco <- raster("mapBuilding/GEBCO_2022_sub_ice_topo.nc")

# Create extent (our map area)
#galap.ex <- extent(-92, -89, -1.55, 0.68)

# Create a crop of the bathymetric data
#gebco.crop <- crop(gebco, galap.ex)

#Save file 
#saveRDS(gebco.crop,file="mapBuilding/GalapBathy.rds")

gebco.crop <- readRDS("mapBuilding/GalapBathy.rds")

#load lat lon
metadata <- read.csv(file = "metadata.site.out.csv")


#function to calculate break points for colours from https://www.benjaminbell.co.uk/2019/08/bathymetric-maps-in-r-colour-palettes.html
## x = raster, b1 & b2 = number of divisions for each sequence, r1 & r2 = rounding value
colbr <- function(x, b1=50, b2=50, r1=-2, r2=-2) {
  # Min/max values of the raster (x)
  mi <- cellStats(x, stat="min")-100
  ma <- cellStats(x, stat="max")+100
  # Create sequences, but only use unique numbers
  s1 <- unique(round(seq(mi, 0, 0-mi/b1),r1))
  s2 <- unique(round(seq(0, ma, ma/b2),r2))
  # Combine sequence for our break points, removing duplicate 0
  s3 <- c(s1, s2[-1])
  # Create a list with the outputs
  # [[1]] = length of the first sequence minus 1 (water)
  # [[2]] = length of the second sequence minus 1 (land)
  # [[3]] = The break points
  x <- list(length(s1)-1, length(s2)-1, s3)
}

galap.br <- colbr(gebco.crop,b2=1)

# Get country shapefiles
eq <- getData("GADM", country="ECU", level=0)

# Colour palette
blue.col <- colorRampPalette(c("darkblue", "lightblue"))

#First plot

pdf(width = 8,height=6.5,file="mapBuilding/test.pdf")

plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend = FALSE)
#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=1)
points(metadata$lon2,metadata$lat2,pch=16, col='white',cex=0.5)

dev.off()

