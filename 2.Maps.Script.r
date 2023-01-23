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

pdf(width = 8,height=6.5,file="mapBuilding/test1.pdf")

plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend = FALSE)
#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=1)
points(metadata$lon2,metadata$lat2,pch=16, col='white',cex=0.5)

dev.off()

#Adding in colours per island
cols <- c("#80B1D3","#FFFFB3","#FFFFB3","#80B1D3","#FB8072","#BEBADA","#FFED6F","#CCEBC5",
          "#80B1D3","#8DD3C7","#FDB462","#FFFFB3","#B3DE69","gray85","#FFFFB3","#80B1D3",
          "#FCCDE5","#80B1D3","#80B1D3","#CCEBC5","#BC80BD","#8DD3C7","#80B1D3")
pdf(width = 8,height=6.5,file="mapBuilding/test2.pdf")

plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend = FALSE)
#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=3)
points(metadata$lon2,metadata$lat2,pch=16, col=cols,cex=2)

dev.off()


#calculate some distances 'as-the-fish-swims' 
library('gdistance') 
library('sf')

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


#trying here to run the whole command sensibly by not looping but it doesnt currently work
test <- costDistance(tr2,cbind(sitePairwiseDist$Flon,sitePairwiseDist$Flat),cbind(sitePairwiseDist$Tlon,sitePairwiseDist$Tlat))

#this works but is stupid and slow
for (row in 1:length(sitePairwiseDist$End)){
  sitePairwiseDist$Calcdistance[row] <-  costDistance(tr2,c(sitePairwiseDist$Flon[row],
                                                            sitePairwiseDist$Flat[row]),
                                                      c(sitePairwiseDist$Tlon[row],
                                                        sitePairwiseDist$Tlat[row]),
                                                      distfun = gdistance::distHaversine)
  print(row)
}

write.csv(sitePairwiseDist,"SiteDistance.csv")

#Another solution

for (row in 1:25){ 

test <- shortestPath(tr2, c(sitePairwiseDist$Flon[row],
                    sitePairwiseDist$Flat[row]),
             c(sitePairwiseDist$Tlon[row],
               sitePairwiseDist$Tlat[row]), output="SpatialLines")
lines(test,col="blue")

}

row <- 4


test <- shortestPath(tr2, c(sitePairwiseDist$Flon[row],
                            sitePairwiseDist$Flat[row]),
                     c(sitePairwiseDist$Tlon[row],
                       sitePairwiseDist$Tlat[row]), output="SpatialLines")

st_line_sample(test, 0.001, type = "regular")
st_line_sample(st_transform(test, 3857), density = c(1/100))

