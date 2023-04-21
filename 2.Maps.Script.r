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

plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend=F)

#-92, -89, -1.55, 0.68)
axis(1,
     at=pretty(c(-92, -89)),
     labels = parse(text=degreeLabelsEW(pretty(c(-92, -89)))))
axis(2,
     at=pretty(c(-1.55, 0.68)),
     labels = parse(text=degreeLabelsNS(pretty(c(-1.55, 0.68)))),las=TRUE)


#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=1)
points(metadata$lon2,metadata$lat2,pch=16, col='white',cex=0.5)

dev.off()

#Adding in colours per island
cols <- c("#80B1D3","#FFFFB3","#FFFFB3","#80B1D3","#FB8072","#BEBADA","#FFED6F","#CCEBC5",
          "#80B1D3","#8DD3C7","#FDB462","#FFFFB3","#B3DE69","gray85","#FFFFB3","#80B1D3",
          "#FCCDE5","#80B1D3","#80B1D3","#CCEBC5","#BC80BD","#8DD3C7","#80B1D3")



##Make a colour index so each ecoregion-island has a different colour

SEasternCols <- colorRampPalette(c("#D55E00","#E69F00","#faf6c1"))
#NorthernCols <-colorRampPalette(c("#0072B2","#56B4E9"))
NorthernCols <-colorRampPalette(c("#003e60","#56B4E9"))

ElizCols <- colorRampPalette(c("#CC79A7","#762d55"))
WesternCols <- colorRampPalette(c("#009E73","#004935"))

AllCols <- c(SEasternCols(8),ElizCols(2),NorthernCols(4),WesternCols(2))
ColIndex <- data.frame("EcoIsland"=sort(unique(paste0(metadat.site$EcoRegion,"-",metadat.site$island))),
                       "Colour"=AllCols)

ColIndex$Colour[match(paste0(metadat.site$EcoRegion,"-",metadat.site$island),ColIndex$EcoIsland)]



pdf(width = 8,height=6.35,file="mapBuilding/mapV1.pdf")

plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend=F)
#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=3)
points(metadata$lon2,metadata$lat2,pch=16, col=ColIndex$Colour[match(paste0(metadat.site$EcoRegion,"-",metadat.site$island),ColIndex$EcoIsland)]
,cex=2)

axis(1,
     at=pretty(c(-92, -89)),
     labels = parse(text=degreeLabelsEW(pretty(c(-92, -89)))),lwd = 0, lwd.ticks = 1)
axis(2,
     at=pretty(c(-1.55, 0.68)),
     labels = parse(text=degreeLabelsNS(pretty(c(-1.55, 0.68)))),las=TRUE,lwd = 0, lwd.ticks = 1)


plot(gebco.crop, legend.only=TRUE, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])),breaks=galap.br[[3]],
     legend.width=1, legend.shrink=0.75,font=2, line=2.5, cex=0.8,
     axis.args=list( at=pretty(galap.br[[3]][1:43]), labels=pretty(galap.br[[3]][1:43])))


dev.off()





### playground basement
plot(gebco.crop, col=c(blue.col(galap.br[[1]]), grey.colors(galap.br[[2]])), breaks=galap.br[[3]],axes = FALSE,box=F,legend = FALSE)
#plot(eq, add=TRUE)
points(metadata$lon2,metadata$lat2,pch=16, col='black',cex=1)
points(metadata$lon2,metadata$lat2,pch=16, col='white',cex=0.5)

text(metadata$lon2,metadata$lat2+0.08,labels=metadata$SiteID,cex=0.5,col="green")


site <- "RED"
site <- "SUAR"

points(pathPointsTable$lon[pathPointsTable$Start==site],
       pathPointsTable$lat[pathPointsTable$Start==site],
       col="lightgreen",
       pch=16,
       cex=0.2)




