#############################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====01.05.2025====####
#############################################

### Script $ - Site <-> site distance and point data generation
### Here we are generating pathways between sites, calculating the distance along the path, and then identifying points along the path to extract oceanogrphic model data from  

# -----------------------------------------------
# 1) Setup: Load packages & set working directory
# -----------------------------------------------
setwd("~/GitHubRepos/GalapeDNA/")
library(gdistance)   # For shortestPath, transition matrices
library(geodata)     # For country admin boundaries (gadm)
library(terra)       # For modern spatial classes
library(raster)      # For raster adjacency functions, etc.

# -----------------------------------------------
# 2) Load data: bathymetry raster & admin polygons
# -----------------------------------------------
# A) Galápagos bathymetry (raster)
gebco.crop  <- readRDS("mapBuilding/GalapBathy.rds")
gebco.crop2 <- gebco.crop  # Make a copy

# B) Ecuador admin boundary (country-level = 0)
ecuador_gadm <- geodata::gadm(country = "ECU", level = 0, path = "data/")

# -----------------------------------------------
# 3) Prepare site coordinates
# -----------------------------------------------
metadatSites <- read.csv("metadata.site.out.csv", row.names = 1)
places       <- cbind(metadatSites$lon2, metadatSites$lat2)
rownames(places) <- metadatSites$SiteID

# -----------------------------------------------
# 4) Convert polygons to 'Spatial' & crop to Galápagos
# -----------------------------------------------
# 'gadm' data often comes as a SpatVector, so we convert:
ecuador_gadm_sp      <- as(ecuador_gadm, "Spatial")
ecuador_gadm_sp_crop <- crop(ecuador_gadm_sp, gebco.crop)

plot(ecuador_gadm_sp_crop, main = "Cropped Ecuador Admin (Galápagos Only)")

# -----------------------------------------------
# 5) Create a water/land raster
# -----------------------------------------------
# We'll rasterize the polygon:
#  - field=1 => land cells = 1
#  - background=0 => water cells = 0
land_raster <- rasterize(
  ecuador_gadm_sp_crop,  # polygon to rasterize
  gebco.crop2,           # template raster (same extent/res)
  field = 1,             
  background = 0
)

# Convert land=1 to NA => not traversable
# Convert water=0 to 1 => traversable
land_raster[land_raster == 1] <- NA  
land_raster[land_raster == 0] <- 1

# -----------------------------------------------
# 6) Build a simple "water-only" transition matrix
# -----------------------------------------------
# We'll call this tr, using mean() as the transition function 
# (for "passable" edges) and 16 directions (8 + diagonals).
tr <- transition(
  x                  = land_raster,
  transitionFunction = mean,
  directions         = 16
)

# Correct for lat/lon geometry (great-circle vs. degrees).
tr_corrected <- geoCorrection(tr, type = "c", multpl = FALSE)

# 'tr_corrected' is a valid TransitionLayer indicating 
# adjacency in water-only cells.

# -----------------------------------------------
# 7) Build a depth-based transition matrix
# -----------------------------------------------
# A) Copy bathymetry raster & set land cells to NA
depth_r <- gebco.crop2
depth_r[is.na(land_raster)] <- NA  # land => NA

# B) Define a function that penalizes big depth changes
alpha   <- 0.001   # tuning parameter
transFun <- function(x) {
  # x has length-2 => c(depthCell_i, depthCell_j)
  diff_depth <- abs(x[2] - x[1])
  exp(-alpha * diff_depth)
}

# C) Build a new transition matrix (depth-based)
tr_depth <- transition(
  x                  = depth_r,
  transitionFunction = transFun,
  directions         = 16
)
# Geo-correct it too:
tr_depth_corrected <- geoCorrection(tr_depth, type = "c", multpl = FALSE)

# -----------------------------------------------
# 8) Snap site coordinates that fall on land to water
# -----------------------------------------------
places2 <- places  # We'll store the "snapped" version here
vals    <- extract(land_raster, places2)  # extract cell values at site coords

# 'land_raster == 1' means water, NA means land or out-of-bounds
bad_idx <- which(is.na(vals))  # indices of sites that are on land

for (i in bad_idx) {
  # Find the cell ID at the site coordinate
  c0 <- cellFromXY(land_raster, places2[i, ])
  
  # Get 8 neighbors (include=TRUE => also the cell itself)
  nbrs  <- adjacent(land_raster, c0, directions = 8, include = TRUE)
  
  # Filter neighbors to find only water cells (land_raster == 1)
  water <- nbrs[ land_raster[nbrs[, 2]] == 1, 2 ]
  
  # If there is at least one water neighbor,
  # pick the one nearest to the original coordinate
  if (length(water) > 0) {
    # Coordinates of water cells
    xy_w <- xyFromCell(land_raster, water)
    
    # Euclidean distance in lat/lon degrees (small offsets => OK).
    d <- sqrt((xy_w[, 1] - places2[i, 1])^2 + (xy_w[, 2] - places2[i, 2])^2)
    
    # Snap to the closest water cell
    places2[i, ] <- xy_w[which.min(d), ]
  }
}

# -----------------------------------------------
# 9) Visualize: Plot paths from site 1 to others
# -----------------------------------------------
plot(gebco.crop2, 
     main   = "Shortest Paths: Simple (red) vs. Depth-based (blue)",
     axes   = FALSE, 
     box    = FALSE,
     legend = FALSE)

# Overlay the land polygon in gray
plot(ecuador_gadm_sp_crop, 
     add    = TRUE, 
     border = NA, 
     col    = "grey34")

# Example: For each site (2..23), draw:
#  - Simple water path in red (solid) & pink (dashed reverse)
#  - Depth-based path in darkblue (solid) & lightblue (dashed reverse)
for (site in 3:6) {
  # Simple water-only path
  lines(shortestPath(tr_corrected, places2[1,], places2[site,],
                     output = "SpatialLines"),
        col = "darkred", 
        lwd = 2)
  
  lines(shortestPath(tr_corrected, places2[site,], places2[1,],
                     output = "SpatialLines"),
        col = "pink",
        lty = 2,
        lwd = 2)
  
  # Depth-based path
  lines(shortestPath(tr_depth_corrected, places2[1,], places2[site,],
                     output = "SpatialLines"),
        col = "darkblue", 
        lwd = 2)
  
  lines(shortestPath(tr_depth_corrected, places2[site,], places2[1,],
                     output = "SpatialLines"),
        col = "lightblue",
        lty = 2,
        lwd = 2)
}


###Now lets generate the distance between sites and points along the path

# 1) Create all pairwise combinations of sites
sitePairwiseDist <- expand.grid(
  Start = metadatSites$SiteID,
  End   = metadatSites$SiteID,
  stringsAsFactors = FALSE
)

# Initialize a column for storing shortest-path distance
sitePairwiseDist$dist <- NA_real_

# Create an empty data.frame to store sampled path points
pathPointsTable <- data.frame()

# 2) Loop over each row in sitePairwiseDist
for (row in seq_len(nrow(sitePairwiseDist))) {
  
  LoopFrom <- sitePairwiseDist$Start[row]
  LoopTo   <- sitePairwiseDist$End[row]
  
  # If the two sites are the same, set dist=0 and skip path creation
  if (LoopFrom == LoopTo) {
    sitePairwiseDist$dist[row] <- 0
    next()
  }
  
  # Extract the (lon, lat) coordinates from 'places2'
  # Assumes rownames(places2) == metadatSites$SiteID
  fromXY <- c(places2[LoopFrom, 1], places2[LoopFrom, 2])
  toXY   <- c(places2[LoopTo,   1], places2[LoopTo,   2])
  
  # Compute the shortest path as a SpatialLines object
  loopPath <- shortestPath(
    x       = tr_depth_corrected, 
    origin  = fromXY,
    goal    = toXY,
    output  = "SpatialLines"
  )
  
  # Measure the length in km (because longlat=TRUE)
  loopLen <- SpatialLinesLengths(loopPath, longlat = TRUE)
  sitePairwiseDist$dist[row] <- loopLen
  
  # 3) Sample points along the path every 1 km
  # 'loopLen / 1' => how many 1-km segments (roughly)
  loopPathPoints <- spsample(loopPath, n = loopLen, type = "regular")
  
  # Build a small table for these sampled points
  nPoints <- length(loopPathPoints)
  if (nPoints > 0) {
    loopPointsTable <- data.frame(
      Order = seq_len(nPoints),
      Start = rep(LoopFrom, nPoints),
      End   = rep(LoopTo,   nPoints),
      lon   = loopPathPoints@coords[, 1],
      lat   = loopPathPoints@coords[, 2]
    )
    
    # Append to the main pathPointsTable
    pathPointsTable <- rbind(pathPointsTable, loopPointsTable)
  }
  
  # (Optional) print status
  print(row)
}


sitePairwiseDist.d <- dcast(sitePairwiseDist,Start~End,value.var = "dist")
sitePairwiseDist.d.out <- as.matrix(sitePairwiseDist.d[,2:24])
rownames(sitePairwiseDist.d.out) <- sitePairwiseDist.d$Start

write.csv(sitePairwiseDist.d.out,"distanceData/SiteDistanceMatrix.csv")
write.csv(sitePairwiseDist,"distanceData/SiteDistance.csv")
write.csv(pathPointsTable,"distanceData/pathPointsTable.csv")


## now we use the path points to extract data from the oceanogrphic model 








