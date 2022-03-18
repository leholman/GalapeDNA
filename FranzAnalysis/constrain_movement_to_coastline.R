## test script for movement along coastline
# 2018-04-27

##
## 1. Set up
##
 ## -- call to core functionality -- ##
  # packages for core functionality
    library(magrittr)
    library(dplyr)
    library(tidyr)

  # functionality for spatial analyses
    library(raster)
    library(rgdal)
    library(sf)
    library(rgeos)

  # call to additional functionality
    library(viridis)
    library(gdistance)

  # set utm details
    utm_details <-
      paste0("+proj=utm +zone=15 +south +datum=WGS84 +units=m",
             " +no_defs +ellps=WGS84") %>% CRS()

  # set data locale     ## -- will need to change for local copy -- ##
    data_locale <- "galcosta/"

  # call to coastline
    galcosta <-
      paste0(data_locale,"galcosta.shp") %>%
      st_read()
      # shapefile()

##
## 2. Subset domain
##
   # set zoom for platform
     pzoom <-
       c(750000,
         820000,
         9935000 - 2.1E4,
         9975000 + 1.1E4) %>%
       extent()

  # constrain the problem
    plataforma <-
      galcosta %>% st_crop(pzoom)

 ## -- have a quick look -- ##
  # open window
    quartz("view platform", 10, 10)

  # visualise
    plataforma %>% dplyr::select(ID) %>% plot()

##
## 3. Create cost layer
##
  # set scale factor
    scale_factor <- 100

  # create empty raster
    empty_raster <-
      raster(resolution = scale_factor,
             xmn  = pzoom[1],
             xmx  = pzoom[2],
             ymn  = pzoom[3],
             ymx  = pzoom[4],
             crs  = utm_details)

  # convert to raster
    galcosta_r <-
      plataforma %>%
      st_geometry() %>%
      as("Spatial") %>%
      rasterize(empty_raster)

  # set arbitrary values for land & sea
    galcosta_r[galcosta_r > 0] <- 2
    galcosta_r[galcosta_r %>% is.na()] <- 1

 ## -- have a quick look -- ##
  # open window
    quartz("view platform", 5, 5)

  # visualise
    galcosta_r %>% plot(col = c("white", "grey70"))

##
## 4. Constrain movement
##
  # create points for movements
    from_here    <- c(810E3, 9920E3)
    to_here      <- c(760E3, 9940E3)
    to_over_here <- c(780E3, 9950E3)

  # set direction cells
    d_cells <- 16

  # create transition layer
    galcosta_t <-
      galcosta_r %>%
      asFactor() %>%
      transition("areas", d_cells)

  # perform geocorrection
    galcosta_t[[1]] %<>% geoCorrection()

  # create cost path
    coastal_path <-
      galcosta_t[[1]] %>%
      shortestPath(from_here,
                   to_here,
                   output = "SpatialLines")

  # test another point
    next_path <-
      galcosta_t[[1]] %>%
      shortestPath(from_here,
                   to_over_here,
                   output = "SpatialLines")

##
## 5. Generate outputs
##
 ## -- visualise -- ##
  # open window
    quartz("test plot", 7, 7)

  # plot points & paths
    galcosta_r   %>% plot(col = c("white", "grey70"))
    from_here[1] %>% points(from_here[2], pch = 19, col = viridis_pal()(5)[4])
    to_here[1]   %>% points(to_here[2],   pch = 19, col = viridis_pal()(5)[2])
    coastal_path %>% lines(col = viridis_pal()(5)[3])
    to_over_here[1] %>%
      points(to_over_here[2], pch = 19,
      col = viridis_pal(opti
                        
                        
                        
                        on = "magma")(7)[5])
    next_path %>% lines(col = viridis_pal(option = "magma")(7)[6])

 ## -- get distances -- ##
  # extract distance
    distance_m <-
      sapply(coastal_path %>% slot("lines"), function(x) LinesLength(x))


    test <- metadat.site
coordinates(test) <- ~lon2+lat2
    

mp <- st_multipoint(paste(metadat.site$lat2,metadat.site$lon2))
    
    
    