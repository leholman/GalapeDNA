################################################
####====== Galapagos eDNA Analysis ======####
####==== Luke E. Holman====01.05.2025====####
#############################################

### Script $ - Site point model data extraction
### Here we are taking the points along the path and then extracting model values for them


###############################################################################
# Point Extraction Script: Mean velocities at surface, top-20m, top-50m for:
#   (1) October subset (days 275 to 305 )
#   (2) Yearly average file (already time-mean)
#
# Then do nearest-neighbor extraction to pathpoints, producing 6 columns each
# for October and Yearly means.
###############################################################################

library(ncdf4)   # For NetCDF file operations
library(RANN)    # For fast nearest-neighbor search
library(dplyr)   # Optional, for data manipulation

# ----------------------------------------------------------------------------
# 1) Define File Paths
# ----------------------------------------------------------------------------
#UNIX - cdo timmean hiRGEMS_vel_2018.nc hiRGEMS_vel_2018_timmean.nc
grid_file       <- "OceanModelDataExtraction/hiRGEMS_grid.nc"
vel_file        <- "OceanModelDataExtraction/hiRGEMS_vel_2018.nc"           # for October subset
vel_file_mean   <- "OceanModelDataExtraction/hiRGEMS_vel_2018_timmean.nc"   # already time-mean (yearly)
hyd_file <- "OceanModelDataExtraction/hiRGEMS_hyd_2018.nc"
hyd_file_mean <-"OceanModelDataExtraction/hiRGEMS_hyd_2018_mean.nc"


# ----------------------------------------------------------------------------
# 2) Open and Read Grid Data
# ----------------------------------------------------------------------------
ncg <- nc_open(grid_file)

lon_grid <- ncvar_get(ncg, "XC")          # Longitudes [630]
lat_grid <- ncvar_get(ncg, "YC")          # Latitudes [768]
maskC    <- ncvar_get(ncg, "maskC")       # [630, 768, 75] (1=water,0=land)
Z        <- ncvar_get(ncg, "Z")           # [75], negative downward
nc_close(ncg)

# 2D surface mask = top cell
mask_surface <- maskC[ , , 1]
mask_flat    <- as.vector(mask_surface)

# ----------------------------------------------------------------------------
# 3) Define a vertical averaging function
# ----------------------------------------------------------------------------
vertical_average_top <- function(U_3d, Z, z_max, maskC = NULL) {
  z_idx <- which(Z >= -z_max)
  if (length(z_idx) == 0) {
    return(array(NA_real_, dim = c(dim(U_3d)[1], dim(U_3d)[2])))
  }
  
  U_subset <- U_3d[ , , z_idx, drop=FALSE]
  
  # Apply 3D land mask if provided
  if (!is.null(maskC)) {
    mask_subset <- maskC[ , , z_idx, drop=FALSE]
    U_subset[mask_subset == 0] <- NA  # mask out land
  }
  
  U_final  <- apply(U_subset, c(1, 2), mean, na.rm = TRUE)
  return(U_final)
}

# ----------------------------------------------------------------------------
# 4) Function to read a time-range from vel_file and average over time
# ----------------------------------------------------------------------------
read_timeavg_3D <- function(filename, varname, day1, day2) {
  nc    <- nc_open(filename)
  t_len <- day2 - day1 + 1
  
  # shape => [630,768,75,t_len]
  U_4d <- ncvar_get(
    nc, varid=varname,
    start=c(1,1,1,day1),
    count=c(-1,-1,-1,t_len)
  )
  nc_close(nc)
  
  # average over time dimension => [630,768,75]
  U_3d_mean <- apply(U_4d, c(1,2,3), mean, na.rm=TRUE)
  return(U_3d_mean)
}

# ----------------------------------------------------------------------------
# 5) Prepare pathpoints DataFrame
# ----------------------------------------------------------------------------
pathpoints <- read.csv("distanceData/pathPointsTable.csv")

# We'll add columns: U_oct_surf, V_oct_surf, U_oct_20, etc.

# ----------------------------------------------------------------------------
# 6) OCTOBER Subset (Days 274..304)
# ----------------------------------------------------------------------------
#### ocean movement
oct_t1 <- 275
oct_t2 <- 305  # about 31 days

# 6a) Read time-averaged U, V for Oct
u_oct_3d <- read_timeavg_3D(vel_file, "UE_VEL_C", oct_t1, oct_t2) # [630,768,75]
v_oct_3d <- read_timeavg_3D(vel_file, "VN_VEL_C", oct_t1, oct_t2)

# 6b) Surface is top layer => z=1
u_oct_surf <- u_oct_3d[ , , 1]
v_oct_surf <- v_oct_3d[ , , 1]

# 6c) Top 20m, Top 50m
u_oct_20 <- vertical_average_top(u_oct_3d, Z, 20)
v_oct_20 <- vertical_average_top(v_oct_3d, Z, 20)

u_oct_50 <- vertical_average_top(u_oct_3d, Z, 50)
v_oct_50 <- vertical_average_top(v_oct_3d, Z, 50)

#### ocean temp
temp_oct_3d <- read_timeavg_3D(hyd_file,"THETA", oct_t1, oct_t2)

temp_oct_surf <- temp_oct_3d[, ,1]
temp_oct_20 <- vertical_average_top(temp_oct_3d, Z, 20, maskC = maskC)
temp_oct_50 <- vertical_average_top(temp_oct_3d, Z, 50, maskC = maskC)

#image(lon_grid, lat_grid, temp_oct_50, main = "Oct Surface Temperature") 

# ----------------------------------------------------------------------------
# 7) YEARLY (Time-Mean) File
# ----------------------------------------------------------------------------
ncy <- nc_open(vel_file_mean)
u_year_3d <- ncvar_get(ncy, "UE_VEL_C", collapse_degen=FALSE) # shape [630,768,75] or [630,768,75,1]
v_year_3d <- ncvar_get(ncy, "VN_VEL_C", collapse_degen=FALSE)
nc_close(ncy)

if (length(dim(u_year_3d)) == 4) u_year_3d <- drop(u_year_3d) # => [630,768,75]
if (length(dim(v_year_3d)) == 4) v_year_3d <- drop(v_year_3d)


# 7a) Surface
u_year_surf <- u_year_3d[ , , 1]
v_year_surf <- v_year_3d[ , , 1]

# 7b) Top 20m, Top 50m
u_year_20 <- vertical_average_top(u_year_3d, Z, 20)
v_year_20 <- vertical_average_top(v_year_3d, Z, 20)

u_year_50 <- vertical_average_top(u_year_3d, Z, 50)
v_year_50 <- vertical_average_top(v_year_3d, Z, 50)

## ocean temp
ncy <- nc_open(hyd_file_mean)
temp_year_3d <-ncvar_get(ncy,"THETA", collapse_degen=FALSE)
nc_close(ncy)
temp_year_3d <-  drop(temp_year_3d)

temp_year_surf <- temp_year_3d[ , , 1]
temp_year_20 <- vertical_average_top(temp_year_3d,Z,20, maskC = maskC)
temp_year_50 <- vertical_average_top(temp_year_3d,Z,50, maskC = maskC)

image(lon_grid, lat_grid, temp_oct_20, main = "Year Surface Temperature") 


# ----------------------------------------------------------------------------
# 8) Nearest-Neighbor Extraction
#    We'll define a small function to flatten & extract to pathpoints
# ----------------------------------------------------------------------------
nn_extract_to_pathpoints <- function(u_2d, v_2d, pathDF, mask_surf, 
                                     lon_grid, lat_grid, colU, colV) {
  # Flatten
  u_flat <- as.vector(u_2d)
  v_flat <- as.vector(v_2d)
  m_flat <- as.vector(mask_surf)
  
  # Create grid_points for nn2
  nx <- length(lon_grid)
  ny <- length(lat_grid)
  grid_lon_2D <- outer(lon_grid, rep(1, ny))
  grid_lat_2D <- outer(rep(1, nx), lat_grid)
  grid_points <- cbind(as.vector(grid_lon_2D), as.vector(grid_lat_2D))
  
  # Query pathpoints
  query_sites <- as.matrix(pathDF[, c("lon", "lat")])
  nn_result   <- nn2(data=grid_points, query=query_sites, k=1)
  nearest_idx <- nn_result$nn.idx[,1]
  
  # Fill pathDF
  pathDF[[colU]] <- NA_real_
  pathDF[[colV]] <- NA_real_
  
  for (i in seq_len(nrow(pathDF))) {
    idx_grid <- nearest_idx[i]
    if (!is.na(m_flat[idx_grid]) && m_flat[idx_grid] == 1) {
      pathDF[[colU]][i] <- u_flat[idx_grid]
      pathDF[[colV]][i] <- v_flat[idx_grid]
    }
  }
  
  return(pathDF)
}

## a second simplified function for temp 

nn_extract_scalar <- function(scalar_2d, lon_grid, lat_grid, query_lon, query_lat, mask2d = NULL) {
  nx <- length(lon_grid)
  ny <- length(lat_grid)
  
  grid_lon_2D <- outer(lon_grid, rep(1, ny))
  grid_lat_2D <- outer(rep(1, nx), lat_grid)
  
  lon_vec <- as.vector(grid_lon_2D)
  lat_vec <- as.vector(grid_lat_2D)
  scalar_flat <- as.vector(scalar_2d)
  
  # Apply mask: only use ocean points for nearest-neighbor search
  if (!is.null(mask2d)) {
    mask_flat <- as.vector(mask2d)
    valid_idx <- which(!is.na(mask_flat) & mask_flat == 1)
    grid_points <- cbind(lon_vec[valid_idx], lat_vec[valid_idx])
    scalar_flat <- scalar_flat[valid_idx]
  } else {
    grid_points <- cbind(lon_vec, lat_vec)
  }
  
  query_points <- cbind(query_lon, query_lat)
  nn_result <- nn2(data=grid_points, query=query_points, k=1)
  nearest_idx <- nn_result$nn.idx[,1]
  
  return(scalar_flat[nearest_idx])
}

# ----------------------------------------------------------------------------
# 9) Add columns to pathpoints
# ----------------------------------------------------------------------------
# (A) October (surface, 20m, 50m)
pathpoints <- nn_extract_to_pathpoints(u_oct_surf, v_oct_surf, pathpoints, 
                                       mask_surface, lon_grid, lat_grid,
                                       "U_oct_surf", "V_oct_surf")
pathpoints <- nn_extract_to_pathpoints(u_oct_20,   v_oct_20,   pathpoints, 
                                       mask_surface, lon_grid, lat_grid,
                                       "U_oct_20", "V_oct_20")
pathpoints <- nn_extract_to_pathpoints(u_oct_50,   v_oct_50,   pathpoints, 
                                       mask_surface, lon_grid, lat_grid,
                                       "U_oct_50", "V_oct_50")

# (B) Yearly means (surface, 20m, 50m)
pathpoints <- nn_extract_to_pathpoints(u_year_surf, v_year_surf, pathpoints,
                                       mask_surface, lon_grid, lat_grid,
                                       "U_year_surf", "V_year_surf")
pathpoints <- nn_extract_to_pathpoints(u_year_20,   v_year_20,   pathpoints,
                                       mask_surface, lon_grid, lat_grid,
                                       "U_year_20", "V_year_20")
pathpoints <- nn_extract_to_pathpoints(u_year_50,   v_year_50,   pathpoints,
                                       mask_surface, lon_grid, lat_grid,
                                       "U_year_50", "V_year_50")

### temp out 

sitemetadata <- read.csv("metadata.site.out.csv")

sitemetadata$temp_year_surf <- nn_extract_scalar(temp_year_surf,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)
sitemetadata$temp_year_20 <- nn_extract_scalar(temp_year_20,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)
sitemetadata$temp_year_50 <- nn_extract_scalar(temp_year_50,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)
sitemetadata$temp_oct_surf <- nn_extract_scalar(temp_oct_surf,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)
sitemetadata$temp_oct_20 <- nn_extract_scalar(temp_oct_20,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)
sitemetadata$temp_oct_50 <- nn_extract_scalar(temp_oct_50,lon_grid, lat_grid,sitemetadata$lon2,sitemetadata$lat2,mask2d = mask_surface)

write.csv(sitemetadata,"distanceData/Temp.csv")

temp <- sitemetadata$temp_year_20
names(temp) <- sitemetadata$SiteID
temp.dist <- outer(temp,temp,"-")

temp.dist2 <- reshape2::melt(temp.dist, varnames = c("Start", "End"), value.name = "TempDiff")

write.csv(temp.dist2,"distanceData/Temp.csv")

# ----------------------------------------------------------------------------
# 10) Inspect final results
# ----------------------------------------------------------------------------
write.csv(pathpoints,"pathPointsTable.out.csv")





plot(pathpoints$V_oct_surf,pathpoints$V_oct_20,cex=0.2,pch=16)
plot(pathpoints$V_oct_surf,pathpoints$V_oct_50,cex=0.2,pch=16)
summary(lm(pathpoints$V_oct_surf~pathpoints$V_oct_20))
summary(lm(pathpoints$V_oct_surf~pathpoints$V_oct_50))


plot(pathpoints$U_oct_surf,pathpoints$U_oct_20,cex=0.2,pch=16)
plot(pathpoints$U_oct_surf,pathpoints$U_oct_50,cex=0.2,pch=16)
summary(lm(pathpoints$U_oct_surf~pathpoints$U_oct_20))
summary(lm(pathpoints$U_oct_surf~pathpoints$U_oct_50))


