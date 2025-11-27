####################################################################################################################
#Data gathering script: in this script, we have extracted environmental data for our specific location and timeframe from larger datasets. 
#Location: Offshore windmill farm Zeeland - Boroselle 1-4 (NL0B and NL0J) - Lat 51.683, Lon 3.0662.
#Timeframe: 1 November 2019 until 31 May 2020
#Environmental data: PAR (Photosynthetically active radiation), DIC (dissolved inorganic carbon), Nitrate,  Temperature

#File created by Luka Biemond and Nubia Middelkoop during nov-dec 2025. 
####################################################################################################################

#Reminder to set working directory to location of data
setwd("D:/R/Wageningen/Seagriculture/Case study")

#Import libraries
library(raster)
library(terra)
library(ncdf4)
library(dplyr)
library(stringr)
library(purrr)
library(suncalc) 
library(lubridate)
library(ggplot2)
library(reshape2)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Irradiance ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Check one day ----
# Open file
IR.dat <- nc_open("IrradianceNovMay/20191101.nc")

# See structure and variables
print(IR.dat)

# List all variables
names(IR.dat$var)
# [1] "dni" "dif" "ghi"
# The meaning of this: shortwave direct normal (to sun) irradiance, diffuse horizontal shortwave radiation, global horizontal shortwave irradiance.

# Extract the global horizontal shortwave downward irradiance
IR.ghi <- ncvar_get(IR.dat, "ghi")

# Close file
nc_close(IR.dat)

# Inspect
dim(IR.ghi) #Equal to the amount of seconds in a day (86400), seconds since 01 01 2020
head(IR.ghi)

plot(IR.ghi, type = "l", 
     main = "Shortwave Downward Irradiance (GHI)",
     xlab = "Time index", 
     ylab = "Irradiance (W/m²)")

# ---- Calculate hour average ----
# Directory with the .nc files
data_dir <- "D:/R/Wageningen/Seagriculture/Case study/IrradianceNovMay"

# List all daily files (assumes format YYYYMMDD.nc)
files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE)

# Function to extract date from filename
get_date <- function(path) {
  fname <- basename(path)
  as.Date(str_sub(fname, 1, 8), format = "%Y%m%d")
}

# Function to compute hourly GHI, DNI and DIF averages for one file
process_file <- function(file) {
  
  date <- get_date(file)
  nc <- nc_open(file)
  ghi <- ncvar_get(nc, "ghi")
  dni <- ncvar_get(nc, "dni")
  dif <- ncvar_get(nc, "dif")
  nc_close(nc)
  
    # Replace negatives with 0 -> BSRN QC procedures typically clamp nighttime values to 0 because irradiance cannot be negative: Pyranometers (like the CMP22 mentioned in our metadata) often report tiny negative values at night.
  ghi[ghi < 0] <- 0
  dni[dni < 0] <- 0
  dif[dif < 0] <- 0
  
  # Number of seconds per hour 
  block_size <- 3600   
  
  # Split into 24 blocks
  blocks.ghi <- split(ghi, rep(1:24, each = block_size))
  blocks.dni <- split(dni, rep(1:24, each = block_size))
  blocks.dif <- split(dif, rep(1:24, each = block_size))
  
  # Compute means
  block.ghi_means <- sapply(blocks.ghi, mean, na.rm = TRUE)
  block.dni_means <- sapply(blocks.dni, mean, na.rm = TRUE)
  block.dif_means <- sapply(blocks.dif, mean, na.rm = TRUE)
  
  # Return as tibble row
  tibble(
    date = date,
    block = 1:24,
    ghi_mean = block.ghi_means,
    dni_mean = block.dni_means,
    dif_mean = block.dif_means
  )
}

# CHECK: Process all files and combine
IR.result <- map_df(files, process_file)
head(IR.result)
range(IR.result$ghi_mean, na.rm = TRUE)
options(scipen = 999, digits = 6) 
View(IR.result)

#Which day has the highest mean ghi value?
IR.result %>% 
  filter(!is.na(ghi_mean)) %>% 
  filter(ghi_mean == max(ghi_mean))
#29-05-2020 has the highest mean ghi value, which makes sense because later in the growing season (november-may), sunlight is more intense.

#upward shortwave flux is the reflected part of sunlight at the surface. It prevents some fraction of the incoming radiation from entering the water column, so it must be subtracted to get the radiation available to S. latissima. So NSW is downward - upward.
#Remember: Lat 51.683, Lon 3.0662
LAT = 51.968063
LON = 4.92774
<<<<<<< Updated upstream
#Workflow
    #Split solar radiation into direct and diffuse components.
    #Convert direct normal irradiance → horizontal.
    #Calculate reflected fraction using Fresnel for direct and a fixed fraction for diffuse.
    #Compute upward reflected radiation.
    #Compute total downward radiation.
    #Subtract reflected upward from incoming → estimate net shortwave radiation.
=======
>>>>>>> Stashed changes

#Workflow
#0. Setting up the Fresnel function
#1. Split solar radiation into direct and diffuse components.
#2. Convert direct normal irradiance → horizontal.
#3. Calculate reflected fraction using Fresnel for direct and a fixed fraction for diffuse.
#4. Compute total downward radiation.
#5. Compute upward reflected radiation.
#6. Subtract reflected upward from incoming → estimate net shortwave radiation (NSW).
#7. Calculate PAR based on estimated NSW

# ---- 0. helper Fresnel ----
#this function gives the fraction of light reflected at the surface (between 0 and 1), given an incidence angle (theta_i_rad) and an refractive index (n). n describes how much light slows down when it enters a material compared to vaccum and how light bends at the interface between two materials. The default from water -> air = 1.33
fresnel_reflectance <- function(theta_i_rad, n = 1.33) {
  # ensure theta_i in [0, pi/2]
  theta_i_rad[theta_i_rad > pi/2] <- pi/2
  theta_t <- asin( pmin(1, sin(theta_i_rad)/n) )
  rs <- (sin(theta_i_rad - theta_t) / sin(theta_i_rad + theta_t))^2    #standard Fresnel formula for perpendicular light.
  rp <- (tan(theta_i_rad - theta_t) / tan(theta_i_rad + theta_t))^2   #standard Fresnel formula for parallel light.
  0.5 * (rs + rp)    #Unpolarized light (such as sunlight) has equal contributions of s- and p-polarization, so the total reflectance is the average.
}

# IR.result must have: date/time, ghi_mean, (optional) dni_mean, dif_mean
# compute sun position
IR.result <- IR.result %>%
  mutate(
    date_time = as.POSIXct(paste(date, block), format="%Y-%m-%d %H", tz="UTC")
  )
pos <- getSunlightPosition(date = IR.result$date_time, lat = LAT, lon = LON)     #computes the position of the Sun for given datetime and location.
alt <- pos$altitude   # Extract the Sun’s altitude, so the angle above the horizon. in radians.
zenith <- pi/2 - alt  # Zenith angle = angle between the Sun and the vertical direction (straight up).
IR.result$zenith_rad <- zenith

# ---- 1.  Split solar radiation into direct and diffuse components ----
#This step is in our case redundant, because the data already disinguishes between direct (dni) and diffuse (dif) components.

# ---- 2.  Convert direct normal irradiance → horizontal. ----
# convert dni (direct normal) to direct horizontal: D_h = dni * cos(zenith)
direct_horiz <- if("dni_mean" %in% names(IR.result)) IR.result$dni_mean * pmax(0, cos(IR.result$zenith_rad)) else NA
#Direct normal irradiance (DNI) is measured perpendicular to the Sun’s rays.
#To get the horizontal component, we multiply by cos(zenith).
#pmax(0, cos(zenith)) ensures we don’t get negative values when the Sun is below the horizon.
#If dni_mean is not available, direct_horiz is set to NA

# ---- 3. Calculate reflected fraction using Fresnel ----
Rdir <- fresnel_reflectance(IR.result$zenith_rad) #So this gives the fraction of direct radiation reflected at the air-water interface. it depends on the zenith angle. More reflection with sharper degrees.
Rdiff <- 0.06   # Diffuse radiation is scattered light from the sky. Typically, surfaces reflect a small fraction of diffuse light. 0.06 is a reasonable default for water or similar surfaces.

up_direct <- if(!all(is.na(direct_horiz))) direct_horiz * Rdir else 0
#Upward radiation from direct sunlight = reflected fraction of horizontal direct radiation. If direct_horiz is missing (NA), assume 0.

up_diffuse <- if("dif_mean" %in% names(IR.result)) IR.result$dif_mean * Rdiff else IR.result$ghi_mean * 0.06
#Upward diffuse radiation = reflected fraction of diffuse radiation. If dif_mean is missing, approximate diffuse as 6% of total global horizontal irradiance (GHI).

# ---- 4. Compute total downward radiation  ----
down_total <- if("ghi_mean" %in% names(IR.result)) IR.result$ghi_mean else direct_horiz + IR.result$dif_mean
#Downward radiation on horizontal plane. use ghi_mean if available (simplest, most direct measure). Otherwise, sum direct horizontal + diffuse.

# ---- 5. Compute upward reflected radiation.  ----
up_total <- up_direct + up_diffuse
 
# ---- 6. Estimate NSW ----
IR.result$NSW_estimate <- down_total - up_total #So this is an estimate of net energy absorbed at the surface.

# ---- 7. Calculate PAR based on estimated NSW ----
# Constants taken from the paper by Venolia.
PARfrac <- 0.43       # fraction of light usable for photosynthesis
C <- 4.56             # μmol γ s⁻¹ W⁻¹
k <- 0.46             # extinction coefficient (m⁻¹)

# NSW from our dataset
NSW <- IR.result$NSW_estimate  # W/m²

# Depths to calculate PAR for
depths <- c(1, 2, 4.5, 7)

# Calculate PAR for each depth and add as new columns
for (z in depths) {
  IR.result[[paste0("PAR_", z, "m")]] <- NSW * PARfrac * C * exp(-k * z) * 3600 / 1000000
}

# View the first few rows
head(IR.result)

#Save the result as a .csv file
write.csv(IR.result, "IR_result_PAR.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Dissolved inorganic carbon (DIC) ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Explore the dissolved carbon data ----
C.dat <- nc_open("TCO2_NNGv2LDEO_climatology.nc")
# See structure and variables
print(C.dat)

# List all variables
names(C.dat$var)
# [1] "depth"          "TCO2_NNGv2LDEO"

# ---- Extracting variables ----
# We want to extract TCO2_NNGv2LDEO because this is the monthly climatology of total dissolved inorganic carbon (TCO2)
TCO2 <- ncvar_get(C.dat, "TCO2_NNGv2LDEO")
dim(TCO2)
# [1]  12 102 168 360
# 12 refers to the numbers of the month, 102 to the amount of depth levels, 168 to latitude range and 360 to longitude range.

#We also want to extract the longitude, latitude, depth and time and store them as numerics. Later we will use it to extract data specific to our location and timeframe.
lon <- ncvar_get(C.dat, "longitude")
lat <- ncvar_get(C.dat, "latitude")
depth <- ncvar_get(C.dat, "depth")
time <- ncvar_get(C.dat, "time")
nc_close(C.dat)

# ---- Visualising a slice of the data ----
#Because TCO2 is a large array, you can not visualize the whole dataset in one plot. Therefore we chose an arbitrary month and depth (both 1) to create and visualize a slice of the data to get an idea about how the data looks like.
depth_index <- 1  # surface
time_index <- 1   # January
tco2_slice <- TCO2[time_index, depth_index, , ]  # dims: lat x lon

# Create a grid
lonlat <- expand.grid(lon=lon, lat=lat)
tco2_df <- cbind(lonlat, TCO2 = as.vector(tco2_slice))

ggplot(tco2_df, aes(x = lon, y = lat, fill = TCO2)) +
  geom_raster() +
  scale_fill_viridis_c(name = "TCO2 (umol/kg)") +
  coord_quickmap() +
  theme_minimal() +
  labs(title = "Surface Total Dissolved Inorganic Carbon",
       subtitle = "Month 1, Surface level")
#This graph shows the amount of total dissolved Inorganic carbon across the world in January.

# ---- Extract DIC values for our location and timeframe ----
# Coordinates: these are the coordinates of the Borselle windfarm.
target_lat <- 51.683
target_lon <- 3.0662

# Find nearest indices: so find the index for which the difference between the latitude in the dataframe and the latitud we target is smallest.
lat_index <- which.min(abs(lat - target_lat))
lon_index <- which.min(abs(lon - target_lon))

lat[lat_index]
lon[lon_index]
#The latitude and longitude that are most close to our location are lat 51.5 and lon 3.5.

#Specifying target depths
target_depths <- c(0, 5, 10) #We take the three most shallow depths, because the deeper you go, the more scarce the light needed for photosynthesis will be.

depth_indices <- sapply(target_depths, function(z) {
  which.min(abs(depth - z))
})

depth_indices

#Selecting the timeframe: we select October till June so that we later can make a subset from the first of November till the 31st of May.
time_indices <- c(10, 11, 12, 1, 2, 3, 4, 5, 6)
months <- c("Oct", "Nov", "Dec", "Jan", "Feb", "Mar", "Apr", "May", "Jun")

# Create an empty output matrix that we later can fill with our extracted values.
tco2_mat <- matrix(NA, nrow = length(time_indices), ncol = length(depth_indices))

for (i in seq_along(time_indices)) {
  for (j in seq_along(depth_indices)) {
    tco2_mat[i, j] <- TCO2[time_indices[i], depth_indices[j], lat_index, lon_index]
  }
}

TCO2_monthly <- data.frame(
  Month = months,
  Depth_0m  = tco2_mat[, 1],
  Depth_5m  = tco2_mat[, 2],
  Depth_10m = tco2_mat[, 3]
)

TCO2_monthly
#Only NAN values because the dataset sees our coordinates as land. So we need to find the nearest grid that does have data.

# ---- Extract DIC values that are not NA ----
# Function to find the nearest non-NaN grid cell around the target
find_nearest_ocean_cell <- function(target_lat, target_lon, lat, lon, TCO2, max_radius = 20) {
  
  # nearest grid indices
  lat_idx0 <- which.min(abs(lat - target_lat))
  lon_idx0 <- which.min(abs(lon - target_lon))
  
  for (r in 0:max_radius) {
    
    lat_candidates <- (lat_idx0 - r):(lat_idx0 + r)
    lon_candidates <- (lon_idx0 - r):(lon_idx0 + r)
    
    # keep valid indices
    lat_candidates <- lat_candidates[lat_candidates >= 1 & lat_candidates <= length(lat)]
    lon_candidates <- lon_candidates[lon_candidates >= 1 & lon_candidates <= length(lon)]
    
    # search through candidates
    for (i in lat_candidates) {
      for (j in lon_candidates) {
        
        # Check if any depth/time layer contains valid ocean values
        if (!all(is.na(TCO2[, , i, j]))) {
          return(list(lat_index = i,
                      lon_index = j,
                      lat = lat[i],
                      lon = lon[j],
                      radius = r))
        }
      }
    }
  }
  
  return(NULL)
}

# Use the function to find the nearest grid cell that does contain data
nearest <- find_nearest_ocean_cell(
  target_lat = 51.683,
  target_lon = 3.0662,
  lat = lat,
  lon = lon,
  TCO2 = TCO2
)

nearest 
# The nearest grid cell that does contain data has the coordinates 51.5 lat, 2.5 lon. Now we can make the dataframe again with the new coordinate indices.
lat_index2 <- nearest$lat_index
lon_index2 <- nearest$lon_index

lat[lat_index2]
lon[lon_index2]
# The correct nearest coordinates were extracted. 

# Create a new empty output matrix
tco2_mat <- matrix(NA, nrow = length(time_indices), ncol = length(depth_indices))

for (i in seq_along(time_indices)) {
  for (j in seq_along(depth_indices)) {
    tco2_mat[i, j] <- TCO2[time_indices[i], depth_indices[j], lat_index2, lon_index2]
  }
}

# We overwrite the former empty TCO2_monthly
TCO2_monthly <- data.frame(
  Month = months,
  Depth_0m  = tco2_mat[, 1],
  Depth_5m  = tco2_mat[, 2],
  Depth_10m = tco2_mat[, 3]
)

TCO2_monthly #The dataframe contains values now.

<<<<<<< Updated upstream
#### Nitrate NOV 2019 MAY 2020 ####
setwd("D:/R/Wageningen/Seagriculture/Case study")
read.table("AL557_tm.tab")
read.table("AL557_tm.tab", sep="\t", header=TRUE, fill=TRUE, comment.char="")


#working directory
setwd("D:/R/Wageningen/Seagriculture/Case study")

#check one day
# Open file
=======
#Save the data as a .csv file.
write.csv(TCO2_monthly, "DICOctJun_Depths0_5_10.csv", row.names = FALSE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Nitrate ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Explore the nitrate data ----
# Check one day, in this case, the first of january in 2022.
>>>>>>> Stashed changes
N.dat <- nc_open("mercatorbiomer4v2r1_global_mean_nut_20220101.nc")

# See structure and variables
print(N.dat)

# List all variables
names(N.dat$var)
# "no3" "po4" "si"  "fe"

#The datafiles have a large size, so we extract it in two parts that we later merge.
# ---- Extract the data for our specific location in November, December and January ----
# Directory containing the .nc files
data_dir2 <- "D:/R/Wageningen/Seagriculture/Case study/NitrateNovMay"

#List all daily nsc files
nc_files <- list.files(data_dir2, pattern = "mercatorbiomer4v2r1_global_mean_nut_.*\\.nc$", full.names = TRUE)

#Extract data for each file
# Initialize empty dataframe
no3_ts <- data.frame(date = as.Date(character()), no3 = numeric())

#A function to extract all the data we want from the individuals nsc files (one for each day)
for (f in nc_files) {
  nc <- nc_open(f)
  
  # Read coordinate grids
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  
  # Find nearest grid cell
  lon_idx <- which.min(abs(lon - target_lon))
  lat_idx <- which.min(abs(lat - target_lat))
  
  # Extract nitrate at surface (depth index = 1)
  no3_val <- ncvar_get(nc, "no3", start = c(lon_idx, lat_idx, 1, 1), count = c(1, 1, 1, 1))
  
  # Parse date from filename (e.g., "20220101.nc")
  date_str <- str_extract(basename(f), "\\d{8}")
  date_val <- as.Date(date_str, format = "%Y%m%d")
  
  # Append to dataframe
  no3_ts <- rbind(no3_ts, data.frame(date = date_val, no3 = no3_val))
  
  nc_close(nc)
}

# Sortthe extracted data chronologically
no3_ts <- arrange(no3_ts, date)

#Save the data for November, December and January as a .csv file.
write.csv(no3_ts, "NitrateNovDecJan.csv", row.names = FALSE)

# ---- Extract the data for our specific location in February, March, April and May ----
#List all daily .nsc files
nc_files2 <- list.files(data_dir2, pattern = "mercatorbiomer4v2r1_global_mean_nut_.*\\.nc$", full.names = TRUE, recursive = TRUE)

#Extract data for each file
# Initialize empty dataframe
no3_ts_FebMay <- data.frame(date = as.Date(character()), no3 = numeric())

#A function to extract all the data we want from the individuals nsc files (one for each day)
for (f in nc_files2) {
  nc <- nc_open(f)
  
  # Read coordinate grids
  lon <- ncvar_get(nc, "longitude")
  lat <- ncvar_get(nc, "latitude")
  
  # Find nearest grid cell
  lon_idx <- which.min(abs(lon - target_lon))
  lat_idx <- which.min(abs(lat - target_lat))
  
  # Extract nitrate at surface (depth index = 1)
  no3_val <- ncvar_get(nc, "no3", start = c(lon_idx, lat_idx, 1, 1), count = c(1, 1, 1, 1))
  
  # Parse date from filename (e.g., "20220101.nc")
  date_str <- str_extract(basename(f), "\\d{8}")
  date_val <- as.Date(date_str, format = "%Y%m%d")
  
  # Append to dataframe
  no3_ts_FebMay <- rbind(no3_ts_FebMay, data.frame(date = date_val, no3 = no3_val))
  
  nc_close(nc)
}

# Sort the extracted data chronologically
no3_ts_FebMay <- arrange(no3_ts_FebMay, date)

#Save the data for February, March, April and May as a .csv file. 
write.csv(no3_ts_FebMay, "NitrateFebMarAprMay.csv", row.names = FALSE)

# ---- Merge the two .csv files ----
# Read the two CSV files
NovDecJan <- read.csv("NitrateNovMay/NitrateNovDecJan.csv")
FebMarAprMay <- read.csv("NitrateNovMay/NitrateFebMarAprMay.csv")

# Merge them
Nitrate_NovMay <- rbind(NovDecJan, FebMarAprMay)

<<<<<<< Updated upstream
#Save as csv
write.csv(Nitrate_NovMay, "Nitrate_NovMay.csv", row.names = FALSE)

=======
#Save the final dataset for nitrate as a .csv file
write.csv(Nitrate_NovMay, "Nitrate_NovMay.csv", row.names = FALSE)
>>>>>>> Stashed changes
