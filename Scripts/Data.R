#### CASE STUDY - DATA ####


#libraries
library(raster)
library(terra)
library(ncdf4)
library(dplyr)
library(stringr)
library(purrr)

#working directory
setwd("D:/R/Wageningen/Seagriculture/Case study/IrradianceNovMay")

#### IRRADIANCE NOV 2019 MAY 2020 ####
#check one day
# Open file
IR.dat <- nc_open("IrradianceNovMay/20191101.nc")

# See structure and variables
print(IR.dat)

# List all variables
names(IR.dat$var)
# [1] "dni" "dif" "ghi"

# Extract the shortwave downward irradiance
IR.ghi <- ncvar_get(IR.dat, "ghi")
range(IR.ghi)
mean(IR.ghi)

# Close file
nc_close(IR.dat)

# Inspect
dim(IR.ghi) #Equal to the amount of seconds in a day (86400), seconds since 01 01 2020
head(IR.ghi)

plot(IR.ghi, type = "l", 
     main = "Shortwave Downward Irradiance (GHI)",
     xlab = "Time index", 
     ylab = "Irradiance (W/m²)")

###Calculate hour average
# Directory with the .nc files
data_dir <- "D:/R/Wageningen/Seagriculture/Case study/IrradianceNovMay"

# List all daily files (assumes format YYYYMMDD.nc)
files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE)

# Function to extract date from filename
get_date <- function(path) {
  fname <- basename(path)
  as.Date(str_sub(fname, 1, 8), format = "%Y%m%d")
}

# Function to compute hourly GHI averages for one file
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

# Process all files and combine
IR.result <- map_df(files, process_file)
head(IR.result)
range(IR.result$ghi_mean, na.rm = TRUE)
options(scipen = 999, digits = 6)
View(IR.result)

#Which day has the highest mean ghi value?
IR.result %>% 
  filter(!is.na(ghi_mean)) %>% 
  filter(ghi_mean == max(ghi_mean))

#upward shortwave flux is the reflected part of sunlight at the surface. It prevents some fraction of the incoming radiation from entering the water column, so it must be subtracted to get the radiation available to S. latissima. So NSW is downward - upward.
#Lon: 4.92774
#Lat: 51.968063

library(suncalc) 
library(lubridate)
LAT = 51.968063
LON = 4.92774
#Workflow
    #Split solar radiation into direct and diffuse components.
    #Convert direct normal irradiance → horizontal.
    #Calculate reflected fraction using Fresnel for direct and a fixed fraction for diffuse.
    #Compute upward reflected radiation.
    #Compute total downward radiation.
    #Subtract reflected upward from incoming → estimate net shortwave radiation.

# ---- helper Fresnel ----
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

# split direct/diffuse if available
# convert dni (direct normal) to direct horizontal: D_h = dni * cos(zenith)
direct_horiz <- if("dni_mean" %in% names(IR.result)) IR.result$dni_mean * pmax(0, cos(IR.result$zenith_rad)) else NA
#Direct normal irradiance (DNI) is measured perpendicular to the Sun’s rays.
#To get the horizontal component, we multiply by cos(zenith).
#pmax(0, cos(zenith)) ensures we don’t get negative values when the Sun is below the horizon.
#If dni_mean is not available, direct_horiz is set to NA

Rdir <- fresnel_reflectance(IR.result$zenith_rad) #So this gives the fraction of direct radiation reflected at the air-water interface. it depends on the zenith angle. More reflection with sharper degrees.
Rdiff <- 0.06   # Diffuse radiation is scattered light from the sky. Typically, surfaces reflect a small fraction of diffuse light. 0.06 is a reasonable default for water or similar surfaces.

up_direct <- if(!all(is.na(direct_horiz))) direct_horiz * Rdir else 0
#Upward radiation from direct sunlight = reflected fraction of horizontal direct radiation. If direct_horiz is missing (NA), assume 0.

up_diffuse <- if("dif_mean" %in% names(IR.result)) IR.result$dif_mean * Rdiff else IR.result$ghi_mean * 0.06
#Upward diffuse radiation = reflected fraction of diffuse radiation. If dif_mean is missing, approximate diffuse as 6% of total global horizontal irradiance (GHI).

down_total <- if("ghi_mean" %in% names(IR.result)) IR.result$ghi_mean else direct_horiz + IR.result$dif_mean
#Downward radiation on horizontal plane. use ghi_mean if available (simplest, most direct measure). Otherwise, sum direct horizontal + diffuse.

up_total <- up_direct + up_diffuse

IR.result$NSW_estimate <- down_total - up_total #So this is an estimate of net energy absorbed at the surface.

###Calculate PAR
# Constants
PARfrac <- 0.43       # fraction of light usable for photosynthesis
C <- 4.56             # μmol γ s⁻¹ W⁻¹
k <- 0.46             # extinction coefficient (m⁻¹)

# NSW from your dataset
NSW <- IR.result$NSW_estimate  # W/m²

# Depths to calculate PAR for
depths <- c(1, 2, 4.5, 7)

# Calculate PAR for each depth and add as new columns
for (z in depths) {
  IR.result[[paste0("PAR_", z, "m")]] <- NSW * PARfrac * C * exp(-k * z) * 3600 / 1000000
}

# View the first few rows
head(IR.result)

#Save
write.csv(IR.result, "IR_result_PAR.csv", row.names = FALSE)

#page 6 - section 2.4

#Calculate 3 hour average
#Use the formula
#Now the unit is right










#### NUTRIENTS ####
setwd("D:/R/Wageningen/Parasite Paradigms/Own model")
library(readr)
data <- read_tsv("AL557_tm - Nitrogen concentrations.tab")

readLines("AL557_tm - Nitrogen concentrations.tab", n = 200) #De eerste 159 zijn niet van belang.

library(readr)

data <- read_tsv("AL557_tm - Nitrogen concentrations.tab", skip = 159)

#DISSOLVED INORGANIC CARBON
setwd("D:/R/Wageningen/Parasite Paradigms/Own model")
library(ncdf4)

# Open file
nc2 <- nc_open("TCO2_NNGv2LDEO_climatology.nc")
# See structure and variables
print(nc2)

# List all variables
names(nc2$var)
# [1] "depth"          "TCO2_NNGv2LDEO"

# 
TCO2 <- ncvar_get(nc2, "TCO2_NNGv2LDEO")

# Close file
nc_close(nc2)

# Inspect
dim(TCO2)  
head(TCO2)

library(ggplot2)
library(reshape2)
library(raster)

#tco2 <- ncvar_get(nc2, "TCO2_NNGv2LDEO")

nc2 <- nc_open("TCO2_NNGv2LDEO_climatology.nc")
lon <- ncvar_get(nc2, "longitude")
lat <- ncvar_get(nc2, "latitude")
depth <- ncvar_get(nc2, "depth")
time <- ncvar_get(nc2, "time")
nc_close(nc2)

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


plot(TCO2, type = "l", 
     main = "Total dissolved inorganic C",
     xlab = "Time index", 
     ylab = "inorgnaic C (W/m²)")

#### Nitrate NOV 2019 MAY 2020 ####
setwd("D:/R/Wageningen/Seagriculture/Case study")
read.table("AL557_tm.tab")
read.table("AL557_tm.tab", sep="\t", header=TRUE, fill=TRUE, comment.char="")


#working directory
setwd("D:/R/Wageningen/Seagriculture/Case study")

#check one day
# Open file
N.dat <- nc_open("mercatorbiomer4v2r1_global_mean_nut_20220101.nc")

# See structure and variables
print(N.dat)

# List all variables
names(N.dat$var)
# "no3" "po4" "si"  "fe"

#Extract for the location with 51,683 lat 3,0662 lon
library(ncdf4)
library(dplyr)
library(stringr)

#NOV DEC JAN
# Directory containing your .nc files
data_dir2 <- "D:/R/Wageningen/Seagriculture/Case study/NitrateNovMay"

# Coordinates of interest
target_lat <- 51.683
target_lon <- 3.0662

#List all daily NetCDF files
nc_files <- list.files(data_dir2, pattern = "mercatorbiomer4v2r1_global_mean_nut_.*\\.nc$", full.names = TRUE)

#Extract data for each file
# Initialize empty dataframe
no3_ts <- data.frame(date = as.Date(character()), no3 = numeric())

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

# Sort chronologically
no3_ts <- arrange(no3_ts, date)

#Save
write.csv(no3_ts, "NitrateNovDecJan.csv", row.names = FALSE)

#FEB MAR APR MAY
# Directory containing your .nc files
data_dir3 <- "D:/R/Wageningen/Seagriculture/Case study/NitrateNovMay"

# Coordinates of interest
target_lat <- 51.683
target_lon <- 3.0662

#List all daily NetCDF files
nc_files2 <- list.files(data_dir3, pattern = "mercatorbiomer4v2r1_global_mean_nut_.*\\.nc$", full.names = TRUE, recursive = TRUE)


#Extract data for each file
# Initialize empty dataframe
no3_ts_FebMay <- data.frame(date = as.Date(character()), no3 = numeric())

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

# Sort chronologically
no3_ts_FebMay <- arrange(no3_ts_FebMay, date)

#Save
write.csv(no3_ts_FebMay, "NitrateFebMarAprMay.csv", row.names = FALSE)

#Merge the two csv files#
# Read the two CSV files
NovDecJan <- read.csv("NitrateNovMay/NitrateNovDecJan.csv")
FebMarAprMay <- read.csv("NitrateNovMay/NitrateFebMarAprMay.csv")

# Merge them
Nitrate_NovMay <- rbind(NovDecJan, FebMarAprMay)

#Save as csv
write.csv(Nitrate_NovMay, "Nitrate_NovMay.csv", row.names = FALSE)

