####################################################################################################################
#Run file of Sugar Kelp model from Venolia et al (in press)
#Site names here begin with names other than those used in the manuscript
#NS_ZL = North Sea, near Zeeland 

#File created by Celeste Venolia in March 2018-December 2019
#Edited by Luka ... and Nubia Middelkoop during nov-dec 2025 
####################################################################################################################

#Reminder to set working directory to location of data
setwd("/Users/nubia/PycharmProjects/seaweedTempsNorthSea/Scripts")

#Import libraries
library(deSolve)
library(tidyverse)
library(lubridate)
library(gridExtra)
library(gdata)
library(Metrics)
library(pracma)
library(glue)


#Required for model runs
source("SolveR_R.R")
source("KelpDEB_model.R")
#Required for Calibration Code
source("N_uptake_Calibration.R")
source("Photosynthesis_Calibration.R")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Minerals and Organics Section #####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Conversion coefficients, organics (n = matrix of chemical indices)
# "food N" "food C" Stucture "N reserves" "C reserves" products
#     X_N   X_C      V    E_N    E_C    P
n_O <- matrix(
  + c(0.00, 1.00, 1.00, 0.00, 1.00, 1.00,  #C/C, equals 1 by definition
      + 0.00, 0.50, 1.33, 0.00, 2.00, 1.80,  #H/C, these values show that we consider dry-mass
      + 3.00, 2.50, 1.00, 2.50, 1.00, 0.50,  #O/C
      + 1.00, 0.00, 0.04, 1.00, 0.00, 0.04), nrow=4, ncol=6, byrow = TRUE) #N/C
#V is the C-mol structure of alginate (Alginic acid: (C6H8O6)n)
#E_N is N03- and N02- averaged
#E_C is glucose C6H12O6 (Laminarin: c18h32o16 and mannitol c6h14o6)
#We aren't using the X_N, X_C, or P collumn here

#Molecular weights
#t() is a matrix transpose function
#organics structure matrix multiplied by the atomic masses (mass in grams of one mole of an element) of C H O N
w_O_step <- t(n_O)*matrix(c(12, 1, 16, 14), nrow=6, ncol=4, byrow= TRUE) #g/mol, molecular weights for organics
w_O <- rowSums(w_O_step) #this provides g/mol of each of the six "pockets of mass" (i.e. X_N, X_C) Remember:  X_N   X_C      V    E_N    E_C    P

#define molecular weights
w_V <- w_O[3]  # g/mol       #molecular weight of structure
w_EN <- w_O[4]  # g/mol      #molecular weight of N reserve
w_EC <- w_O[5]  #g/mol       #molecular weight of C reserve
w_O2 <- 32 #g/mol

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Parameters compiled #####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
params_NS <- c(#maximum volume-specific assimilation rate of N before temperature correction
  JENAM = 1.5e-4, #mol N / molV / h
  #half saturation constant of N uptake
  K_N = 2.5e-6, #molNO3 and NO2/L
  #max volume-specific carbon dioxide assimilation rate
  JCO2M = 0.0075, #mol DIC/molV/h
  #half saturation constant of C uptake
  K_C = 4e-7, #mol DIC/L
  #maximum volume-specific carbon assimilation rate
  JECAM = 0.282, #molC/molV/h
  #Photosynthetic unit density
  rho_PSU = 0.5, #mol PSU/ mol V
  #binding probability of photons to a Light SU
  b_I = 0.5, #dimensionless
  #Specific photon arrival cross section
  alpha = 1, #m^2 mol PSU–1
  #dissociation rate
  k_I = 0.075, #molγ molPS–1 h–1
  #Yield factor of C reserve to photon
  y_I_C = 10, #mol γ mol C-1
  #Yield factor of C reserve to DIC
  y_CO2_C = 1, #mol DIC mol C-1
  #Yield factor of photon to O2
  y_LO2 = 0.125, #molO2 molγ –1
  #reserve turnover
  kE_C = 0.02, #0.05, #1/h
  kE_N = 0.04, #0.01, #1/h
  #fraction of rejection flux from growth SU incorporated back into i-reserve
  kappa_Ei = 0.9, #dimensionless
  #yield of structure on N reserve (percent of N in structure)
  y_EN_V = 0.04, #mol N/mol V
  #yield of structure on C reserve (percent of C in structure)
  y_EC_V = 1, #mol C/mol V
  #specific maintenance costs requiring N before temp correction
  JENM = 4*10^-6, #4e-6, #mol N/molM_V/h
  #specific maintenance costs requiring C before temp correction
  JECM = 1*10^-6, #1e-6, #mol C/molM_V/h
  #Arrhenius temperature
  T_A = 6314.3, # K
  #Upper boundary of temperature tolerance
  T_H = 13.386 + 273.15, # K
  #Lower boundary of temperature tolerance
  T_L = 273.15, # K
  #Arrhenius temperature outside T_H
  T_AH = 18702, #K
  #Arrhenius temperature outside T_L
  T_AL = 4391.9, #K
  #temperature at which rate parameters are given
  T_0 = 20 + 273.15) # K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
####### State Initial conditions ############
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Initial conditions of state variables
#these values are not coming from any field data or literature information, estimated
state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

state_Johansson <- c(m_EC = 0.3, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                     m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                     M_V = 0.05/(w_V+0.01*w_EN+0.3*w_EC)) #molM_V #initial mass of structure

# Our own state initial conditions that we calculated from European data (Broch at all (2012)) 
state_Clemente <- c(m_EC = 0.321, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                    m_EN = 0.017, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                    M_V = 0.05/(w_V+0.01*w_EN+0.3*w_EC)) #molM_V #initial mass of structure

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#######Time step to run the model on#######
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#(First number of time step, last number of time step, interval to step)
times_NS <- seq(1, 5112, 1) #213 days stepped hourly

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Irradiance
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Irradiance <- read.csv("IR_result_PAR.csv", header = TRUE)

Irradiance$date_time[Irradiance$block == 24] <- paste0(Irradiance$date_time[Irradiance$block == 24], " 00:00:00")

Irradiance$date_time <- ymd_hms(Irradiance$date_time, tz = "UTC") #NOAA data in UTC (5 hours ahead)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Nitrate 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nitrate <- read.csv("Nitrate_NovMay.csv") 

nitrate_hourly <- nitrate %>%
  # create 24 rows per date (0–23 hours)
  mutate(hour = list(0:23)) %>%
  tidyr::unnest(hour) %>%          # specify tidyr:: just to be safe
  dplyr::mutate(
    datetime = as.POSIXct(date) + lubridate::hours(hour)
  ) %>%
  dplyr::select(datetime, no3) %>%
  dplyr::arrange(datetime)

dupe_rows <- nitrate_hourly[2857:2880, ]

# insert them after 2880
nitrate_hourly <- rbind(
  nitrate_hourly[1:2880, ],      # everything up to row 2880
  dupe_rows,         # duplicated rows for Feb 29
  nitrate_hourly[2881:nrow(nitrate_hourly), ]  # the rest of the original data
)

# N unit is now in: millimol/m^3, we convert it to mol/L 
nitrate_hourly$no3 <- nitrate_hourly$no3/1000000

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Sea surface temperature 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
temp <- read.csv("temperatue_20192020.csv")
temp <- temp[7273:12385,]
temp <- temp %>%
  mutate(
    # parse date
    dates = ymd(dates),
    # ensure hour is numeric
    hour = as.numeric(hour),
    # create datetime
    date_time = dates + hours(hour),
    TZ = TZ/10
  )
temp$TZ_K <- temp$TZ+273.15 #kelvin
temp <- temp[-nrow(temp), ]
#T_field <- make_function(temp$TZ_K)
TZ_K <- temp$TZ_K

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### DIC forcing ###
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Dit mist de maand mei, omdat als je lineair interpoleerd doe je tussen twee punten, nu pakt hij de eerste van de maand als de punt dat hij ingelezen heeft. Maar dan is er niks aan t einde van mei. Dus moeten we ook de maand Juni in lezen anders is er in de maand mei niks. En eigen lijk ook kijken waneer de meting is genomen in de maand. 

TCO2_monthly <- read.csv("DICOctJun_Depths0_5_10.csv")
density_seawater <- 1.026 #kilo/liter
TCO2_monthly$Depth_0m <- TCO2_monthly$Depth_0m/1000000 * density_seawater
TCO2_monthly$Depth_5m <- TCO2_monthly$Depth_5m/1000000 * density_seawater
TCO2_monthly$Depth_10m <- TCO2_monthly$Depth_10m/1000000 * density_seawater
# Monthly timestamps (start of each month)
month_dates <- as.POSIXct(c(
  "2019-10-15 00:00",
  "2019-11-15 00:00",
  "2019-12-15 00:00",
  "2020-01-15 00:00",
  "2020-02-15 00:00",
  "2020-03-15 00:00",
  "2020-04-15 00:00",
  "2020-05-15 00:00",
  "2020-06-15 00:00"
), tz = "UTC")

# Hourly sequence from 1 Nov 2019 (new = 15 Oct 2019) to 31 May 2020 (nes = 15 of june) (inclusive)
hourly_time <- seq(
  from = as.POSIXct("2019-10-15 00:00", tz="UTC"),
  to   = as.POSIXct("2020-06-15 00:00", tz="UTC"),
  by   = "1 hour"
)

length(hourly_time)


t_months_num  <- as.numeric(month_dates)
t_hours_num   <- as.numeric(hourly_time)


f_0m  <- approxfun(t_months_num, TCO2_monthly$Depth_0m,  method = "linear")
f_5m  <- approxfun(t_months_num, TCO2_monthly$Depth_5m,  method = "linear")
f_10m <- approxfun(t_months_num, TCO2_monthly$Depth_10m, method = "linear")


TCO2_hourly <- data.frame(
  Datetime   = hourly_time,
  Depth_0m   = f_0m(t_hours_num),
  Depth_5m   = f_5m(t_hours_num),
  Depth_10m  = f_10m(t_hours_num)
)

plotplot_TCO2 <- ggplot() + 
  geom_line(data = TCO2_hourly, aes(Datetime, Depth_0m), color = "gray0") 

plotplot_TCO2
str(TCO2_hourly$Datetime)
TCO2_hourly$Datetime <- as.POSIXct(TCO2_hourly$Datetime, tz = "UTC")

#Making a subset that includes all the datapoints between first of november and 31 of may.
length(TCO2_hourly$Datetime)
#This is 5880
start_time <- as.POSIXct("2019-11-01 01:00", tz="UTC")
end_time   <- as.POSIXct("2020-05-31 24:00", tz="UTC")
TCO2_hourly <- subset(TCO2_hourly,
                      Datetime >= start_time & Datetime <= end_time)
length(TCO2_hourly$Datetime)


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rates_NS over time 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#inital biomass for conversions (cannot put in initial conditions)
W <- 0.05 # Chosen by origonal paper 


T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                     y = TCO2_hourly$Depth_0m,
                     rule = 2)
#-------------------------------------------------------------------------------------------------------
# Start MODEL ode 
#-------------------------------------------------------------------------------------------------------
# MODEL North Sea (region Zeeland)
#sol_NS_ZL <- ode(y= state_Clemente, t = times_NS, func = rates_NS, parms = params_NS)
#conversions to dataframes
#sol_NS_ZL <- as.data.frame(sol_NS_ZL)
#addition of a date variable
#sol_NS_ZL$Date <- seq(as_datetime("2019-11-1 01:00:00"), as_datetime("2020-05-31 24:00:00"), by="hour")
#conversion back to Celsius from Kelvin
#sol_NS_ZL$TZ_C <- TZ_K - 273.15
#create source collumn to prepare for binding all these dataframes together
#sol_NS_ZL$source <- "North Sea, just of the coast of Zeeland"
#combine all Y1 field data into one dataframe
#sol_all <- rbind(sol_NS_ZL)

# Function to run the model and create sol_all (the dataframe with the simulation)
run_model <- function(state, times, rates_func, params, TZ_K) {
  # Solve ODE system
  sol <- ode(y = state, times = times, func = rates_func, parms = params)
  
  # Convert to data frame
  sol_df <- as.data.frame(sol)
  
  # Add date sequence
  sol_df$Date <- seq(
    from = as_datetime("2019-11-01 01:00:00"),
    to   = as_datetime("2020-05-31 24:00:00"),
    by   = "hour"
  )
  
  # Convert Kelvin to Celsius
  sol_df$TZ_C <- TZ_K - 273.15
  
  # Add source column
  sol_df$source <- "North Sea, just off the coast of Zeeland"
  
  # Return the processed data frame
  return(sol_df)
}

sol_all <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K)

#-------------------------------------------------------------------------------------------------------
# Model plots 
#-------------------------------------------------------------------------------------------------------
## Irradiance plot ###
plot_I <- ggplot() + 
  geom_line(data = sol_all, aes(Date, I), color = "gray0") +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)")) +
  ggtitle("a) Irradiance") +
  theme_minimal()


## Temperature plot ##
plot_T <- ggplot(data = sol_all, aes(Date, TZ_C)) + 
  geom_line() +
  ylim(5, 12.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = "Temperature (°C)") +
  ggtitle("b) Sea surface temperature") + 
  theme(legend.position = "none") + 
  theme_minimal()

## Nitrate plot ##
plot_N <- ggplot() + 
  geom_line(data = sol_all, aes(Date, N), size = 1) +
  #geom_point(data = sol_NS_ZL, aes(Date, N)) +
  #theme_bw() +
  #theme(legend.position="none") +
  #theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  #xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  #ylim(0, 10) +
  labs(x= "Date (2019-2020)", y = bquote('N concentration (mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'L'^"-1"*')')) +
  ggtitle("c) Nitrate concentration") + 
  theme_minimal()

## DIC plot ##
plot_DIC <- ggplot() + 
    geom_line(data = sol_all, aes(Date, C), size = 1) +
    #geom_point(data = sol_NS_ZL, aes(Date, N)) +
    #theme_bw() +
    #theme(legend.position="none") +
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
    #xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
    #ylim(0, 10) +
    labs(x= "Date (2019-2020)", y = bquote('DIC concentration (mol'~ 'L'^"-1"*')')) +
    ggtitle("d) DIC concentration") + 
    theme_minimal()

grid.arrange(plot_I, plot_T, plot_N, plot_DIC, ncol=2) #gridded plot


## Structure, C and N reserves ##
(plot_reserves <- ggplot() +
  geom_line(data = sol_all, aes(Date, M_V,  color = "Mol structure"),  size = 1) +
  geom_line(data = sol_all, aes(Date, m_EC, color = "Mol C reserve"), size = 1) +
  geom_line(data = sol_all, aes(Date, m_EN, color = "Mol N reserve"), size = 1) +
  scale_color_manual(values = c("Mol structure" = "blue",
                                "Mol C reserve" = "red",
                                "Mol N reserve" = "green"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Mol/L')) +
  ggtitle("Structure and reserve density (of C and N)") +
  theme_minimal())

## Whole blade dry weight in grams ##
plot_mass <- ggplot() +
  geom_line(data = sol_all, aes(x = Date, y = W), color = "orange", size = 1) +
  labs(x = "Date (2019-2020)",
       y = bquote('Blade dry weight (gram per blade)')) +
  ggtitle("Whole blade dry weight") +
  theme_minimal()
## Blade length in cm ##
# This is the allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
plot_length <- ggplot() +
  geom_line(data = sol_all, aes(x = Date, y = L_allometric), color = "darkgreen", size = 1) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length") +
  theme_minimal()
grid.arrange(plot_mass, plot_length, ncol=2) #grided plot


#Bar plot 
sol_last <- sol_all %>% filter(Date == max(Date, na.rm = TRUE)) #take the last value of the df 
barplot_mass <- ggplot() + # plot last value of df 
  geom_col(data = sol_last, aes(x = factor(Date), y = W), fill = "orange", width = 0.3) +
  labs(x = "Date (2019-2020)",
       y = bquote('Blade dry weight (gram per blade)')) +
  ggtitle("Whole blade dry weight") +
  theme_minimal()
## Blade length in cm ##
# This is the allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
barplot_length <- ggplot() +
  geom_col(data = sol_last, aes(x = factor(Date), y = L_allometric), fill = "darkgreen", width = 0.3) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length") +
  theme_minimal()
grid.arrange(barplot_mass, barplot_length, ncol=2) #grided plot

#Growth rate plot 
(plot_growthrate <- ggplot() +
  geom_line(data = sol_all, aes(x = Date, y = r*24), color = "darkgreen", size = 1) +
  labs(x = "Date (2019-2020)",
       y = bquote('per day')) +
  ggtitle("Growth rate") +
  theme_minimal())


#-------------------------------------------------------------------------------------------------------
# Calulations 
#-------------------------------------------------------------------------------------------------------
#Mass structure created in Mol 
mass_in_mol <- tail(sol_all$M_V, 1) - sol_all$M_V[1]
mass_created <- mass_in_mol* w_V #w_V = molecular weight of structure(g/mol)
print(glue("The total growth in blade structure is {mass_in_mol} Mol and {mass_created} grams"))

# Dry body weight (of the blades) created in grams 
weight = tail(sol_all$W, 1) - sol_all$W[1]
print(glue("The total growth in Dry body weight (of the blades) is {weight} grams "))

# Length grown (of the blades) in cm 
blade_growth <- tail(sol_all$L_allometric, 1) - sol_all$L_allometric[1]
print(glue("The total growth in blade length is {blade_growth} cm"))


#-------------------------------------------------------------------------------------------------------
# Depth analysis 
#-------------------------------------------------------------------------------------------------------
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)


## Sea surface ##  
depth_irradiance <- Irradiance$PAR_1m
I_field <- approxfun(x = seq(0, length(depth_irradiance) - 1),
                     y = depth_irradiance,
                     rule = 2)
depth_CO_2 <- TCO2_hourly$Depth_0m
CO_2 <- approxfun(x = seq(0, length(depth_CO_2) - 1),
                  y = depth_CO_2,
                  rule = 2)
sea_surface_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K)

## Medium depth ##
depth_irradiance <- Irradiance$PAR_4.5m
I_field <- approxfun(x = seq(0, length(depth_irradiance) - 1),
                     y = depth_irradiance,
                     rule = 2)
depth_CO_2 <- TCO2_hourly$Depth_5m
CO_2 <- approxfun(x = seq(0, length(depth_CO_2) - 1),
                  y = depth_CO_2,
                  rule = 2)
medium_depth_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K)

## Large depth ## 
depth_irradiance <- Irradiance$PAR_7m
I_field <- approxfun(x = seq(0, length(depth_irradiance) - 1),
                     y = depth_irradiance,
                     rule = 2)
depth_CO_2 <- TCO2_hourly$Depth_10m
CO_2 <- approxfun(x = seq(0, length(depth_CO_2) - 1),
                  y = depth_CO_2,
                  rule = 2)
large_depth_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


## Irradiance plot ###
plot_depth_I <- ggplot() + 
  geom_line(data = sea_surface_df, aes(Date, I, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, I, color = "Medium depth"), size = 1) +
  geom_line(data = large_depth_df, aes(Date, I, color = "Large depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Large depth" = "darkblue"),
                     name = "Variable") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  #xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1"*')')) +
  theme_minimal() + 
  theme(legend.position = "none") +
  ggtitle("a) Irradiance at varing depths")

## DIC plot ##
plot_depth_DIC <- ggplot() + 
    geom_line(data = sea_surface_df, aes(Date, C, color = "Sea surface"), size = 1) +
    geom_line(data = medium_depth_df, aes(Date, C, color = "Medium depth"), size = 1) +
    geom_line(data = large_depth_df, aes(Date, C, color = "Large depth"), size = 1) +
    scale_color_manual(values = c("Sea surface" = "lightblue",
                                  "Medium depth" = "blue",
                                  "Large depth" = "darkblue"),
                       name = "Variable") +
    #xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
    labs(x= "Date (2019-2020)", y = bquote('DIC concentration (mol L'^"-1"*')')) +
    theme_minimal() + 
    theme(legend.position = "none") +
    ggtitle("b) DIC concentration at varing depths")

grid.arrange(plot_depth_I, plot_depth_DIC, ncol=2) #gridded plot



## Whole blade dry weight in grams ##
plot_depth_mass <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, W, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, W, color = "Medium depth"), size = 1) +
  geom_line(data = large_depth_df, aes(Date, W, color = "Large depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Large depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Blade dry weight \n(g per blade)')) +
  ggtitle("c) Whole blade dry weight") +
  theme_minimal() + 
  theme(legend.position = "none") 

## Blade length in cm ##
# This is the allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
plot_depth_length <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, L_allometric, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, L_allometric, color = "Medium depth"), size = 1) +
  geom_line(data = large_depth_df, aes(Date, L_allometric, color = "Large depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Large depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length \n(cm per blade)')) +
  ggtitle("d) Blade length") +
  theme_minimal() + 
  theme(legend.position = "none") 

grid.arrange(plot_depth_mass, plot_depth_length, ncol=2) #gridded plot

#BAR PLOT 
sea_surface_last <- sea_surface_df %>% filter(Date == max(Date, na.rm = TRUE))
medium_depth_last <- medium_depth_df %>% filter(Date == max(Date, na.rm = TRUE))
large_depth_last <- large_depth_df %>% filter(Date == max(Date, na.rm = TRUE))

# Add 'Depth' label to each
sea_surface_last$Depth <- "Sea surface"
medium_depth_last$Depth <- "Medium depth"
large_depth_last$Depth <- "Large depth"

# Combine into one data frame
df_last <- bind_rows(sea_surface_last, medium_depth_last, large_depth_last)

plot_mass_bar <- ggplot(df_last, aes(x = Depth, y = W, fill = Depth)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("Sea surface" = "lightblue",
                               "Medium depth" = "blue",
                               "Large depth" = "darkblue")) +
  labs(x = "Depth",
       y = bquote('Blade dry weight \n(g per blade)')) +
  ggtitle("e) Final blade dry weight") +
  theme_minimal() + 
  theme(legend.position = "none") 

plot_length_bar <- ggplot(df_last, aes(x = Depth, y = L_allometric, fill = Depth)) +
  geom_col(width = 0.5) +
  scale_fill_manual(values = c("Sea surface" = "lightblue",
                               "Medium depth" = "blue",
                               "Large depth" = "darkblue")) +
  labs(x = "Depth",
       y = bquote('Blade length \n(cm per blade)')) +
  ggtitle("f) Final blade length") +
  theme_minimal() + 
  theme(legend.position = "none") 
grid.arrange(plot_depth_I, plot_depth_DIC, plot_depth_mass, plot_depth_length, plot_mass_bar, plot_length_bar, ncol = 2)

#-------------------------------------------------------------------------------------------------------
#Sensitivity analysis 
#-------------------------------------------------------------------------------------------------------
#Before starting the Sensitivity analysis, make sure the baseline stats are set to sea surface depth
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                  y = TCO2_hourly$Depth_0m,
                  rule = 2)
#----------------------------------
## Temperature ##
#----------------------------------
temp_minus_20 <- temp$TZ_K * 0.8 
T_field <- approxfun(x = seq(0, length(temp_minus_20) - 1), y = temp_minus_20, rule = 2) 
# Run model again 
temp_minus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp_minus_20) 

temp_minus_10 <- temp$TZ_K * 0.9 
T_field <- approxfun(x = seq(0, length(temp_minus_10) - 1), y = temp_minus_10, rule = 2)
# Run model again 
temp_minus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp_minus_10) 

temp_nul <- temp$TZ_K
T_field <- approxfun(x = seq(0, length(temp_nul) - 1), y = temp_nul, rule = 2)
# Run model again 
temp_nul_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp_nul) 

temp_plus_10 <- temp$TZ_K * 1.1 
T_field <- approxfun(x = seq(0, length(temp_plus_10) - 1), y = temp_plus_10, rule = 2)
# Run model again 
temp_plus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp_plus_10) 

temp_plus_20 <- temp$TZ_K * 1.2 
T_field <- approxfun(x = seq(0, length(temp_plus_20) - 1), y = temp_plus_20, rule = 2)
# Run model again 
temp_plus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp_plus_20) 


plot_temp_sens <- ggplot() + 
  geom_line(data = temp_plus_20_df, aes(Date, TZ_C, color = "+ 20 %"), size = 1) +
  geom_line(data = temp_plus_10_df, aes(Date, TZ_C, color = "+ 10 %"), size = 1) +
  geom_line(data = temp_nul_df, aes(Date, TZ_C, color = "0 %"), size = 1) +
  geom_line(data = temp_minus_10_df, aes(Date, TZ_C, color = "- 10 %"), size = 1) +
  geom_line(data = temp_minus_20_df, aes(Date, TZ_C, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Temperature in gedrees Celcius')) +
  ggtitle("Temperature sensitivity analysis")

plot_biomass_temp_sens <- ggplot() + 
  geom_line(data = temp_plus_20_df, aes(Date, W, color = "+ 20 %"), size = 1) +
  geom_line(data = temp_plus_10_df, aes(Date, W, color = "+ 10 %"), size = 1) +
  geom_line(data = temp_nul_df, aes(Date, W, color = "0 %"), size = 1) +
  geom_line(data = temp_minus_10_df, aes(Date, W, color = "- 10 %"), size = 1) +
  geom_line(data = temp_minus_20_df, aes(Date, W, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
                       "+ 20 %" = "red",
                       "+ 10 %" = "orange",
                       "0 %"     = "black",
                       "- 10 %" = "purple",
                       "- 20 %" = "lightblue"),
                     breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Blade dry weight in grams per blade')) +
  ggtitle("Biomass growth in varing temperatures")

plot_length_temp_sens <- ggplot() + 
  geom_line(data = temp_plus_20_df, aes(Date, L_allometric, color = "+ 20 %"), size = 1) +
  geom_line(data = temp_plus_10_df, aes(Date, L_allometric, color = "+ 10 %"), size = 1) +
  geom_line(data = temp_nul_df, aes(Date, L_allometric, color = "0 %"), size = 1) +
  geom_line(data = temp_minus_10_df, aes(Date, L_allometric, color = "- 10 %"), size = 1) +
  geom_line(data = temp_minus_20_df, aes(Date, L_allometric, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length growth in varing temperatures")

grid.arrange(plot_temp_sens, plot_length_temp_sens, plot_biomass_temp_sens, ncol=3) #gridded plot



#--------------------------------------------------------------------
# Nitrate sensitivity 
#--------------------------------------------------------------------
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                  y = TCO2_hourly$Depth_0m,
                  rule = 2)


N_minus_20 <- nitrate_hourly$no3 * 0.8
N_field <- approxfun(x = seq(0, length(N_minus_20) - 1), y = N_minus_20, rule = 2) 
# Run model again 
N_minus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

N_minus_10 <- nitrate_hourly$no3 * 0.9 
N_field <- approxfun(x = seq(0, length(N_minus_10) - 1), y = N_minus_10, rule = 2) 
# Run model again 
N_minus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

N_nul <- nitrate_hourly$no3 
N_field <- approxfun(x = seq(0, length(N_nul) - 1), y = N_nul,rule = 2) 
# Run model again 
N_nul_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


N_plus_10 <- nitrate_hourly$no3 * 1.1 
N_field <- approxfun(x = seq(0, length(N_plus_10) - 1), y = N_plus_10,rule = 2) 
# Run model again 
N_plus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

N_plus_20 <- nitrate_hourly$no3 * 1.2 
N_field <- approxfun(x = seq(0, length(N_plus_20) - 1), y = N_plus_20,rule = 2) 
# Run model again 
N_plus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K)  

plot_N_sens <- ggplot() + 
  geom_line(data = N_plus_20_df, aes(Date, N, color = "+ 20 %"), size = 1) +
  geom_line(data = N_plus_10_df, aes(Date, N, color = "+ 10 %"), size = 0.5) +
  geom_line(data = N_nul_df, aes(Date, N, color = "0 %"), size = 0.5) +
  geom_line(data = N_minus_10_df, aes(Date, N, color = "- 10 %"), size = 1) +
  geom_line(data = N_minus_20_df, aes(Date, N, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('N concentration mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'L'^"-1")) +
  ggtitle("Nitrate sensitivity analysis")


plot_biomass_N_sens <- ggplot() + 
  geom_line(data = N_plus_20_df, aes(Date, W, color = "+ 20 %"), size = 1) +
  geom_line(data = N_plus_10_df, aes(Date, W, color = "+ 10 %"), size = 0.5) +
  geom_line(data = N_nul_df, aes(Date, W, color = "0 %"), size = 0.5) +
  geom_line(data = N_minus_10_df, aes(Date, W, color = "- 10 %"), size = 1) +
  geom_line(data = N_minus_20_df, aes(Date, W, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Blade dry weight in grams per blade')) +
  ggtitle("Biomass growth per blade")

plot_length_N_sens <- ggplot() + 
  geom_line(data = N_plus_20_df, aes(Date, L_allometric, color = "+ 20 %"), size = 1) +
  geom_line(data = N_plus_10_df, aes(Date, L_allometric, color = "+ 10 %"), size = 1) +
  geom_line(data = N_nul_df, aes(Date, L_allometric, color = "0 %"), size = 0.5) +
  geom_line(data = N_minus_10_df, aes(Date, L_allometric, color = "- 10 %"), size = 1) +
  geom_line(data = N_minus_20_df, aes(Date, L_allometric, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length growth") 


grid.arrange(plot_N_sens, plot_length_N_sens, plot_biomass_N_sens, ncol=3) #gridded plot


#--------------------------------------------------------------------
# Irradiance sensitivity 
#--------------------------------------------------------------------
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                  y = TCO2_hourly$Depth_0m,
                  rule = 2)

# The change is in Par, not in Iraadiance. The conversion formula is: 
# Irradiance$PAR_1m = Irradiance$NSW_estimate * PARfrac * C * exp(-k * z) * 3600 / 1000000
# It does not matter if you multiply Irradiance$NSW_estimate or Irradiance$PAR_1m by the multiplyer. Result is the same. 

irradiance_minus_20 <- Irradiance$PAR_1m * 0.8 
I_field <- approxfun(x = seq(0, length(irradiance_minus_20) - 1), y = irradiance_minus_20, rule = 2)
irradiance_minus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

irradiance_minus_10 <- Irradiance$PAR_1m * 0.9 
I_field <- approxfun(x = seq(0, length(irradiance_minus_10) - 1),y = irradiance_minus_10,rule = 2)
irradiance_minus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

irradiance_nul <- Irradiance$PAR_1m  
I_field <- approxfun(x = seq(0, length(irradiance_nul) - 1), y = irradiance_nul,rule = 2)
irradiance_nul_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

irradiance_plus_10 <- Irradiance$PAR_1m * 1.1
I_field <- approxfun(x = seq(0, length(irradiance_plus_10) - 1), y = irradiance_plus_10,rule = 2)
irradiance_plus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

irradiance_plus_20 <- Irradiance$PAR_1m * 1.2 
I_field <- approxfun(x = seq(0, length(irradiance_plus_20) - 1), y = irradiance_plus_20,rule = 2)
irradiance_plus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

plot_I_sens <- ggplot() + 
  geom_line(data = irradiance_plus_20_df, aes(Date, I, color = "+ 20 %"), size = 1) +
  geom_line(data = irradiance_plus_10_df, aes(Date, I, color = "+ 10 %"), size = 1) +
  geom_line(data = irradiance_nul_df, aes(Date, I, color = "0 %"), size = 1) +
  geom_line(data = irradiance_minus_10_df, aes(Date, I, color = "- 10 %"), size = 1) +
  geom_line(data = irradiance_minus_20_df, aes(Date, I, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)")) +
  ggtitle("Irradiance sensitivity analysis")


plot_biomass_I_sens <- ggplot() + 
  geom_line(data = irradiance_plus_20_df, aes(Date, W, color = "+ 20 %"), size = 1) +
  geom_line(data = irradiance_plus_10_df, aes(Date, W, color = "+ 10 %"), size = 1) +
  geom_line(data = irradiance_nul_df, aes(Date, W, color = "0 %"), size = 1) +
  geom_line(data = irradiance_minus_10_df, aes(Date, W, color = "- 10 %"), size = 1) +
  geom_line(data = irradiance_minus_20_df, aes(Date, W, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Blade dry weight in grams per blade')) +
  ggtitle("Biomass growth per blade")

plot_length_I_sens <- ggplot() + 
  geom_line(data = irradiance_plus_20_df, aes(Date, L_allometric, color = "+ 20 %"), size = 1) +
  geom_line(data = irradiance_plus_10_df, aes(Date, L_allometric, color = "+ 10 %"), size = 1) +
  geom_line(data = irradiance_nul_df, aes(Date, L_allometric, color = "0 %"), size = 0.5) +
  geom_line(data = irradiance_minus_10_df, aes(Date, L_allometric, color = "- 10 %"), size = 1) +
  geom_line(data = irradiance_minus_20_df, aes(Date, L_allometric, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length growth") 


grid.arrange(plot_I_sens, plot_length_I_sens, plot_biomass_I_sens, ncol=3) #gridded plot



#--------------------------------------------------------------------
# DIC sensitivity 
#--------------------------------------------------------------------
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                  y = TCO2_hourly$Depth_0m,
                  rule = 2)

DIC_minus_20 <- TCO2_hourly$Depth_0m * 0.8 
CO_2 <- approxfun(x = seq(0, length(DIC_minus_20) - 1), y = DIC_minus_20, rule = 2)
DIC_minus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

DIC_minus_10 <- TCO2_hourly$Depth_0m * 0.9 
CO_2 <- approxfun(x = seq(0, length(DIC_minus_10) - 1), y = DIC_minus_10, rule = 2)
DIC_minus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

DIC_nul <- TCO2_hourly$Depth_0m 
CO_2 <- approxfun(x = seq(0, length(DIC_nul) - 1), y = DIC_nul, rule = 2)
DIC_nul_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

DIC_plus_10 <- TCO2_hourly$Depth_0m * 1.1
CO_2 <- approxfun(x = seq(0, length(DIC_plus_10) - 1), y = DIC_plus_10, rule = 2)
DIC_plus_10_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

DIC_plus_20 <- TCO2_hourly$Depth_0m * 1.1
CO_2 <- approxfun(x = seq(0, length(DIC_plus_20) - 1), y = DIC_plus_20, rule = 2)
DIC_plus_20_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


plot_DIC_sens <- ggplot() + 
  geom_line(data = DIC_plus_20_df, aes(Date, C, color = "+ 20 %"), size = 1) +
  geom_line(data = DIC_plus_10_df, aes(Date, C, color = "+ 10 %"), size = 1) +
  geom_line(data = DIC_nul_df, aes(Date, C, color = "0 %"), size = 1) +
  geom_line(data = DIC_minus_10_df, aes(Date, C, color = "- 10 %"), size = 1) +
  geom_line(data = DIC_minus_20_df, aes(Date, C, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('DIC concentration mol'~ 'L'^"-1")) +
  ggtitle("DIC sensitivity analysis")


plot_biomass_DIC_sens <- ggplot() + 
  geom_line(data = DIC_plus_20_df, aes(Date, W, color = "+ 20 %"), size = 1) +
  geom_line(data = DIC_plus_10_df, aes(Date, W, color = "+ 10 %"), size = 1) +
  geom_line(data = DIC_nul_df, aes(Date, W, color = "0 %"), size = 1) +
  geom_line(data = DIC_minus_10_df, aes(Date, W, color = "- 10 %"), size = 1) +
  geom_line(data = DIC_minus_20_df, aes(Date, W, color = "- 20 %"), size = 1) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Blade dry weight in grams per blade')) +
  ggtitle("Biomass growth per blade")

plot_length_DIC_sens <- ggplot() + 
  geom_line(data = DIC_plus_20_df, aes(Date, L_allometric, color = "+ 20 %"), size = 1) +
  geom_line(data = DIC_plus_10_df, aes(Date, L_allometric, color = "+ 10 %"), size = 0.5) +
  geom_line(data = DIC_nul_df, aes(Date, L_allometric, color = "0 %"), size = 0.5) +
  geom_line(data = DIC_minus_10_df, aes(Date, L_allometric, color = "- 10 %"), size = 0.5, linetype = 2) +
  geom_line(data = DIC_minus_20_df, aes(Date, L_allometric, color = "- 20 %"), size = 0.5) +
  scale_color_manual(values = c(
    "+ 20 %" = "red",
    "+ 10 %" = "orange",
    "0 %"     = "black",
    "- 10 %" = "purple",
    "- 20 %" = "lightblue"),
    breaks = c("+ 20 %", "+ 10 %", "0 %", "- 10 %", "- 20 %")) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length growth") 


grid.arrange(plot_DIC_sens, plot_length_DIC_sens, plot_biomass_DIC_sens, ncol=3) #gridded plot


#----------------------------------------------------------------------------------------
# Plot sensitivity analysis 
#----------------------------------------------------------------------------------------
generate_sensitivity_df <- function(base_vector,       # e.g. TCO2_hourly$Depth_0m
                                    param_name,        # e.g. "DIC concentration"
                                    variable_name,     # e.g. "CO_2" or "I_field"
                                    multiplier_seq = seq(0.2, 1.8, by = 0.1),
                                    var_prefix = "env") {
  results_list <- list()
  
  for (m in multiplier_seq) {
    # Create percentage label for plotting
    percent_change <- round((m - 1) * 100)
    label <- ifelse(percent_change > 0,
                    paste0("+ ", percent_change),
                    ifelse(percent_change == 0, "0", paste0("- ", abs(percent_change))))
    
    # Create the scaled input vector
    scaled_input <- base_vector * m
   
     # Reset baseline variables before starting analysis
    T_field <<- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
    N_field <<- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1), y = nitrate_hourly$no3, rule = 2)
    I_field <<- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1), y = Irradiance$PAR_1m, rule = 2)
    CO_2    <<- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1), y = TCO2_hourly$Depth_0m, rule = 2)
    
    # Dynamically assign the environmental driver function
    assign(variable_name,
           eval(parse(text = paste0("approxfun(x = seq(0, length(scaled_input) - 1), y = scaled_input, rule = 2)"))),
           envir = .GlobalEnv)
    
    # Run the model
    df_result <- run_model(state = state_Clemente,
                           times = times_NS,
                           rates_func = rates_NS,
                           params = params_NS,
                           TZ_K = temp$TZ_K)
    
    # Get the last timestep
    last_row <- df_result %>%
      filter(Date == max(Date, na.rm = TRUE)) %>%
      select(W, L_allometric) %>%
      mutate(percentage = label,
             env = param_name)
    
    # Store
    results_list[[label]] <- last_row
  }
  
  # Combine all into one data.frame
  bind_rows(results_list)
}


sensitivity_analysis_DIC <- generate_sensitivity_df(
  base_vector = TCO2_hourly$Depth_0m,
  param_name = "DIC concentration",
  variable_name = "CO_2",
)


sensitivity_analysis_I <- generate_sensitivity_df(
  base_vector = Irradiance$PAR_1m,
  param_name = "Irradiance",
  variable_name = "I_field",
)

sensitivity_analysis_N <- generate_sensitivity_df(
  base_vector = nitrate_hourly$no3,
  param_name = "N concentration",
  variable_name = "N_field",
)

sensitivity_analysis_T <- generate_sensitivity_df(
  base_vector = temp$TZ_K,
  param_name = "Temperature",
  variable_name = "T_field",
)


# Combine into one data frame
sensitivity_analysis <- bind_rows(
  sensitivity_analysis_DIC,
  sensitivity_analysis_I,
  sensitivity_analysis_N,
  sensitivity_analysis_T
)


plot_sensitivity_mass_line <- ggplot(
  sensitivity_analysis,
  aes(x = factor(percentage, levels = c("- 80", "- 70", "- 60", "- 50", "- 40", "- 30", "- 20", "- 10", "0", "+ 10", "+ 20", "+ 30", "+ 40", "+ 50", "+ 60", "+ 70", "+ 80")),
      y = W,
      color = env,
      group = env)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Temperature" = "darkred",
    "N concentration" = "black",
    "Irradiance" = "blue",
    "DIC concentration" = 'darkgreen'
  )) +
  labs(
    x = "Change (%)",
    y = bquote('Blade dry weight \n(g per blade)'),
    color = "Environmental factor"
  ) +
  ggtitle("Final blade dry weight") +
  theme_minimal()


plot_sensitivity_length_line <- ggplot(
  sensitivity_analysis,
  aes(x = factor(percentage, levels = c("- 80", "- 70", "- 60", "- 50", "- 40", "- 30", "- 20", "- 10", "0", "+ 10", "+ 20", "+ 30", "+ 40", "+ 50", "+ 60", "+ 70", "+ 80")),
      y = L_allometric,
      color = env,
      group = env)
) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Temperature" = "darkred",
    "N concentration" = "black",
    "Irradiance" = "blue",
    "DIC concentration" = 'darkgreen'
  )) +
  labs(
    x = "Change (%)",
    y = bquote('Blade length \n(cm per blade)'),
    color = "Environmental factor"
  ) +
  ggtitle("Final blade length") +
  theme_minimal()

# Combine plots side-by-side
grid.arrange(plot_sensitivity_mass_line, plot_sensitivity_length_line, ncol = 2)



#---------------------------------------------------------------------------------------------------
# Climate change scenario's 
#---------------------------------------------------------------------------------------------------
#Reset orginal values 
T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(Irradiance$PAR_1m) - 1),
                     y = Irradiance$PAR_1m,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(TCO2_hourly$Depth_0m) - 1),
                  y = TCO2_hourly$Depth_0m,
                  rule = 2)

# To add Irradiance (in W/m²) to Irradiance$NSW_estimate. We need to update the PAR formula. 
# formula = Irradiance$NSW_estimate * PARfrac * C * exp(-k * z) * 3600 / 1000000
# It does not matter if you multiply Irradiance$NSW_estimate or Irradiance$PAR_1m by the multiplyer. Result is the same. 
#Adding to Kelvin is the same ass adding to degrees celcius, so no conversion is needed.

# Current situation 
current_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 

avg_IR2019 = mean(Irradiance$NSW_estimate, na.rm = TRUE) # 20/11/2019 19:00:00 till 23/11/2019 23:00:00 there is missing data # W/m^2 
avg_temp2019 = mean(temp$TZ_K)  # Kelvin



# RCP 2.6 Irradiance + 0.3, Temp + 0.5 
avg_IR2100 = avg_IR2019 + 0.3 
avg_temp2100 = avg_temp2019 + 0.5
irradiance_2_6 <- Irradiance$PAR_1m * (avg_IR2100 / avg_IR2019)
temp_2_6 <- temp$TZ_K * (avg_temp2100 / avg_temp2019)
#irradiance_2_6 <-  (ifelse(Irradiance$NSW_estimate > 0, Irradiance$NSW_estimate + 0.3, Irradiance$NSW_estimate)) * PARfrac * C * exp(-k * z) * 3600 / 1000000
#temp_2_6 <- temp$TZ_K + 0.5 
I_field <- approxfun(x = seq(0, length(irradiance_2_6) - 1), y = irradiance_2_6, rule = 2)
T_field <- approxfun(x = seq(0, length(temp_2_6) - 1), y = temp_2_6, rule = 2)
RCP2.6_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


# RCP 4.5 Irradiance + 2.2, Temp + 1 
avg_IR2100 = avg_IR2019 + 2.2 
avg_temp2100 = avg_temp2019 + 1
temp_4_5 <- temp$TZ_K * (avg_temp2100 / avg_temp2019)
irradiance_4_5 <- Irradiance$PAR_1m * (avg_IR2100 / avg_IR2019)
#irradiance_4_5 <- (ifelse(Irradiance$NSW_estimate > 0, Irradiance$NSW_estimate + 2.2, Irradiance$NSW_estimate)) * PARfrac * C * exp(-k * z) * 3600 / 1000000
#temp_4_5 <- temp$TZ_K + 1 # plus 1 degree temp
I_field <- approxfun(x = seq(0, length(irradiance_4_5) - 1), y = irradiance_4_5, rule = 2)
T_field <- approxfun(x = seq(0, length(temp_4_5) - 1), y = temp_4_5, rule = 2)
RCP4.5_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


# RCP 6 Irradiance + 3.7 , Temp + 1.5 
avg_IR2100 = avg_IR2019 + 3.7  
avg_temp2100 = avg_temp2019 + 1.5
temp_6 <- temp$TZ_K * (avg_temp2100 / avg_temp2019)
irradiance_6 <- Irradiance$PAR_1m * (avg_IR2100 / avg_IR2019)
#irradiance_6 <- (ifelse(Irradiance$NSW_estimate > 0, Irradiance$NSW_estimate + 3.7, Irradiance$NSW_estimate)) * PARfrac * C * exp(-k * z) * 3600 / 1000000
#temp_6 <- temp$TZ_K + 1.5 # plus 1.5 degree temp
I_field <- approxfun(x = seq(0, length(irradiance_6) - 1), y = irradiance_6, rule = 2)
T_field <- approxfun(x = seq(0, length(temp_6) - 1), y = temp_6, rule = 2)
RCP6_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 


# RCP 8.5 Irradiance + 6.2, Temp + 3.4
avg_IR2100 = avg_IR2019 + 6.2  
avg_temp2100 = avg_temp2019 + 3.4
irradiance_8_5 <- Irradiance$PAR_1m * (avg_IR2100 / avg_IR2019)
temp_8_5 <- temp$TZ_K + (avg_temp2100 / avg_temp2019)
#irradiance_8_5 <- (ifelse(Irradiance$NSW_estimate > 0, Irradiance$NSW_estimate + 6.2, Irradiance$NSW_estimate)) * PARfrac * C * exp(-k * z) * 3600 / 1000000
#temp_8_5 <- temp$TZ_K + 3.4 # plus 3.4 degree temp
I_field <- approxfun(x = seq(0, length(irradiance_8_5) - 1), y = irradiance_8_5, rule = 2)
T_field <- approxfun(x = seq(0, length(temp_8_5) - 1), y = temp_8_5, rule = 2)
RCP8.5_df <- run_model(state = state_Clemente, times = times_NS, rates_func = rates_NS, params = params_NS, TZ_K = temp$TZ_K) 



plot_RCP_biomass_sens <- ggplot() + 
  geom_line(data = current_df, aes(Date, W, color = "Now"), size = 1) +
  geom_line(data = RCP2.6_df, aes(Date, W, color = "RCP 2.6"), size = 1) +
  geom_line(data = RCP4.5_df, aes(Date, W, color = "RCP 4.5"), size = 1) +
  geom_line(data = RCP6_df, aes(Date, W, color = "RCP 6"), size = 1) +
  geom_line(data = RCP8.5_df, aes(Date, W, color = "RCP 8.5"), size = 1) +
  scale_color_manual(values = c(
    "Now" = "grey",
    "RCP 2.6" = "black",
    "RCP 4.5" = "darkblue",
    "RCP 6"     = "darkgreen",
    "RCP 8.5" = "red"))+
    #breaks = c("RCP 2.6", "RCP 4.5", "RCP 6", "RCP 8.5")) +
  #xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  labs(x= "Date (2019-2020)", y = bquote('Blade dry weight in grams per blade')) +
  ggtitle("Biomass growth per blade") +
  theme_minimal()

plot_RCP_length_sens <- ggplot() + 
  geom_line(data = current_df, aes(Date, L_allometric, color = "Now"), size = 1) +
  geom_line(data = RCP2.6_df, aes(Date, L_allometric, color = "RCP 2.6"), size = 1) +
  geom_line(data = RCP4.5_df, aes(Date, L_allometric, color = "RCP 4.5"), size = 1) +
  geom_line(data = RCP6_df, aes(Date, L_allometric, color = "RCP 6"), size = 1) +
  geom_line(data = RCP8.5_df, aes(Date, L_allometric, color = "RCP 8.5"), size = 1) +
  scale_color_manual(values = c(
    "Now" = "grey",
    "RCP 2.6" = "black",
    "RCP 4.5" = "darkblue",
    "RCP 6"     = "darkgreen",
    "RCP 8.5" = "red"))+
  #breaks = c("RCP 2.6", "RCP 4.5", "RCP 6", "RCP 8.5")) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm per blade)')) +
  ggtitle("Blade length growth") + 
  theme_minimal()

grid.arrange(plot_RCP_biomass_sens, plot_RCP_length_sens, ncol = 2)



# BAR PLOT 

current_last <- current_df %>% filter(Date == max(Date, na.rm = TRUE))
RCP2.6_last <- RCP2.6_df %>% filter(Date == max(Date, na.rm = TRUE))
RCP4.5_last <- RCP4.5_df %>% filter(Date == max(Date, na.rm = TRUE))
RCP6_last <- RCP6_df %>% filter(Date == max(Date, na.rm = TRUE))
RCP8.5_last <- RCP8.5_df %>% filter(Date == max(Date, na.rm = TRUE))

# Add 'Depth' label to each
current_last$RCP <- "2019"
RCP2.6_last$RCP <- "RCP 2.6"
RCP4.5_last$RCP <- "RCP 4.5"
RCP6_last$RCP <- "RCP 6"
RCP8.5_last$RCP <- "RCP 8.5"

# Combine into one data frame
RCPdf_last <- bind_rows(current_last, RCP2.6_last, RCP4.5_last, RCP6_last, RCP8.5_last)

RCPplot_mass_bar <- ggplot(RCPdf_last, aes(x = RCP, y = W)) +
  geom_col(width = 0.5) +
  labs(x = "RCP",
       y = bquote('Blade dry weight \n(g per blade)')) +
  ggtitle("Final blade dry weight per RCP scenario") +
  theme_minimal() 
RCPplot_length_bar <- ggplot(RCPdf_last, aes(x = RCP, y = L_allometric)) +
  geom_col(width = 0.5) +
  labs(x = "RCP scenario",
       y = bquote('Blade length \n(cm per blade)')) +
  ggtitle("Final blade length per RCP scenario") +
  theme_minimal() 
grid.arrange(RCPplot_mass_bar, RCPplot_length_bar, ncol = 2)

