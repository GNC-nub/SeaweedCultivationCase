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
#state_Lo <- c(m_EC = 0.002, #0.1, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
#              m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
#              M_V = 0.05/(w_V+0.01*w_EN+0.002*w_EC)) #molM_V #initial mass of structure

state_LoY2 <- c(m_EC = 0.01, #0.9 #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                m_EN = 0.09, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                M_V = 0.05/(w_V+0.09*w_EN+0.01*w_EC)) #molM_V #initial mass of structure

state_Johansson <- c(m_EC = 0.3, #mol C/molM_V  #Reserve density of C reserve (initial mass of C reserve per intital mass of structure)
                     m_EN = 0.01, #mol N/molM_V #Reserve density of N reserve (initial mass of N reserve per intital mass of structure)
                     M_V = 0.005/(w_V+0.01*w_EN+0.3*w_EC)) #molM_V #initial mass of structure

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
TCO2_monthly <- read.csv("DICNovDecJan_Depths0_5_10.csv")
density_seawater <- 1.026 #kilo/liter
TCO2_monthly$Depth_0m <- TCO2_monthly$Depth_0m/1000000 * density_seawater
TCO2_monthly$Depth_5m <- TCO2_monthly$Depth_5m/1000000 * density_seawater
TCO2_monthly$Depth_10m <- TCO2_monthly$Depth_10m/1000000 * density_seawater
# Monthly timestamps (start of each month)
month_dates <- as.POSIXct(c(
  "2019-11-01 00:00",
  "2019-12-01 00:00",
  "2020-01-01 00:00",
  "2020-02-01 00:00",
  "2020-03-01 00:00",
  "2020-04-01 00:00",
  "2020-05-01 00:00"
), tz = "UTC")

# Hourly sequence from 1 Nov 2019 to 31 May 2020 (inclusive)
hourly_time <- seq(
  from = as.POSIXct("2019-11-01 00:00", tz="UTC"),
  to   = as.POSIXct("2020-05-31 23:00", tz="UTC"),
  by   = "1 hour"
)

length(hourly_time)
# should be 5112

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Rates_NS over time 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#inital biomass for conversions (cannot put in initial conditions)
W <- 0.05 # Chosen by origonal paper 

# Choose the right Irradiance depth here 
depth_irradiance <- Irradiance$PAR_1m # 1m
depth_irradiance <- Irradiance$PAR_2m # 2m 
depth_irradiance <- Irradiance$PAR_4.5m # 4.5 m
depth_irradiance <- Irradiance$PAR_7m # 7m

# Choose the right DIC depth here 
depth_CO_2 <- TCO2_hourly$Depth_0m # 0m
depth_CO_2 <- TCO2_hourly$Depth_5m # 5m
depth_CO_2 <- TCO2_hourly$Depth_10m # 10m


T_field <- approxfun(x = seq(0, length(temp$TZ_K) - 1), y = temp$TZ_K, rule = 2)
N_field <- approxfun(x = seq(0, length(nitrate_hourly$no3) - 1),
                     y = nitrate_hourly$no3,
                     rule = 2)
I_field <- approxfun(x = seq(0, length(depth_irradiance) - 1),
                     y = depth_irradiance,
                     rule = 2)
CO_2 <- approxfun(x = seq(0, length(depth_CO_2) - 1),
                     y = depth_CO_2,
                     rule = 2)
#-------------------------------------------------------------------------------------------------------
# Start MODEL ode 
#-------------------------------------------------------------------------------------------------------
# MODEL North Sea (region Zeeland)
sol_NS_ZL <- ode(y= state_Johansson, t = times_NS, func = rates_NS, parms = params_NS)

###### Convert DeSolve solutions into data frame for broader plotting use ####
#conversions to dataframes
sol_NS_ZL <- as.data.frame(sol_NS_ZL)

#addition of a date variable
sol_NS_ZL$Date <- seq(as_datetime("2019-11-1 01:00:00"), as_datetime("2020-05-31 24:00:00"), by="hour")

#conversion back to Celsius from Kelvin
sol_NS_ZL$TZ_C <- TZ_K - 273.15

#create source collumn to prepare for binding all these dataframes together
sol_NS_ZL$source <- "North Sea, just of the coast of Zeeland"

#combine all Y1 field data into one dataframe
sol_all <- rbind(sol_NS_ZL)

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
  ggtitle("a)")

## Temperature plot ##
plot_T <- ggplot(data = sol_all, aes(Date, TZ_C, color = source)) + 
  geom_line() +
  scale_color_manual(values = c("blue", "blueviolet", "cyan", "coral", "darkgoldenrod1", "firebrick", "black")) +
  ylim(5, 12.5) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = "Temperature (°C)") +
  ggtitle("Figure Sea surface Temperature in Zeeland")

## Nitrate plot ##
(plot_N <- ggplot() + 
  geom_line(data = sol_all, aes(Date, N), size = 1) +
  #geom_point(data = sol_NS_ZL, aes(Date, N)) +
  #theme_bw() +
  #theme(legend.position="none") +
  #theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
  #ylim(0, 10) +
  labs(x= "Date (2019-2020)", y = bquote('N concentration mol' ~NO[3]^{"-"}~ 'and' ~NO[2]^{"-"}~ 'L'^"-1")) +
  ggtitle("Nitrate concentration in Zeeland"))

## CO2 plot ##
(plot_CO2 <- ggplot() + 
    geom_line(data = sol_all, aes(Date, C), size = 1) +
    #geom_point(data = sol_NS_ZL, aes(Date, N)) +
    #theme_bw() +
    #theme(legend.position="none") +
    #theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
    xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
    #ylim(0, 10) +
    labs(x= "Date (2019-2020)", y = bquote('CO2 concentration mol'~ 'L'^"-1")) +
    ggtitle("CO2 concentration in Zeeland"))

grid.arrange(plot_I, plot_T, plot_N, plot_CO2, ncol=2) #gridded plot


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
       y = bquote('Mass in Mol')) +
  ggtitle("Mass of structure and reserve density \n(of C and N) for kelp growth in the \nNorth Sea near Zeeland") +
  theme_minimal())

## Whole blade dry weight in grams ##
plot_mass <- ggplot() +
  geom_line(data = sol_all, aes(x = Date, y = W), color = "orange", size = 1) +
  labs(x = "Date (2019-2020)",
       y = bquote('Blade dry weight (g)')) +
  ggtitle("Whole *S. latissima* blade dry weight in the\nNorth Sea near Zeeland") +
  theme_minimal()


## Blade length in cm ##
# This is the allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
plot_length <- ggplot() +
  geom_line(data = sol_all, aes(x = Date, y = L_allometric), color = "darkgreen", size = 1) +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm)')) +
  ggtitle("S. latissima blade length in the \nNorth Sea near Zeeland") +
  theme_minimal()
grid.arrange(plot_mass, plot_length, ncol=2) #gridded plot


#-------------------------------------------------------------------------------------------------------
# Caluclations 
#-------------------------------------------------------------------------------------------------------
#Mass strucute created in Moll 
mass_in_mol <- tail(sol_all$M_V, 1) - sol_all$M_V[1]
mass_created <- mass_in_mol* w_V #w_V = molecular weight of structure(g/mol)
mass_created # in grams

# Dry body weight (of the blades) created in grams 
weight = tail(sol_all$W, 1) - sol_all$W[1]
weight # in grams 

# Length grown (of the blades) in cm 
blade_growth <- tail(sol_all$L_allometric, 1) - sol_all$L_allometric[1]
blade_growth #in cm


#-------------------------------------------------------------------------------------------------------
# Depth analysis 
#-------------------------------------------------------------------------------------------------------
# Options for Irradiance at 1m, 2m, 4.5 and 7m. 
# Options for CO2 at 0m, 5m and 10m 


## Sea surface ##  
#Choose option: CO2 = 0m and Irradiance = 1m (at line 250) 
sea_surface_df <- sol_all

## Medium depth ##
#Choose option: CO2 = 5m and Irradiance = 4.5m (at line 250) 
medium_depth_df <- sol_all

## Deep depth ## 
#Choose option: CO2 = 10m and Irradiance = 7m (at line 250) 
deep_depth_df <- sol_all 



## Irradiance plot ###
plot_I <- ggplot() + 
  geom_line(data = sea_surface_df, aes(Date, I, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, I, color = "Medium depth"), size = 1) +
  geom_line(data = deep_depth_df, aes(Date, I, color = "Deep depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Deep depth" = "darkblue"),
                     name = "Variable") +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)")) +
  ggtitle("a)")

## CO2 plot ##
plot_CO2 <- ggplot() + 
    geom_line(data = sea_surface_df, aes(Date, C, color = "Sea surface"), size = 1) +
    geom_line(data = medium_depth_df, aes(Date, C, color = "Medium depth"), size = 1) +
    geom_line(data = deep_depth_df, aes(Date, C, color = "Deep depth"), size = 1) +
    scale_color_manual(values = c("Sea surface" = "lightblue",
                                  "Medium depth" = "blue",
                                  "Deep depth" = "darkblue"),
                       name = "Variable") +
    xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-05-31 24:00:00"))) +
    labs(x= "Date (2019-2020)", y = bquote('CO2 concentration mol'~ 'L'^"-1")) +
    ggtitle("CO2 concentration in Zeeland")

grid.arrange(plot_I, plot_CO2, ncol=2) #gridded plot


## Structure, C and N reserves ##
plot_structure_reserves <- ggplot() +
    geom_line(data = sea_surface_df, aes(Date, M_V, color = "Sea surface"), size = 1) +
    geom_line(data = medium_depth_df, aes(Date, M_V, color = "Medium depth"), size = 1) +
    geom_line(data = deep_depth_df, aes(Date, M_V, color = "Deep depth"), size = 1) +
    scale_color_manual(values = c("Sea surface" = "lightblue",
                                  "Medium depth" = "blue",
                                  "Deep depth" = "darkblue"),
                       name = "Variable") +
    labs(x = "Date (2019-2020)",
         y = bquote('Mol structure')) +
    ggtitle("Structure in Mol") +
    theme_minimal()

plot_N_reserves <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, m_EN, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, m_EN, color = "Medium depth"), size = 1) +
  geom_line(data = deep_depth_df, aes(Date, m_EN, color = "Deep depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Deep depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Mol N reserve')) +
  ggtitle("Reserve density of Nitrate") +
  theme_minimal()

plot_C_reserves <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, m_EC, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, m_EC, color = "Medium depth"), size = 1) +
  geom_line(data = deep_depth_df, aes(Date, m_EC, color = "Deep depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Deep depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Mol C reserve')) +
  ggtitle("Reserve density of Carbon") +
  theme_minimal()

## Whole blade dry weight in grams ##
plot_mass <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, W, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, W, color = "Medium depth"), size = 1) +
  geom_line(data = deep_depth_df, aes(Date, W, color = "Deep depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Deep depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Blade dry weight')) +
  ggtitle("Whole S. latissima blade dry weight") +
  theme_minimal()

## Blade length in cm ##
# This is the allometic relationship between length (cm) and dry weight (g) from Gevaert (2001)
plot_length <- ggplot() +
  geom_line(data = sea_surface_df, aes(Date, L_allometric, color = "Sea surface"), size = 1) +
  geom_line(data = medium_depth_df, aes(Date, L_allometric, color = "Medium depth"), size = 1) +
  geom_line(data = deep_depth_df, aes(Date, L_allometric, color = "Deep depth"), size = 1) +
  scale_color_manual(values = c("Sea surface" = "lightblue",
                                "Medium depth" = "blue",
                                "Deep depth" = "darkblue"),
                     name = "Variable") +
  labs(x = "Date (2019-2020)",
       y = bquote('Physical length (cm)')) +
  ggtitle("S. latissima blade length") +
  theme_minimal()

grid.arrange(plot_structure_reserves, plot_N_reserves, plot_C_reserves, plot_mass, plot_length, ncol=2) #gridded plot


#-------------------------------------------------------------------------------------------------------
#Sensitivity analysis 
#-------------------------------------------------------------------------------------------------------






#-------------------------------------------------------------------------------------------------------
# From here on useless
#-------------------------------------------------------------------------------------------------------
#IGNORE
#read in GSO N data
#GSO_N <- read.csv("T98BayNitrate.csv", header = TRUE, fileEncoding="UTF-8-BOM") #Import water quality data
#GSO_N$Date <- mdy(GSO_N$Date) #convert dates

#ggplot() + 
#  geom_line(data = GSO_N, aes(Date, NO3NO2)) +
#  geom_line(data = Wickford_WSA, aes(Date, NitrateNitrite_uM), color ="red") +                    
#  geom_line(data = RomePt_WSA, aes(Date, NitrateNitrite_uM), color ="blue") +
#  geom_line(data = RomePt_WSA2, aes(Date, NO3NO2_µM), color ="blue") +
#  geom_line(data = Wickford_WSA2, aes(Date, NO3NO2_µM*1000000), color ="red") +  
  #scale_color_manual(values = c("gray77", "gray60", "gray60", "gray30", "gray30", "gray0", "gray0")) +
#  theme_bw() +
#  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  #theme(legend.position="none") +
#  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  #ylim(0, 10) +
#  labs(x= "Date", y = bquote('N concentration μmol' ~NO[3]^{"-"}~ 'and'~NO[2]^{"-"}~ 'L'^"-1"))



# Rejected C and R (Figure 8????) 
plot_J_EC_R_PJ <- ggplot() +
  geom_line(data = sol_all, aes(Date, J_EC_R, color = source)) +
  #geom_line(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EC_R, color = source)) +
  #geom_line(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EC_R, color = source)) +
  scale_color_grey() +
  #xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  #theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2019-2020)", y = bquote('Rejected C (mol C mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("A)")

plot_J_EN_R_PJ <- ggplot() +
  geom_line(data = sol_all, aes(Date, J_EN_R, color = source)) +
  #geom_line(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, J_EN_R, color = source)) +
  #geom_line(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, J_EN_R, color = source)) +
  scale_color_grey() +
  #xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-06-01 23:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2019-2020)", y = bquote('Rejected N (mol N mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("B)")

grid.arrange(plot_J_EC_R_PJ, plot_J_EN_R_PJ, ncol=2)


# Relaxation rate --> Figure 9
plot_T_PJ <- ggplot() + 
  geom_line(data = sol_all, aes(Date, C_T)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-06-01 00:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.title = element_blank()) +
  theme(legend.position="none") + 
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = "Temperature correction factor") +
  ggtitle("A)")
plot_I_PJ <- ggplot() + 
  geom_line(data = sol_all, aes(Date, I)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-06-01 00:00:00"))) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= "Date (2019-2020)", y = bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)")) +
  ggtitle("B)")
plot_J_I_PJ <- ggplot() + #???? look at the time (2 h 06-01 extra??) 
  geom_line(data = sol_all, aes(Date, J_I, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-06-01 02:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2019-2020)", y = bquote('Specific relaxation rate (mol γ mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("C)")
plot_J_EC_A_PJ <- ggplot() +
  geom_line(data = sol_all, aes(Date, J_EC_A, color = source)) +
  scale_color_grey() +
  xlim(as.POSIXct(c("2019-11-01 00:00:00", "2020-06-01 02:00:00"))) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") + 
  labs(x= "Date (2019-2020)", y = bquote('Specific C assimilation (mol C mol V'^"-1"*' h'^"-1"*')')) +
  ggtitle("D)")
grid.arrange(plot_T_PJ, plot_I_PJ, plot_J_I_PJ, plot_J_EC_A_PJ, ncol=2)


# ----------------------------------------------------------------------
# Now its only field data and Literature data for comparison/Calibration
#----------------------------------------------------------------------


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#END FIELD DATA, START LITERATURE DATA FOR CALIBRATION
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #9C
###### N forcing set-up##############
Nmax <- 73.1221719/1000000 #M

###### Temperature forcing set-Up #############
#9 C is 282.15 K
T_dat <- 9 #C (conversion in Nuptake function to K)
###################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_9 <- Nuptake(params_NS, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_9 <- as.data.frame(sol_EspinozaChapman1983_N_9) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Setting up the forcing functions for Espinoza and Chapman (1983) nitrogen uptake #18C
###### N forcing set-up##############
Nmax <- 76.9543147/1000000 #M

###### Temperature forcing set-Up #############
#18 C is 282.15 K
T_dat <- 18 #C (conversion in Nuptake function to K)
#################################
#Model run (the differential equation solver)
sol_EspinozaChapman1983_N_18 <- Nuptake(params_NS, T_dat, Nmax, w_EN) #function from N_uptake_Calibration.R code
sol_EspinozaChapman1983_N_18 <- as.data.frame(sol_EspinozaChapman1983_N_18) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Photosynthesis model calibration
##### Temp forcing set-up ####
T_dat <- 14 #C (maintained for entire experiment
###### I forcing set-up #####
I_max <- 3233205 #micromol photons m-2 h-1
##############
sol_Johansson2002 <- Photosynthesis(params_NS, state_Johansson, w_V, w_EN, w_EC, w_O2, T_dat, I_max) #function from Photosynthesis_Calibration.R
sol_Johansson2002 <- as.data.frame(sol_Johansson2002) #conversion to dataframe for later use
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


##### Kelp Field Data Comparison plot (Figure 7) ####
#import field data
KelpY1 <- read.csv("Year1kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY1)[2] <- "Site"
KelpY1 <- filter(KelpY1, Site != "Fox Island")
KelpY1$Date <- mdy(KelpY1$SamplingDate)
KelpY1$SiteLine <- paste(KelpY1$Site, KelpY1$Line)
KelpY1 <- filter(KelpY1, SiteLine != "Narragansett Bay N 2")
KelpY2 <- read.csv("Year2kelpdata.csv", header = TRUE, fileEncoding="UTF-8-BOM")
names(KelpY2)[2] <- "Site"
KelpY2 <- filter(KelpY2, Site != "Fox Island")
KelpY2$Date <- mdy(KelpY2$SamplingDate)
KelpY2$SiteLine <- paste(KelpY2$Site, KelpY2$Line)
KelpY2 <- filter(KelpY2, SiteLine != "Narragansett Bay N 2")

NBN1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length))
NBN1_meandat$Date <- as.POSIXct(NBN1_meandat$Date)
NBN1_meandat_sub <- NBN1_meandat[2:6,]
erNBN1 <- merge(NBN1_meandat_sub, sol_all[sol_all$source == "Narragansett Bay N 1",], all.x = TRUE)
NBN1_rmse <- rmse(erNBN1$mean_length, erNBN1$L_allometric)
NBN1_rmse <- round(NBN1_rmse, 2)()

NBN1_meandat$Date <- as.POSIXct(NBN1_meandat$Date)
NBN1 <- ggplot() + 
  geom_point(data = NBN1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(NBN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBN1_rmse)) +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay N") +
  theme(legend.position="none")

NBS1_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS1_meandat$Date <- as.POSIXct(NBS1_meandat$Date)
NBS1_meandat_sub <- NBS1_meandat[2:6,]
erNBS1 <- merge(NBS1_meandat_sub, sol_all[sol_all$source == "Narragansett Bay S 1",], all.x = TRUE)
NBS1_rmse <- rmse(erNBS1$mean_length, erNBS1$L_allometric)

NBS2_meandat <- KelpY1[KelpY1$SiteLine == "Narragansett Bay S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS2_meandat$Date <- as.POSIXct(NBS2_meandat$Date)
NBS2_meandat_sub <- NBS2_meandat[2:6,]
erNBS2 <- merge(NBS2_meandat_sub, sol_all[sol_all$source == "Narragansett Bay S 2",], all.x = TRUE)
NBS2_rmse <- rmse(erNBS2$mean_length, erNBS2$L_allometric)

NBS1_2 <- ggplot() + 
  geom_point(data = NBS1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(NBS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBS1_rmse)) +
  geom_point(data = NBS2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(NBS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Narragansett Bay S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", NBS2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay S") +
  theme(legend.position="none")

PJS1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_meandat$Date <- as.POSIXct(PJS1_meandat$Date)
PJS1_meandat_sub <- PJS1_meandat[2:6,]
erPJS1 <- merge(PJS1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 1",], all.x = TRUE)
PJS1_rmse <- rmse(erPJS1$mean_length, erPJS1$L_allometric)

PJS2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_meandat$Date <- as.POSIXct(PJS2_meandat$Date)
PJS2_meandat_sub <- PJS2_meandat[2:6,]
erPJS2 <- merge(PJS2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond S 2",], all.x = TRUE)
PJS2_rmse <- rmse(erPJS2$mean_length, erPJS2$L_allometric)

PJS1_2 <- ggplot() + 
  geom_point(data = PJS1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJS1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_rmse)) +
  geom_point(data = PJS2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(PJS2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-01 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJS2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond S") +
  theme(legend.position="none")

PJN1_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_meandat$Date <- as.POSIXct(PJN1_meandat$Date)
PJN1_meandat_sub <- PJN1_meandat[2:6,]
erPJN1 <- merge(PJN1_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_rmse <- rmse(erPJN1$mean_length, erPJN1$L_allometric)

PJN2_meandat <- KelpY1[KelpY1$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_meandat$Date <- as.POSIXct(PJN2_meandat$Date)
PJN2_meandat_sub <- PJN2_meandat[2:6,]
erPJN2 <- merge(PJN2_meandat_sub, sol_all[sol_all$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_rmse <- rmse(erPJN2$mean_length, erPJN2$L_allometric)


PJN1_2 <- ggplot() + 
  geom_point(data = PJN1_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJN1_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_rmse)) +
  geom_point(data = PJN2_meandat, aes(Date, mean_length), color ="gray50", size = 3) +
  geom_errorbar(PJN2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all[sol_all$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2018-04-24 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJN2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2017-10-30 23:00:00", "2018-04-24 23:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2017-2018)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond N") +
  theme(legend.position="none")

NBN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBN1_Y2_meandat$Date <- as.POSIXct(NBN1_Y2_meandat$Date)
NBN1_Y2_meandat_sub <- NBN1_Y2_meandat[2:5,]
erNBN1_Y2 <- merge(NBN1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], all.x = TRUE)
NBN1_Y2_rmse <- rmse(erNBN1_Y2$mean_length, erNBN1_Y2$L_allometric)

NBN1_Y2 <- ggplot() + 
  geom_point(data = NBN1_Y2_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(NBN1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBN1_Y2_rmse)) +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay N") +
  theme(legend.position="none")

NBS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS1_Y2_meandat$Date <- as.POSIXct(NBS1_Y2_meandat$Date)
NBS1_Y2_meandat_sub <- NBS1_Y2_meandat[2:4,]
erNBS1_Y2 <- merge(NBS1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], all.x = TRUE)
NBS1_Y2_rmse <- rmse(erNBS1_Y2$mean_length, erNBS1_Y2$L_allometric)

NBS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Narragansett Bay S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
NBS2_Y2_meandat$Date <- as.POSIXct(NBS2_Y2_meandat$Date)
NBS2_Y2_meandat_sub <- NBS2_Y2_meandat[2:3,]
erNBS2_Y2 <- merge(NBS2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 2",], all.x = TRUE)
NBS2_Y2_rmse <- rmse(erNBS2_Y2$mean_length, erNBS2_Y2$L_allometric)

NBS1_2_Y2 <- ggplot() + 
  geom_point(data = NBS1_Y2_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(NBS1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", NBS1_Y2_rmse)) +
  geom_point(data = NBS2_Y2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(NBS2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Narragansett Bay S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", NBS2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Narragansett Bay S") +
  theme(legend.position="none")

PJS1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS1_Y2_meandat$Date <- as.POSIXct(PJS1_Y2_meandat$Date)
PJS1_Y2_meandat_sub <- PJS1_Y2_meandat[2:6,]
erPJS1_Y2 <- merge(PJS1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], all.x = TRUE)
PJS1_Y2_rmse <- rmse(erPJS1_Y2$mean_length, erPJS1_Y2$L_allometric)

PJS2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond S 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJS2_Y2_meandat$Date <- as.POSIXct(PJS2_Y2_meandat$Date)
PJS2_Y2_meandat_sub <- PJS2_Y2_meandat[2:4,]
erPJS2_Y2 <- merge(PJS2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 2",], all.x = TRUE)
PJS2_Y2_rmse <- rmse(erPJS2_Y2$mean_length, erPJS2_Y2$L_allometric)

PJS1_2_Y2 <- ggplot() + 
  geom_point(data = PJS1_Y2_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJS1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJS1_Y2_rmse)) +
  geom_point(data = PJS2_Y2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(PJS2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond S 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJS2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) + 
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond S") +
  theme(legend.position="none")

PJN1_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 1",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN1_Y2_meandat$Date <- as.POSIXct(PJN1_Y2_meandat$Date)
PJN1_Y2_meandat_sub <- PJN1_Y2_meandat[2:6,]
erPJN1_Y2 <- merge(PJN1_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], all.x = TRUE)
PJN1_Y2_rmse <- rmse(erPJN1_Y2$mean_length, erPJN1_Y2$L_allometric)

PJN2_Y2_meandat <- KelpY2[KelpY2$SiteLine == "Point Judith Pond N 2",] %>%
  group_by(Date) %>%
  summarize(mean_length = mean(Length, na.rm = TRUE), sd_length = sd(Length, na.rm = TRUE))
PJN2_Y2_meandat$Date <- as.POSIXct(PJN2_Y2_meandat$Date)
PJN2_Y2_meandat_sub <- PJN2_Y2_meandat[2:4,]
erPJN2_Y2 <- merge(PJN2_Y2_meandat_sub, sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 2",], all.x = TRUE)
PJN2_Y2_rmse <- rmse(erPJN2_Y2$mean_length, erPJN2_Y2$L_allometric)

PJN1_2_Y2 <- ggplot() + 
  geom_point(data = PJN1_Y2_meandat, aes(Date, mean_length), shape = 'diamond', size = 3) +
  geom_errorbar(PJN1_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1) +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 1",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 250, hjust = 1), label = sprintf("RMSE: %f", PJN1_Y2_rmse)) +
  geom_point(data = PJN2_Y2_meandat, aes(Date, mean_length), color = "gray50", size = 3) +
  geom_errorbar(PJN2_Y2_meandat, mapping = aes(x = Date, ymin = mean_length-sd_length, ymax = mean_length+sd_length), width = 1, color = "gray50") +
  geom_smooth(data = sol_all_Y2[sol_all_Y2$source == "Point Judith Pond N 2",], aes(Date, L_allometric, color = source)) +
  geom_label(aes(as.POSIXct("2019-05-30 01:00:00"), 230, hjust = 1), label = sprintf("RMSE: %f", PJN2_Y2_rmse), color = "gray50") +
  ylim(0, 258) +
  xlim(as.POSIXct(c("2018-12-01 12:00:00", "2019-05-30 12:00:00"))) +
  scale_color_manual(values = c("gray0", "gray50")) +
  theme_bw() +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  labs(x= "Date (2018-2019)", y = "Blade length (cm)") +
  ggtitle("Point Judith Pond N") +
  theme(legend.position="none")

grid.arrange(NBN1, NBS1_2, PJN1_2, PJS1_2, NBN1_Y2, NBS1_2_Y2, PJN1_2_Y2, PJS1_2_Y2, ncol=4)




########
##### Literature data for comparison/Calibration ####
######## Nitrate uptake ####
#Espinoza and Chapman (1983) and Ahn et al. (1998)
EC1983_9C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_9C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")
EC1983_18C_Nuptake_StM <- read.csv("EspinozaChapman1983_Nuptake_18C_StMargaretsBay.csv", header = TRUE, fileEncoding="UTF-8-BOM")

#conversions 9C
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$ResidualNitrateConcentration
EC1983_9C_Nuptake_StM$N <- round(EC1983_9C_Nuptake_StM$N, digits = 2)
EC1983_9C_Nuptake_StM$N <- EC1983_9C_Nuptake_StM$N/1000000 #microM to M
EC1983_9C_Nuptake_StM$NuptakeRate <- EC1983_9C_Nuptake_StM$NuptakeRate/1000000/w_EN #convert micro g N gDW–1 h–1 to mol N gDW–1 h–1
#conversions 18C
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$ResidualNitrateConcentration
EC1983_18C_Nuptake_StM$N <- round(EC1983_18C_Nuptake_StM$N, digits = 2)
EC1983_18C_Nuptake_StM$N <- EC1983_18C_Nuptake_StM$N/1000000 #microM to M
EC1983_18C_Nuptake_StM$NuptakeRate <- EC1983_18C_Nuptake_StM$NuptakeRate/1000000/w_EN

#testing rounding
sol_EspinozaChapman1983_N_9$N <- round(sol_EspinozaChapman1983_N_9$N*1000000, digits = 3)/1000000
sol_EspinozaChapman1983_N_18$N <- round(sol_EspinozaChapman1983_N_18$N*1000000, digits = 3)/1000000

N_calibration <- ggplot() +
  geom_line(data = sol_EspinozaChapman1983_N_9, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 9°C")) +
  geom_line(data = sol_EspinozaChapman1983_N_18, mapping = aes(x = N*1000000, y = J_EN_A*1000000, color = "Model of Espinoza and Chapman (1983) at 18°C")) +
  geom_point(data = EC1983_9C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 9°C"), size=3) +
  geom_point(data = EC1983_18C_Nuptake_StM, mapping = aes(x = N*1000000, y = NuptakeRate*1000000, color="Espinoza and Chapman (1983), St. Margaret's Bay, 18°C"), shape = 23, fill = 'grey', size=3) +
  xlim(0, 80) +
  scale_color_manual(values = c("gray60", "gray0", "gray60", "gray0")) +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(legend.position="none") +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('N Concentration (μmol' ~NO[3]^{"-"}~ 'L'^"-1"*')'), y = bquote('N uptake (μmol' ~NO[3]^{"-"}~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('a)')

#Error calculations
er9 <- merge(EC1983_9C_Nuptake_StM, sol_EspinozaChapman1983_N_9, all.x = TRUE)
rmse(er9$NuptakeRate, er9$J_EN_A) #3.683799e-07

er18 <- merge(EC1983_18C_Nuptake_StM, sol_EspinozaChapman1983_N_18, all.x = TRUE)
rmse(er18$NuptakeRate, er18$J_EN_A) #2.606024e-07

######## Photosynthesis related ####
#Johansson2002
Johansson2002 <- read.csv("Johansson2002.csv", header = TRUE, fileEncoding="UTF-8-BOM")
#conversions
Johansson2002$Irradiance <- Johansson2002$Irradiance*3600*1e-6 #micromol photons m-2 s-1 to mol photons m-2 h-1
Johansson2002$O2production <- Johansson2002$O2production/1e+6*32/1000*3600 #micromol O2 kg DW-1 s-1 to g O2/g/h
Johansson2002$O2productionSHIFT <- Johansson2002$O2production + 0.001720976 #from net to gross

Photosynthesis_calibration <- ggplot(data = Johansson2002) +
  geom_line(data = sol_Johansson2002, mapping = aes(x = I, y = J_O*1000)) +
  geom_point(mapping = aes(x = Irradiance, y = O2productionSHIFT*1000), size = 3) +
  scale_color_grey() +
  theme_bw() +
  theme(axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank()) +
  theme(axis.text=element_text(size=12), axis.title=element_text(size=16)) +
  labs(x= bquote('Irradiance (mol γ m'^"-2"*' h'^"-1)"), y = bquote('Oxygen production (mg' ~O[2]~ 'g DW'^"-1"*' h'^"-1"*')')) +
  ggtitle('b)')

#error calculations
Johansson2002$I <- round(Johansson2002$Irradiance, digits = 6)
sol_Johansson2002$I <- round(sol_Johansson2002$I, digits = 6)
erPhoto <- merge(Johansson2002, sol_Johansson2002, all.x = TRUE)
rmse(erPhoto$O2productionSHIFT, erPhoto$J_O)

######## Combine calibration plot (Figure 5) #######
#Figure 5
grid.arrange(N_calibration, Photosynthesis_calibration, ncol=2)
