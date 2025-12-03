# Sugar Kelp DEB Model (R) — Scripts Overview and Usage

Inside `Scripts/`:
- `Data.R`
- `MODEL.R`
- `KelpDEB_model.R`
- `SolveR_R.R`

## Quick start overview
1. Load environmental data with `Data.R` first.
  - Update the line: `setwd("/path/to/your/SeaweedCultivationCase/Scripts")`.
  - It extracts and processes environmental datasets (PAR/irradiance, DIC, nitrate, temperature) and saves hourly results to `IR_result_PAR.csv`,`DICOctJun_Depths0_5_10.csv`, `Nitrate_NovMay.csv` and `temperature_20192020.csv`.
2. Open `MODEL.R`, set your working directory, and run it.
   - Update the line: `setwd("/path/to/your/SeaweedCultivationCase/Scripts")`.
   - Make sure all the environmental data is set up correctly
   - Then run `MODEL.R` to execute the DEB model using the prepared data.

## Script Details

### Data.R
- Purpose: Gather and preprocess environmental data for the specific location and timeframe (North Sea, Borssele farm; Nov 1, 2019 to May 31, 2020).
- Run this file first 
- What it does:
  - Reads daily irradiance files, computes hourly data, estimates net shortwave radiation (NSW), and converts to PAR at various depths using Fresnel reflection and attenuation.
  - Reads global DIC climatology (`TCO2_NNGv2LDEO_climatology.nc`), selects nearest grid points to the target coordinates, and prepares monthly values for depths 0/5/10 m.
  - Reads daily nitrate files (`mercatorbiomer4v2r1_global_mean_nut_20220101.nc`), takes grid cell closest to our study location, computes hourly data from November to May. 
  - Reads file `uurgeg_321_2011-2020.txt` from Europlatform (near our study location), takes the TZ (Sea surface temperature), takes the data from 2019-2020. 
  - Outputs `IR_result_PAR.csv`,`DICOctJun_Depths0_5_10.csv`, `Nitrate_NovMay.csv` and `temperature_20192020.csv` files with data on a hourly basis in 2019-2020.


### MODEL.R
- Purpose: Main run file for the sugar kelp DEB model.
- What it does:
  - Sets working directory. Update this to your local path.
  - Loads libraries and sources model components: `SolveR_R.R` and `KelpDEB_model.R`.
  - Loads and prepares environmental inputs:
    - Irradiance: reads `IR_result_PAR.csv` produced by `Data.R`, parses timestamp blocks (UTC), and uses `PAR_1m` as a forcing.
    - Nitrate: reads `Nitrate_NovMay.csv`, interpolates to hourly, fixes a small gap by duplicating rows, converts to mol/L.
    - Temperature: reads `temperature_20192020.csv`, subsets the date range, converts TZ to Kelvin.
    - DIC: reads `DICOctJun_Depths0_5_10.csv`, converts units to mol/L with seawater density, interpolates monthly values to hourly along the Nov–May window.
  - Constructs hourly forcing functions (`approxfun`): `T_field`, `N_field`, `I_field`, `CO_2`.
  - Defines model parameters, initial states, and time vector (`times_NS` hourly sequence).
  - Runs the ODE model via `deSolve::ode` through a helper `run_model`, producing `sol_all` with date, temperature in °C, and model outputs.
  - Generates diagnostic plots for irradiance, temperature, DIC, and nitrate.
  - Generates depth analysis
  - Generates sensitivity analysis
  - Generates a climate change growth prediction model for RCP 2.6, RCP 4.5, RCP 6 and RCP 8.5 


### KelpDEB_model.R
- Purpose: Defines the rate function `rates_NS` for the DEB model (North Sea context).
- What it does:
  - Computes temperature correction and resource-specific assimilation/uptake rates for nitrogen and carbon, including light-driven terms.
  - Calls `SolveR_R` to determine the specific growth rate `r` and associated fluxes through an iterative Newton–Raphson procedure.
  - Returns derivatives for state variables (`dm_ECdt`, `dm_ENdt`, `dM_Vdt`) and a set of auxiliary outputs (growth rate, biomass W, oxygen production, and internal fluxes).

### SolveR_R.R
- Purpose: Internal solver used by `KelpDEB_model.R` to compute growth rate `r` and flux partitioning.
- What it does:
  - Implements an iteration over reserve densities and maintenance costs to find `r` using a norm-based Newton–Raphson update.
  - Returns `r` and flux components (catabolic, maintenance, structural maintenance, rejections) with a convergence flag.
  - `MODEL.R` sources this file, and `rates_NS` calls `SolveR_R` for each time step.

## Run Order and Notes
- Order:
  1. Open `Data.R`, set your working directory to the `Scripts` folder.
  2. Run `Data.R` to generate `IR_result_PAR.csv`,`DICOctJun_Depths0_5_10.csv`, `Nitrate_NovMay.csv` and `temperature_20192020.csv`.
  3. Open `MODEL.R`, set your working directory to the `Scripts` folder.
  4. Run `MODEL.R` to build forcing functions and execute the DEB model.

- File expectations:
  - `IR_result_PAR.csv`,`DICOctJun_Depths0_5_10.csv`, `Nitrate_NovMay.csv` and `temperature_20192020.csv` must exist in `Scripts/`.
- Optional calibration:
  - Use `N_uptake_Calibration.R` and `Photosynthesis_Calibration.R` separately to fit parameters; they do not affect the basic `MODEL.R` run unless you integrate calibrated parameters. (We did not do this, only done by original authors)

## Tips
- Time zones: `MODEL.R` uses UTC for irradiance and DIC timestamps; ensure consistency across datasets.
- Verify input lengths match `times_NS` or adjust the sequence as needed.
- Gaps: `MODEL.R` includes a small nitrate gap fix; review if your nitrate inputs differ.

