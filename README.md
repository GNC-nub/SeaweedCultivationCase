# Seaweed Temps North Sea

Utilities to load and process KNMI station time-series and satellite SST data for the North Sea.

Default data directory:
- `/Users/nubia/Desktop/Seeweed/Case/data`

Key station files (examples):
- `uurgeg_214_2011-2020.txt`
- `uurgeg_252_2011-2020.txt`
- `uurgeg_321_2011-2020.txt`

Optional satellite file example:
- `NOAA_1_AVHRR_SST_NZ7_20251110.h5`

## Install dependencies

```bash
pip install -r requirements.txt
```

Packages: pandas, numpy, matplotlib, requests, h5py

## Load KNMI station text files

Module: `tempsloading.py`

- Header handling: the header line starts with `#` and fields are comma separated. Lines before the first `#` are ignored.
- Parsing: `build_dataframe(path, min_date=None)` returns a DataFrame with columns:
	- `dates` (string YYYYMMDD)
	- `hour` (string/integer-like hour 0–23)
	- `TZ` (last field of each line)
	- `T` (8th field of each line)
- Missing values: whitespace-only fields are converted to `None` and become NaN in pandas.
- Filtering: pass `min_date="YYYYMMDD"` to keep only rows with `dates > min_date`.

Quick usage:

```bash
python3 main.py
```

Or in Python:

```python
from tempsloading import build_dataframe, FILE_321
df = build_dataframe(FILE_321, min_date="20190101")
print(df.head())
```

## Export station data to CSV

Script: `export_to_csv.py`

Writes CSV files next to their source `.txt` files (same data directory).

```bash
python3 export_to_csv.py
```

To apply a date filter, set `min_date` in `export_to_csv.py`.

## TZ vs T analysis and plots

Script: `diff_TZ_T.py`

- Computes per-timestamp difference `TZ_minus_T = TZ - T`
- Prints preview, mean, and standard deviation
- Plots:
	- Boxplot of `TZ_minus_T` with mean and std annotated
	- Line plot of `TZ` and `T` over time using a combined datetime from `dates` + `hour`

Run:

```bash
python3 diff_TZ_T.py
```

## Satellite SST visualization

- Reader: either `tempsloading.py::read_satellite_h5` or your `satalite/tempsloading_satalite.py` (depending on your local setup) reads SST data plus optional lat/lon grids.
- Plotter: `satalite/viz_satellite.py` creates a map-like visualization.

Run:

```bash
python3 satalite/viz_satellite.py
```

Behavior:
- Uses `pcolormesh` when lat/lon grids are present and shape-compatible.
- Falls back to `imshow` if coordinates are unavailable; axes will show pixel indices.
- Optionally transposes data to match desired orientation.

## KNMI Open Data API helper (satellite dataset)

Script: `satellite/API_loader.py`

- Lists available files and can download a selected one to the data directory.
- Auth: uses the provided token by default; you can override with env var `KNMI_API_KEY`.

Run:

```bash
python3 satellite/API_loader.py
```

## Project tips

- Data directory is hard-coded as `/Users/nubia/Desktop/Seeweed/Case/data`. Adjust in `tempsloading.py` if needed.
- If station files have different field order, update the extraction indexes in `build_dataframe` accordingly.
- For satellite data, dataset names can vary; update candidate keys in the reader if necessary.
