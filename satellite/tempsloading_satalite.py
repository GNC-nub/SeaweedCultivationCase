
from pathlib import Path
import pandas as pd 


DATA_DIR = Path("/Users/nubia/Desktop/Seeweed/Case/data")
FILENAMES = [
	"NOAA_1_AVHRR_SST_NZ7_20251110.h5",
]

FILE_satalite = DATA_DIR / "NOAA_1_AVHRR_SST_NZ7_20251110.h5"

def read_satellite_h5(path: Path):
	"""Read satellite HDF5 file returning a dict with data arrays.

	Attempts to locate common dataset names for sea surface temperature (SST) along with latitude and longitude grids.
	Returns: {"data": <2D array>, "lat": <2D or 1D array>, "lon": <2D or 1D array>, "attributes": <file attrs>}.

	Requires h5py. Raises FileNotFoundError if the file doesn't exist.
	"""
	if not path.exists():
		raise FileNotFoundError(f"Satellite file not found: {path}")
	try:
		import h5py  # type: ignore
	except ImportError as e:
		raise ImportError("h5py is required to read satellite HDF5 files. Install with 'pip install h5py'.") from e

	candidate_data_keys = ["sst", "SST", "sea_surface_temperature", "data", "SST_Data"]
	candidate_lat_keys = ["lat", "latitude", "Latitude", "nav_lat"]
	candidate_lon_keys = ["lon", "longitude", "Longitude", "nav_lon"]

	result = {"data": None, "lat": None, "lon": None, "attributes": {}}
	with h5py.File(path, "r") as h5:
		# Capture global attrs
		for k, v in h5.attrs.items():
			result["attributes"][k] = v
		# Breadth-first search through groups
		def find_dataset(keys):
			for name in h5:
				obj = h5[name]
				if hasattr(obj, 'keys'):
					# group: search inside
					for sub in obj.keys():
						full = f"{name}/{sub}"
						for candidate in keys:
							if candidate.lower() in sub.lower():
								return h5[full][()]
				elif hasattr(obj, 'shape'):
					for candidate in keys:
						if candidate.lower() in name.lower():
							return obj[()]
			return None

		data = find_dataset(candidate_data_keys)
		lat = find_dataset(candidate_lat_keys)
		lon = find_dataset(candidate_lon_keys)

		result["data"] = data
		result["lat"] = lat
		result["lon"] = lon

	if result["data"] is None:
		raise ValueError("Could not locate SST data in HDF5 file using known candidate keys.")
	return result



