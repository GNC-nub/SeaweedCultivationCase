
from pathlib import Path
import pandas as pd 
DATA_DIR = Path("/Users/nubia/Desktop/Seeweed/Case/data")
FILENAMES = [
	"uurgeg_214_2011-2020.txt",
	"uurgeg_321_2011-2020.txt",
	"uurgeg_201_2011-2020.txt",
	"uurgeg_252_2011-2020.txt"
]

FILE_214 = DATA_DIR / "uurgeg_214_2011-2020.txt"
FILE_321 = DATA_DIR / "uurgeg_321_2011-2020.txt"
FILE_201 = DATA_DIR / "uurgeg_201_2011-2020.txt"
FILE_252 = DATA_DIR / "uurgeg_252_2011-2020.txt"


def extract_headers(path: Path) -> list:
	"""Extract comma-separated headers from the first line that starts with '#'.

	Returns a list of header names with surrounding whitespace stripped.
	Raises FileNotFoundError if path missing and ValueError if no header line found.
	"""

	with path.open("r", encoding="utf-8", errors="ignore") as f:
		for line in f:
			# Allow leading whitespace before '#'
			stripped = line.lstrip()
			if stripped.startswith('#'):
				raw = stripped[1:].strip()  # remove leading '#'
				headers = [h.strip() for h in raw.split(',') if h.strip()]
				if not headers:
					raise ValueError("Header line found but no comma-separated fields detected")
				return headers
	raise ValueError("No header line starting with '#' found in file")


def load_dataframe_with_headers(path: Path) -> pd.DataFrame:
	"""Load a single file into a DataFrame using extracted '#' comma headers.

	The header line is treated as comments (skipped) when reading the rest.
	"""
	headers = extract_headers(path)
	# pandas will skip commented header line; provide names since header is commented
	df = pd.read_csv(path, comment='#', header=None, names=headers)
	return df


def build_dataframe(path: Path, min_date: str | None = None) -> pd.DataFrame:
	"""Create a dataframe with columns 'dates', 'hour', 'TZ', 'T' from raw file.

	File format assumptions:
	  * First header line starts with '#'
	  * Subsequent lines are comma-separated values
	  * Second field = date, third field = hour, last field = TZ
	"""
	rows = []
	seen_header = False
	with path.open('r', encoding='utf-8', errors='ignore') as f:
		for line in f:
			if not line.strip():
				continue
			stripped = line.lstrip()
			if not seen_header:
				# Skip all lines until we meet the first header/comment line
				if stripped.startswith('#'):
					seen_header = True
				continue
			# After the header line, skip any further comment lines
			if stripped.startswith('#'):
				continue
			raw_parts = line.rstrip('\n').split(',')
			# Normalize fields: strip whitespace; map empty strings to None
			parts = [p.strip() if p is not None else None for p in raw_parts]
			parts = [val if val != '' else None for val in parts]
			if len(parts) < 8:
				# skip malformed line
				continue
			date_val = parts[1]
			hour_val = parts[2]
			tz_val = parts[-1]
			t_val = parts[7]
			if min_date is not None:
				# Keep rows strictly greater than min_date
				# Only compare when date_val present
				if date_val is None or date_val <= str(min_date):
					continue
			rows.append((date_val, hour_val, tz_val, t_val))
	df = pd.DataFrame(rows, columns=['dates', 'hour', 'TZ', 'T'])
	return df



