from __future__ import annotations
import os
from pathlib import Path
from typing import List, Dict, Any, Optional

import requests


DATASET_URL = "https://api.dataplatform.knmi.nl/open-data/v1/datasets/nl_sat_data_knmi_osi_7d_nz/versions/1"
FILES_ENDPOINT = f"{DATASET_URL}/files"
DOWNLOAD_ENDPOINT = f"{DATASET_URL}/files/{{filename}}/url"

# Default download directory (same as project data folder if present)
DEFAULT_DOWNLOAD_DIR = Path("/Users/nubia/Desktop/Seeweed/Case/data")

# Provided API key; in practice prefer environment variable KNMI_API_KEY
DEFAULT_API_KEY = "eyJvcmciOiI1ZTU1NGUxOTI3NGE5NjAwMDEyYTNlYjEiLCJpZCI6ImVlNDFjMWI0MjlkODQ2MThiNWI4ZDViZDAyMTM2YTM3IiwiaCI6Im11cm11cjEyOCJ9"


def _auth_headers(api_key: Optional[str] = None) -> Dict[str, str]:
	token = api_key or os.getenv("KNMI_API_KEY") or DEFAULT_API_KEY
	return {"Authorization": token}


def list_files(limit: int = 100, api_key: Optional[str] = None) -> List[Dict[str, Any]]:
	"""List available files in the KNMI dataset.

	Returns a list of dicts with at least a 'filename' field.
	"""
	params = {"maxKeys": str(limit)}
	resp = requests.get(FILES_ENDPOINT, headers=_auth_headers(api_key), params=params, timeout=30)
	resp.raise_for_status()
	data = resp.json()
	return data.get("files", data)


def get_download_url(filename: str, api_key: Optional[str] = None) -> str:
	"""Get a pre-signed download URL for a specific filename."""
	url = DOWNLOAD_ENDPOINT.format(filename=filename)
	resp = requests.get(url, headers=_auth_headers(api_key), timeout=30)
	resp.raise_for_status()
	data = resp.json()
	# API returns {'temporaryDownloadUrl': 'https://...'}
	return data.get("temporaryDownloadUrl") or data.get("url") or data[list(data.keys())[0]]


def download_file(filename: str, dest_dir: Path = DEFAULT_DOWNLOAD_DIR, api_key: Optional[str] = None) -> Path:
	"""Download a file to dest_dir using the temporary download URL."""
	dest_dir.mkdir(parents=True, exist_ok=True)
	download_url = get_download_url(filename, api_key=api_key)
	dest_path = dest_dir / filename
	with requests.get(download_url, stream=True, timeout=300) as r:
		r.raise_for_status()
		with open(dest_path, "wb") as f:
			for chunk in r.iter_content(chunk_size=1024 * 1024):
				if chunk:
					f.write(chunk)
	print(f"Downloaded: {dest_path}")
	return dest_path


if __name__ == "__main__":
	# Example usage: list files and download the first one
	files = list_files(limit=400)
	print(f"Found {len(files)} files")
	for i, item in enumerate(files[210:313]):
		print(f"{i+1}. {item.get('filename')}")

	if files:
		first = files[210].get("filename")
		if first:
			print(f"Downloading: {first}")
			download_file(first)




