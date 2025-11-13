from tempsloading import (
	build_dataframe,
	FILE_321,
	FILE_252
)
import pandas as pd
import matplotlib.pyplot as plt
from datetime import datetime


def compute_diff_and_stats(df: pd.DataFrame) -> dict:
	"""Add TZ-T difference per timestamp and compute its standard deviation.

	Returns a dict with summary stats.
	"""
	# Ensure numeric
	df = df.copy()
	df["TZ"] = pd.to_numeric(df["TZ"], errors="coerce")
	df["T"] = pd.to_numeric(df["T"], errors="coerce")

	# Row-wise difference for each timestamp
	df["TZ_minus_T"] = df["TZ"] - df["T"]

	std_diff = float(df["TZ_minus_T"].std(skipna=True))
	mean_diff = float(df["TZ_minus_T"].mean(skipna=True))

	# Print a short preview
	print("Preview with difference (first 5 rows):")
	print(df[["dates", "hour", "TZ", "T", "TZ_minus_T"]].head(5))
	print(f"mean diff: {mean_diff}")
	print(f"Standard deviation of (TZ - T): {std_diff}")

	return {
		"std": std_diff,
		"mean": mean_diff,
		"series": df["TZ_minus_T"],
	}


def plot_diff_boxplot(series: pd.Series, mean_val: float, std_val: float, title: str = "TZ - T boxplot") -> None:
	clean = series.dropna()
	fig, ax = plt.subplots(figsize=(6, 5), constrained_layout=True)
	bp = ax.boxplot(clean, vert=True, patch_artist=True, showmeans=True,
					meanprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black'),
					boxprops=dict(facecolor='#87CEFA'))
	ax.set_title(title)
	ax.set_ylabel("TZ - T")

	# Annotate mean and std
	y_min, y_max = ax.get_ylim()
	text_y = y_max - 0.05 * (y_max - y_min)
	ax.text(1.15, text_y, f"mean = {mean_val:.3f}\nstd = {std_val:.3f}",
			transform=ax.get_yaxis_transform(which='grid'),
			va='top', ha='left', bbox=dict(boxstyle='round', facecolor='white', alpha=0.7))

	plt.show()


def plot_time_series(df: pd.DataFrame, title: str = "TZ and T over time") -> None:
	"""Plot TZ and T as a function of time using a combined datetime index.

	Assumes 'dates' is in YYYYMMDD format and 'hour' is integer-like [0..23].
	"""
	# Build datetime column
	out = df.copy()
	out["dates"] = out["dates"].astype(str)
	# Zero-pad hour and build timestamp 'YYYYMMDD HH:00'
	hour_str = pd.to_numeric(out["hour"], errors="coerce").fillna(0).astype(int).astype(str).str.zfill(2)
	ts_str = out["dates"] + " " + hour_str + ":00"
	out["timestamp"] = pd.to_datetime(ts_str, format="%Y%m%d %H:%M", errors="coerce")

	# Numeric TZ and T
	out["TZ"] = pd.to_numeric(out["TZ"], errors="coerce")
	out["T"] = pd.to_numeric(out["T"], errors="coerce")

	# Sort by time
	out = out.sort_values("timestamp")

	fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
	ax.plot(out["timestamp"], out["TZ"], label="TZ", color="#1f77b4", linewidth=1)
	ax.plot(out["timestamp"], out["T"], label="T", color="#ff7f0e", linewidth=1)
	ax.set_title(title)
	ax.set_xlabel("Time")
	ax.set_ylabel("Value")
	ax.legend()
	ax.grid(True, linestyle=":", alpha=0.5)
	plt.show()


if __name__ == "__main__":
	# Build DataFrame for the requested file, filtering dates > 20190101
	df_321 = build_dataframe(FILE_321, min_date="20190101")
	df_252 = build_dataframe(FILE_252, min_date="20190101")
	stats = compute_diff_and_stats(df_321)
	plot_diff_boxplot(stats["series"], stats["mean"], stats["std"], title="TZ - T (FILE_321)")
	plot_time_series(df_321, title="TZ and T over time (FILE_321)")
	plot_time_series(df_252, title="TZ and T over time (FILE_252)")




