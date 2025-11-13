from __future__ import annotations
import numpy as np
import matplotlib.pyplot as plt

from tempsloading_satalite import FILE_satalite, read_satellite_h5


def plot_satellite_map(sat: dict, title: str | None = None, save_path: str | None = "satellite/satellite_map.png") -> None:
    """Render a simple map-like plot from satellite data.

    - If lat/lon 2D grids are available, uses pcolormesh.
    - If 1D lat/lon vectors are available, builds a meshgrid.
    - Otherwise, falls back to imshow.
    """
    data = sat.get("data")
    lat = sat.get("lat")
    lon = sat.get("lon")

    if data is None:
        raise ValueError("sat['data'] is required")

    data = np.array(data)

    # Compute reasonable vmin/vmax ignoring NaNs
    vmin = np.nanpercentile(data, 5) if np.isfinite(data).any() else None
    vmax = np.nanpercentile(data, 95) if np.isfinite(data).any() else None

    fig, ax = plt.subplots(figsize=(8, 6), constrained_layout=True)

    plotted = False
    if lat is not None and lon is not None:
        lat = np.array(lat)
        lon = np.array(lon)
        # If 1D vectors, make grid
        if lat.ndim == 1 and lon.ndim == 1 and data.ndim == 2:
            Lon, Lat = np.meshgrid(lon, lat)
            h = ax.pcolormesh(Lon, Lat, data, shading='auto', cmap='turbo', vmin=vmin, vmax=vmax)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            plotted = True
        elif lat.ndim == 2 and lon.ndim == 2 and data.shape == lat.shape == lon.shape:
            h = ax.pcolormesh(lon, lat, data, shading='auto', cmap='turbo', vmin=vmin, vmax=vmax)
            ax.set_xlabel('Longitude')
            ax.set_ylabel('Latitude')
            plotted = True

    if not plotted:
        # Fallback to imshow without geographic axes
        h = ax.imshow(data, cmap='turbo', origin='lower', vmin=vmin, vmax=vmax)
        ax.set_xlabel('X index')
        ax.set_ylabel('Y index')

    cbar = plt.colorbar(h, ax=ax, label='Value')

    if title is None:
        attrs = sat.get('attributes', {}) or {}
        title = str(attrs.get('title', 'Satellite SST'))
    ax.set_title(title)

    if save_path:
        plt.savefig(save_path, dpi=150)
        print(f"Saved map to {save_path}")

    # Also show interactively if running locally
    try:
        plt.show()
    except Exception:
        pass


if __name__ == "__main__":
    # Minimal CLI-like usage: read FILE_satalite and plot
    sat = read_satellite_h5(FILE_satalite)
    plot_satellite_map(sat)
