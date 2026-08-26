#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import cartopy.crs as ccrs
import cartopy.feature as cfeature


def plot_mesh(grid_file, output_file="mesh_plot.png", extent=None):
    ds = Dataset(grid_file, 'r')

    lat_cell = np.degrees(ds.variables['latCell'][:])
    lon_cell = np.degrees(ds.variables['lonCell'][:])
    lat_vtx = np.degrees(ds.variables['latVertex'][:])
    lon_vtx = np.degrees(ds.variables['lonVertex'][:])
    n_edges_on_cell = ds.variables['nEdgesOnCell'][:]
    vertices_on_cell = ds.variables['verticesOnCell'][:]
    dc_edge = ds.variables['dcEdge'][:]

    # Shift longitudes to [-180, 180]
    lon_cell = np.where(lon_cell > 180, lon_cell - 360, lon_cell)
    lon_vtx = np.where(lon_vtx > 180, lon_vtx - 360, lon_vtx)

    r_earth = 6371.229
    dc_edge_km = dc_edge * r_earth

    ds.close()

    if extent is None:
        extent = [-180, 180, -90, 90]

    # Filter cells within the plot extent (with margin)
    margin = 5
    mask = (
        (lon_cell >= extent[0] - margin) & (lon_cell <= extent[1] + margin) &
        (lat_cell >= extent[2] - margin) & (lat_cell <= extent[3] + margin)
    )
    cell_indices = np.where(mask)[0]

    proj = ccrs.PlateCarree()
    fig, ax = plt.subplots(figsize=(14, 10), subplot_kw={'projection': proj})
    ax.set_extent(extent, crs=proj)

    # Build polygons for each cell
    polygons = []
    for idx in cell_indices:
        n = n_edges_on_cell[idx]
        vtx_ids = vertices_on_cell[idx, :n] - 1
        lons = lon_vtx[vtx_ids]
        lats = lat_vtx[vtx_ids]

        # Skip cells that wrap around the date line
        if np.max(lons) - np.min(lons) > 180:
            continue

        poly = np.column_stack((lons, lats))
        polygons.append(poly)

    pc = PolyCollection(
        polygons, closed=True, transform=proj,
        linewidths=0.25, edgecolors='#999999', facecolors='white',
        zorder=2
    )
    ax.add_collection(pc)

    # Draw coastlines and borders on top of the mesh
    ax.add_feature(cfeature.COASTLINE, linewidth=1.8, color='darkblue', zorder=3)
    ax.add_feature(cfeature.BORDERS, linewidth=0.8, color='darkblue',
                   linestyle='--', zorder=3)
    ax.gridlines(draw_labels=True, linewidth=0.3, alpha=0.5)

    ax.set_title(f'MPAS Mesh ({len(cell_indices)} cells shown)\n'
                 f'dcEdge range: {dc_edge_km.min():.1f} - {dc_edge_km.max():.1f} km')

    plt.savefig(output_file, dpi=200, bbox_inches='tight')
    print(f"Saved plot to {output_file}")


if __name__ == '__main__':
    grid_file = sys.argv[1] if len(sys.argv) > 1 else 'grid.nc'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'mesh_plot.png'

    # Optional: set extent as lon_min lon_max lat_min lat_max
    extent = None
    if len(sys.argv) > 5:
        extent = [float(sys.argv[i]) for i in range(3, 7)]

    plot_mesh(grid_file, output_file, extent)
