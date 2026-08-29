#!/usr/bin/env python

import sys
import numpy as np
from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import cartopy.io.shapereader as shpreader


def latlon_to_xyz(lat_deg, lon_deg, r=1.0):
    lat = np.radians(lat_deg)
    lon = np.radians(lon_deg)
    return (r * np.cos(lon) * np.cos(lat),
            r * np.sin(lon) * np.cos(lat),
            r * np.sin(lat))


def get_coastlines():
    """Load Natural Earth coastlines as 3D line segments on the unit sphere."""
    coast_file = shpreader.natural_earth(
        resolution='110m', category='physical', name='coastline'
    )
    reader = shpreader.Reader(coast_file)
    lines_3d = []
    for geom in reader.geometries():
        if geom.geom_type == 'MultiLineString':
            parts = list(geom.geoms)
        else:
            parts = [geom]
        for part in parts:
            coords = np.array(part.coords)
            x, y, z = latlon_to_xyz(coords[:, 1], coords[:, 0], r=1.001)
            seg = np.column_stack((x, y, z))
            lines_3d.append(seg)
    return lines_3d


def plot_mesh_sphere(grid_file, output_file="mesh_sphere.png",
                     center_lat=None, center_lon=None, zoom=1.0):
    ds = Dataset(grid_file, 'r')

    x_cell = ds.variables['xCell'][:]
    y_cell = ds.variables['yCell'][:]
    z_cell = ds.variables['zCell'][:]
    x_vtx = ds.variables['xVertex'][:]
    y_vtx = ds.variables['yVertex'][:]
    z_vtx = ds.variables['zVertex'][:]
    n_edges_on_cell = ds.variables['nEdgesOnCell'][:]
    vertices_on_cell = ds.variables['verticesOnCell'][:]

    ds.close()

    # Set view direction
    if center_lat is not None and center_lon is not None:
        view_x, view_y, view_z = latlon_to_xyz(center_lat, center_lon)
    else:
        view_x, view_y, view_z = 1.0, 0.0, 0.0

    # Only draw cells on the visible hemisphere
    dot = x_cell * view_x + y_cell * view_y + z_cell * view_z
    cell_indices = np.where(dot > 0.05)[0]

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Draw mesh polygons (light gray edges, white fill)
    polygons = []
    for idx in cell_indices:
        n = n_edges_on_cell[idx]
        vtx_ids = vertices_on_cell[idx, :n] - 1
        verts = np.column_stack((x_vtx[vtx_ids], y_vtx[vtx_ids], z_vtx[vtx_ids]))
        polygons.append(verts)

    pc = Poly3DCollection(
        polygons, edgecolors='#999999', facecolors='white',
        linewidths=0.25, alpha=0.95
    )
    ax.add_collection3d(pc)

    # Draw coastlines (thick dark blue, rendered on top)
    coastlines = get_coastlines()
    for seg in coastlines:
        dots = seg[:, 0] * view_x + seg[:, 1] * view_y + seg[:, 2] * view_z
        mask = dots > 0.0
        splits = np.where(np.diff(mask.astype(int)) != 0)[0] + 1
        for sub in np.split(np.arange(len(seg)), splits):
            if len(sub) > 1 and mask[sub[0]]:
                s = seg[sub]
                ax.plot3D(s[:, 0], s[:, 1], s[:, 2],
                          color='darkblue', linewidth=1.8, zorder=10)

    # Draw refinement region boundaries and transition zones
    try:
        import hfun
        for bd in hfun.get_region_boundaries(mesh_file=grid_file):
            if bd['kind'] == 'boundary':
                lw, ls = 1.5, '-'
            else:
                lw, ls = 0.8, '--'
            bx, by, bz = latlon_to_xyz(bd['lats_deg'], bd['lons_deg'], r=1.002)
            pts = np.column_stack((bx, by, bz))
            dots = pts[:, 0] * view_x + pts[:, 1] * view_y + pts[:, 2] * view_z
            mask = dots > 0.0
            splits = np.where(np.diff(mask.astype(int)) != 0)[0] + 1
            for sub in np.split(np.arange(len(pts)), splits):
                if len(sub) > 1 and mask[sub[0]]:
                    s = pts[sub]
                    ax.plot3D(s[:, 0], s[:, 1], s[:, 2],
                              color='red', linewidth=lw, linestyle=ls, zorder=10)
    except Exception:
        pass

    # Zoom: smaller limits = more zoomed in
    lim = 1.0 / zoom
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_aspect('equal')
    ax.axis('off')

    elev = np.degrees(np.arcsin(view_z))
    azim = np.degrees(np.arctan2(view_y, view_x))
    ax.view_init(elev=elev, azim=azim)

    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {output_file}")


if __name__ == '__main__':
    grid_file = sys.argv[1] if len(sys.argv) > 1 else 'grid.nc'
    output_file = sys.argv[2] if len(sys.argv) > 2 else 'mesh_sphere.png'
    center_lat = float(sys.argv[3]) if len(sys.argv) > 3 else None
    center_lon = float(sys.argv[4]) if len(sys.argv) > 4 else None
    zoom = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0

    plot_mesh_sphere(grid_file, output_file, center_lat, center_lon, zoom)
