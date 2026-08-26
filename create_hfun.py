#!/usr/bin/env python

import numpy as np
import hfun

hfun_min = hfun.finest_resolution()

#
# Determine resolution of output lat-lon grid based on hfun_min. The approach
# below is probably quite conservative if we assume that the mesh refinement
# region(s) and resolution transition zones are much larger in scale than
# hfun_min.
#
r_earth = 6371.229
deg_to_km = 2.0 * np.pi * r_earth / 360.0
nlat = int(180.0 * deg_to_km / hfun_min) + 1

#
# Generate 2-d lat-lon meshgrid (radians)
#
lats = np.linspace(-0.5 * np.pi, 0.5 * np.pi, num=nlat, endpoint=True)
lons = np.linspace(-np.pi, np.pi, num=2 * nlat, endpoint=True)
latgrid, longrid = np.meshgrid(lats, lons)

nlats = lats.size
nlons = lons.size
npts = nlats * nlons

#
# Invoke code in hfun.py to return an array of grid distances at the coordinates
# provided by the longrid and latgrid arrays
#
distance = hfun.get_hfun(longrid, latgrid)

#
# Write grid distances out to a JIGSAW-compatible file
#
with open('HFUN.msh', 'w') as f:
    f.write('MSHID=3;ellipsoid-grid\n')
    f.write('NDIMS=2\n')
    f.write(f'COORD=1;{nlons}\n')
    for lon in lons:
        f.write(f'{lon}\n')
    f.write(f'COORD=2;{nlats}\n')
    for lat in lats:
        f.write(f'{lat}\n')

    f.write(f'VALUE={npts}; 1\n')
    for d in distance.flatten():
        f.write(f'{d}\n')
