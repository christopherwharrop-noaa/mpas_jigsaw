#!/usr/bin/env python

import subprocess
import numpy as np
from mpi4py import MPI
import hfun

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

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
# Generate 1-d lat and lon arrays (radians)
#
lats = np.linspace(-0.5 * np.pi, 0.5 * np.pi, num=nlat, endpoint=True)
lons = np.linspace(-np.pi, np.pi, num=2 * nlat, endpoint=True)

nlats = lats.size
nlons = lons.size

#
# Split longitude columns across MPI ranks
#
counts = np.array([(nlons // nprocs) + (1 if i < nlons % nprocs else 0)
                   for i in range(nprocs)])
offsets = np.cumsum(counts) - counts
my_lons = lons[offsets[rank]:offsets[rank] + counts[rank]]

#
# Build local meshgrid and compute hfun for this rank's columns
#
local_latgrid, local_longrid = np.meshgrid(lats, my_lons)
local_distance = hfun.get_hfun(local_longrid, local_latgrid)

#
# Gather results to rank 0
#
if rank == 0:
    distance = np.empty((nlons, nlats), dtype=np.float64)
else:
    distance = None

sendcounts = counts * nlats
displacements = offsets * nlats
comm.Gatherv(np.ascontiguousarray(local_distance),
             [distance, sendcounts, displacements, MPI.DOUBLE],
             root=0)

#
# Rank 0 pipes binary arrays to write_hfun which writes HFUN.msh
#
if rank == 0:
    proc = subprocess.Popen(
        ['./bin/write_hfun', str(nlons), str(nlats)],
        stdin=subprocess.PIPE,
    )
    proc.stdin.write(lons.tobytes())
    proc.stdin.write(lats.tobytes())
    proc.stdin.write(distance.tobytes())
    proc.stdin.close()
    rc = proc.wait()
    if rc != 0:
        raise RuntimeError(f'write_hfun exited with code {rc}')
