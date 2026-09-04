#!/usr/bin/env python

from pathlib import Path
import shutil

import numpy as np
import hfun

hfun_min = hfun.finest_resolution()


def cart_to_geo(x, y, z):
    from numpy import arctan2, arcsin
    lam = arctan2(y, x)
    phi = arcsin(z)
    return (lam, phi)

coords = np.loadtxt('SaveVertices')
coords = coords / 6371.229

longitude, latitude = cart_to_geo(coords[:,0], coords[:,1], coords[:,2])

dx = hfun.get_hfun(longitude, latitude)

density = (1.0 / (dx / hfun_min))**4

with open('SaveDensity', 'w') as f:
    for d in density:
        f.write(f'{d}\n')

shutil.copyfile(Path(__file__).with_name('hfun.py'), 'SaveCode')
shutil.copyfile('mesh.yaml', 'SaveConfig')
