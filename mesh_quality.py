#!/usr/bin/env python

from netCDF4 import Dataset
import numpy as np


r_earth = 6371229.0    # MPAS-Atmosphere's assumed Earth radius (m)


def vec_dot(A, B):
    """ For two (n x 3) matrices, A and B, compute the inner product
    of row i of A with row i of B, for i =1..n, and return the result
    as an array of length n.
    """
    return np.einsum('ij,ij->i', A, B)


def plane_distance(a, b):
    """ Compute the Euclidean distance between two points, a and b, given
    as arrays of coordinates.
    """
    d = b - a
    return np.sqrt(np.sum(d * d, axis=1))


def sphere_distance(a, b):
    """ Compute the great-circle arc distance between two points, a and b,
    both of which lie on the surface of the same sphere.
    """
    c = plane_distance(a, b)
    return 2.0 * np.arcsin(c / 2.0)


def sphere_angle(A, B, C):
    """: Computes the angle between arcs AB and AC, given points A, B, and C
    Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
    """
    a = sphere_distance(B, C)
    b = sphere_distance(A, C)
    c = sphere_distance(A, B)

    AB = B - A
    AC = C - A

    D = np.cross(AB, AC)

    s = 0.5 * (a + b + c)
    sin_angle = np.sqrt((np.sin(s-b)*np.sin(s-c))/(np.sin(b)*np.sin(c))) # Eqn. (28)

    return np.where(vec_dot(D, A) >= 0.0,
                    2.0 * np.arcsin(sin_angle),
                   -2.0 * np.arcsin(sin_angle))


def check_distances(dcEdge, dvEdge):
    print('')
    print('Min/max dcEdge (m):', r_earth * np.min(dcEdge), r_earth * np.max(dcEdge))
    print('Min/max dvEdge (m):', r_earth * np.min(dvEdge), r_earth * np.max(dvEdge))


def check_areas(areaCell, areaTriangle, kiteAreasOnVertex):
    print('')
    print('sum(areaCell):         ', np.sum(areaCell))
    print('sum(areaTriangle):     ', np.sum(areaTriangle))
    print('sum(kiteAreasOnVertex):', np.sum(kiteAreasOnVertex))
    print('4 * Pi:                ', 4.0 * np.pi)


def check_obtuse_triangles(nVertices, vertexDegree, xCell, yCell, zCell, xVertex, yVertex, zVertex, cellsOnVertex):
    locCell = np.column_stack((xCell, yCell, zCell))
    locVertex = np.column_stack((xVertex, yVertex, zVertex))

    sum_angles = np.zeros((nVertices))
    for k in range(vertexDegree):
        sum_angles += sphere_angle(locVertex, locCell[cellsOnVertex[:,k]], locCell[cellsOnVertex[:,(k+1)%vertexDegree]])

    obtuse = np.ma.array(np.arange(0,nVertices,1,dtype=np.int32), mask=np.logical_not(sum_angles < np.pi))
    print('')
    if np.any(sum_angles < np.pi):
        print('Obtuse triangles:', obtuse.compressed())
    else:
        print('No obtuse triangles!')


def check_resolution_gradient(nominalMinDc, meshDensity, nEdgesOnCell, edgesOnCell, cellsOnEdge, dcEdge):
    nominalDx = r_earth * nominalMinDc * np.power(1.0 / meshDensity, 0.25)
    gradient = np.abs(nominalDx[cellsOnEdge[:,0]] - nominalDx[cellsOnEdge[:,1]]) / dcEdge / r_earth

    print('')
    print('Min nominal cell size gradient:', np.min(gradient))
    print('Max nominal cell size gradient:', np.max(gradient))


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_file', help='the name of the netCDF file with mesh fields')
    args = parser.parse_args()

    f = Dataset(args.mesh_file)

    if f.sphere_radius != 1.0:
        print('Error: Argument must be a mesh defined on the unit sphere.')
        print('       Global attribute sphere_radius is set to', f.sphere_radius)
        sys.exit(1)

    nCells = f.dimensions['nCells'].size
    nVertices = f.dimensions['nVertices'].size
    nEdges = f.dimensions['nEdges'].size
    maxEdges = f.dimensions['maxEdges'].size
    vertexDegree = f.dimensions['vertexDegree'].size
    nEdgesOnCell = f.variables['nEdgesOnCell'][:] - 1
    cellsOnVertex = f.variables['cellsOnVertex'][:] - 1
    cellsOnEdge = f.variables['cellsOnEdge'][:] - 1
    edgesOnCell = f.variables['edgesOnCell'][:] - 1
    xCell = f.variables['xCell'][:]
    yCell = f.variables['yCell'][:]
    zCell = f.variables['zCell'][:]
    xVertex = f.variables['xVertex'][:]
    yVertex = f.variables['yVertex'][:]
    zVertex = f.variables['zVertex'][:]
    dcEdge = f.variables['dcEdge'][:]
    dvEdge = f.variables['dvEdge'][:]
    areaCell = f.variables['areaCell'][:]
    areaTriangle = f.variables['areaTriangle'][:]
    kiteAreasOnVertex = f.variables['kiteAreasOnVertex'][:]
    nominalMinDc = f.variables['nominalMinDc'][:]
    meshDensity = f.variables['meshDensity'][:]

    check_distances(dcEdge, dvEdge)

    check_areas(areaCell, areaTriangle, kiteAreasOnVertex)

    check_obtuse_triangles(nVertices, vertexDegree, xCell, yCell, zCell, xVertex, yVertex, zVertex, cellsOnVertex)

    check_resolution_gradient(nominalMinDc, meshDensity, nEdgesOnCell, edgesOnCell, cellsOnEdge, dcEdge)

    f.close()
