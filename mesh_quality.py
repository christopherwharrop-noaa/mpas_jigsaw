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


def check_cell_types(nEdgesOnCell):
    names = {3: 'triangles', 4: 'quadrilaterals', 5: 'pentagons',
             6: 'hexagons', 7: 'heptagons', 8: 'octagons',
             9: 'nonagons', 10: 'decagons'}
    counts = dict(zip(*np.unique(nEdgesOnCell, return_counts=True)))
    nCells = len(nEdgesOnCell)
    n_hex = counts.get(6, 0)
    print('')
    print(f'Cell type census ({nCells} cells):')
    for nsides in sorted(counts):
        label = names.get(nsides, f'{nsides}-gons')
        pct = 100.0 * counts[nsides] / nCells
        marker = '' if nsides == 6 else '  *'
        print(f'  {label:>16s}: {counts[nsides]:>8d}  ({pct:5.2f}%){marker}')
    print(f'  {"non-hexagons":>16s}: {nCells - n_hex:>8d}  ({100.0*(nCells - n_hex)/nCells:5.2f}%)')


def nice_bin_width(bin_width_km):
    """Round a positive bin width to a human-friendly 1-2-2.5-5 sequence."""
    scale = 10.0 ** np.floor(np.log10(bin_width_km))
    choices = scale * np.array([1.0, 2.0, 2.5, 5.0, 10.0])
    return choices[np.argmin(np.abs(choices - bin_width_km))]


def cell_size_bin_edges(cell_size_km, bin_width_km=None,
                        bin_center_km=None):
    """Return fixed or Freedman-Diaconis histogram bin edges."""
    min_size = cell_size_km.min()
    max_size = cell_size_km.max()

    if bin_width_km is None:
        q25, q75 = np.percentile(cell_size_km, [25, 75])
        automatic = True
        bin_width_km = 2.0 * (q75 - q25) / cell_size_km.size ** (1.0 / 3.0)
    else:
        automatic = False

    if not np.isfinite(bin_width_km) or bin_width_km <= 0.0:
        bin_width_km = max_size - min_size

    if bin_width_km == 0.0:
        return np.array([min_size - 0.5, max_size + 0.5])

    if automatic:
        bin_width_km = max(bin_width_km, (max_size - min_size) / 50.0)
        bin_width_km = nice_bin_width(bin_width_km)

    if bin_center_km is not None:
        first_center = (bin_center_km + bin_width_km *
                        np.floor((min_size - bin_center_km) / bin_width_km))
        last_center = (bin_center_km + bin_width_km *
                       np.ceil((max_size - bin_center_km) / bin_width_km))
        return np.arange(first_center - 0.5 * bin_width_km,
                         last_center + 1.5 * bin_width_km, bin_width_km)

    lower = np.floor(min_size / bin_width_km) * bin_width_km
    upper = np.ceil(max_size / bin_width_km) * bin_width_km
    return np.arange(lower, upper + bin_width_km, bin_width_km)


def check_cell_size_distribution(areaCell, nominalMinDc, meshDensity,
                                 bin_width_km=None):
    """Report equivalent regular-hexagon spacing derived from dual-cell area."""
    cell_area_m2 = np.asarray(areaCell) * r_earth ** 2
    cell_size_km = np.sqrt(2.0 * cell_area_m2 / np.sqrt(3.0)) / 1000.0
    percentiles = np.percentile(cell_size_km, [1, 5, 25, 50, 75, 95, 99])
    target_size_km = (r_earth * nominalMinDc *
                      np.power(1.0 / meshDensity, 0.25) / 1000.0)
    is_quasi_uniform = np.allclose(target_size_km, target_size_km[0],
                                   rtol=1.0e-10, atol=1.0e-10)

    if is_quasi_uniform:
        if bin_width_km is None:
            bin_width_km = nice_bin_width(0.1 * target_size_km[0])
        bin_edges = cell_size_bin_edges(
            cell_size_km, bin_width_km, bin_center_km=target_size_km[0])
    else:
        bin_edges = cell_size_bin_edges(cell_size_km, bin_width_km)

    counts, bin_edges = np.histogram(cell_size_km, bins=bin_edges)
    width = bin_edges[1] - bin_edges[0]

    print('')
    print('Equivalent-hexagon cell spacing (km):')
    if is_quasi_uniform:
        print(f'  Target spacing: {target_size_km[0]:.3f} (quasi-uniform)')
    else:
        print(f'  Target range: {target_size_km.min():.3f} / '
              f'{target_size_km.max():.3f} (variable-resolution)')
    print(f'  min/max: {cell_size_km.min():.3f} / {cell_size_km.max():.3f}')
    print(f'  mean:    {cell_size_km.mean():.3f}')
    print(f'  p01/p05: {percentiles[0]:.3f} / {percentiles[1]:.3f}')
    print(f'  p25/p50: {percentiles[2]:.3f} / {percentiles[3]:.3f}')
    print(f'  p75/p95: {percentiles[4]:.3f} / {percentiles[5]:.3f}')
    print(f'  p99:     {percentiles[6]:.3f}')
    print(f'  Distribution ({width:.3g} km bins):')
    for start, end, count in zip(bin_edges[:-1], bin_edges[1:], counts):
        pct = 100.0 * count / cell_size_km.size
        print(f'    [{start:7.3f}, {end:7.3f}): {count:>10d}  ({pct:6.2f}%)')


def check_dc_edge_distribution(dcEdge, nominalMinDc, meshDensity,
                               bin_width_km=None):
    """Report cell-center-to-cell-center spacing from dcEdge."""
    edge_size_km = np.asarray(dcEdge) * r_earth / 1000.0
    percentiles = np.percentile(edge_size_km, [1, 5, 25, 50, 75, 95, 99])
    target_size_km = (r_earth * nominalMinDc *
                      np.power(1.0 / meshDensity, 0.25) / 1000.0)
    is_quasi_uniform = np.allclose(target_size_km, target_size_km[0],
                                   rtol=1.0e-10, atol=1.0e-10)

    if is_quasi_uniform:
        if bin_width_km is None:
            bin_width_km = nice_bin_width(0.1 * target_size_km[0])
        bin_edges = cell_size_bin_edges(
            edge_size_km, bin_width_km, bin_center_km=target_size_km[0])
    else:
        bin_edges = cell_size_bin_edges(edge_size_km, bin_width_km)

    counts, bin_edges = np.histogram(edge_size_km, bins=bin_edges)
    width = bin_edges[1] - bin_edges[0]

    print('')
    print('Cell-center spacing, dcEdge (km):')
    print(f'  min/max: {edge_size_km.min():.3f} / {edge_size_km.max():.3f}')
    print(f'  mean:    {edge_size_km.mean():.3f}')
    print(f'  p01/p05: {percentiles[0]:.3f} / {percentiles[1]:.3f}')
    print(f'  p25/p50: {percentiles[2]:.3f} / {percentiles[3]:.3f}')
    print(f'  p75/p95: {percentiles[4]:.3f} / {percentiles[5]:.3f}')
    print(f'  p99:     {percentiles[6]:.3f}')
    print(f'  Distribution ({width:.3g} km bins):')
    for start, end, count in zip(bin_edges[:-1], bin_edges[1:], counts):
        pct = 100.0 * count / edge_size_km.size
        print(f'    [{start:7.3f}, {end:7.3f}): {count:>10d}  ({pct:6.2f}%)')


if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser()
    parser.add_argument('mesh_file', help='the name of the netCDF file with mesh fields')
    parser.add_argument('--cell-size-bin-width', type=float,
                        help='fixed cell-size histogram bin width in km')
    args = parser.parse_args()

    if (args.cell_size_bin_width is not None and
            args.cell_size_bin_width <= 0.0):
        parser.error('--cell-size-bin-width must be positive')

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
    nEdgesOnCell = f.variables['nEdgesOnCell'][:]
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

    check_cell_size_distribution(areaCell, nominalMinDc, meshDensity,
                                 args.cell_size_bin_width)

    check_dc_edge_distribution(dcEdge, nominalMinDc, meshDensity,
                               args.cell_size_bin_width)

    check_cell_types(nEdgesOnCell)

    f.close()
