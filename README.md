## A Simple JIGSAW-based Workflow for MPAS-Atmosphere

This repository accompanies the "Generating meshes for MPAS-Atmosphere" mini-tutorial, presented at the 2026 MPAS/WRF Users Workshop.

### Prerequisites

Before generating meshes with the workflow provided by this repository, you'll first need to install [JIGSAW]([url](https://github.com/dengwirda/jigsaw)). The JIGSAW `README.md` file provides installation guidance, though the following should generally be sufficient:
```
git clone https://github.com/dengwirda/jigsaw.git
cd jigsaw && mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=<WHERE_TO_INSTALL_JIGSAW>
make -j 4 install
```
(Be sure to set the installation directory by replacing `WHERE_TO_INSTALL_JIGSAW` in the `cmake` command above.)

Additionally, you'll also need a Python environment with at least the `numpy` and `netCDF4` packages. If you don't already have a suitable environemnt, setting up a Python virtual environment may be an easiest:
```
python -m venv mpas_mesh
source mpas_mesh/bin/activate
pip install --upgrade pip
pip install numpy netCDF4
```

Finally, in order to compile the `mkgrid` program, which derives MPAS's required mesh geometry and connectivity information from generating points and their triangulation, you will need an MPI implementation (OpenMPI and MPICH are common options) plus the [PnetCDF library]([url](https://parallel-netcdf.github.io/)). After ensuring that the `mpicc` command is in your `PATH`, and the `PNETCDF` environment variable points to the installation path of the PnetCDF library, you can simply run
```
make
```
to build the `mkgrid` program from `mkgrid.c`.

### Generating your first variable-resolution mesh with JIGSAW

Generating a variable-resolution mesh with a 12-km circular refinement region centered at 38 N, 95 W, relaxing to 60-km grid spacing over a distance of 1600 km is accomplished with the following steps.

The `radius` of a refinement region defines its constant, fine-resolution core. Outside that core, the default `transition_profile: linear` gives the most compact transition zone. Use `transition_profile: smoothstep` when zero slope and zero curvature at both transition endpoints are more important than keeping the transition compact.

```yaml
mesh:
  transition_profile: smoothstep
```

Use `transition_scale` to widen either transition profile when a broader graded region is desired. The transition width is chosen to preserve the configured maximum resolution gradient.

1. Run the `create_hfun.py` script to generate an `HFUN.msh` file
2. Run `jigsaw`, specifying `MESH.jig` as its command-line argument, to produce a `MESH.msh` file
3. Run `convert_jigsaw.py` to produce `SaveVertices` and `SaveTriangles` files from the `MESH.msh` file
4. Run `create_density.py` to produce `SaveDensity`, `SaveCode`, and `SaveConfig` files
5. Run `mkgrid`, specifying `12000.0` as its command-line argument, to produce `grid.nc` and `graph.info` files. The generated grid includes the density-function source and YAML mesh configuration as provenance. `mesh_quality.py` uses that embedded configuration to report cell shapes by refinement zone.
