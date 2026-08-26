/*******************************************************************************
 * Copyright 2026 The University Corporation for Atmospheric Research (UCAR)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 ******************************************************************************/

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <strings.h>
#include "pnetcdf.h"
#include "mpi.h"


#define min(A,B) ((A) < (B) ? (A) : (B))
#define max(A,B) ((A) > (B) ? (A) : (B))

struct var_atts {
	const char *varname;
	const char *units;
	const char *long_name;
};

struct var_atts attlist[] = {
	{ "indexToCellID", "-", "Mapping from local array index to global cell ID" },
	{ "latCell", "rad", "Latitude of cells" },
	{ "lonCell", "rad", "Longitude of cells" },
	{ "xCell", "m", "Cartesian x-coordinate of cells" },
	{ "yCell", "m", "Cartesian y-coordinate of cells" },
	{ "zCell", "m", "Cartesian z-coordinate of cells" },
	{ "indexToVertexID", "-", "Mapping from local array index to global vertex ID" },
	{ "latVertex", "rad", "Latitude of vertices" },
	{ "lonVertex", "rad", "Longitude of vertices" },
	{ "xVertex", "m", "Cartesian x-coordinate of vertices" },
	{ "yVertex", "m", "Cartesian y-coordinate of vertices" },
	{ "zVertex", "m", "Cartesian z-coordinate of vertices" },
	{ "indexToEdgeID", "-", "Mapping from local array index to global edge ID" },
	{ "latEdge", "rad", "Latitude of edges" },
	{ "lonEdge", "rad", "Longitude of edges" },
	{ "xEdge", "m", "Cartesian x-coordinate of edges" },
	{ "yEdge", "m", "Cartesian y-coordinate of edges" },
	{ "zEdge", "m", "Cartesian z-coordinate of edges" },
	{ "nEdgesOnCell", "-", "Number of edges forming the boundary of a cell" },
	{ "cellsOnCell", "-", "IDs of cells neighboring a cell" },
	{ "verticesOnCell", "-", "IDs of vertices (corner points) of a cell" },
	{ "edgesOnCell", "-", "IDs of edges forming the boundary of a cell" },
	{ "cellsOnVertex", "-", "IDs of the cells that meet at a vertex" },
	{ "edgesOnVertex", "-", "IDs of the edges that meet at a vertex" },
	{ "cellsOnEdge", "-", "IDs of cells divided by an edge" },
	{ "verticesOnEdge", "-", "IDs of the two vertex endpoints of an edge" },
	{ "nEdgesOnEdge", "-", "Number of edges involved in reconstruction of tangential velocity for an edge" },
	{ "edgesOnEdge", "-", "IDs of edges involved in reconstruction of tangential velocity for an edge" },
	{ "areaCell", "m^2", "Spherical area of a Voronoi cell" },
	{ "areaTriangle", "m^2", "Spherical area of a Delaunay triangle" },
	{ "kiteAreasOnVertex", "m^2", "Intersection areas between primal (Voronoi) and dual (triangular) mesh cells" },
	{ "dvEdge", "m", "Spherical distance between vertex endpoints of an edge" },
	{ "dcEdge", "m", "Spherical distance between cells separated by an edge" },
	{ "angleEdge", "rad", "Angle between local north and the positive tangential direction of an edge" },
	{ "weightsOnEdge", "-", "Weights used in reconstruction of tangential velocity for an edge" },
	{ "meshDensity", "unitless", "Mesh density function (used when generating the mesh) evaluated at a cell" },
	{ "nominalMinDc", "m", "Nominal minimum dcEdge value where meshDensity == 1.0" },
	{ "densityFunctionCode", "-", "Source code for mesh density function" }
};

size_t n_attlist = sizeof(attlist) / sizeof(struct var_atts);

/*
 * Returns the number of lines in the named file with a specified format.
 * The format for lines in the file is given by the format argument, and the
 * number of fields to be read via the format is given by the nfields argument.
 * The remaining arguments must be nfields in number and correspond to variables
 * into which the fields can be read for each line in the file.
 */
size_t formatted_line_count(const char *filename, const char *format, int nfields, ...)
{
	FILE *fd = NULL;
	size_t nlines = 0;
	int nconv = 0;
	va_list args;

	va_start(args, nfields);

	fd = fopen(filename, "r");
	if (fd == NULL) return 0;

	while (nconv != EOF) {
		nconv = vfscanf(fd, format, args);
		if (nconv == nfields) {
			nlines++;
		} else if (nconv != EOF) {
			fprintf(stderr, "Invalid input on line %zi of %s\n", (nlines + 1), filename);
			fclose(fd);
			va_end(args);
			return 0;
		}
	}

	fclose(fd);
	va_end(args);

	return nlines;
}


/*
 * Returns the number of cell locations in the named file.
 * The file must be a plain text file with one cell location per line, given as
 * the Cartesian coordinates of the cell location.
 */
size_t cell_count(const char *filename)
{
	double x, y, z;

	return formatted_line_count(filename, "%lf %lf %lf\n", 3, &x, &y, &z);
}


/*
 * Returns the number of triangles in the named file.
 * The file must be a plain text file with one triangle per line, given as the
 * triplet of cell indices that form the corner points of the triangle.
 */
size_t vertex_count(const char *filename)
{
	long i, j, k;

	return formatted_line_count(filename, "%li %li %li\n", 3, &i, &j, &k);
}


struct point {
	double x;
	double y;
	double z;
};

/*
 * Normalize a vector in R^3 given by its components, returning the magnitude of
 * the vector
 */
double normalize_vect_r3(double *x, double *y, double *z)
{
	double m = sqrt((*x) * (*x) + (*y) * (*y) + (*z) * (*z));

	(*x) /= m;
	(*y) /= m;
	(*z) /= m;

	return m;
}


/*
 * Returns the result of (a) - (b) as a point structure
 */
struct point subtract(struct point a, struct point b)
{
	struct point c;

	c.x = a.x - b.x;
	c.y = a.y - b.y;
	c.z = a.z - b.z;

	return c;
}


/*
 * Returns the cross product of (u) and (v) as a point structure
 * Here, (u) and (v) are assumed to be vectors from the origin to
 * the respective point.
 * If the 'normalize' argument is set to 1, the result is made
 * to have a magnitude of 1.0.
 */
struct point cross_product(struct point u, struct point v, int normalize)
{
	struct point w;
	long double u1, u2, u3;
	long double v1, v2, v3;
	long double w1, w2, w3;
	long double mag;

	u1 = (long double)u.x;
	u2 = (long double)u.y;
	u3 = (long double)u.z;

	v1 = (long double)v.x;
	v2 = (long double)v.y;
	v3 = (long double)v.z;

	w1 = u2*v3 - u3*v2;
	w2 = u3*v1 - u1*v3;
	w3 = u1*v2 - u2*v1;

	if (normalize == 1)  {
		mag = sqrtl(w1*w1 + w2*w2 + w3*w3);
		w1 /= mag;
		w2 /= mag;
		w3 /= mag;
	}

	w.x = (double)w1;
	w.y = (double)w2;
	w.z = (double)w3;

	return w;
}


/*
 * Returns the dot product of (a) and (b). Here, (a) and (b) are
 * assumed to be vectors from the origin to the respective point.
 */
double dot_product(struct point a, struct point b)
{
	return (a.x * b.x + a.y * b.y + a.z * b.z);
}


/*
 * Returns the 'magnitude' of a point, treating that point as a vector
 * from the origin to the coordinates stored in the point.
 */
double magnitude(struct point a)
{
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}


/*
 * Computes the circumcenter of the spherical triangle with corners at a, b, and c.
 * It is assumed that a, b, and c are in CCW order.
 */
struct point circumcenter(struct point a, struct point b, struct point c)
{
	struct point u;
	struct point v;

	u = subtract(b, a);
	v = subtract(c, a);

	return cross_product(u, v, 1);
}


/*
 * Swaps two points
 */
void swap_points(struct point *a, struct point *b)
{
	double temp;

	temp = a->x;
	a->x = b->x;
	b->x = temp;

	temp = a->y;
	a->y = b->y;
	b->y = temp;

	temp = a->z;
	a->z = b->z;
	b->z = temp;
}



/*
 * Permutes three input points so that they are in CCW order. If the input
 * points were not in CCW order, this function returns 1, and otherwise it
 * returns 0.
 */
int ccw_triangle(struct point *a, struct point *b, struct point *c)
{
	struct point u;
	struct point v;
	struct point w;

	u = subtract(*b, *a);
	v = subtract(*c, *a);

	w = cross_product(u, v, 1);

	if (dot_product(*a, w) < (double)0.0) {
		swap_points(b, c);
		return 1;
	}
	return 0;
}


/*
 * Returns the distance between two points
 */
long double plane_distance(struct point a, struct point b)
{
	long double dx, dy, dz;

	dx = (long double)b.x - (long double)a.x;
	dy = (long double)b.y - (long double)a.y;
	dz = (long double)b.z - (long double)a.z;

	return sqrtl(dx * dx + dy * dy + dz * dz);
}


/*
 * Returns the squared distance between two points
 */
long double plane_distance_sq(struct point a, struct point b)
{
	long double dx, dy, dz;

	dx = (long double)b.x - (long double)a.x;
	dy = (long double)b.y - (long double)a.y;
	dz = (long double)b.z - (long double)a.z;

	return dx * dx + dy * dy + dz * dz;
}


/*
 * Returns the length of the great circle arc from (a) to (b) on the unit sphere.
 */
long double sphere_distance(struct point a, struct point b)
{
	long double r, c;

	r = (long double)1.0;   /* Assume a unit sphere */

	c = plane_distance(a, b);
	return (long double)2.0 * r * asinl(c / ((long double)2.0 * r));
}


/*
 * Computes the angle between arcs AB and AC, given points A, B, and C
 * Equation numbers w.r.t. http://mathworld.wolfram.com/SphericalTrigonometry.html
 */
double sphere_angle(struct point A, struct point B, struct point C)
{
	long double a, b, c;		/* Side lengths of spherical triangle ABC */

	struct point AB, AC;		/* The vectors AB and AC */

	struct point D;			/* The cross product AB x AC */

	long double s;			/* Semiperimeter of the triangle */
	long double sin_angle;

	a = sphere_distance(B, C);
	b = sphere_distance(A, C);
	c = sphere_distance(A, B);

	AB = subtract(B, A);
	AC = subtract(C, A);

	D = cross_product(AB, AC, 0);

	s = (long double)0.5 * (a + b + c);
	sin_angle = sqrt(min((long double)1.0,max((long double)0.0,(sinl(s-b)*sinl(s-c))/(sinl(b)*sinl(c))))); /* Eqn. (28) */

	if (dot_product(D, A) >= (double)0.0) {
		return (double)2.0 * (double)asinl(max(min(sin_angle,(long double)1.0),(long double)-1.0));
	}
	else {
		return (double)-2.0 * (double)asinl(max(min(sin_angle,(long double)1.0),(long double)-1.0));
	}
}


/*
 * Given the coordinates of the corners of a triangle, returns the area of the triangle
 */
double triangle_area(struct point A, struct point B, struct point C)
{
	double a, b, c, s, tanqe;

	/* Compute lengths of sides of triangle */
	a = (double)sphere_distance(B, C);
	b = (double)sphere_distance(A, C);
	c = (double)sphere_distance(A, B);

	/* Compute semi-perimeter */
	s = (double)0.5 * (a+b+c);

	tanqe = sqrt(tan((double)0.5*s)*tan((double)0.5*(s-a))*tan((double)0.5*(s-b))*tan((double)0.5*(s-c)));

	return (double)4.0 * atan(tanqe);
}


/*
 * Returns the angle between vectors AB and AC in the range [0,2*pi).
 */
double plane_angle(struct point a, struct point b, struct point c)
{
	long double ab2, ac2, bc2, cosc, theta;

	ab2 = plane_distance_sq(a, b);
	ac2 = plane_distance_sq(a, c);
	bc2 = plane_distance_sq(b, c);

	cosc = bc2 - ac2 - ab2 + (long double)2.0 * sqrtl(ab2) * sqrtl(ac2);

	theta = acosl(cosc);

	return (double)theta;
}


int get_var_piecewise(int ncid, int varid, int ndims, int typesize, MPI_Offset *start, MPI_Offset *count, void *buf)
{
	int ncerr;
	int i;
	MPI_Offset *chunk_start;
	MPI_Offset *chunk_count;
	MPI_Offset chunk_read = (MPI_Offset)0;
	size_t outer_dim_lens = 1;
	const MPI_Offset chunk_size = (MPI_Offset)1000000;

	chunk_start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);
	chunk_count = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);

	for (i = 1; i < ndims; i++) {
		chunk_start[i] = start[i];
		chunk_count[i] = count[i];
		outer_dim_lens *= (size_t)chunk_count[i];
	}
	chunk_start[0] = start[0];
	chunk_count[0] = min(count[0], chunk_size);

	while (chunk_read < count[0]) {
		ncerr = ncmpi_get_vara(ncid, varid, chunk_start, chunk_count, (uint8_t *)buf + ((size_t)typesize * (size_t)chunk_read*outer_dim_lens), chunk_count[0]*(MPI_Offset)outer_dim_lens, MPI_DATATYPE_NULL);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading varid %i: %s\n", varid, ncmpi_strerror(ncerr));
			return ncerr;
		}
		chunk_read += chunk_count[0];
		chunk_start[0] = start[0] + chunk_read;
		chunk_count[0] = min(count[0]-chunk_read, chunk_size);
	}

	free(chunk_start);
	free(chunk_count);

	return ncerr;
}


int put_var_piecewise(int ncid, int varid, int ndims, int typesize, MPI_Offset *start, MPI_Offset *count, void *buf)
{
	int ncerr;
	int i;
	MPI_Offset *chunk_start;
	MPI_Offset *chunk_count;
	MPI_Offset chunk_written = (MPI_Offset)0;
	size_t outer_dim_lens = 1;
	const MPI_Offset chunk_size = (MPI_Offset)1000000;

	chunk_start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);
	chunk_count = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);

	for (i = 1; i < ndims; i++) {
		chunk_start[i] = start[i];
		chunk_count[i] = count[i];
		outer_dim_lens *= (size_t)chunk_count[i];
	}
	chunk_start[0] = start[0];
	chunk_count[0] = min(count[0], chunk_size);

	while (chunk_written < count[0]) {
		ncerr = ncmpi_put_vara(ncid, varid, chunk_start, chunk_count, (uint8_t *)buf + ((size_t)typesize * (size_t)chunk_written*outer_dim_lens), chunk_count[0]*(MPI_Offset)outer_dim_lens, MPI_DATATYPE_NULL);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing varid %i: %s\n", varid, ncmpi_strerror(ncerr));
			return ncerr;
		}
		chunk_written += chunk_count[0];
		chunk_start[0] = start[0] + chunk_written;
		chunk_count[0] = min(count[0]-chunk_written, chunk_size);
	}

	free(chunk_start);
	free(chunk_count);

	return ncerr;
}


int get_var(int ncid, const char *varname, void **buf, MPI_Offset *start_in, MPI_Offset *count_in)
{
	int ncerr;
	int varid;
	int ndims;
	nc_type vartype;
	int typesize;
	int *dimids;
	MPI_Offset *start;
	MPI_Offset *count;
	MPI_Offset total_size;
	int i;

	ncerr = ncmpi_inq_varid(ncid, varname, &varid);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error inquiring variable ID for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}


	/*
	 * Inquire about the variable dimensions
	 */
	ncerr = ncmpi_inq_varndims(ncid, varid, &ndims);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error inquiring number of dims for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	dimids = (int *)malloc(sizeof(int) * (size_t)ndims);
	start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);
	count = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);

	ncerr = ncmpi_inq_vardimid(ncid, varid, dimids);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		fprintf(stderr, "Error inquiring dim IDs for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	total_size = (MPI_Offset)1;
	for (i=0; i<ndims; i++) {
		if (count_in == NULL) {
			ncerr = ncmpi_inq_dimlen(ncid, dimids[i], &count[i]);
			if (ncerr != NC_NOERR) {
				free(dimids);
				free(start);
				free(count);
				fprintf(stderr, "Error inquiring dim length: %s\n", ncmpi_strerror(ncerr));
				return ncerr;
			}
		}
		else {
			count[i] = count_in[i];
		}

		if (start_in == NULL) {
			start[i] = (MPI_Offset)0;
		}
		else {
			start[i] = start_in[i];
		}
		total_size *= count[i];
	}


	/*
	 * Inquire about the variable type
	 */
	ncerr = ncmpi_inq_vartype(ncid, varid, &vartype);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		fprintf(stderr, "Error inquiring variable type for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	if (vartype == NC_INT) {
		typesize = sizeof(int);
	}
	else if (vartype == NC_FLOAT) {
		typesize = sizeof(float);
	}
	else if (vartype == NC_DOUBLE) {
		typesize = sizeof(double);
	}
	else if (vartype == NC_CHAR) {
		typesize = sizeof(char);
	}
	else {
		fprintf(stderr, "Unsupported variable type in get_var\n");
		free(dimids);
		free(start);
		free(count);
		return -1;
	}


	/*
	 * Allocate storage for the variable
	 */
	*buf = malloc(typesize * (size_t)total_size);

	ncerr = get_var_piecewise(ncid, varid, ndims, typesize, start, count, *buf);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		free(*buf);
		return ncerr;
	}

	free(dimids);
	free(start);
	free(count);

	return 0;
}


int put_var(int ncid, const char *varname, void *buf, MPI_Offset *start_in, MPI_Offset *count_in)
{
	int ncerr;
	int varid;
	int ndims;
	nc_type vartype;
	int typesize;
	int *dimids;
	MPI_Offset *start;
	MPI_Offset *count;
	int i;

	ncerr = ncmpi_inq_varid(ncid, varname, &varid);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error inquiring variable ID for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}


	/*
	 * Inquire about the variable dimensions
	 */
	ncerr = ncmpi_inq_varndims(ncid, varid, &ndims);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error inquiring number of dims for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	dimids = (int *)malloc(sizeof(int) * (size_t)ndims);
	start = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);
	count = (MPI_Offset *)malloc(sizeof(MPI_Offset) * (size_t)ndims);

	ncerr = ncmpi_inq_vardimid(ncid, varid, dimids);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		fprintf(stderr, "Error inquiring dim IDs for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	for (i=0; i<ndims; i++) {
		if (count_in == NULL) {
			ncerr = ncmpi_inq_dimlen(ncid, dimids[i], &count[i]);
			if (ncerr != NC_NOERR) {
				free(dimids);
				free(start);
				free(count);
				fprintf(stderr, "Error inquiring dim length: %s\n", ncmpi_strerror(ncerr));
				return ncerr;
			}
		}
		else {
			count[i] = count_in[i];
		}

		if (start_in == NULL) {
			start[i] = (MPI_Offset)0;
		}
		else {
			start[i] = start_in[i];
		}
	}


	/*
	 * Inquire about the variable type
	 */
	ncerr = ncmpi_inq_vartype(ncid, varid, &vartype);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		fprintf(stderr, "Error inquiring variable type for %s: %s\n", varname, ncmpi_strerror(ncerr));
		return ncerr;
	}

	if (vartype == NC_INT) {
		typesize = sizeof(int);
	}
	else if (vartype == NC_FLOAT) {
		typesize = sizeof(float);
	}
	else if (vartype == NC_DOUBLE) {
		typesize = sizeof(double);
	}
	else if (vartype == NC_CHAR) {
		typesize = sizeof(char);
	}
	else {
		fprintf(stderr, "Unsupported variable type in get_var\n");
		free(dimids);
		free(start);
		free(count);
		return -1;
	}


	ncerr = put_var_piecewise(ncid, varid, ndims, typesize, start, count, buf);
	if (ncerr != NC_NOERR) {
		free(dimids);
		free(start);
		free(count);
		free(buf);
		return ncerr;
	}

	free(dimids);
	free(start);
	free(count);

	return 0;
}


/*
 * Find index of neighbor in the cellsOnCell list of a cell.
 * Note that cellsOnCell is a 1-d array dimensioned by maxEdges, and
 * represents the neighbors for just the cell in question. The neighbor
 * index, as well as the returned index, are both 1-based.
 */
int find_neighbor(int *cellsOnCell, MPI_Offset maxEdges, int neighbor)
{
	int i;

	for (i = 0; i < maxEdges; i++) {
		if (cellsOnCell[i] == neighbor) {
			return i+1;
		}
	}
	return 0;
}


int found_cell(int *cellsOnCell, MPI_Offset maxEdges, int cell)
{
	int i;

	for (i = 0; i < maxEdges; i++) {
		if (cellsOnCell[i] == cell) {
			return 1;
		}
	}
	return 0;
}


void handle_netcdf_error(int ncerr, const char *errstr)
{
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "%s: %s\n", errstr, ncmpi_strerror(ncerr));
	}
}


void add_variable_attributes(int ncid, int varid, const char *varname)
{
	int ncerr;

	for(size_t i = 0; i < n_attlist; i++) {
		if (!strcmp(attlist[i].varname, varname)) {
			ncerr = ncmpi_put_att_text(ncid, varid, "units", strlen(attlist[i].units), attlist[i].units);
			if (ncerr != NC_NOERR) {
				fprintf(stderr, "Error writing units attribute for variable %s\n", varname);
			}
			ncerr = ncmpi_put_att_text(ncid, varid, "long_name", strlen(attlist[i].long_name), attlist[i].long_name);
			if (ncerr != NC_NOERR) {
				fprintf(stderr, "Error writing long_name attribute for variable %s\n", varname);
			}
		}
	}
}


void write_netcdf(MPI_Offset nCells, MPI_Offset nEdges, MPI_Offset nVertices, MPI_Offset maxEdges, MPI_Offset vertexDegree,
                  double * xCell, double * yCell, double * zCell,
                  double * xVertex, double * yVertex, double * zVertex,
                  int * nEdgesOnCell, int * cellsOnCell, int * verticesOnCell, int * cellsOnVertex,
                  long density_len, char *density_code,
                  double min_dc_m
                 )
{
	int ncerr;
	int ncid;
	int dimIDnCells, dimIDnEdges, dimIDnVertices, dimIDmaxEdges, dimIDmaxEdges2, dimIDvertexDegree, dimIDTWO, dimIDlen_code;
	int varIDindexToCellID;
	int varIDxCell, varIDyCell, varIDzCell;
	int varIDlatCell, varIDlonCell;
	int varIDindexToEdgeID;
	int varIDxEdge, varIDyEdge, varIDzEdge;
	int varIDlatEdge, varIDlonEdge;
	int varIDindexToVertexID;
	int varIDxVertex, varIDyVertex, varIDzVertex;
	int varIDlatVertex, varIDlonVertex;
	int varIDnEdgesOnCell, varIDnEdgesOnEdge;
	int varIDcellsOnCell, varIDedgesOnCell, varIDverticesOnCell;
	int varIDcellsOnEdge, varIDverticesOnEdge, varIDedgesOnEdge, varIDweightsOnEdge;
	int varIDedgesOnVertex, varIDcellsOnVertex, varIDkiteAreasOnVertex;
	int varIDdvEdge, varIDdcEdge, varIDareaCell, varIDareaTriangle, varIDangleEdge;
	int varIDmeshDensity;
	int varIDnominalMinDc;
	int varIDdensity_code;

	int dimids1[1];
	int dimids2[2];
	MPI_Offset start1[1];
	MPI_Offset start2[2];
	MPI_Offset count1[1];
	MPI_Offset count2[2];

	size_t max_var_size;
	int format;

	double sphere_radius = 1.0;


	start1[0] = (MPI_Offset)0;
	start2[0] = (MPI_Offset)0;
	start2[1] = (MPI_Offset)0;

	/*
	 * The largest variable to be written is currently weightsOnEdge, and so
	 * the product of its dimension sizes can be used to determine whether
	 * the mesh should be written to a CDF-5 or a CDF-2 file.
	 */
	max_var_size = sizeof(double) * (size_t)nEdges * (size_t)maxEdges * 2ul;
	fprintf(stderr, "Maximum variable size is %lu bytes\n", max_var_size);
	if (max_var_size > INT_MAX) {
		format = NC_64BIT_DATA;
		fprintf(stderr, "Creating grid.nc in CDF-5 format\n");
	} else {
		format = NC_64BIT_OFFSET;
		fprintf(stderr, "Creating grid.nc in CDF-2 format\n");
	}


	ncerr = ncmpi_create(MPI_COMM_WORLD, "grid.nc", NC_WRITE|format, MPI_INFO_NULL, &ncid);
	handle_netcdf_error(ncerr, "Error creating grid.nc");

	/*
	 * Dimensions
	 */
	ncerr = ncmpi_def_dim(ncid, "nCells", nCells, &dimIDnCells);
	handle_netcdf_error(ncerr, "Error defining dimension nCells");
	ncerr = ncmpi_def_dim(ncid, "nVertices", nVertices, &dimIDnVertices);
	handle_netcdf_error(ncerr, "Error defining dimension nVertices");
	ncerr = ncmpi_def_dim(ncid, "nEdges", nEdges, &dimIDnEdges);
	handle_netcdf_error(ncerr, "Error defining dimension nEdges");
	ncerr = ncmpi_def_dim(ncid, "maxEdges", maxEdges, &dimIDmaxEdges);
	handle_netcdf_error(ncerr, "Error defining dimension maxEdges");
	ncerr = ncmpi_def_dim(ncid, "maxEdges2", ((MPI_Offset)2*maxEdges), &dimIDmaxEdges2);
	handle_netcdf_error(ncerr, "Error defining dimension maxEdges2");
	ncerr = ncmpi_def_dim(ncid, "vertexDegree", (MPI_Offset)vertexDegree, &dimIDvertexDegree);
	handle_netcdf_error(ncerr, "Error defining dimension vertexDegree");
	ncerr = ncmpi_def_dim(ncid, "TWO", (MPI_Offset)2, &dimIDTWO);
	handle_netcdf_error(ncerr, "Error defining dimension TWO");
	if (density_len > 0 && density_code != NULL) {
		ncerr = ncmpi_def_dim(ncid, "len_code", (MPI_Offset)density_len, &dimIDlen_code);
		handle_netcdf_error(ncerr, "Error defining dimension len_code");
	}

	/*
	 * Variables
	 */
	dimids1[0] = dimIDnCells;
	ncerr = ncmpi_def_var(ncid, "indexToCellID", NC_INT, 1, dimids1, &varIDindexToCellID);
	handle_netcdf_error(ncerr, "Error defining variable indexToCellID");
	add_variable_attributes(ncid, varIDindexToCellID, "indexToCellID");
	ncerr = ncmpi_def_var(ncid, "xCell", NC_DOUBLE, 1, dimids1, &varIDxCell);
	handle_netcdf_error(ncerr, "Error defining variable xCell");
	add_variable_attributes(ncid, varIDxCell, "xCell");
	ncerr = ncmpi_def_var(ncid, "yCell", NC_DOUBLE, 1, dimids1, &varIDyCell);
	handle_netcdf_error(ncerr, "Error defining variable yCell");
	add_variable_attributes(ncid, varIDyCell, "yCell");
	ncerr = ncmpi_def_var(ncid, "zCell", NC_DOUBLE, 1, dimids1, &varIDzCell);
	handle_netcdf_error(ncerr, "Error defining variable zCell");
	add_variable_attributes(ncid, varIDzCell, "zCell");
	ncerr = ncmpi_def_var(ncid, "latCell", NC_DOUBLE, 1, dimids1, &varIDlatCell);
	handle_netcdf_error(ncerr, "Error defining variable latCell");
	add_variable_attributes(ncid, varIDlatCell, "latCell");
	ncerr = ncmpi_def_var(ncid, "lonCell", NC_DOUBLE, 1, dimids1, &varIDlonCell);
	handle_netcdf_error(ncerr, "Error defining variable lonCell");
	add_variable_attributes(ncid, varIDlonCell, "lonCell");
	ncerr = ncmpi_def_var(ncid, "nEdgesOnCell", NC_INT, 1, dimids1, &varIDnEdgesOnCell);
	handle_netcdf_error(ncerr, "Error defining variable nEdgesOnCell");
	add_variable_attributes(ncid, varIDnEdgesOnCell, "nEdgesOnCell");
	ncerr = ncmpi_def_var(ncid, "areaCell", NC_DOUBLE, 1, dimids1, &varIDareaCell);
	handle_netcdf_error(ncerr, "Error defining variable areaCell");
	add_variable_attributes(ncid, varIDareaCell, "areaCell");
	dimids1[0] = dimIDnEdges;
	ncerr = ncmpi_def_var(ncid, "indexToEdgeID", NC_INT, 1, dimids1, &varIDindexToEdgeID);
	handle_netcdf_error(ncerr, "Error defining variable indexToEdgeID");
	add_variable_attributes(ncid, varIDindexToEdgeID, "indexToEdgeID");
	ncerr = ncmpi_def_var(ncid, "xEdge", NC_DOUBLE, 1, dimids1, &varIDxEdge);
	handle_netcdf_error(ncerr, "Error defining variable xEdge");
	add_variable_attributes(ncid, varIDxEdge, "xEdge");
	ncerr = ncmpi_def_var(ncid, "yEdge", NC_DOUBLE, 1, dimids1, &varIDyEdge);
	handle_netcdf_error(ncerr, "Error defining variable yEdge");
	add_variable_attributes(ncid, varIDyEdge, "yEdge");
	ncerr = ncmpi_def_var(ncid, "zEdge", NC_DOUBLE, 1, dimids1, &varIDzEdge);
	handle_netcdf_error(ncerr, "Error defining variable zEdge");
	add_variable_attributes(ncid, varIDzEdge, "zEdge");
	ncerr = ncmpi_def_var(ncid, "latEdge", NC_DOUBLE, 1, dimids1, &varIDlatEdge);
	handle_netcdf_error(ncerr, "Error defining variable latEdge");
	add_variable_attributes(ncid, varIDlatEdge, "latEdge");
	ncerr = ncmpi_def_var(ncid, "lonEdge", NC_DOUBLE, 1, dimids1, &varIDlonEdge);
	handle_netcdf_error(ncerr, "Error defining variable lonEdge");
	add_variable_attributes(ncid, varIDlonEdge, "lonEdge");
	ncerr = ncmpi_def_var(ncid, "nEdgesOnEdge", NC_INT, 1, dimids1, &varIDnEdgesOnEdge);
	handle_netcdf_error(ncerr, "Error defining variable nEdgesOnEdge");
	add_variable_attributes(ncid, varIDnEdgesOnEdge, "nEdgesOnEdge");
	ncerr = ncmpi_def_var(ncid, "dvEdge", NC_DOUBLE, 1, dimids1, &varIDdvEdge);
	handle_netcdf_error(ncerr, "Error defining variable dvEdge");
	add_variable_attributes(ncid, varIDdvEdge, "dvEdge");
	ncerr = ncmpi_def_var(ncid, "dcEdge", NC_DOUBLE, 1, dimids1, &varIDdcEdge);
	handle_netcdf_error(ncerr, "Error defining variable dcEdge");
	add_variable_attributes(ncid, varIDdcEdge, "dcEdge");
	ncerr = ncmpi_def_var(ncid, "angleEdge", NC_DOUBLE, 1, dimids1, &varIDangleEdge);
	handle_netcdf_error(ncerr, "Error defining variable angleEdge");
	add_variable_attributes(ncid, varIDangleEdge, "angleEdge");
	dimids1[0] = dimIDnVertices;
	ncerr = ncmpi_def_var(ncid, "indexToVertexID", NC_INT, 1, dimids1, &varIDindexToVertexID);
	handle_netcdf_error(ncerr, "Error defining variable indexToVertexID");
	add_variable_attributes(ncid, varIDindexToVertexID, "indexToVertexID");
	ncerr = ncmpi_def_var(ncid, "xVertex", NC_DOUBLE, 1, dimids1, &varIDxVertex);
	handle_netcdf_error(ncerr, "Error defining variable xVertex");
	add_variable_attributes(ncid, varIDxVertex, "xVertex");
	ncerr = ncmpi_def_var(ncid, "yVertex", NC_DOUBLE, 1, dimids1, &varIDyVertex);
	handle_netcdf_error(ncerr, "Error defining variable yVertex");
	add_variable_attributes(ncid, varIDyVertex, "yVertex");
	ncerr = ncmpi_def_var(ncid, "zVertex", NC_DOUBLE, 1, dimids1, &varIDzVertex);
	handle_netcdf_error(ncerr, "Error defining variable zVertex");
	add_variable_attributes(ncid, varIDzVertex, "zVertex");
	ncerr = ncmpi_def_var(ncid, "latVertex", NC_DOUBLE, 1, dimids1, &varIDlatVertex);
	handle_netcdf_error(ncerr, "Error defining variable latVertex");
	add_variable_attributes(ncid, varIDlatVertex, "latVertex");
	ncerr = ncmpi_def_var(ncid, "lonVertex", NC_DOUBLE, 1, dimids1, &varIDlonVertex);
	handle_netcdf_error(ncerr, "Error defining variable lonVertex");
	add_variable_attributes(ncid, varIDlonVertex, "lonVertex");
	ncerr = ncmpi_def_var(ncid, "areaTriangle", NC_DOUBLE, 1, dimids1, &varIDareaTriangle);
	handle_netcdf_error(ncerr, "Error defining variable areaTriangle");
	add_variable_attributes(ncid, varIDareaTriangle, "areaTriangle");
	dimids2[0] = dimIDnCells;
	dimids2[1] = dimIDmaxEdges;
	ncerr = ncmpi_def_var(ncid, "cellsOnCell", NC_INT, 2, dimids2, &varIDcellsOnCell);
	handle_netcdf_error(ncerr, "Error defining variable cellsOnCell");
	add_variable_attributes(ncid, varIDcellsOnCell, "cellsOnCell");
	ncerr = ncmpi_def_var(ncid, "edgesOnCell", NC_INT, 2, dimids2, &varIDedgesOnCell);
	handle_netcdf_error(ncerr, "Error defining variable edgesOnCell");
	add_variable_attributes(ncid, varIDedgesOnCell, "edgesOnCell");
	ncerr = ncmpi_def_var(ncid, "verticesOnCell", NC_INT, 2, dimids2, &varIDverticesOnCell);
	handle_netcdf_error(ncerr, "Error defining variable verticesOnCell");
	add_variable_attributes(ncid, varIDverticesOnCell, "verticesOnCell");
	dimids2[0] = dimIDnEdges;
	dimids2[1] = dimIDTWO;
	ncerr = ncmpi_def_var(ncid, "cellsOnEdge", NC_INT, 2, dimids2, &varIDcellsOnEdge);
	handle_netcdf_error(ncerr, "Error defining variable cellsOnEdge");
	add_variable_attributes(ncid, varIDcellsOnEdge, "cellsOnEdge");
	ncerr = ncmpi_def_var(ncid, "verticesOnEdge", NC_INT, 2, dimids2, &varIDverticesOnEdge);
	handle_netcdf_error(ncerr, "Error defining variable verticesOnEdge");
	add_variable_attributes(ncid, varIDverticesOnEdge, "verticesOnEdge");
	dimids2[0] = dimIDnEdges;
	dimids2[1] = dimIDmaxEdges2;
	ncerr = ncmpi_def_var(ncid, "edgesOnEdge", NC_INT, 2, dimids2, &varIDedgesOnEdge);
	handle_netcdf_error(ncerr, "Error defining variable edgesOnEdge");
	add_variable_attributes(ncid, varIDedgesOnEdge, "edgesOnEdge");
	ncerr = ncmpi_def_var(ncid, "weightsOnEdge", NC_DOUBLE, 2, dimids2, &varIDweightsOnEdge);
	handle_netcdf_error(ncerr, "Error defining variable weightsOnEdge");
	add_variable_attributes(ncid, varIDweightsOnEdge, "weightsOnEdge");
	dimids2[0] = dimIDnVertices;
	dimids2[1] = dimIDvertexDegree;
	ncerr = ncmpi_def_var(ncid, "edgesOnVertex", NC_INT, 2, dimids2, &varIDedgesOnVertex);
	handle_netcdf_error(ncerr, "Error defining variable edgesOnVertex");
	add_variable_attributes(ncid, varIDedgesOnVertex, "edgesOnVertex");
	ncerr = ncmpi_def_var(ncid, "cellsOnVertex", NC_INT, 2, dimids2, &varIDcellsOnVertex);
	handle_netcdf_error(ncerr, "Error defining variable cellsOnVertex");
	add_variable_attributes(ncid, varIDcellsOnVertex, "cellsOnVertex");
	ncerr = ncmpi_def_var(ncid, "kiteAreasOnVertex", NC_DOUBLE, 2, dimids2, &varIDkiteAreasOnVertex);
	handle_netcdf_error(ncerr, "Error defining variable kiteAreasOnVertex");
	add_variable_attributes(ncid, varIDkiteAreasOnVertex, "kiteAreasOnVertex");
	dimids1[0] = dimIDnCells;
	ncerr = ncmpi_def_var(ncid, "meshDensity", NC_DOUBLE, 1, dimids1, &varIDmeshDensity);
	handle_netcdf_error(ncerr, "Error defining variable meshDensity");
	add_variable_attributes(ncid, varIDmeshDensity, "meshDensity");
	ncerr = ncmpi_def_var(ncid, "nominalMinDc", NC_DOUBLE, 0, NULL, &varIDnominalMinDc);
	handle_netcdf_error(ncerr, "Error defining variable nominalMinDc");
	add_variable_attributes(ncid, varIDnominalMinDc, "nominalMinDc");
	if (density_len > 0 && density_code != NULL) {
		dimids1[0] = dimIDlen_code;
		ncerr = ncmpi_def_var(ncid, "densityFunctionCode", NC_CHAR, 1, dimids1, &varIDdensity_code);
		handle_netcdf_error(ncerr, "Error defining variable densityFunctionCode");
		add_variable_attributes(ncid, varIDdensity_code, "densityFunctionCode");
	}

	/*
	 * Global attributes
	 */
	ncerr = ncmpi_put_att_text(ncid, NC_GLOBAL, "on_a_sphere", 16, "YES             ");
	handle_netcdf_error(ncerr, "Error defining attribute on_a_sphere");
	ncerr = ncmpi_put_att_double(ncid, NC_GLOBAL, "sphere_radius", NC_DOUBLE, 1, &sphere_radius);
	handle_netcdf_error(ncerr, "Error defining attribute sphere_radius");


	ncerr = ncmpi_enddef(ncid);
	handle_netcdf_error(ncerr, "Error exiting define mode");
	ncerr = ncmpi_begin_indep_data(ncid);
	handle_netcdf_error(ncerr, "Error beginning independent data mode");


	count1[0] = nCells;
	ncerr = put_var_piecewise(ncid, varIDxCell, 1, sizeof(double), start1, count1, (void *)xCell);
	handle_netcdf_error(ncerr, "Error writing variable xCell");
	ncerr = put_var_piecewise(ncid, varIDyCell, 1, sizeof(double), start1, count1, (void *)yCell);
	handle_netcdf_error(ncerr, "Error writing variable yCell");
	ncerr = put_var_piecewise(ncid, varIDzCell, 1, sizeof(double), start1, count1, (void *)zCell);
	handle_netcdf_error(ncerr, "Error writing variable zCell");
	ncerr = put_var_piecewise(ncid, varIDnEdgesOnCell, 1, sizeof(int), start1, count1, (void *)nEdgesOnCell);
	handle_netcdf_error(ncerr, "Error writing variable nEdgesOnCell");

	count1[0] = nVertices;
	ncerr = put_var_piecewise(ncid, varIDxVertex, 1, sizeof(double), start1, count1, (void *)xVertex);
	handle_netcdf_error(ncerr, "Error writing variable xVertex");
	ncerr = put_var_piecewise(ncid, varIDyVertex, 1, sizeof(double), start1, count1, (void *)yVertex);
	handle_netcdf_error(ncerr, "Error writing variable yVertex");
	ncerr = put_var_piecewise(ncid, varIDzVertex, 1, sizeof(double), start1, count1, (void *)zVertex);
	handle_netcdf_error(ncerr, "Error writing variable zVertex");

	count2[0] = nCells;
	count2[1] = maxEdges;
	ncerr = put_var_piecewise(ncid, varIDcellsOnCell, 2, sizeof(int), start2, count2, (void *)cellsOnCell);
	handle_netcdf_error(ncerr, "Error writing variable cellsOnCell");
	ncerr = put_var_piecewise(ncid, varIDverticesOnCell, 2, sizeof(int), start2, count2, (void *)verticesOnCell);
	handle_netcdf_error(ncerr, "Error writing variable verticesOnCell");

	count2[0] = nVertices;
	count2[1] = vertexDegree;
	ncerr = put_var_piecewise(ncid, varIDcellsOnVertex, 2, sizeof(int), start2, count2, (void *)cellsOnVertex);
	handle_netcdf_error(ncerr, "Error writing variable cellsOnVertex");

	ncerr = ncmpi_put_vara(ncid, varIDnominalMinDc, NULL, NULL,
	                       &min_dc_m, (MPI_Offset)1, MPI_DATATYPE_NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing nominalMinDc: %s\n", ncmpi_strerror(ncerr));
	}

	if (density_len > 0 && density_code != NULL) {
		start1[0] = (MPI_Offset)0;
		count1[0] = (MPI_Offset)density_len;
		ncerr = ncmpi_put_vara(ncid, varIDdensity_code, start1, count1,
		                       density_code, (MPI_Offset)density_len, MPI_DATATYPE_NULL);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing density_function_code: %s\n", ncmpi_strerror(ncerr));
		}
	}

	ncerr = ncmpi_close(ncid);
}


/*
 * Find the point of intersection between arcs (a1,a2) and (b1,b2)
 */
struct point arc_intersect(struct point a1, struct point a2, struct point b1, struct point b2)
{
	struct point n;
	struct point m;
	struct point s;
	double dot;

	n = cross_product(a1, a2, 1);
	m = cross_product(b1, b2, 1);

	s = cross_product(n, m, 1);

	dot = dot_product(a1, s);

	if (dot < (double)0.0) {
		s.x = -s.x;
		s.y = -s.y;
		s.z = -s.z;
	}

	return s;
}


/*
 * Convert (x, y, z) to a (lat, lon) location on a unit sphere
 */
void xyz_to_latlon(double x, double y, double z, double *lat, double *lon)
{
	const double eps = (double)1.0e-10;

	double pii;

	pii = (double)2.0 * asin((double)1.0);

	*lat = asin(z);

	/* check for being close to either pole */
	if (fabs(x) > eps) {

		if (fabs(y) > eps) {
			*lon = atan(fabs(y/x));

			if ((x <= (double)0.0) && (y >= (double)0.0)) {
				*lon = pii - *lon;
			}
			else if ((x <= (double)0.0) && (y < (double)0.0)) {
				*lon = *lon + pii;
			}
			else if ((x >= (double)0.0) && (y <= (double)0.0)) {
				*lon = (double)2.0 * pii - *lon;
			}
		}
		else { /* we're either on longitude 0 or 180 */
			if (x > (double)0.0) {
				*lon = (double)0.0;
			}
			else {
				*lon = pii;
			}
		}

	}
	else if (fabs(y) > eps) {
		if (y > (double)0.0) {
			*lon = (double)0.5 * pii;
		}
		else {
			*lon = (double)1.5 * pii;
		}
	}
	else {  /* we are at a pole */
		*lon = (double)0.0;
	}
}


/*
 * Returns 1.0 if edgeID is an outgoing edge for cellID, and -1.0 if edgeID is an incoming edge
 * Note: Both cellID and edgeID are 1-based IDs.
 */
double sign_for_edge(int cellID, int edgeID, int *cellsOnEdge)
{
	return (cellID == cellsOnEdge[(size_t)(edgeID-1) * (size_t)2 + (size_t)0]) ? (double)1.0 : (double)(-1.0);
}


/*
 * Rotate the kiteAreasOnCell list for a cell so that it begins with the kite immediately following
 * the specified edge in CCW order. The rotated output list, rotated_kiteAreasOnCell, must be allocated
 * with size at least nEdgesOnCell.
 * Note: edgeID is a 1-based ID.
 */
void rotate_kiteAreasOnCell(double *kiteAreasOnCell, int edgeID, int *edgesOnCell, int nEdgesOnCell, double *rotated_kiteAreasOnCell)
{
	int i, j, istart;

	/*
	 * Find starting edge in the edgesOnCell list
	 * NB: verticesOnCell (and therefore, kiteAreasOnCell), "lead" edgesOnCell, so we actually
	 *     need to start with the vertex whose index is one more than the index of the edgeID
	 */
	for (i=0; i<nEdgesOnCell; i++) {
		if (edgesOnCell[i] == edgeID) {
			istart = (i+1) % nEdgesOnCell;
			break;
		}
	}
	if (i == nEdgesOnCell) {
		fprintf(stderr, "Could not find the edge we need...\n");
	}


	/*
	 * Copy kiteAreasOnCell
	 */
	j = 0;
	for (i=istart; i<nEdgesOnCell; i++) {
		rotated_kiteAreasOnCell[j++] = kiteAreasOnCell[i];
	}
	for (i=0; i<istart; i++) {
		rotated_kiteAreasOnCell[j++] = kiteAreasOnCell[i];
	}
}


int main(int argc, char **argv)
{
	int stat;
	FILE *fd;
	size_t k, kk;
	size_t i, j1, j2, j3;
	int j, ii;
	int cell1, cell2;
	int vtx1, vtx2;
	int edge_id;
	size_t ichunk, chunksize;
	MPI_Offset nCells, nEdges, nVertices, maxEdges, maxEdges2, vertexDegree;
	MPI_Offset start[2];
	MPI_Offset count[2];
	struct point a, b, c, vtx;
	struct point c1, c2, v1, v2, e;
	double *angles;
	double minangle;
	int *elements;
	int minelement;
	double v;
	double sum_r;
	double norm_min, norm_max;
	int tri_min, tri_max;

	int *indexToCellID;
	int *indexToEdgeID;
	int *indexToVertexID;
	double *xCell;
	double *yCell;
	double *zCell;
	double *latCell;
	double *lonCell;
	double *xEdge;
	double *yEdge;
	double *zEdge;
	double *latEdge;
	double *lonEdge;
	double *xVertex;
	double *yVertex;
	double *zVertex;
	double *latVertex;
	double *lonVertex;
	int *nEdgesOnCell;
	int *nEdgesOnEdge;
	int *cellsOnCell;
	int *edgesOnCell;
	int *verticesOnCell;
	int *cellsOnVertex;
	int *cellsOnEdge;
	int *verticesOnEdge;
	int *edgesOnEdge;
	int *edgesOnVertex;
	double *kiteAreasOnVertex;
	double *kiteAreasOnCell;  /* Never actually written to grid.nc */
	double *rotated_kiteAreasOnCell;  /* Used when computing weightsOnEdge */
	double *dvEdge;
	double *dcEdge;
	double *areaCell;
	double *areaTriangle;
	double *angleEdge;
	double *weightsOnEdge;
	double *meshDensity;

	int ncid;
	int ncerr;
	int ncells_id;
	int nedges_id;
	int nvertices_id;
	int maxedges_id;
	int maxedges2_id;
	int vertexdegree_id;

	char *density_code;
	long density_len;
	FILE *density_file; 
	double min_dc_m = 0.0;

	double twopi = (double)4.0 * asin((double)1.0);
	const double r_earth = 6371229.0;


	stat = MPI_Init(&argc, &argv);
	if (stat != MPI_SUCCESS) {
		fprintf(stderr, "Error: failed call to MPI_Init\n");
		return 1;
	}

	/*
	 * No command-line arguments are expected.
	 */
	if (argc != 2) {
		fprintf(stderr, "\nUsage: mkgrid <nominalMinDc_meters>\n\n");
		fprintf(stderr, "    where nominalMinDc_meters is the nominal minimum grid distance\n");
		fprintf(stderr, "    assuming a sphere radius of %lf m.\n\n", r_earth);
		return 0;
	}

	/*
	 * Set nominal minimum grid distance assuming a sphere radius of r_earth
	 */
	min_dc_m = strtod(argv[1], NULL);
	fprintf(stderr, "Nominal minimum grid distance (m): %lf\n", min_dc_m);
	min_dc_m /= r_earth;

	nCells = (MPI_Offset)cell_count("SaveVertices");
	nVertices = (MPI_Offset)vertex_count("SaveTriangles");

	if (nVertices != (nCells - (MPI_Offset)2) * (MPI_Offset)2) {
		fprintf(stderr, "The number of cells and vertices are invalid for a spherical mesh.\n");
		fprintf(stderr, "  %li cells found in SaveVertices\n", (long)nCells);
		fprintf(stderr, "  %li vertices found in SaveTriangles\n", (long)nVertices);
		return 1;
	}

	/*
	 * For spherical meshes, the number of edges and vertices are calculated
	 * assuming an Euler characteristic of 2
	 */
	nEdges = (nCells - (MPI_Offset)2) * (MPI_Offset)3;
	vertexDegree = (MPI_Offset)3;

	if (nCells == (MPI_Offset)0) {
		fprintf(stderr, "Error: invalid value for nCells\n");
		return 1;
	}
	else {
		fprintf(stderr, "nCells = %lu\n", (unsigned long)nCells);
	}

	/*
	 * Read SaveCode contents into a buffer
	 */
	density_file = fopen("SaveCode", "r");
	if (density_file != NULL) {
		fprintf(stderr, "Reading SaveCode\n");
		fseek(density_file, 0L, SEEK_END);	
		density_len = ftell(density_file);
		fseek(density_file, 0L, SEEK_SET);	
		density_code = malloc((size_t)density_len * sizeof(char));
		if (fread(density_code, sizeof(char), (size_t)density_len, density_file) != (size_t)density_len) {
			fprintf(stderr, "Error: short read from SaveCode\n");
			return 1;
		}
		fclose(density_file);
	} else {
		fprintf(stderr, "Unable to open SaveCode file; density_function_code will not be written to grid.nc\n");
		density_len = 0L;
		density_code = NULL;
	}


	xCell = (double *)malloc(sizeof(double) * (size_t)nCells);
	yCell = (double *)malloc(sizeof(double) * (size_t)nCells);
	zCell = (double *)malloc(sizeof(double) * (size_t)nCells);

	fprintf(stderr, "Reading SaveVertices\n");
	fd = fopen("SaveVertices", "r");
	norm_min = DBL_MAX;
	norm_max = -DBL_MAX;
	for (i=0; i<nCells; i++) {
		if (fscanf(fd, "%lf %lf %lf\n", &xCell[i], &yCell[i], &zCell[i]) != 3) {
			fprintf(stderr, "Error: failed to read cell %zu from SaveVertices\n", i);
			return 1;
		}
		v = normalize_vect_r3(&xCell[i], &yCell[i], &zCell[i]);
		norm_min = min(v, norm_min);
		norm_max = max(v, norm_max);
	}
	fclose(fd);
	fprintf(stderr, "  found sphere radii between %f and %f\n", norm_min, norm_max);

	cellsOnVertex = (int *)malloc(sizeof(int) * (size_t)nVertices * (size_t)vertexDegree);

	xVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);
	yVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);
	zVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);

	fd = fopen("SaveTriangles", "r");
	fprintf(stderr, "Reading SaveTriangles\n");
	tri_min = INT_MAX;
	tri_max = INT_MIN;
	for (i=0; i<(size_t)nVertices; i++) {
		if (fscanf(fd, "%d %d %d\n", &cellsOnVertex[(size_t)vertexDegree*i], &cellsOnVertex[(size_t)vertexDegree*i+1], &cellsOnVertex[(size_t)vertexDegree*i+2]) != 3) {
			fprintf(stderr, "Error: failed to read triangle %zu from SaveTriangles\n", i);
			return 1;
		}
		for (j=0; j<(int)vertexDegree; j++) {
			tri_min = min(cellsOnVertex[(size_t)vertexDegree*i+(size_t)j], tri_min);
			tri_max = max(cellsOnVertex[(size_t)vertexDegree*i+(size_t)j], tri_max);
		}
	}
	fclose(fd);
	if (tri_min == 1 && tri_max == nCells) {
		fprintf(stderr, "  detected 1-based indices\n");
	} else if (tri_min == 0 && tri_max == (nCells - 1)) {
		fprintf(stderr, "  detected 0-based indices\n");
		for (i=0; i<(size_t)nVertices; i++) {
			for (j=0; j<(int)vertexDegree; j++) {
				cellsOnVertex[(size_t)vertexDegree*i+(size_t)j] += 1;
			}
		}
	} else {
		fprintf(stderr, "  invalid index range %i to %i for triangle vertices\n", tri_min, tri_max);
		return 1;
	}

	/*
	 * Determine the maximum number of neighbors that any cell has
	 */
	nEdgesOnCell = (int *)malloc(sizeof(int) * (size_t)nCells);
	bzero((void *)nEdgesOnCell, sizeof(int) * (size_t)nCells);

	for (i=0; i<(size_t)nVertices; i++) {
		for (j=0; j<(int)vertexDegree; j++) {
			nEdgesOnCell[cellsOnVertex[(size_t)vertexDegree*i+(size_t)j]-1]++;
		}
	}

	maxEdges = 0;
	for (i=0; i<(size_t)nCells; i++) {
		if (nEdgesOnCell[i] > (int)maxEdges) {
			maxEdges = nEdgesOnCell[i];
		}
	}
	fprintf(stderr, "maxEdges = %i\n", (int)maxEdges);

	free(nEdgesOnCell);
	nEdgesOnCell = NULL;


	fprintf(stderr, "Computing {x,y,z}Vertex\n");
	for (i=0; i<(size_t)nVertices; i++) {
		a.x = xCell[cellsOnVertex[(size_t)vertexDegree*i]-1];
		a.y = yCell[cellsOnVertex[(size_t)vertexDegree*i]-1];
		a.z = zCell[cellsOnVertex[(size_t)vertexDegree*i]-1];

		b.x = xCell[cellsOnVertex[(size_t)vertexDegree*i+1]-1];
		b.y = yCell[cellsOnVertex[(size_t)vertexDegree*i+1]-1];
		b.z = zCell[cellsOnVertex[(size_t)vertexDegree*i+1]-1];

		c.x = xCell[cellsOnVertex[(size_t)vertexDegree*i+2]-1];
		c.y = yCell[cellsOnVertex[(size_t)vertexDegree*i+2]-1];
		c.z = zCell[cellsOnVertex[(size_t)vertexDegree*i+2]-1];

		/*
		 * If points are not in CCW order, swap the second and third vertices
		 */
		if (ccw_triangle(&a, &b, &c)) {
			k = cellsOnVertex[(size_t)vertexDegree*i+1];
			cellsOnVertex[(size_t)vertexDegree*i+1] = cellsOnVertex[(size_t)vertexDegree*i+2];
			cellsOnVertex[(size_t)vertexDegree*i+2] = k;
		}

		vtx = circumcenter(a, b, c);

		xVertex[i] = vtx.x;
		yVertex[i] = vtx.y;
		zVertex[i] = vtx.z;
	}

	nEdgesOnCell = (int *)malloc(sizeof(int) * (size_t)nCells);
	bzero((void *)nEdgesOnCell, (size_t)nCells * sizeof(int));

	fprintf(stderr, "Computing cellsOnCell\n");
	cellsOnCell = (int *)malloc(sizeof(int) * (size_t)nCells * (size_t)maxEdges);
	for (i=0; i<(size_t)nVertices; i++) {
		j1 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i]-1);
		j2 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i+1]-1);
		j3 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i+2]-1);

		/*
		 *  Add j2 and j3 as neighbors of j1
		 */
		if (!found_cell(&cellsOnCell[j1*(size_t)maxEdges], nEdgesOnCell[j1], (j2+1))) {
			cellsOnCell[j1*(size_t)maxEdges + nEdgesOnCell[j1]] = (j2+1);
			nEdgesOnCell[j1]++;
		}
		if (!found_cell(&cellsOnCell[j1*(size_t)maxEdges], nEdgesOnCell[j1], (j3+1))) {
			cellsOnCell[j1*(size_t)maxEdges + nEdgesOnCell[j1]] = (j3+1);
			nEdgesOnCell[j1]++;
		}

		/*
		 *  Add j1 and j3 as neighbors of j2
		 */
		if (!found_cell(&cellsOnCell[j2*(size_t)maxEdges], nEdgesOnCell[j2], (j1+1))) {
			cellsOnCell[j2*(size_t)maxEdges + nEdgesOnCell[j2]] = (j1+1);
			nEdgesOnCell[j2]++;
		}
		if (!found_cell(&cellsOnCell[j2*(size_t)maxEdges], nEdgesOnCell[j2], (j3+1))) {
			cellsOnCell[j2*(size_t)maxEdges + nEdgesOnCell[j2]] = (j3+1);
			nEdgesOnCell[j2]++;
		}

		/*
		 *  Add j1 and j2 as neighbors of j3
		 */
		if (!found_cell(&cellsOnCell[j3*(size_t)maxEdges], nEdgesOnCell[j3], (j1+1))) {
			cellsOnCell[j3*(size_t)maxEdges + nEdgesOnCell[j3]] = (j1+1);
			nEdgesOnCell[j3]++;
		}
		if (!found_cell(&cellsOnCell[j3*(size_t)maxEdges], nEdgesOnCell[j3], (j2+1))) {
			cellsOnCell[j3*(size_t)maxEdges + nEdgesOnCell[j3]] = (j2+1);
			nEdgesOnCell[j3]++;
		}
	}


	/*
	 * Check that all nEdgesOnCell values are at most maxEdges
	 */
	for (i=0; i<(size_t)nCells; i++) {
		if (nEdgesOnCell[i] > (int)maxEdges) {
			fprintf(stderr, "Error: cell %i (1-based) has %i neighbors, which is more than maxEdges=%i\n",
				(int)(i+1), nEdgesOnCell[i], (int)maxEdges);
			return 1;
		}
	}


	/*
	 * Place cellsOnCell in CCW order. This is done by arbitrarily chosing
	 * the first cellOnCell to be at the beginning of the order, and computing
	 * the angle at the center of the cell (a) from the start of the order (b)
	 * to each of the other cells (c). The cellsOnCell are then rewritten in
	 * increasing order of angle.
	 */
	fprintf(stderr, "Ordering cellsOnCell\n");
	angles = (double *)malloc(sizeof(double) * (size_t)maxEdges);
	elements = (int *)malloc(sizeof(int) * (size_t)maxEdges);
	for (i=0; i<(size_t)nCells; i++) {
		a.x = xCell[i];
		a.y = yCell[i];
		a.z = zCell[i];

		b.x = xCell[cellsOnCell[(size_t)maxEdges*i]-1];
		b.y = yCell[cellsOnCell[(size_t)maxEdges*i]-1];
		b.z = zCell[cellsOnCell[(size_t)maxEdges*i]-1];
		angles[0] = (double)0.0;
		elements[0] = cellsOnCell[(size_t)maxEdges*i];

		/*
		 * Compute angles for all cellsOnCell
		 */
		for (k = 1; k < (size_t)nEdgesOnCell[i]; k++) {
			c.x = xCell[cellsOnCell[(size_t)maxEdges*i+k]-1];
			c.y = yCell[cellsOnCell[(size_t)maxEdges*i+k]-1];
			c.z = zCell[cellsOnCell[(size_t)maxEdges*i+k]-1];
			angles[k] = sphere_angle(a, b, c);
			if (angles[k] < (double)0.0) angles[k] += twopi;
			elements[k] = cellsOnCell[(size_t)maxEdges*i+k];
		}

		/*
		 * Sort angles by storing their sorted order in the elements array
		 * Since the list to be sorted will only have 6-7 elements, insertion sort is fine...
		 */
		for (k = 1; k < (size_t)nEdgesOnCell[i]; k++) {
			kk = k;
			while (kk > 0 && angles[kk-1] > angles[kk]) {
				minangle = angles[kk];
				angles[kk] = angles[kk-1];
				angles[kk-1] = minangle;

				minelement = elements[kk];
				elements[kk] = elements[kk-1];
				elements[kk-1] = minelement;

				kk--;
			}
		}

		/*
		 * Set cellsOnCell to sorted index order
		 */
		for (k = 0; k < (size_t)nEdgesOnCell[i]; k++) {
			cellsOnCell[(size_t)maxEdges*i+k] = elements[k];
		}

		/*
		 * Set cellsOnCell to zero for non-existent neighbors
		 */
		for (k = (size_t)nEdgesOnCell[i]; k < (size_t)maxEdges; k++) {
			cellsOnCell[(size_t)maxEdges*i+k] = 0;
		}
	}
	free(angles);
	free(elements);


	/*
	 * Write a graph.info file
	 */
	fprintf(stderr, "Writing graph.info\n");
	fd = fopen("graph.info", "w");
	fprintf(fd, "%i %i\n", (int)nCells, (int)nEdges);
	for (i=0; i<nCells; i++) {
		for (k=0; k<nEdgesOnCell[i]; k++) {
			fprintf(fd, "%i ", cellsOnCell[(size_t)maxEdges*i+k]);
		}
		fprintf(fd, "\n");
	}
	fclose(fd);


	fprintf(stderr, "Computing verticesOnCell\n");
	verticesOnCell = (int *)malloc(sizeof(int) * (size_t)nCells * (size_t)maxEdges);

	/* Zero-out the nEdgesOnCell array for use as an index of the last-used entry
	 * in the verticesOnCell arrays for each cell, below.
	 */
	bzero((void *)nEdgesOnCell, (size_t)nCells * sizeof(int));

	for (i=0; i<(size_t)nVertices; i++) {
		j1 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i]-1);
		j2 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i+1]-1);
		j3 = (MPI_Offset)(cellsOnVertex[(size_t)vertexDegree*i+2]-1);

		/*
		 *  Add i as neighbor of j1
		 */
		verticesOnCell[j1*(size_t)maxEdges + nEdgesOnCell[j1]] = (i+1);
		nEdgesOnCell[j1]++;

		/*
		 *  Add i as neighbor of j2
		 */
		verticesOnCell[j2*(size_t)maxEdges + nEdgesOnCell[j2]] = (i+1);
		nEdgesOnCell[j2]++;

		/*
		 *  Add i as neighbor of j3
		 */
		verticesOnCell[j3*(size_t)maxEdges + nEdgesOnCell[j3]] = (i+1);
		nEdgesOnCell[j3]++;
	}


	/*
	 * Place verticesOnCell in CCW order. This is done by arbitrarily chosing
	 * the last cellOnCell to be at the beginning of the order, and computing
	 * the angle at the center of the cell (a) from the start of the order (b)
	 * to each of the vertices (c). The verticesOnCell are then rewritten in
	 * increasing order of angle.
	 * The end result should be that both cellsOnCell and verticesOnCell are in
	 * CCW order, but verticesOnCell "leads" cellsOnCell.
	 */
	fprintf(stderr, "Ordering verticesOnCell\n");
	angles = (double *)malloc(sizeof(double) * (size_t)maxEdges);
	elements = (int *)malloc(sizeof(int) * (size_t)maxEdges);
	for (i=0; i<(size_t)nCells; i++) {
		a.x = xCell[i];
		a.y = yCell[i];
		a.z = zCell[i];

		/* Last cellOnCell */
		b.x = xCell[cellsOnCell[(size_t)maxEdges*i + (size_t)(nEdgesOnCell[i]-1)]-1];
		b.y = yCell[cellsOnCell[(size_t)maxEdges*i + (size_t)(nEdgesOnCell[i]-1)]-1];
		b.z = zCell[cellsOnCell[(size_t)maxEdges*i + (size_t)(nEdgesOnCell[i]-1)]-1];

		/*
		 * Compute angles for all verticesOnCell
		 */
		for (k = 0; k < (size_t)nEdgesOnCell[i]; k++) {
			c.x = xVertex[verticesOnCell[(size_t)maxEdges*i+k]-1];
			c.y = yVertex[verticesOnCell[(size_t)maxEdges*i+k]-1];
			c.z = zVertex[verticesOnCell[(size_t)maxEdges*i+k]-1];
			angles[k] = sphere_angle(a, b, c);
			if (angles[k] < (double)0.0) angles[k] += twopi;
			elements[k] = verticesOnCell[(size_t)maxEdges*i+k];
		}

		/*
		 * Sort angles by storing their sorted order in the elements array
		 * Since the list to be sorted will only have 6-7 elements, insertion sort is fine...
		 */
		for (k = 1; k < (size_t)nEdgesOnCell[i]; k++) {
			kk = k;
			while (kk > 0 && angles[kk-1] > angles[kk]) {
				minangle = angles[kk];
				angles[kk] = angles[kk-1];
				angles[kk-1] = minangle;

				minelement = elements[kk];
				elements[kk] = elements[kk-1];
				elements[kk-1] = minelement;

				kk--;
			}
		}

		/*
		 * Set verticesOnCell to sorted index order
		 */
		for (k = 0; k < (size_t)nEdgesOnCell[i]; k++) {
			verticesOnCell[(size_t)maxEdges*i+k] = elements[k];
		}

		/*
		 * Set verticesOnCell to zero for non-existent corners
		 */
		for (k = (size_t)nEdgesOnCell[i]; k < (size_t)maxEdges; k++) {
			verticesOnCell[(size_t)maxEdges*i+k] = 0;
		}
	}
	free(angles);
	free(elements);

	fprintf(stderr, "Writing grid.nc\n");
	write_netcdf(nCells, nEdges, nVertices, maxEdges, vertexDegree,
	             xCell, yCell, zCell,
	             xVertex, yVertex, zVertex,
	             nEdgesOnCell,
	             cellsOnCell, verticesOnCell,
	             cellsOnVertex,
	             density_len, density_code,
	             min_dc_m
	            );

	free(density_code);

	free(xCell);
	free(yCell);
	free(zCell);
	free(cellsOnVertex);
	free(xVertex);
	free(yVertex);
	free(zVertex);
	free(nEdgesOnCell);
	free(cellsOnCell);
	free(verticesOnCell);


/* Begin part 2 */

	ncerr = ncmpi_open(MPI_COMM_WORLD, "grid.nc", NC_WRITE, MPI_INFO_NULL, &ncid);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error opening grid.nc: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	ncerr = ncmpi_inq_dimid(ncid, "nCells", &ncells_id);
	ncerr = ncmpi_inq_dimlen(ncid, ncells_id, &nCells);

	ncerr = ncmpi_inq_dimid(ncid, "nEdges", &nedges_id);
	ncerr = ncmpi_inq_dimlen(ncid, nedges_id, &nEdges);

	ncerr = ncmpi_inq_dimid(ncid, "nVertices", &nvertices_id);
	ncerr = ncmpi_inq_dimlen(ncid, nvertices_id, &nVertices);

	ncerr = ncmpi_inq_dimid(ncid, "maxEdges", &maxedges_id);
	ncerr = ncmpi_inq_dimlen(ncid, maxedges_id, &maxEdges);

	ncerr = ncmpi_inq_dimid(ncid, "maxEdges2", &maxedges2_id);
	ncerr = ncmpi_inq_dimlen(ncid, maxedges2_id, &maxEdges2);

	ncerr = ncmpi_inq_dimid(ncid, "vertexDegree", &vertexdegree_id);
	ncerr = ncmpi_inq_dimlen(ncid, vertexdegree_id, &vertexDegree);

	ncerr = ncmpi_begin_indep_data(ncid);

	fprintf(stderr, "nCells = %ld\n", (long)nCells);
	fprintf(stderr, "nEdges = %ld\n", (long)nEdges);
	fprintf(stderr, "nVertices = %ld\n", (long)nVertices);
	fprintf(stderr, "maxEdges = %ld\n", (long)maxEdges);
	fprintf(stderr, "maxEdges2 = %ld\n", (long)maxEdges2);
	fprintf(stderr, "vertexDegree = %ld\n", (long)vertexDegree);

	ncerr = get_var(ncid, "nEdgesOnCell", (void **)&nEdgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading nEdgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "cellsOnCell", (void **)&cellsOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading cellsOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "edgesOnCell", (void **)&edgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading edgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "cellsOnEdge", (void **)&cellsOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading cellsOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	/*
	 * Compute edgesOnCell and cellsOnEdge
	 */
	fprintf(stderr, "Computing edgesOnCell and cellsOnEdge\n");
	edge_id = 1;
	for (i=0; i<nCells; i++) {
		for (j=0; j<nEdgesOnCell[i]; j++) {

			/*
			 * Only "create" a new edge (by using and incrementing edge_id)
			 * if the edge is "owned" by this cell (i.e., this cell has a lower
			 * index than the cellOnCell on the other side of the edge).
			 */
			ii = cellsOnCell[i*maxEdges + (size_t)j] - 1;    /* 0-based neighbor index */
			if (i < (size_t)ii) {

				/* Add edge for ourself */
				edgesOnCell[i*maxEdges + (size_t)j] = edge_id;

				/* Add edge for our neighbor */
				if ((k = find_neighbor(&cellsOnCell[(size_t)ii*maxEdges], maxEdges, i+1)) > 0) {
					edgesOnCell[(size_t)ii*maxEdges + (size_t)(k-1)] = edge_id;
				}
				else {
					fprintf(stderr, "Error: Did not find ourself in our neighbor\'s neighbor list\n");
					fprintf(stderr, "       Cell ID (1-based): %i\n", (int)(i+1));
					fprintf(stderr, "       Neighbor number (1-based): %i\n", (j+1));
					fprintf(stderr, "       Neighbor ID (1-based): %i\n", (ii+1));
					return 1;
				}

				/* Add i and ii to cellsOnEdge */
				cellsOnEdge[(size_t)(edge_id-1)*(size_t)2 + (size_t)0] = i+1;
				cellsOnEdge[(size_t)(edge_id-1)*(size_t)2 + (size_t)1] = ii+1;

				edge_id++;
			}
		}
	}

	ncerr = put_var(ncid, "edgesOnCell", (void *)edgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing edgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = put_var(ncid, "cellsOnEdge", (void *)cellsOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing cellsOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(edgesOnCell);
	free(cellsOnEdge);

	ncerr = get_var(ncid, "verticesOnCell", (void **)&verticesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "verticesOnEdge", (void **)&verticesOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	/*
	 * Compute verticesOnEdge
	 */
	fprintf(stderr, "Computing verticesOnEdge\n");
	edge_id = 1;
	for (i=0; i<nCells; i++) {
		for (j=0; j<nEdgesOnCell[i]; j++) {

			ii = cellsOnCell[i*maxEdges + (size_t)j] - 1;    /* 0-based neighbor index */
			if (i < (size_t)ii) {

				/* Since the verticesOnCell ordering "leads" the cellsOnCell ordering, the two vertices
				 * that form the endpoints for the jth edgeOnCell should be verticesOnCell[i][j] and
				 * verticesOnCell[i][j+1] (with some modular arithmetic).
				 */
				verticesOnEdge[(size_t)(edge_id-1)*(size_t)2 + (size_t)0] = verticesOnCell[i*maxEdges + (size_t)j];
				verticesOnEdge[(size_t)(edge_id-1)*(size_t)2 + (size_t)1] = verticesOnCell[i*maxEdges + (size_t)((j+1) % nEdgesOnCell[i])];

				edge_id++;
			}
		}
	}

	ncerr = put_var(ncid, "verticesOnEdge", (void *)verticesOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing verticesOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(verticesOnCell);
	free(nEdgesOnCell);
	free(cellsOnCell);
	free(verticesOnEdge);


	/*
	 * Compute xEdge, yEdge, zEdge
	 * In order to limit memory usage, this is done for chunks of the edges at a time
	 */
	fprintf(stderr, "Computing {x,y,z}Edge\n");
	ncerr = get_var(ncid, "xCell", (void **)&xCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yCell", (void **)&yCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zCell", (void **)&zCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "xVertex", (void **)&xVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yVertex", (void **)&yVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zVertex", (void **)&zVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	chunksize = (size_t)10000000;

	xEdge = (double *)malloc(sizeof(double) * min(chunksize,(size_t)nEdges));
	yEdge = (double *)malloc(sizeof(double) * min(chunksize,(size_t)nEdges));
	zEdge = (double *)malloc(sizeof(double) * min(chunksize,(size_t)nEdges));

	start[1] = (MPI_Offset)0;
	count[1] = (MPI_Offset)2;
	for (ichunk=0; ichunk<(size_t)nEdges; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nEdges - ichunk));

		ncerr = get_var(ncid, "verticesOnEdge", (void **)&verticesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading verticesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "cellsOnEdge", (void **)&cellsOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading cellsOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {
			c1.x = xCell[cellsOnEdge[(size_t)2 * i + (size_t)0]-1];
			c1.y = yCell[cellsOnEdge[(size_t)2 * i + (size_t)0]-1];
			c1.z = zCell[cellsOnEdge[(size_t)2 * i + (size_t)0]-1];

			c2.x = xCell[cellsOnEdge[(size_t)2 * i + (size_t)1]-1];
			c2.y = yCell[cellsOnEdge[(size_t)2 * i + (size_t)1]-1];
			c2.z = zCell[cellsOnEdge[(size_t)2 * i + (size_t)1]-1];

			v1.x = xVertex[verticesOnEdge[(size_t)2 * i + (size_t)0]-1];
			v1.y = yVertex[verticesOnEdge[(size_t)2 * i + (size_t)0]-1];
			v1.z = zVertex[verticesOnEdge[(size_t)2 * i + (size_t)0]-1];

			v2.x = xVertex[verticesOnEdge[(size_t)2 * i + (size_t)1]-1];
			v2.y = yVertex[verticesOnEdge[(size_t)2 * i + (size_t)1]-1];
			v2.z = zVertex[verticesOnEdge[(size_t)2 * i + (size_t)1]-1];

			e = arc_intersect(c1, c2, v1, v2);

			xEdge[i] = e.x;
			yEdge[i] = e.y;
			zEdge[i] = e.z;
		}

		free(verticesOnEdge);
		free(cellsOnEdge);

		ncerr = put_var(ncid, "xEdge", (void *)xEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing xEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = put_var(ncid, "yEdge", (void *)yEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing yEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = put_var(ncid, "zEdge", (void *)zEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing zEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
	}

	free(xCell);
	free(yCell);
	free(zCell);
	free(xVertex);
	free(yVertex);
	free(zVertex);
	free(xEdge);
	free(yEdge);
	free(zEdge);


	/*
	 * Compute edgesOnVertex
	 */
	fprintf(stderr, "Computing edgesOnVertex\n");
	ncerr = get_var(ncid, "verticesOnEdge", (void **)&verticesOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "edgesOnVertex", (void **)&edgesOnVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading edgesOnVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	memset((void *)edgesOnVertex, 0, sizeof(int) * (size_t)nVertices * (size_t)vertexDegree);
	for (i=0; i<(size_t)nEdges; i++) {

		/* Add this edge to edgesOnVertex for the first vertex on this edge */
		ii = verticesOnEdge[i * (size_t)2 + (size_t)0];
		for (j=0; j<vertexDegree; j++) {
			if (edgesOnVertex[(size_t)(ii-1) * (size_t)vertexDegree + (size_t)j] == 0) {
				edgesOnVertex[(size_t)(ii-1) * (size_t)vertexDegree + (size_t)j] = (i+1);
				break;
			}
		}
		if (j > vertexDegree) {
			fprintf(stderr, "We found too many edges incident with vertex %i (1-based)\n", ii);
			return 1;
		}

		/* Add this edge to edgesOnVertex for the second vertex on this edge */
		ii = verticesOnEdge[i * (size_t)2 + (size_t)1];
		for (j=0; j<vertexDegree; j++) {
			if (edgesOnVertex[(size_t)(ii-1) * (size_t)vertexDegree + (size_t)j] == 0) {
				edgesOnVertex[(size_t)(ii-1) * (size_t)vertexDegree + (size_t)j] = (i+1);
				break;
			}
		}
		if (j > vertexDegree) {
			fprintf(stderr, "We found too many edges incident with vertex %i (1-based)\n", ii);
			return 1;
		}
	}

	ncerr = put_var(ncid, "edgesOnVertex", (void *)edgesOnVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing edgesOnVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(verticesOnEdge);
	free(edgesOnVertex);


	/*
	 * Place edgesOnVertex in CCW order, and starting at the correct position w.r.t. cellsOnVertex
	 */
	fprintf(stderr, "Ordering edgesOnVertex\n");

	ncerr = get_var(ncid, "xCell", (void **)&xCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yCell", (void **)&yCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zCell", (void **)&zCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "xEdge", (void **)&xEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yEdge", (void **)&yEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zEdge", (void **)&zEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	angles = (double *)malloc(sizeof(double) * (size_t)vertexDegree);
	elements = (int *)malloc(sizeof(int) * (size_t)vertexDegree);
	start[1] = (MPI_Offset)0;
	count[1] = vertexDegree;
	for (ichunk=0; ichunk<(size_t)nVertices; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nVertices - ichunk));

		ncerr = get_var(ncid, "xVertex", (void **)&xVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading xVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "yVertex", (void **)&yVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading yVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "zVertex", (void **)&zVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading zVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "edgesOnVertex", (void **)&edgesOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading edgesOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "cellsOnVertex", (void **)&cellsOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading cellsOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {
			a.x = xVertex[i];
			a.y = yVertex[i];
			a.z = zVertex[i];

			/* All angles w.r.t. last cellOnVertex so that edgesOnVertex will "lead" cellsOnVertex */
			b.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];
			b.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];
			b.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];

			for (j=0; j<(int)vertexDegree; j++) {
if ((edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1) < 0 || (edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1) > nEdges) {
fprintf(stderr, "Bad access of edge %i\n", edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1);
fprintf(stderr, "Vertex = %li\n", i);
fprintf(stderr, "EdgeOnVertex = %i\n", j);
}
				c.x = xEdge[edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1];
				c.y = yEdge[edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1];
				c.z = zEdge[edgesOnVertex[i*(size_t)vertexDegree + (size_t)j]-1];

				angles[j] = sphere_angle(a, b, c);
				if (angles[j] < (double)0.0) angles[j] += twopi;
				elements[j] = edgesOnVertex[i*(size_t)vertexDegree + (size_t)j];
			}

			/*
			 * Sort angles by storing their sorted order in the elements array
			 * Since the list to be sorted will only have 3-4 elements, insertion sort is fine...
			 */
			for (k = 1; k < (size_t)vertexDegree; k++) {
				kk = k;
				while (kk > 0 && angles[kk-1] > angles[kk]) {
					minangle = angles[kk];
					angles[kk] = angles[kk-1];
					angles[kk-1] = minangle;

					minelement = elements[kk];
					elements[kk] = elements[kk-1];
					elements[kk-1] = minelement;

					kk--;
				}
			}

			/*
			 * Set edgesOnVertex to sorted index order
			 */
			for (k = 0; k < (size_t)vertexDegree; k++) {
				edgesOnVertex[i*(size_t)vertexDegree + k] = elements[k];
			}
		}

		free(xVertex);
		free(yVertex);
		free(zVertex);
		free(cellsOnVertex);

		ncerr = put_var(ncid, "edgesOnVertex", (void *)edgesOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing edgesOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		free(edgesOnVertex);
	}
	free(angles);
	free(elements);

	free(xCell);
	free(yCell);
	free(zCell);
	free(xEdge);
	free(yEdge);
	free(zEdge);


	/*
	 * Compute nEdgesOnEdge and edgesOnEdge
	 */
	fprintf(stderr, "Computing nEdgeOnEdge and edgesOnEdge\n");

	ncerr = get_var(ncid, "nEdgesOnCell", (void **)&nEdgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading nEdgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "edgesOnCell", (void **)&edgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading edgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	start[1] = (MPI_Offset)0;
	for (ichunk=0; ichunk<(size_t)nEdges; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nEdges - ichunk));

		count[1] = (MPI_Offset)2;
		ncerr = get_var(ncid, "cellsOnEdge", (void **)&cellsOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading cellsOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		count[1] = (MPI_Offset)1;  /* Does not really matter for nEdgesOnEdge... */
		ncerr = get_var(ncid, "nEdgesOnEdge", (void **)&nEdgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading nEdgesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		count[1] = maxEdges2;
		ncerr = get_var(ncid, "edgesOnEdge", (void **)&edgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading edgesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {
			cell1 = cellsOnEdge[i * (size_t)2 + (size_t)0];
			cell2 = cellsOnEdge[i * (size_t)2 + (size_t)1];

			nEdgesOnEdge[i] = nEdgesOnCell[cell1-1] + nEdgesOnCell[cell2-1] - 2;

			/* kk holds the next position in the edgesOnEdge array to fill */
			kk = 0;


			/* Find our index in the edgesOnCell list of cell1 */
			if ((k = find_neighbor(&edgesOnCell[(size_t)(cell1-1)*(size_t)maxEdges], maxEdges, (int)(i + start[0])+1)) > 0) {
				for (j=(int)(k-1)+1; j<(int)nEdgesOnCell[(size_t)(cell1-1)]; j++) {
					edgesOnEdge[i*(size_t)maxEdges2 + (kk++)] = edgesOnCell[(size_t)(cell1-1)*(size_t)maxEdges + (size_t)j];
				}
				for (j=0; j<(int)(k-1); j++) {
					edgesOnEdge[i*(size_t)maxEdges2 + (kk++)] = edgesOnCell[(size_t)(cell1-1)*(size_t)maxEdges + (size_t)j];
				}
			}
			else {
				fprintf(stderr, "Could not locate edge %i in edgesOnCell field for cell %i (1-based indices)\n", (int)(i+1), cell1);
				fprintf(stderr, "although cell %i is the first cellOnEdge for edge %i\n", cell1, (int)(i+1));
				return 1;
			}


			/* Find our index in the edgesOnCell list of cell2 */
			if ((k = find_neighbor(&edgesOnCell[(size_t)(cell2-1)*(size_t)maxEdges], maxEdges, (int)(i + start[0])+1)) > 0) {
				for (j=(int)(k-1)+1; j<(int)nEdgesOnCell[(size_t)(cell2-1)]; j++) {
					edgesOnEdge[i*(size_t)maxEdges2 + (kk++)] = edgesOnCell[(size_t)(cell2-1)*(size_t)maxEdges + (size_t)j];
				}
				for (j=0; j<(int)(k-1); j++) {
					edgesOnEdge[i*(size_t)maxEdges2 + (kk++)] = edgesOnCell[(size_t)(cell2-1)*(size_t)maxEdges + (size_t)j];
				}
			}
			else {
				fprintf(stderr, "Could not locate edge %i in edgesOnCell field for cell %i (1-based indices)\n", (int)(i+1), cell2);
				fprintf(stderr, "although cell %i is the second cellOnEdge for edge %i\n", cell2, (int)(i+1));
				return 1;
			}

			if ((int)kk != nEdgesOnEdge[i]) {
				fprintf(stderr, "The actual number of edgesOnEdge for edge %i (found to be %i) does not match the expected number (%i)\n", (int)(i+1), (int)kk, nEdgesOnEdge[i]);
				return 1;
			}

		}

		free(cellsOnEdge);

		count[1] = (MPI_Offset)1;  /* Does not really matter for nEdgesOnEdge... */
		ncerr = put_var(ncid, "nEdgesOnEdge", (void *)nEdgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing nEdgesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		count[1] = maxEdges2;
		ncerr = put_var(ncid, "edgesOnEdge", (void *)edgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing edgesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		free(nEdgesOnEdge);
		free(edgesOnEdge);
	}

	free(nEdgesOnCell);
	free(edgesOnCell);


	/*
	 * Compute areaTriangle
	 */
	fprintf(stderr, "Computing areaTriangle\n");

	ncerr = get_var(ncid, "xCell", (void **)&xCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yCell", (void **)&yCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zCell", (void **)&zCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	areaTriangle = (double *)malloc(sizeof(double) * min(chunksize,(size_t)nVertices));

	start[1] = (MPI_Offset)0;
	for (ichunk=0; ichunk<(size_t)nVertices; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nVertices - ichunk));

		count[1] = vertexDegree;
		ncerr = get_var(ncid, "cellsOnVertex", (void **)&cellsOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading cellsOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {
			a.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)0]-1];
			a.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)0]-1];
			a.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)0]-1];

			b.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)1]-1];
			b.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)1]-1];
			b.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)1]-1];

			c.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];
			c.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];
			c.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + (size_t)2]-1];

			/* TODO: here and elsewhere, we assume that vertexDegree == 3 */

			areaTriangle[i] = triangle_area(a, b, c);
		}

		free(cellsOnVertex);

		count[1] = (size_t)1; /* Does not really matter for areaTriangle */
		ncerr = put_var(ncid, "areaTriangle", (void *)areaTriangle, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing areaTriangle beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

	}

	free(areaTriangle);

	free(xCell);
	free(yCell);
	free(zCell);


	/*
	 * Compute areaCell
	 */
	fprintf(stderr, "Computing areaCell\n");

	ncerr = get_var(ncid, "xVertex", (void **)&xVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yVertex", (void **)&yVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zVertex", (void **)&zVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	areaCell = (double *)malloc(sizeof(double) * min(chunksize,(size_t)nCells));

	start[1] = (MPI_Offset)0;
	for (ichunk=0; ichunk<(size_t)nCells; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nCells - ichunk));

		count[1] = maxEdges;
		ncerr = get_var(ncid, "verticesOnCell", (void **)&verticesOnCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading verticesOnCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		count[1] = (size_t)1; /* Does not really matter for {x,y,z}Cell */
		ncerr = get_var(ncid, "xCell", (void **)&xCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading xCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "yCell", (void **)&yCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading yCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "zCell", (void **)&zCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading zCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		ncerr = get_var(ncid, "nEdgesOnCell", (void **)&nEdgesOnCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading nEdgesOnCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {

			areaCell[i] = (double)0.0;

			a.x = xCell[i];
			a.y = yCell[i];
			a.z = zCell[i];

			for (j=0; j<nEdgesOnCell[i]; j++) {
				b.x = xVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)j]-1];
				b.y = yVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)j]-1];
				b.z = zVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)j]-1];

				c.x = xVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)((j+1)%nEdgesOnCell[i])]-1];
				c.y = yVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)((j+1)%nEdgesOnCell[i])]-1];
				c.z = zVertex[verticesOnCell[i*(size_t)maxEdges + (size_t)((j+1)%nEdgesOnCell[i])]-1];

				areaCell[i] += triangle_area(a, b, c);
			}
		}

		free(verticesOnCell);
		free(nEdgesOnCell);
		free(xCell);
		free(yCell);
		free(zCell);

		count[1] = (size_t)1; /* Does not really matter for areaCell */
		ncerr = put_var(ncid, "areaCell", (void *)areaCell, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing areaCell beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

	}

	free(areaCell);

	free(xVertex);
	free(yVertex);
	free(zVertex);


	/*
	 * Compute kiteAreasOnVertex
	 */
	fprintf(stderr, "Computing kiteAreasOnVertex\n");

	ncerr = get_var(ncid, "xCell", (void **)&xCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yCell", (void **)&yCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zCell", (void **)&zCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	ncerr = get_var(ncid, "xEdge", (void **)&xEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yEdge", (void **)&yEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zEdge", (void **)&zEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	kiteAreasOnVertex = (double *)malloc(sizeof(double) * (size_t)vertexDegree * min(chunksize,(size_t)nVertices));

	start[1] = (MPI_Offset)0;
	for (ichunk=0; ichunk<(size_t)nVertices; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nVertices - ichunk));

		count[1] = vertexDegree;
		ncerr = get_var(ncid, "edgesOnVertex", (void **)&edgesOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading edgesOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "cellsOnVertex", (void **)&cellsOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading cellsOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		count[1] = (size_t)1; /* Does not matter for {x,y,z}Vertex */
		ncerr = get_var(ncid, "xVertex", (void **)&xVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading xVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "yVertex", (void **)&yVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading yVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
		ncerr = get_var(ncid, "zVertex", (void **)&zVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading zVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		for (i=0; i<(size_t)count[0]; i++) {

			a.x = xVertex[i];
			a.y = yVertex[i];
			a.z = zVertex[i];

			for (j=0; j<(int)vertexDegree; j++) {
				j1 = (size_t)j;
				j2 = (size_t)((j+1) % (int)vertexDegree);

				b.x = xEdge[edgesOnVertex[i*(size_t)vertexDegree + j1]-1];
				b.y = yEdge[edgesOnVertex[i*(size_t)vertexDegree + j1]-1];
				b.z = zEdge[edgesOnVertex[i*(size_t)vertexDegree + j1]-1];

				c.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];
				c.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];
				c.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];

				kiteAreasOnVertex[i * (size_t)vertexDegree + (size_t)j] = triangle_area(a, b, c);


				b.x = xCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];
				b.y = yCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];
				b.z = zCell[cellsOnVertex[i*(size_t)vertexDegree + j1]-1];

				c.x = xEdge[edgesOnVertex[i*(size_t)vertexDegree + j2]-1];
				c.y = yEdge[edgesOnVertex[i*(size_t)vertexDegree + j2]-1];
				c.z = zEdge[edgesOnVertex[i*(size_t)vertexDegree + j2]-1];

				kiteAreasOnVertex[i * (size_t)vertexDegree + (size_t)j] += triangle_area(a, b, c);
			}

		}

		free(edgesOnVertex);
		free(cellsOnVertex);

		free(xVertex);
		free(yVertex);
		free(zVertex);

		count[1] = (size_t)vertexDegree;
		ncerr = put_var(ncid, "kiteAreasOnVertex", (void *)kiteAreasOnVertex, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing kiteAreasOnVertex beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
	}

	free(kiteAreasOnVertex);

/* The following are re-used in the computation of dcEdge, below
	free(xCell);
	free(yCell);
	free(zCell);
*/

	free(xEdge);
	free(yEdge);
	free(zEdge);


	/*
	 * Compute dcEdge
	 */
	fprintf(stderr, "Computing dcEdge\n");

	ncerr = get_var(ncid, "cellsOnEdge", (void **)&cellsOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading cellsOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	dcEdge = (double *)malloc(sizeof(double) * (size_t)nEdges);

	for (i=0; i<(size_t)nEdges; i++) {
		cell1 = cellsOnEdge[i*(size_t)2 + (size_t)0] - 1;
		cell2 = cellsOnEdge[i*(size_t)2 + (size_t)1] - 1;

		a.x = xCell[cell1];
		a.y = yCell[cell1];
		a.z = zCell[cell1];

		b.x = xCell[cell2];
		b.y = yCell[cell2];
		b.z = zCell[cell2];

		dcEdge[i] = (double)sphere_distance(a, b);
	}

	free(xCell);
	free(yCell);
	free(zCell);

	free(cellsOnEdge);

	ncerr = put_var(ncid, "dcEdge", (void *)dcEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing dcEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(dcEdge);


	/*
	 * Compute dvEdge
	 */
	fprintf(stderr, "Computing dvEdge\n");

	ncerr = get_var(ncid, "verticesOnEdge", (void **)&verticesOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "xVertex", (void **)&xVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yVertex", (void **)&yVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zVertex", (void **)&zVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	dvEdge = (double *)malloc(sizeof(double) * (size_t)nEdges);

	for (i=0; i<(size_t)nEdges; i++) {
		vtx1 = verticesOnEdge[i*(size_t)2 + (size_t)0] - 1;
		vtx2 = verticesOnEdge[i*(size_t)2 + (size_t)1] - 1;

		a.x = xVertex[vtx1];
		a.y = yVertex[vtx1];
		a.z = zVertex[vtx1];

		b.x = xVertex[vtx2];
		b.y = yVertex[vtx2];
		b.z = zVertex[vtx2];

		dvEdge[i] = (double)sphere_distance(a, b);
	}

/* The following are re-used below in the computation of latVertex and lonVertex
	free(xVertex);
	free(yVertex);
	free(zVertex);
*/

	free(verticesOnEdge);

	ncerr = put_var(ncid, "dvEdge", (void *)dvEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing dvEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(dvEdge);


	/*
	 * Compute latVertex and lonVertex
	 */
	fprintf(stderr, "Computing latVertex and lonVertex\n");

/*
 * The xVertex, yVertex, and zVertex fields were already read in the code above
	ncerr = get_var(ncid, "xVertex", (void **)&xVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yVertex", (void **)&yVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zVertex", (void **)&zVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
*/

	latVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);
	lonVertex = (double *)malloc(sizeof(double) * (size_t)nVertices);

	for (i=0; i<(size_t)nVertices; i++) {
		xyz_to_latlon(xVertex[i], yVertex[i], zVertex[i], &latVertex[i], &lonVertex[i]);
	}

	free(xVertex);
	free(yVertex);
	free(zVertex);

	ncerr = put_var(ncid, "latVertex", (void *)latVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing latVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = put_var(ncid, "lonVertex", (void *)lonVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing lonVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(latVertex);
	free(lonVertex);


	/*
	 * Compute latCell and lonCell
	 */
	fprintf(stderr, "Computing latCell and lonCell\n");

	ncerr = get_var(ncid, "xCell", (void **)&xCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yCell", (void **)&yCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zCell", (void **)&zCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	latCell = (double *)malloc(sizeof(double) * (size_t)nCells);
	lonCell = (double *)malloc(sizeof(double) * (size_t)nCells);

	for (i=0; i<(size_t)nCells; i++) {
		xyz_to_latlon(xCell[i], yCell[i], zCell[i], &latCell[i], &lonCell[i]);
	}

	free(xCell);
	free(yCell);
	free(zCell);

	ncerr = put_var(ncid, "latCell", (void *)latCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing latCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = put_var(ncid, "lonCell", (void *)lonCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing lonCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(latCell);
	free(lonCell);


	/*
	 * Compute latEdge and lonEdge
	 */
	fprintf(stderr, "Computing latEdge and lonEdge\n");

	ncerr = get_var(ncid, "xEdge", (void **)&xEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yEdge", (void **)&yEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zEdge", (void **)&zEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	latEdge = (double *)malloc(sizeof(double) * (size_t)nEdges);
	lonEdge = (double *)malloc(sizeof(double) * (size_t)nEdges);

	for (i=0; i<(size_t)nEdges; i++) {
		xyz_to_latlon(xEdge[i], yEdge[i], zEdge[i], &latEdge[i], &lonEdge[i]);
	}

	free(xEdge);
	free(yEdge);
	free(zEdge);

	ncerr = put_var(ncid, "latEdge", (void *)latEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing latEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = put_var(ncid, "lonEdge", (void *)lonEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing lonEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(latEdge);
	free(lonEdge);


	/*
	 * Set meshDensity
	 */
	fprintf(stderr, "Adding meshDensity\n");

	meshDensity = (double *)malloc(sizeof(double) * (size_t)nCells);
        if ((fd = fopen("SaveDensity", "r")) != NULL) {
		fprintf(stderr, "Reading meshDensity from SaveDensity file\n");
		for (i=0; i<nCells; i++) {
			if (fscanf(fd, "%lf\n", &meshDensity[i]) != 1) {
				fprintf(stderr, "Error: failed to read meshDensity %zu from SaveDensity\n", i);
				return 1;
			}
		}
		fclose(fd);
	}
	else {
		fprintf(stderr, "Setting meshDensity to 1.0 everywhere\n");
		for (i=0; i<(size_t)nCells; i++) {
			meshDensity[i] = (double)1.0;
		}
	}

	ncerr = put_var(ncid, "meshDensity", (void *)meshDensity, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing meshDensity: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(meshDensity);


	/*
	 * Set indexToCellID
	 */
	fprintf(stderr, "Setting indexToCellID\n");

	indexToCellID = (int *)malloc(sizeof(int) * (size_t)nCells);

	for (i=0; i<(size_t)nCells; i++) {
		indexToCellID[i] = (i+1);
	}

	ncerr = put_var(ncid, "indexToCellID", (void *)indexToCellID, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing indexToCellID: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(indexToCellID);


	/*
	 * Set indexToVertexID
	 */
	fprintf(stderr, "Setting indexToVertexID\n");

	indexToVertexID = (int *)malloc(sizeof(int) * (size_t)nVertices);

	for (i=0; i<(size_t)nVertices; i++) {
		indexToVertexID[i] = (i+1);
	}

	ncerr = put_var(ncid, "indexToVertexID", (void *)indexToVertexID, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing indexToVertexID: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(indexToVertexID);


	/*
	 * Set indexToEdgeID
	 */
	fprintf(stderr, "Setting indexToEdgeID\n");

	indexToEdgeID = (int *)malloc(sizeof(int) * (size_t)nEdges);

	for (i=0; i<(size_t)nEdges; i++) {
		indexToEdgeID[i] = (i+1);
	}

	ncerr = put_var(ncid, "indexToEdgeID", (void *)indexToEdgeID, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing indexToEdgeID: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(indexToEdgeID);


	/*
	 * Compute angleEdge
	 */
	fprintf(stderr, "Computing angleEdge\n");

	ncerr = get_var(ncid, "verticesOnEdge", (void **)&verticesOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "xVertex", (void **)&xVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading xVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "yVertex", (void **)&yVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading yVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "zVertex", (void **)&zVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading zVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "latEdge", (void **)&latEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading latEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "lonEdge", (void **)&lonEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading lonEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	angleEdge = (double *)malloc(sizeof(double) * (size_t)nEdges);

	for (i=0; i<(size_t)nEdges; i++) {
		a.x = xVertex[verticesOnEdge[i * (size_t)2 + (size_t)0]-1];
		a.y = yVertex[verticesOnEdge[i * (size_t)2 + (size_t)0]-1];
		a.z = zVertex[verticesOnEdge[i * (size_t)2 + (size_t)0]-1];

		b.x = a.x - cos(lonEdge[i]) * sin(latEdge[i]);
		b.y = a.y - sin(lonEdge[i]) * sin(latEdge[i]);
		b.z = a.z + cos(latEdge[i]);

		v = sqrt(dot_product(b, b));
		b.x = b.x / v;
		b.y = b.y / v;
		b.z = b.z / v;

		c.x = xVertex[verticesOnEdge[i * (size_t)2 + (size_t)1]-1];
		c.y = yVertex[verticesOnEdge[i * (size_t)2 + (size_t)1]-1];
		c.z = zVertex[verticesOnEdge[i * (size_t)2 + (size_t)1]-1];

		angleEdge[i] = sphere_angle(a, b, c);
	}

	free(verticesOnEdge);
	free(xVertex);
	free(yVertex);
	free(zVertex);
	free(latEdge);
	free(lonEdge);

	ncerr = put_var(ncid, "angleEdge", (void *)angleEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error writing angleEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	free(angleEdge);


	/*
	 * Compute weightsOnEdge
	 */
	fprintf(stderr, "Computing weightsOnEdge\n");

	ncerr = get_var(ncid, "kiteAreasOnVertex", (void **)&kiteAreasOnVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading kiteAreasOnVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "cellsOnVertex", (void **)&cellsOnVertex, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading cellsOnVertex: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "verticesOnCell", (void **)&verticesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading verticesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "nEdgesOnCell", (void **)&nEdgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading nEdgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	kiteAreasOnCell = (double *)malloc(sizeof(double) * (size_t)nCells * (size_t)maxEdges);

	for (i=0; i<(size_t)nCells; i++) {
		for (j=0; j<nEdgesOnCell[i]; j++) {
			ii = verticesOnCell[i*(size_t)maxEdges + (size_t)j]-1;
			k = find_neighbor(&cellsOnVertex[(size_t)ii * (size_t)vertexDegree], vertexDegree, (int)(i+1));
			if (k > 0) {
				kiteAreasOnCell[i * (size_t)maxEdges + (size_t)j] = kiteAreasOnVertex[(size_t)ii * (size_t)vertexDegree + (k-1)];
			}
			else {
				fprintf(stderr, "Error finding cell %i in cellsOnVertex list for vertex %i (1-based)\n", (int)(i+1), (ii+1));
				return 1;
			}
		}
		for (j=nEdgesOnCell[i]; j<(int)maxEdges; j++) {
			kiteAreasOnCell[i * (size_t)maxEdges + (size_t)j] = (double)0.0;
		}
	}

	free(kiteAreasOnVertex);
	free(cellsOnVertex);
	free(verticesOnCell);

	ncerr = get_var(ncid, "areaCell", (void **)&areaCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading areaCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	for (i=0; i<(size_t)nCells; i++) {
		v = (double)1.0 / areaCell[i];
		for (j=0; j<nEdgesOnCell[i]; j++) {
			kiteAreasOnCell[i * (size_t)maxEdges + (size_t)j] *= v;
		}
	}

	free(areaCell);

	ncerr = get_var(ncid, "edgesOnCell", (void **)&edgesOnCell, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading edgesOnCell: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "cellsOnEdge", (void **)&cellsOnEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading cellsOnEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}
	ncerr = get_var(ncid, "dvEdge", (void **)&dvEdge, NULL, NULL);
	if (ncerr != NC_NOERR) {
		fprintf(stderr, "Error reading dvEdge: %s\n", ncmpi_strerror(ncerr));
		return 1;
	}

	weightsOnEdge = (double *)malloc(sizeof(double) * (size_t)maxEdges2 * min(chunksize,(size_t)nEdges));
	rotated_kiteAreasOnCell = (double *)malloc(sizeof(double) * (size_t)maxEdges);


	start[1] = (MPI_Offset)0;
	for (ichunk=0; ichunk<(size_t)nEdges; ichunk+=chunksize) {
		start[0] = (MPI_Offset)ichunk;
		count[0] = min(chunksize,((size_t)nEdges - ichunk));

		count[1] = (size_t)1; /* Does not really matter for nEdgesOnEdge */
		ncerr = get_var(ncid, "nEdgesOnEdge", (void **)&nEdgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading nEdgesOnEdge beginning at %lu (count %lu): %s\n", (unsigned long)start[0], (unsigned long)count[0], ncmpi_strerror(ncerr));
			fprintf(stderr, "ncerr = %i\n", ncerr);
			return 1;
		}

		count[1] = (size_t)1; /* Does not really matter for dcEdge */
		ncerr = get_var(ncid, "dcEdge", (void **)&dcEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading dcEdge beginning at %lu (count %lu): %s\n", (unsigned long)start[0], (unsigned long)count[0], ncmpi_strerror(ncerr));
			fprintf(stderr, "ncerr = %i\n", ncerr);
			return 1;
		}

		count[1] = (size_t)maxEdges2;
		ncerr = get_var(ncid, "edgesOnEdge", (void **)&edgesOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error reading edgesOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}

		/* Initialize weightsOnEdge to zero */
		memset((void *)weightsOnEdge, 0, sizeof(double) * (size_t)count[0] * (size_t)maxEdges2);

		for (i=0; i<(size_t)count[0]; i++) {
			j = 0;

			/* Handle weights from first cell on this edge */
			edge_id = (int)(start[0] + i);
			cell1 = cellsOnEdge[edge_id * (size_t)2 + (size_t)0] - 1;
			sum_r = (double)0.0;
			rotate_kiteAreasOnCell(&kiteAreasOnCell[(size_t)cell1 * (size_t)maxEdges],
						(edge_id+1),
						&edgesOnCell[(size_t)cell1 * (size_t)maxEdges],
						nEdgesOnCell[cell1],
						rotated_kiteAreasOnCell);

			for (ii=0; ii<nEdgesOnCell[cell1]-1; ii++) {
				sum_r += rotated_kiteAreasOnCell[ii];
				edge_id = edgesOnEdge[i*(size_t)maxEdges2 + (size_t)j] - 1;
				weightsOnEdge[i*(size_t)maxEdges2 + (size_t)j] = sign_for_edge((cell1+1), (edge_id+1), cellsOnEdge) * ((double)0.5 - sum_r) * dvEdge[edge_id] / dcEdge[i];
				j++;
			}

			/* Handle weights from second cell on this edge */
			edge_id = (int)(start[0] + i);
			cell2 = cellsOnEdge[edge_id * (size_t)2 + (size_t)1] - 1;
			sum_r = (double)0.0;
			rotate_kiteAreasOnCell(&kiteAreasOnCell[(size_t)cell2 * (size_t)maxEdges],
						(edge_id+1),
						&edgesOnCell[(size_t)cell2 * (size_t)maxEdges],
						nEdgesOnCell[cell2],
						rotated_kiteAreasOnCell);

			for (ii=0; ii<nEdgesOnCell[cell2]-1; ii++) {
				sum_r += rotated_kiteAreasOnCell[ii];
				edge_id = edgesOnEdge[i*(size_t)maxEdges2 + (size_t)j] - 1;
				weightsOnEdge[i*(size_t)maxEdges2 + (size_t)j] = (double)(-1.0) * sign_for_edge((cell2+1), (edge_id+1), cellsOnEdge) * ((double)0.5 - sum_r) * dvEdge[edge_id] / dcEdge[i];
				j++;
			}

			if (j != nEdgesOnEdge[i]) {
				fprintf(stderr, "nEdgesOnEdge[%i] does not agree with nEdgesOnCell values for cells %i and %i\n", (int)(start[0]+i+1), (cell1+1), (cell2+1));
				return 1;
			}
		}

		free(nEdgesOnEdge);
		free(dcEdge);
		free(edgesOnEdge);

		count[1] = (size_t)maxEdges2;
		ncerr = put_var(ncid, "weightsOnEdge", (void *)weightsOnEdge, start, count);
		if (ncerr != NC_NOERR) {
			fprintf(stderr, "Error writing weightsOnEdge beginning at %lu: %s\n", (unsigned long)start[0], ncmpi_strerror(ncerr));
			return 1;
		}
	}

	free(kiteAreasOnCell);
	free(nEdgesOnCell);
	free(edgesOnCell);
	free(cellsOnEdge);
	free(dvEdge);
	free(weightsOnEdge);
	free(rotated_kiteAreasOnCell);

	stat = MPI_Finalize();
	if (stat != MPI_SUCCESS) {
		fprintf(stderr, "Error: failed call to MPI_Finalize\n");
		return 1;
	}

	fprintf(stderr, "All done.\n");

	return 0;
}
