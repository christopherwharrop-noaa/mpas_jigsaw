/* Reads binary doubles from stdin and writes HFUN.msh as text.
   Usage: write_hfun <nlons> <nlats>
   Stdin: nlons doubles (lons), nlats doubles (lats), nlons*nlats doubles (values) */

#include <stdio.h>
#include <stdlib.h>

static double *read_doubles(size_t n)
{
    double *buf = malloc(n * sizeof(double));
    if (!buf) {
        fprintf(stderr, "write_hfun: malloc failed for %zu doubles\n", n);
        exit(1);
    }
    if (fread(buf, sizeof(double), n, stdin) != n) {
        fprintf(stderr, "write_hfun: short read, expected %zu doubles\n", n);
        exit(1);
    }
    return buf;
}

int main(int argc, char *argv[])
{
    if (argc != 3) {
        fprintf(stderr, "Usage: write_hfun <nlons> <nlats>\n");
        return 1;
    }

    size_t nlons = (size_t)atol(argv[1]);
    size_t nlats = (size_t)atol(argv[2]);
    size_t npts  = nlons * nlats;

    double *lons = read_doubles(nlons);
    double *lats = read_doubles(nlats);
    double *vals = read_doubles(npts);

    FILE *fp = fopen("HFUN.msh", "w");
    if (!fp) {
        perror("write_hfun: fopen");
        return 1;
    }

    /* Use a large buffer for faster writes */
    setvbuf(fp, NULL, _IOFBF, 4 * 1024 * 1024);

    fprintf(fp, "MSHID=3;ellipsoid-grid\n");
    fprintf(fp, "NDIMS=2\n");
    fprintf(fp, "COORD=1;%zu\n", nlons);
    for (size_t i = 0; i < nlons; i++)
        fprintf(fp, "%.17g\n", lons[i]);

    fprintf(fp, "COORD=2;%zu\n", nlats);
    for (size_t i = 0; i < nlats; i++)
        fprintf(fp, "%.17g\n", lats[i]);

    fprintf(fp, "VALUE=%zu; 1\n", npts);
    for (size_t i = 0; i < npts; i++)
        fprintf(fp, "%.17g\n", vals[i]);

    fclose(fp);
    free(lons);
    free(lats);
    free(vals);
    return 0;
}
