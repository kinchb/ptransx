#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NX 10000
#define ELX 1.0
#define EUX 1.0e6
#define NL 500000

typedef int bool;
enum {false, true};

extern void xstarsub_(double *x_energies, double *x_mi, double *x_intmi, double *density, double *Fe_abund, double *heat, double *temp, int *niter, double *xee, double *new_temp, double *x_c_emis, double *x_absorp, double *x_elines, double *x_linemis, double *heat_vals);
extern void fbg2_(double *u, double *gam, double *fbg);

double interpolate(double *x_vals, double *y_vals, int num, double x);
int xstarcomm(double *energies, double *mi, double density, double Fe_abund, double heat, double temp, int niter, int N, double *emis, double *c_emis, double *absorp, double *heat_vals, double *xee, double *new_temp);

