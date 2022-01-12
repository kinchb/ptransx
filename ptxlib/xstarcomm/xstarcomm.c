#include "xstarcomm.h"

double interpolate(double *x_vals, double *y_vals, int num, double x) {
    int i;

    // find i such that x_vals[i] <= x < x_vals[i+1]
    i = 0;
    while (i < num-1) {
        if ((x_vals[i] <= x) && (x < x_vals[i+1]))
            break;
        else
            i++;
    }

    // do some bounds checking
    if (i == num-1) {
        if (x < x_vals[0])
            return y_vals[0];
        if (x >= x_vals[num-1])
            return 0.0;
    }

    // interpolate
    return y_vals[i] + ((y_vals[i+1] - y_vals[i])/(x_vals[i+1] - x_vals[i]))*(x - x_vals[i]);
}

int xstarcomm(double *energies, double *mi, double density, double Fe_abund, double heat, double temp, int niter, int N, double *emis, double *c_emis, double *absorp, double *heat_vals, double *xee, double *new_temp) {
    // more variables we need for sending to XSTAR
    double *energies_c = NULL;

    double *x_energies = NULL;
    double *x_mi       = NULL;
    double *x_intmi    = NULL;

    // and for receiving from XSTAR
    double *x_c_emis = NULL;
    double *x_absorp = NULL;

    double *x_elines  = NULL;
    double *x_linemis = NULL;

    double mult = pow(10.0, log10(EUX/ELX)/(NX-1));

    double norm1, norm2, tmpsum;

    int i, j;

    setbuffer(stdout, NULL, 0);
    setbuffer(stderr, NULL, 0);

    while (energies_c == NULL)
        energies_c = malloc((N-1) * sizeof(*energies_c));

    while (x_energies == NULL)
        x_energies = malloc(NX * sizeof(*x_energies));
    while (x_mi == NULL)
        x_mi = malloc(NX * sizeof(*x_mi));
    while (x_intmi == NULL)
        x_intmi = malloc(NX * sizeof(*x_intmi));

    while (x_c_emis == NULL)
        x_c_emis = malloc(NX * sizeof(*x_c_emis));
    while (x_absorp == NULL)
        x_absorp = malloc(NX * sizeof(*x_absorp));

    while (x_elines == NULL)
        x_elines = malloc(NL * sizeof(*x_elines));
    while (x_linemis == NULL)
        x_linemis = malloc(NL * sizeof(*x_linemis));

    // make energies_c
    for (i = 0; i < N-1; i++)
        energies_c[i] = 0.5 * (energies[i+1] + energies[i]); 

    // populate the XSTAR energy grid
    x_energies[0] = ELX;
    for (i = 1; i < NX-1; i++)
        x_energies[i] = mult * x_energies[i-1];
    x_energies[NX-1] = EUX;

    // map mi onto x_mi
    for (i = 0; i < NX; i++)
        x_mi[i] = interpolate(energies_c, mi, N-1, x_energies[i]);

    // make x_intmi;
    x_intmi[NX-1] = 0.0;
    for (i = NX-2; i >= 0; i--)
        x_intmi[i] = x_intmi[i+1] + 0.5 * 1.602177e-12 * (x_mi[i] + x_mi[i+1]) * (x_energies[i+1] - x_energies[i]);

    xstarsub_(x_energies, x_mi, x_intmi, &density, &Fe_abund, &heat, &temp, &niter, xee, new_temp, x_c_emis, x_absorp, x_elines, x_linemis, heat_vals);

    if (*xee < 1.0)
        printf("xee is small: %e\n", *xee);

    /*
    for (i = 0; i < NX; i++) {
        if (x_c_emis[i] != x_c_emis[i]) {
            fprintf(stdout, "x_c_emis nan %d\n", i); fflush(stdout);
        }
        if (x_absorp[i] != x_absorp[i]) {
            fprintf(stdout, "x_absorp nan %d\n", i); fflush(stdout);
        }
    }
    */

    // map x_c_emis onto c_emis
    for (i = 0; i < N-1; i++)
        c_emis[i] = interpolate(x_energies, x_c_emis, NX, energies_c[i]);

    // map x_absorp onto absorp
    for (i = 0; i < N-1; i++) {
        absorp[i] = interpolate(x_energies, x_absorp, NX, energies_c[i]);
//      if (absorp[i] != absorp[i]) {
//          fprintf(stdout, "absorp nan %d\n", i); fflush(stdout);
//      }
    }

    // bin lines into continuum
    for (i = 0; i < N-1; i++) {
        tmpsum = 0.0;
        for (j = 0; j < NL; j++) {
            if ((x_elines[j] >= energies[i]) && (x_elines[j] < energies[i+1]))
                tmpsum += x_linemis[j];
        }
        tmpsum /= 1.602177e-12 * (energies[i+1] - energies[i]) * 2.0 * 3.14159;
        emis[i] = c_emis[i] + tmpsum;
    }

    // turn into photon- (from energy-) intensities
    for (i = 0; i < N-1; i++) {
        c_emis[i] /= energies_c[i];
        emis[i]   /= energies_c[i];
    }

    // free what isn't passed up
    free(energies_c);
    free(x_energies);
    free(x_mi);
    free(x_intmi);
    free(x_c_emis);
    free(x_absorp);
    free(x_elines);
    free(x_linemis);

    return 0;
}

