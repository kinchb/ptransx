#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define K_TO_THETA (1.68637005e-10) // T_e (in Kelvin) * K_TO_THETA = k_B T_e / m_e c^2 

#define MC_N (long)1e7 // the number of photons to simulate n-scattering when building up redistribution function

#define N         8000   // the number of grid points in eta to use when calculating E_0 -> eta * E_0 probability
#define ETA_BND_L 1.0e-4 // the smallest resolved amplification factor eta
#define ETA_BND_U 1.0e4  // the largest resolved amplification

#define GAM_M_N  (long)1e6 // number of gamma grid-points used to determine gam_max 
#define GAM_N    (long)1e5 // number of gamma grid-points used to construct MJ CDF for direct inverse CDF method
#define GAM_N_SA (long)1e5 // number of gamma grid-points used for integration over MJ for semi-analytic mean amplification calculation
#define MU_I_N   128       // number of supplied Gauss-Legendre quadrature points used to perform the fluid frame integral over angle for SA mean amplification calc

// below is the collection of functions used to properly sample an electron pre-scatter velocity accounting for
// all relativistic effects and the Klein-Nishina cross section;
// some are also used for sampling the electron rest frame photon scattering angle
double sig_adjust(double e);
double kn_cdf(double M, double e, double u);
double brent_kn(double a, double b, double enrg, double u);
double draw_kn(double e);
double G_cdf(double gam, double Theta, double u);
double brent_mj(double a, double b, double Theta, double u);
double draw_gam_rej(double Theta);
double find_gam_max(double T);
void make_gam_cdf(double T, double *gam_grid, double *gam_cdf);
double draw_gam_cdf(double *gam_grid, double *gam_cdf);
double draw_B(double T);
double draw_mu(double beta);
void sample_electron(double T, double E, double *beta, double *gam, double *mu);

// functions for directly simulating series of Compton scattering events off thermal electrons populations
void boost(double gam, double B, double A[4], double A_p[4]);
void rotate(double th, double A[4], double A_p[4]);
double mc_scatters(int n, double T, double E_0, double *E_record, double *mu_record);

// functions for calculating the fluid frame mean amplification ratio
void make_gam_pdf(double T, double *gam, double *dgam, double *gam_pdf, double *B);
double find_B_max(double T);
void make_B_pdf(double T, double *B, double *dB, double *B_pdf, double *gam);
double mf1(double e);
double mf2(double e);
double sa_calc_ratio(double T, double E_i, double *mu_i, double *w);

// functions for calculating the fluid frame total (outscattering) opacity
double sigma(double e);
double sa_calc_sigma(double T, double E_i, double *mu_i, double *w);

// simulate a single scattering event
void single_scatter(double T, double E_i, double *E_f, double *mu);
