#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <sys/time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define K_TO_THETA (1.68637005e-10)

#define ETA_BND_L 1.0e-4
#define ETA_BND_U 1.0e4

#define MC_N (long)1e6

#define GAM_M_N  (long)1e6
#define GAM_N    (long)1e5
#define GAM_N_SA (long)1e5
#define MU_I_N   128

#define TOTAL_PHOTONS   (long)1e5
#define TOTAL_THICKNESS 100.0

double sig_adjust(double e);
double kn_cdf(double M, double e, double u);
double brent(double a, double b, double enrg, double u);
double draw_kn(double e);
double G_cdf(double gam, double Theta, double u);
double brent_mj(double a, double b, double Theta, double u);
double draw_gam_rej(double Theta);
double draw_mu(double beta);
double find_gam_max(double T);
void make_gam_cdf(double T, double *gam_grid, double *gam_cdf);
double draw_gam(double *gam_grid, double *gam_cdf);
double draw_B(double T);
void boost(double gam, double B, double A[4], double A_p[4]);
void rotate(double th, double A[4], double A_p[4]);
void sample_electron(double T, double E, double *beta, double *gam, double *mu);
double mc_scatters(double T, double E_i, double *E_record, double *mu_record);
void make_gam_pdf(double T, double *gam, double *dgam, double *gam_pdf, double *B);
double find_B_max(double T);
void make_B_pdf(double T, double *B, double *dB, double *B_pdf, double *gam);
double mf1(double e);
double mf2(double e);
double sa_calc_ratio(double T, double E_i, double *mu_i, double *w);
double sigma(double e);
double sa_calc_sigma(double T, double E_i, double *mu_i, double *w);
void single_scatter(double T, double E_i, double *E_f, double *mu);
// void simulate_scatter(double E_0, double th_0, double phi_0, double T, double *gam_grid, double *gam_cdf, double *E_f, double *th_f, double *phi_f);
void simulate_scatter(double E_0, double th_0, double phi_0, double T, double *E_f, double *th_f, double *phi_f);
void simulate_trajectory(double E_0, double T, double *E_record, double *mu_record);
double draw_E_pl(double Gamma, double E_min, double E_max);
void simulate_trajectory_pl(double Gamma, double E_min, double E_max, double T, double *E_record, double *mu_record);
void simulate_reflection(double E_0, double T, double *E_record, double *mu_record, double *scatter_record);
int new_z(double E, double *z_pos, double dtau, double *z, double *e, double *tau_grids, int D, int N, int *ndx);
int response(int top_or_bot, double E_0, double *z, double *e, double *T, double *tau_grids, int D, int N, double *cdf);
