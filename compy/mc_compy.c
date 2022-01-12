#include "mc_compy.h"

// compute the ratio of Thomson cross section to Klein-Nishina cross section, in the electron rest frame (erf),
// for a photon of energy e [ev]
double sig_adjust(double e) {
    double eps = e/511.0e3;

    if (e < 100.0) {
        return (8.0/3.0)/(8.0/3.0 - (16.0/3.0)*eps + (208.0/15.0)*eps*eps - (522.0/15.0)*eps*eps*eps);
    } else {
        return (8.0/3.0)/(((2.0*eps*(2.0 + eps*(1.0 + eps)*(8.0 + eps)))/((1.0 + 2.0*eps)*(1.0 + 2.0*eps)) + (-2.0 + (-2.0 + eps)*eps)*log(1.0 + 2.0*eps))/(eps*eps*eps));
    }
}

// compute the probability P for a photon of dimensionless pre-scatter energy (i.e., energy / m_e c^2) to scatter,
// in the electron rest frame (erf), through an angle M (in cosine); return P - u
double kn_cdf(double M, double e, double u) {
    double tmp;

    if (e < 0.001)
        return 0.125*(4 + 3*M + M*M*M) - u;
    tmp = (-((2 + e*(6 + e - 8*e*e))/pow(1 + 2*e, 2)) + 2*e*M + e*e/pow(1 + e - e*M, 2) + (2 + 4*e)/(1 + e - e*M) + 2*e*e*log(1 + 2*e) - 2*e*e*log(1 + e - e*M) + 4*(1 + e)*log((1 + e - e*M)/(1 + 2*e)))/(2.*((2*e*(2 + e*(1 + e)*(8 + e)))/pow(1 + 2*e, 2) + (-2 + (-2 + e)*e)*log(1 + 2*e)));
    if (tmp < 0.0) {
        return -u;
    } else if (tmp > 1.0) {
        return 1.0 - u;
    } else {
        return tmp - u;
    }
}

// use Brent's method to solve the inverse CDF problem: KN_CDF(mu, enrg) = u;
// where mu is the scattering angle to be solved for given a supplied u on [0, 1),
// enrg is the dimensionless pre-scatter photon energy (i.e., energy / m_e c^2),
// and KN_CDF is the electron rest frame (erf) integral of the
// Klein-Nishina partial differential cross section (i.e., kn_cdf);
// solution must be bound by a, b (typically -1 and 1)
double brent_kn(double a, double b, double enrg, double u) {
    double machep, t, c, d, e, fa, fb, fc, m, p, q, r, s, sa, sb, tol;

    machep = DBL_EPSILON;
    t      = DBL_EPSILON;

    sa = a;
    sb = b;
    fa = kn_cdf(sa, enrg, u);
    fb = kn_cdf(sb, enrg, u);

    c  = sa;
    fc = fa;
    e  = sb - sa;
    d  = e;

    for (;;) {
        if (fabs(fc) < fabs (fb)) {
            sa = sb;
            sb = c;
            c  = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 2.0 * machep * fabs (sb) + t;
        m   = 0.5 * (c - sb);

        if (fabs(m) <= tol || fb == 0.0)
            break;

        if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
            e = m;
            d = e;
        } else {
            s = fb / fa;

            if (sa == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (0.0 < p) {
                q = -q;
            } else {
                p = -p;
            }

            s = e;
            e = d;

            if (2.0 * p < 3.0 * m * q - fabs(tol * q) && p < fabs(0.5 * s * q)) {
                d = p / q;
            } else {
                e = m;
                d = e;
            }
        }

        sa = sb;
        fa = fb;

        if (tol < fabs(d)) {
            sb = sb + d;
        } else if (0.0 < m) {
            sb = sb + tol;
        } else {
            sb = sb - tol;
        }

        fb = kn_cdf(sb, enrg, u);

        if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
            c  = sa;
            fc = fa;
            e  = sb - sa;
            d  = e;
        }
    }

    return sb;
}

// for a pre-scatter photon of energy e [eV], properly samples (via inverse CDF) the scattering angle (in cosine)
// in the electron rest frame (erf)
double draw_kn(double e) {
    return brent_kn(-1.0, 1.0, e/511.0e3, (double)rand()/(double)RAND_MAX);
}

// compute the probability P for an electron sampled from the *auxiliary* pdf to Maxwell-Juttner to have a gamma factor
// less than the supplied value (gam), given a dimensionless electron temperature Theta;
// return P - u
double G_cdf(double gam, double Theta, double u) {
    double tmp;

    tmp = 1.0 - exp((1.0-gam)/Theta)*(gam*gam + 2.0*gam*Theta + 2.0*Theta*Theta)/(2.0*Theta*Theta + 2.0*Theta + 1.0);
    if (tmp < 0.0) {
        return -u;
    } else if (tmp > 1.0) {
        return 1.0 - tmp;
    } else {
        return tmp - u;
    }
}

// use Brent's method to solve the inverse CDF problem: G_cdf(gam, Theta) = u;
// where gam is the electron gamma factor to be solved for given a supplied u on [0, 1),
// Theta is the dimensionless electron temperature,
// and G_cdf is integral of the *auxiliary* pdf to Maxwell-Juttner;
// solution must be bound by a, b (typically 1 and some very large number proportional to Theta)
double brent_mj(double a, double b, double Theta, double u) {
    double machep, t, c, d, e, fa, fb, fc, m, p, q, r, s, sa, sb, tol;

    machep = DBL_EPSILON;
    t      = DBL_EPSILON;

    sa = a;
    sb = b;
    fa = G_cdf(sa, Theta, u);
    fb = G_cdf(sb, Theta, u);

    c  = sa;
    fc = fa;
    e  = sb - sa;
    d  = e;

    for (;;) {
        if (fabs(fc) < fabs (fb)) {
            sa = sb;
            sb = c;
            c  = sa;
            fa = fb;
            fb = fc;
            fc = fa;
        }

        tol = 2.0 * machep * fabs (sb) + t;
        m   = 0.5 * (c - sb);

        if (fabs(m) <= tol || fb == 0.0)
            break;

        if (fabs(e) < tol || fabs(fa) <= fabs(fb)) {
            e = m;
            d = e;
        } else {
            s = fb / fa;

            if (sa == c) {
                p = 2.0 * m * s;
                q = 1.0 - s;
            } else {
                q = fa / fc;
                r = fb / fc;
                p = s * (2.0 * m * q * (q - r) - (sb - sa) * (r - 1.0));
                q = (q - 1.0) * (r - 1.0) * (s - 1.0);
            }

            if (0.0 < p) {
                q = -q;
            } else {
                p = -p;
            }

            s = e;
            e = d;

            if (2.0 * p < 3.0 * m * q - fabs(tol * q) && p < fabs(0.5 * s * q)) {
                d = p / q;
            } else {
                e = m;
                d = e;
            }
        }

        sa = sb;
        fa = fb;

        if (tol < fabs(d)) {
            sb = sb + d;
        } else if (0.0 < m) {
            sb = sb + tol;
        } else {
            sb = sb - tol;
        }

        fb = G_cdf(sb, Theta, u);

        if ((0.0 < fb && 0.0 < fc) || (fb <= 0.0 && fc <= 0.0)) {
            c  = sa;
            fc = fa;
            e  = sb - sa;
            d  = e;
        }
    }

    return sb;
}

// sample a gamma factor from the relativistic Maxwell-Juttner distribution
// given a supplied dimensionless electron temperature, using rejection and numerical inverse CDF methods
double draw_gam_rej(double Theta) {
    double u, gam;

    for (;;) {
        u = (double)rand()/(double)RAND_MAX;

        gam = brent_mj(1.0, 1.0 + 1000.0*Theta, Theta, u);

        u = (double)rand()/(double)RAND_MAX;

        if (u < sqrt(1.0 - (1.0/(gam*gam))))
            return gam;
    }
}

// find gamma factor such that only 10^-9 of the total electrons are faster, according to Maxwell-Juttner distribution, given T (in Kelvin)
double find_gam_max(double T) {
    double *gam_grid, *gam_pdf;
    double th, norm, cuml, gam_max;

    long i;

    gam_grid = NULL;
    while (gam_grid == NULL)
        gam_grid = malloc(GAM_M_N * sizeof(*gam_grid));
    gam_pdf = NULL;
    while (gam_pdf == NULL)
        gam_pdf = malloc(GAM_M_N * sizeof(*gam_pdf));

    th = K_TO_THETA * T;

    for (i = 0; i < GAM_M_N; i++)
        gam_grid[i] = exp(i*(log(1.0e5)/(GAM_M_N-1)));

    for (i = 0; i < GAM_M_N; i++) {
        gam_pdf[i] = gam_grid[i]*gam_grid[i] * sqrt(1.0 - (1.0/(gam_grid[i]*gam_grid[i]))) * exp(-gam_grid[i]/th);
        if (gam_pdf[i] != gam_pdf[i])
            gam_pdf[i] = 0.0;
    }

    norm = 0.0;
    for (i = 0; i < GAM_M_N-1; i++)
        norm += 0.5*(gam_pdf[i+1] + gam_pdf[i])*(gam_grid[i+1] - gam_grid[i]);

    cuml = 0.0;
    for (i = 1; i < GAM_M_N-1; i++) {
        cuml += 0.5*(gam_pdf[i-1] + gam_pdf[i])*(gam_grid[i] - gam_grid[i-1])/norm;
        if (cuml > 1.0 - 1.0e-9)
            break;
    }

    gam_max = gam_grid[i];

    free(gam_grid);
    free(gam_pdf);

    return gam_max;
}

// construct Maxwell-Juttner CDF for given T (Kelvin), return logarithmically-spaced grid points and CDF
void make_gam_cdf(double T, double *gam_grid, double *gam_cdf) {
    double *gam_pdf;
    double th, gam_max, norm;

    long i;

    gam_pdf = NULL;
    while (gam_pdf == NULL)
        gam_pdf = malloc(GAM_N * sizeof(*gam_pdf));

    th = K_TO_THETA * T;

    gam_max = find_gam_max(T);
    for (i = 0; i < GAM_N; i++)
        gam_grid[i] = exp(i*(log(gam_max)/(GAM_N-1)));

    for (i = 0; i < GAM_N; i++) {
        gam_pdf[i] = gam_grid[i]*gam_grid[i] * sqrt(1.0 - (1.0/(gam_grid[i]*gam_grid[i]))) * exp(-gam_grid[i]/th);
        if (gam_pdf[i] != gam_pdf[i])
            gam_pdf[i] = 0.0;
    }

    norm = 0.0;
    for (i = 0; i < GAM_N-1; i++)
        norm += 0.5*(gam_pdf[i+1] + gam_pdf[i])*(gam_grid[i+1] - gam_grid[i]);

    for (i = 0; i < GAM_N; i++)
        gam_pdf[i] /= norm;

    gam_cdf[0] = 0.0;
    for (i = 1; i < GAM_N; i++) {
        gam_cdf[i] = gam_cdf[i-1] + 0.5*(gam_pdf[i] + gam_pdf[i-1])*(gam_grid[i] - gam_grid[i-1]);
        if (gam_cdf[i] > 1.0)
            gam_cdf[i] = 1.0;
    }
    gam_cdf[GAM_N-1] = 1.0;

    free(gam_pdf);
}

// directly invert the gridded, supplied Maxwell-Juttner CDF via linear interpolation in order to sample gamma
// NOTE: currently unused in favor of rejection sampling
double draw_gam_cdf(double *gam_grid, double *gam_cdf) {
    double u;

    long i, i_lb, i_ub;

    u = (double)rand()/(double)RAND_MAX;

    if (u < gam_cdf[GAM_N/2]) {
        i_lb = 0;
        i_ub = GAM_N/2;
    } else {
        i_lb = GAM_N/2;
        i_ub = GAM_N;
    }
    for (;;) {
        if (i_ub - i_lb == 1) {
            i = i_lb;
            break;
        } else {
            if (u < gam_cdf[(i_lb + i_ub)/2]) {
                i_ub = (i_lb + i_ub)/2;
            } else {
                i_lb = (i_lb + i_ub)/2;
            }
        }
    }

    return ((gam_grid[i+1] - gam_grid[i])/(gam_cdf[i+1] - gam_cdf[i]))*(u - gam_cdf[i]) + gam_grid[i];
}

// sample the nonrelativistic Maxwell-Boltzmann distribution for B = v/c, given T (in Kelvin),
// using the Box-Muller transform
double draw_B(double T) {
    double a, u1, u2, B, Bx, By, Bz;

    a = sqrt(K_TO_THETA * T);

    B = 1.0;
    while (B >= 1.0) {
        u1 = (double)rand()/(double)RAND_MAX;
        u2 = (double)rand()/(double)RAND_MAX;
        Bx = a*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);

        u1 = (double)rand()/(double)RAND_MAX;
        u2 = (double)rand()/(double)RAND_MAX;
        By = a*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);

        u1 = (double)rand()/(double)RAND_MAX;
        u2 = (double)rand()/(double)RAND_MAX;
        Bz = a*sqrt(-2.0*log(u1))*cos(2.0*M_PI*u2);

        B = sqrt(Bx*Bx + By*By + Bz*Bz);
    }

    return B;
}

// sample the *auxiliary* scattering angle mu (in cosine), *given* a properly-sampled electron velocity beta = v/c;
// the pre-scatter photon energy enters via the rejection critera, not here
double draw_mu(double beta) {
    double u, mu;

    u = (double)rand()/(double)RAND_MAX;

    if (beta < 0.001) {
        mu = (2.0*u - 1.0) + 2.0*beta*u*(u - 1.0);
    } else {
        mu = (1.0 - sqrt(beta*beta - 4.0*beta*u + 2.0*beta + 1.0))/beta;
    }
    if (mu < -1.0) {
        return -1.0;
    } else if (mu > 1.0) {
        return 1.0;
    } else {
        return mu;
    }
}

// properly sample an electron's pre-scatter velocity and angle with respect to the pre-scatter photon's trajectory:
// input: electron temperature T (in Kelvin), photon energy E (in eV) --- in the fluid rest frame
// output: sampled, pre-scatter electron beta (= v/c), gam (= gamma factor), and mu (= cosine of angle w.r.t. pre-scatter photon's trajectory)
// as described in accompanying literature, this method accurately accounts for the foreshortening of the apparent electron cross section
// due to special relativistic effects
void sample_electron(double T, double E, double *beta, double *gam, double *mu) {
    double u;

    for (;;) {
        // draw an electron speed
        if (T > 1.0e7) {
            *gam  = draw_gam_rej(K_TO_THETA * T);
            *beta = sqrt(1.0 - (1.0/((*gam)*(*gam))));
        } else {
            *beta = draw_B(T);
            *gam = 1.0/sqrt(1.0 - (*beta)*(*beta));
        }

        // draw angle between initial photon direction and initial electron direction
        *mu = draw_mu(*beta);

        // check rejection criteion
        u = (double)rand()/(double)RAND_MAX;
        if (u < (1.0/sig_adjust((*gam)*(1.0 - (*beta)*(*mu))*E))) {
            return;
        }
    }
}

// perform Lorentz boost into the sampled electron's rest frame
void boost(double gam, double B, double A[4], double A_p[4]) {
    A_p[0] = gam*A[0] - B*gam*A[3];
    A_p[1] = A[1];
    A_p[2] = A[2];
    A_p[3] = -B*gam*A[0] + gam*A[3];
}

// align pre-scatter photon's trajectory (in erf) with local z-axis
void rotate(double th, double A[4], double A_p[4]) {
    A_p[0] = A[0];
    A_p[1] = cos(th)*A[1] - sin(th)*A[3];
    A_p[2] = 0.0;
    A_p[3] = sin(th)*A[1] + cos(th)*A[3];
}

// simulates MC_N n-scatter events for an electron temperature T (in Kelvin) and pre-scatter photon energy E_0 (in eV);
// post-scatter energy and fluid rest frame scattering angle are stored in E_record and mu_record, respectively;
// returns also the statistically-determined n-scatter mean amplification ratio
double mc_scatters(int n, double T, double E_0, double *E_record ,double *mu_record) {
    double k_0[4], k_i[4], k_i_erf[4], k_i_erf_rot[4], k_f_erf_rot[4], k_f_erf[4], k_f[4];
    double E_i, mean, mu_i, th_i, gam, B, E_i_erf, th_i_erf, mu_s, sin_th_s, phi_s, E_s, E_f;

    long i, j;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    mean = 0.0;
    for (i = 0; i < MC_N; i++) {
        E_i = E_0;
        for (j = 0; j < n; j++) {
            sample_electron(T, E_i, &B, &gam, &mu_i);

            th_i = acos(mu_i);

            // initial photon 4-vector (initial electron direction set to z-axis)
            k_i[0] = E_i;
            k_i[1] = E_i * sin(th_i);
            k_i[2] = 0.0;
            k_i[3] = E_i * cos(th_i);

            // save first set for later
            if (j == 0) {
                k_0[0] = k_i[0];
                k_0[1] = k_i[1];
                k_0[2] = k_i[2];
                k_0[3] = k_i[3];
            }

            // boost to electron rest frame (e.r.f.)
            boost(gam, B, k_i, k_i_erf);

            // photon energy and angle (with respect to local z-axis, i.e., the direction of boost)
            E_i_erf  = k_i_erf[0];
            th_i_erf = acos(k_i_erf[3]/k_i_erf[0]);

            // rotate axes to align photon 4-vector with z-axis
            rotate(th_i_erf, k_i_erf, k_i_erf_rot);

            // draw scattering angles
            mu_s  = draw_kn(E_i_erf);
            phi_s = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

            // post-scatter photon energy in the e.r.f.
            E_s = E_i_erf/(1.0 + (E_i_erf/511.0e3)*(1.0 - mu_s));

            // post-scatter photon 4-vector in the rotated e.r.f.
            sin_th_s = sqrt(1.0 - mu_s*mu_s);
            k_f_erf_rot[0] = E_s;
            k_f_erf_rot[1] = E_s * sin_th_s * cos(phi_s);
            k_f_erf_rot[2] = E_s * sin_th_s * sin(phi_s);
            k_f_erf_rot[3] = E_s * mu_s;

            // un-rotate the axes
            rotate(-th_i_erf, k_f_erf_rot, k_f_erf);

            // reverse boost
            boost(gam, -B, k_f_erf, k_f);

            // final post-scatter photon energy in lab frame
            E_f = k_f[0];

            E_i = E_f;
        }
        // record result
        E_record[i]  = E_f;
        mu_record[i] = (1.0/(E_0 * E_f)) * (k_0[1] * k_f[1] + k_0[2] * k_f[2] + k_0[3] * k_f[3]);

        mean += E_f;
    }
    mean = (mean/MC_N)/E_0;

    return mean;
}

// construct Maxwell-Juttner PDF for given T (Kelvin), return logarithmically-spaced grid points and PDF
void make_gam_pdf(double T, double *gam, double *dgam, double *gam_pdf, double *B) {
    double th, gam_max, norm;

    long i;

    th = (8.617e-5 * T)/511.0e3;

    gam_max = find_gam_max(T);
    for (i = 0; i < GAM_N_SA; i++)
        gam[i] = exp(i*(log(gam_max)/(GAM_N_SA-1)));

    for (i = 1; i < GAM_N_SA-1; i++)
        dgam[i] = 0.5*(gam[i+1] - gam[i-1]);
    dgam[0] = 0.5*(gam[1] - gam[0]);
    dgam[GAM_N_SA-1] = 0.5*(gam[GAM_N_SA-1] - gam[GAM_N_SA-2]);

    for (i = 0; i < GAM_N_SA; i++) {
        gam_pdf[i] = gam[i]*gam[i] * sqrt(1.0 - (1.0/(gam[i]*gam[i]))) * exp(-gam[i]/th);
        if (gam_pdf[i] != gam_pdf[i])
            gam_pdf[i] = 0.0;
    }

    norm = 0.0;
    for (i = 0; i < GAM_N_SA; i++)
        norm += gam_pdf[i] * dgam[i];

    for (i = 0; i < GAM_N_SA; i++)
        gam_pdf[i] /= norm;

    for (i = 0; i < GAM_N_SA; i++) {
        B[i] = sqrt(1.0 - (1.0/(gam[i]*gam[i])));
        if (B[i] != B[i])
            B[i] = 0.0;
    }
}

// find beta = v/c such that only 10^-9 of the total electrons are faster, according to Maxwell-Boltzmann distribution, given T (in Kelvin)
double find_B_max(double T) {
    double *B_grid, *B_pdf;
    double th, norm, cuml, B_max;

    long i;

    B_grid = NULL;
    while (B_grid == NULL)
        B_grid = malloc(GAM_M_N * sizeof(*B_grid));
    B_pdf = NULL;
    while (B_pdf == NULL)
        B_pdf = malloc(GAM_M_N * sizeof(*B_pdf));

    th = (8.617e-5 * T)/511.0e3;

    for (i = 0; i < GAM_M_N; i++)
        B_grid[i] = exp(i*(log(1.0/1.0e-6)/(GAM_M_N-1)) + log(1.0e-6));
    B_grid[0] = 0.0;

    for (i = 0; i < GAM_M_N; i++) {
        B_pdf[i] = B_grid[i]*B_grid[i] * exp(-B_grid[i]*B_grid[i]/(2.0*th));
        if (B_pdf[i] != B_pdf[i])
            B_pdf[i] = 0.0;
    }

    norm = 0.0;
    for (i = 0; i < GAM_M_N-1; i++)
        norm += 0.5*(B_pdf[i+1] + B_pdf[i])*(B_grid[i+1] - B_grid[i]);

    cuml = 0.0;
    for (i = 1; i < GAM_M_N-1; i++) {
        cuml += 0.5*(B_pdf[i-1] + B_pdf[i])*(B_grid[i] - B_grid[i-1])/norm;
        if (cuml > 1.0 - 1.0e-9)
            break;
    }

    B_max = B_grid[i];

    free(B_grid);
    free(B_pdf);

    return B_max;
}

// construct Maxwell-Juttner PDF for given T (Kelvin), return logarithmically-spaced grid points and PDF
void make_B_pdf(double T, double *B, double *dB, double *B_pdf, double *gam) {
    double th, B_max, norm;

    long i;

    th = (8.617e-5 * T)/511.0e3;

    B_max = find_B_max(T);
    for (i = 0; i < GAM_N_SA; i++)
        B[i] = i*(B_max/(GAM_N_SA-1));

    for (i = 1; i < GAM_N_SA-1; i++)
        dB[i] = 0.5*(B[i+1] - B[i-1]);
    dB[0] = 0.5*(B[1] - B[0]);
    dB[GAM_N_SA-1] = 0.5*(B[GAM_N_SA-1] - B[GAM_N_SA-2]);

    for (i = 0; i < GAM_N_SA; i++) {
        B_pdf[i] = B[i]*B[i] * exp(-B[i]*B[i]/(2.0*th));
        if (B_pdf[i] != B_pdf[i])
            B_pdf[i] = 0.0;
    }

    norm = 0.0;
    for (i = 0; i < GAM_N_SA; i++)
        norm += B_pdf[i] * dB[i];

    for (i = 0; i < GAM_N_SA; i++)
        B_pdf[i] /= norm;

    for (i = 0; i < GAM_N_SA; i++) {
        gam[i] = 1.0/sqrt(1.0 - B[i]*B[i]);
        if (gam[i] != gam[i])
            gam[i] = 1.0;
    }
}

// first moment function in e.r.f. amplification integral
double mf1(double e) {
    if (e < 0.001)
        return 1.0 - e + 2.2*e*e;
    return ((2*e*(-3 + e*(-15 + 2*e*(-9 + e*(3 + 8*e)))))/pow(1 + 2*e,3) + 3*log(1 + 2*e))/(3.*((2*e*(2 + e*(1 + e)*(8 + e)))/pow(1 + 2*e,2) + (-2 + (-2 + e)*e)*log(1 + 2*e)));
}

// second moment function in e.r.f. amplification integral
double mf2(double e) {
    if (e < 0.001)
        return 1.2*e - 2.9*e*e;
    return ((2*e*(-3 + e*(1 + e)*(-9 + 4*e))*(3 + e*(9 + 4*e)))/pow(1 + 2*e,3) + (9 - 3*(-3 + e)*e)*log(1 + 2*e))/(3.*e*((2*e*(2 + e*(1 + e)*(8 + e)))/pow(1 + 2*e,2) + (-2 + (-2 + e)*e)*log(1 + 2*e)));
}

// "semi-analytic" calculation of the thermally- and angle-averaged mean amplification ratio
double sa_calc_ratio(double T, double E_i, double *mu_i, double *w) {
    double gam[GAM_N_SA], dgam[GAM_N_SA], B[GAM_N_SA], dB[GAM_N_SA], gam_pdf[GAM_N_SA], B_pdf[GAM_N_SA], th_i[MU_I_N];
    double k_i[4], k_i_erf[4];
    double E_i_erf, mu_i_erf, E_avg;
    double num, den, wght;

    long i, j;

    if (T > 1.0e7) {
        make_gam_pdf(T, gam, dgam, gam_pdf, B);
    } else {
        make_B_pdf(T, B, dB, B_pdf, gam);
    }

    for (i = 0; i < MU_I_N; i++) {
        th_i[i] = acos(mu_i[i]);
    }

    num = 0.0;
    den = 0.0;
    for (i = 0; i < GAM_N_SA; i++) {
        for (j = 0; j < MU_I_N; j++) {
            // initial photon 4-vector (initial electron direction set to z-axis)
            k_i[0] = E_i;
            k_i[1] = E_i * sin(th_i[j]);
            k_i[2] = 0.0;
            k_i[3] = E_i * cos(th_i[j]);

            // boost to electron rest frame (e.r.f.)
            boost(gam[i], B[i], k_i, k_i_erf);

            // photon energy and angle (with respect to local z-axis, i.e., the direction of boost)
            E_i_erf  = k_i_erf[0];
            mu_i_erf = k_i_erf[3]/k_i_erf[0];

            // properly averaged post-scatter photon energy in the e.r.f.
            E_avg = E_i_erf * gam[i] * (mf1(E_i_erf/511.0e3) + B[i]*mu_i_erf*mf2(E_i_erf/511.0e3));

            // contribute to average
            if (T > 1.0e7) {
                wght = w[j] * gam_pdf[i] * dgam[i] * (1.0 - B[i] * mu_i[j]) * (1.0/sig_adjust((gam[i])*(1.0 - (B[i])*(mu_i[j]))*E_i));
            } else {
                wght = w[j] * B_pdf[i] * dB[i] * (1.0 - B[i] * mu_i[j]) * (1.0/sig_adjust((gam[i])*(1.0 - (B[i])*(mu_i[j]))*E_i));
            }
            num += wght * E_avg;
            den += wght;
        }
    }

    return (num/den)/E_i;
}

// the e.r.f. value of the scattering angle-integrated KN cross section, in units of the "classical electron cross section," i.e., pi * r_e^2
double sigma(double e) {
    double eps = e/511.0e3;

    if (e < 100.0) {
        return (8.0/3.0 - (16.0/3.0)*eps + (208.0/15.0)*eps*eps - (522.0/15.0)*eps*eps*eps);
    } else {
        return (((2.0*eps*(2.0 + eps*(1.0 + eps)*(8.0 + eps)))/((1.0 + 2.0*eps)*(1.0 + 2.0*eps)) + (-2.0 + (-2.0 + eps)*eps)*log(1.0 + 2.0*eps))/(eps*eps*eps));
    }
}

// "semi-analytic" calculation of the thermally- and angle-averaged total scattering opacity, in units of the Thomson cross section
double sa_calc_sigma(double T, double E_i, double *mu_i, double *w) {
    double gam[GAM_N_SA], dgam[GAM_N_SA], B[GAM_N_SA], dB[GAM_N_SA], gam_pdf[GAM_N_SA], B_pdf[GAM_N_SA], th_i[MU_I_N];
    double num, den, integrand, wght;

    long i, j;

    if (T > 1.0e7) {
        make_gam_pdf(T, gam, dgam, gam_pdf, B);
    } else {
        make_B_pdf(T, B, dB, B_pdf, gam);
    }

    for (i = 0; i < MU_I_N; i++) {
        th_i[i] = acos(mu_i[i]);
    }

    num = 0.0;
    den = 0.0;
    for (i = 0; i < GAM_N_SA; i++) {
        for (j = 0; j < MU_I_N; j++) {
            // contribute to integral
            integrand = 0.5 * (1.0 - B[i] * mu_i[j]) * sigma((gam[i])*(1.0 - (B[i])*(mu_i[j]))*E_i);
            if (T > 1.0e7) {
                wght = w[j] * gam_pdf[i] * dgam[i];
            } else {
                wght = w[j] * B_pdf[i] * dB[i];
            }
            num += wght * integrand;
            den += wght;
        }
    }

    return num/(8.0/3.0);
}

// simulate a single scattering event for an electron of initial (fluid frame) energy E_i (eV)
// from an electron at temperature T (K); store final energy in E_f (eV) and
// store (fluid frame) (cosine) scattering angle in mu
void single_scatter(double T, double E_i, double *E_f, double *mu) {
    double k_i[4], k_i_erf[4], k_i_erf_rot[4], k_f_erf_rot[4], k_f_erf[4], k_f[4];
    double mu_i, th_i, gam, B, E_i_erf, th_i_erf, mu_s, sin_th_s, phi_s, E_s;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    sample_electron(T, E_i, &B, &gam, &mu_i);

    th_i = acos(mu_i);

    // initial photon 4-vector (initial electron direction set to z-axis)
    k_i[0] = E_i;
    k_i[1] = E_i * sin(th_i);
    k_i[2] = 0.0;
    k_i[3] = E_i * cos(th_i);

    // boost to electron rest frame (e.r.f.)
    boost(gam, B, k_i, k_i_erf);

    // photon energy and angle (with respect to local z-axis, i.e., the direction of boost)
    E_i_erf  = k_i_erf[0];
    th_i_erf = acos(k_i_erf[3]/k_i_erf[0]);

    // rotate axes to align photon 4-vector with z-axis
    rotate(th_i_erf, k_i_erf, k_i_erf_rot);

    // draw scattering angles
    mu_s  = draw_kn(E_i_erf);
    phi_s = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

    // post-scatter photon energy in the e.r.f.
    E_s = E_i_erf/(1.0 + (E_i_erf/511.0e3)*(1.0 - mu_s));

    // post-scatter photon 4-vector in the rotated e.r.f.
    sin_th_s = sqrt(1.0 - mu_s*mu_s);
    k_f_erf_rot[0] = E_s;
    k_f_erf_rot[1] = E_s * sin_th_s * cos(phi_s);
    k_f_erf_rot[2] = E_s * sin_th_s * sin(phi_s);
    k_f_erf_rot[3] = E_s * mu_s;

    // un-rotate the axes
    rotate(-th_i_erf, k_f_erf_rot, k_f_erf);

    // reverse boost
    boost(gam, -B, k_f_erf, k_f);

    // final post-scatter photon energy in lab frame
    *E_f = k_f[0];

    // record result
    *mu = (1.0/(E_i * k_f[0])) * (k_i[1] * k_f[1] + k_i[2] * k_f[2] + k_i[3] * k_f[3]);
}
