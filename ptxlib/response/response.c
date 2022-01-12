#include "response.h"

double sig_adjust(double e) {
    double eps = e/511.0e3;

    if (e < 100.0) {
        return (8.0/3.0)/(8.0/3.0 - (16.0/3.0)*eps + (208.0/15.0)*eps*eps - (522.0/15.0)*eps*eps*eps);
    } else {
        return (8.0/3.0)/(((2.0*eps*(2.0 + eps*(1.0 + eps)*(8.0 + eps)))/((1.0 + 2.0*eps)*(1.0 + 2.0*eps)) + (-2.0 + (-2.0 + eps)*eps)*log(1.0 + 2.0*eps))/(eps*eps*eps));
    }
}

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

double draw_kn(double e) {
//  return acos(brent(-1, 1, e/511.0e3, (double)rand()/(double)RAND_MAX));
    return brent_kn(-1.0, 1.0, e/511.0e3, (double)rand()/(double)RAND_MAX);
}

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

    th = (8.617e-5 * T)/511.0e3;

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

void make_gam_cdf(double T, double *gam_grid, double *gam_cdf) {
    double *gam_pdf;
    double th, gam_max, norm;

    long i;

    gam_pdf = NULL;
    while (gam_pdf == NULL)
        gam_pdf = malloc(GAM_N * sizeof(*gam_pdf));

    th = (8.617e-5 * T)/511.0e3;

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

double draw_gam(double *gam_grid, double *gam_cdf) {
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

double draw_B(double T) {
    double a, u1, u2, B, Bx, By, Bz;

    a = sqrt((8.617e-5 * T)/511.0e3);

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

void boost(double gam, double B, double A[4], double A_p[4]) {
    A_p[0] = gam*A[0] - B*gam*A[3];
    A_p[1] = A[1];
    A_p[2] = A[2];
    A_p[3] = -B*gam*A[0] + gam*A[3];
}

void rotate(double th, double A[4], double A_p[4]) {
    A_p[0] = A[0];
    A_p[1] = cos(th)*A[1] - sin(th)*A[3];
    A_p[2] = 0.0;
    A_p[3] = sin(th)*A[1] + cos(th)*A[3];
}

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

double mc_scatters(double T, double E_i, double *E_record ,double *mu_record) {
//  double *gam_grid, *gam_cdf;
    double k_i[4], k_i_erf[4], k_i_erf_rot[4], k_f_erf_rot[4], k_f_erf[4], k_f[4];
    double mean, wght, den, mu_i, th_i, gam, B, E_i_erf, th_i_erf, mu_s, sin_th_s, phi_s, E_s, E_f;

    long i, ndx;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    if (T > 1.0e7) {
        gam_grid = NULL;
        while (gam_grid == NULL)
            gam_grid = malloc(GAM_N * sizeof(*gam_grid));
        gam_cdf = NULL;
        while (gam_cdf == NULL)
            gam_cdf = malloc(GAM_N * sizeof(*gam_cdf));
        make_gam_cdf(T, gam_grid, gam_cdf);
    }
    */

    mean = 0.0;
    for (i = 0; i < MC_N; i++) {
//      if (i % 1000 == 0)
//          printf("%ld\n", i);

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
//      th_s  = draw_kn(E_i_erf);
        mu_s  = draw_kn(E_i_erf);
        phi_s = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

        // post-scatter photon energy in the e.r.f.
        E_s = E_i_erf/(1.0 + (E_i_erf/511.0e3)*(1.0 - mu_s));

        // post-scatter photon 4-vector in the rotated e.r.f.
        sin_th_s = sqrt(1.0 - mu_s*mu_s);
        k_f_erf_rot[0] = E_s;
        k_f_erf_rot[1] = E_s * sin_th_s * cos(phi_s);
        k_f_erf_rot[2] = E_s * sin_th_s * sin(phi_s);
//      k_f_erf_rot[3] = E_s * cos(th_s);
        k_f_erf_rot[3] = E_s * mu_s;

        // un-rotate the axes
        rotate(-th_i_erf, k_f_erf_rot, k_f_erf);

        // reverse boost
        boost(gam, -B, k_f_erf, k_f);

        // final post-scatter photon energy in lab frame
        E_f = k_f[0];

        // record result
        E_record[i]  = E_f;
        mu_record[i] = (1.0/(E_i * E_f)) * (k_i[1] * k_f[1] + k_i[2] * k_f[2] + k_i[3] * k_f[3]);

        mean += E_f;
    }
    mean = (mean/MC_N)/E_i;

    /*
    if (T > 1.0e7) {
        free(gam_grid);
        free(gam_cdf);
    }
    */

    return mean;
}

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
        if (cuml > 1.0 - 1.0e-8)
            break;
    }

    B_max = B_grid[i];

    free(B_grid);
    free(B_pdf);

    return B_max;
}

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

double mf1(double e) {
    if (e < 0.001)
        return 1.0 - e + 2.2*e*e;
    return ((2*e*(-3 + e*(-15 + 2*e*(-9 + e*(3 + 8*e)))))/pow(1 + 2*e,3) + 3*log(1 + 2*e))/(3.*((2*e*(2 + e*(1 + e)*(8 + e)))/pow(1 + 2*e,2) + (-2 + (-2 + e)*e)*log(1 + 2*e)));
}

double mf2(double e) {
    if (e < 0.001)
        return 1.2*e - 2.9*e*e;
    return ((2*e*(-3 + e*(1 + e)*(-9 + 4*e))*(3 + e*(9 + 4*e)))/pow(1 + 2*e,3) + (9 - 3*(-3 + e)*e)*log(1 + 2*e))/(3.*e*((2*e*(2 + e*(1 + e)*(8 + e)))/pow(1 + 2*e,2) + (-2 + (-2 + e)*e)*log(1 + 2*e)));
}

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

double sigma(double e) {
    double eps = e/511.0e3;

    if (e < 100.0) {
        return (8.0/3.0 - (16.0/3.0)*eps + (208.0/15.0)*eps*eps - (522.0/15.0)*eps*eps*eps);
    } else {
        return (((2.0*eps*(2.0 + eps*(1.0 + eps)*(8.0 + eps)))/((1.0 + 2.0*eps)*(1.0 + 2.0*eps)) + (-2.0 + (-2.0 + eps)*eps)*log(1.0 + 2.0*eps))/(eps*eps*eps));
    }
}

double sa_calc_sigma(double T, double E_i, double *mu_i, double *w) {
    double gam[GAM_N_SA], dgam[GAM_N_SA], B[GAM_N_SA], dB[GAM_N_SA], gam_pdf[GAM_N_SA], B_pdf[GAM_N_SA], th_i[MU_I_N];
    double k_i[4], k_i_erf[4];
    double E_i_erf, mu_i_erf, E_avg;
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

void single_scatter(double T, double E_i, double *E_f, double *mu) {
//  double *gam_grid, *gam_cdf;
    double k_i[4], k_i_erf[4], k_i_erf_rot[4], k_f_erf_rot[4], k_f_erf[4], k_f[4];
    double mu_i, th_i, gam, B, E_i_erf, th_i_erf, mu_s, sin_th_s, phi_s, E_s;

    long i, ndx;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    if (T > 1.0e7) {
        gam_grid = NULL;
        while (gam_grid == NULL)
            gam_grid = malloc(GAM_N * sizeof(*gam_grid));
        gam_cdf = NULL;
        while (gam_cdf == NULL)
            gam_cdf = malloc(GAM_N * sizeof(*gam_cdf));
        make_gam_cdf(T, gam_grid, gam_cdf);
    }
    */

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
//  th_s  = draw_kn(E_i_erf);
    mu_s  = draw_kn(E_i_erf);
    phi_s = 2.0*M_PI*((double)rand()/(double)RAND_MAX);

    // post-scatter photon energy in the e.r.f.
    E_s = E_i_erf/(1.0 + (E_i_erf/511.0e3)*(1.0 - mu_s));

    // post-scatter photon 4-vector in the rotated e.r.f.
    sin_th_s = sqrt(1.0 - mu_s*mu_s);
    k_f_erf_rot[0] = E_s;
    k_f_erf_rot[1] = E_s * sin_th_s * cos(phi_s);
    k_f_erf_rot[2] = E_s * sin_th_s * sin(phi_s);
//  k_f_erf_rot[3] = E_s * cos(th_s);
    k_f_erf_rot[3] = E_s * mu_s;

    // un-rotate the axes
    rotate(-th_i_erf, k_f_erf_rot, k_f_erf);

    // reverse boost
    boost(gam, -B, k_f_erf, k_f);

    // final post-scatter photon energy in lab frame
    *E_f = k_f[0];

    // record result
    *mu = (1.0/(E_i * k_f[0])) * (k_i[1] * k_f[1] + k_i[2] * k_f[2] + k_i[3] * k_f[3]);

    /*
    if (T > 1.0e7) {
        free(gam_grid);
        free(gam_cdf);
    }
    */
}

/*
int main() {
    double cdf[N];
    double T = 1.0e13;
    double E = 1.0e6;

    printf("%e\n", mc_scatters(T, E, cdf));
    printf("%e\n", sa_calc_ratio(T, E));

    return 0;
}
*/

// void simulate_scatter(double E_0, double th_0, double phi_0, double T, double *gam_grid, double *gam_cdf, double *E_f, double *th_f, double *phi_f) {
void simulate_scatter(double E_0, double th_0, double phi_0, double T, double *E_f, double *th_f, double *phi_f) {
    double gam, B, mu_e, th_e, phi_e, kt, kx, ky, kz, th, phi, th_s, phi_s, kst, ksx, ksy, ksz, kft, kfx, kfy, kfz;

    /*
    if (T > 1.0e7) {
//      gam = draw_gam(gam_grid, gam_cdf);
        gam = draw_gam_rej(K_TO_THETA * T);
        B   = sqrt(1.0 - (1.0/(gam*gam)));
    } else {
        B = draw_B(T);
    }
    */

    sample_electron(T, E_0, &B, &gam, &mu_e);

    th_e  = acos(mu_e);
    phi_e = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;

//  th_e  = acos(2.0 * (double)rand()/(double)RAND_MAX - 1.0); // the electron's direction in the "lab" frame is drawn from an isotropic distribution
//  phi_e = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;

    // in the rest frame of the electron (ERF), the components of the pre-scatter photon's 4-momentum are (with 1/c factored out, so kt is just the energy in this frame)
    kt = (E_0/sqrt(1.0 - B*B))*(1.0 - B*cos(th_0)*cos(th_e) - B*cos(phi_0-phi_e)*sin(th_0)*sin(th_e));
    kx = (E_0/sqrt(1.0 - B*B))*(cos(th_0)*cos(th_e) + cos(phi_0-phi_e)*sin(th_0)*sin(th_e) - B);
    ky = E_0*sin(th_0)*sin(phi_0-phi_e);
    kz = E_0*(cos(th_0)*sin(th_e) - cos(th_e)*cos(phi_0-phi_e)*sin(th_0));

    // so, in the ERF, the pre-scatter photon has angle theta
    th = acos(kz/sqrt(kx*kx + ky*ky + kz*kz));
    // and angle phi (ccw from +x axis)
    if ((kx > 0.0) && (ky > 0.0)) {
        // QI
        phi = atan(ky/kx);
    }
    else if ((kx < 0.0) && (ky > 0.0)) {
        // QII
        phi = M_PI - atan(ky/fabs(kx));
    }
    else if ((kx < 0.0) && (ky < 0.0)) {
        // QIII
        phi = M_PI + atan(ky/kx);
    }
    else {
        // QIV
        phi = 2.0*M_PI - atan(fabs(ky)/kx);
    }

    // the electron undergoes scattering with angles drawn from the KN cross section
    th_s  = acos(draw_kn(kt));
    phi_s = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;

    // the scattering angle th_s determines the post-scatter energy, kst, in the ERF
    // the spatial components are derived from careful rotations (see scatter.nb)
    kst = kt/(1.0 + (1.956951e-6 * kt)*(1.0 - cos(th_s)));

    // REMOVE THIS LINE!
//  th_s = acos(2.0*((double)rand()/(double)RAND_MAX) - 1.0);

    ksx = kst*(cos(th_s)*cos(phi)*sin(th) + sin(th_s)*(cos(th)*cos(phi)*cos(phi_s) - sin(phi)*sin(phi_s)));
    ksy = kst*(cos(th_s)*sin(th)*sin(phi) + sin(th_s)*(cos(th)*cos(phi_s)*sin(phi) + cos(phi)*sin(phi_s)));
    ksz = kst*(cos(th)*cos(th_s) - cos(phi_s)*sin(th)*sin(th_s));

    // boost and re-rotate back to the lab frame to get the components of the post-scatter photon's 4-momentum
    kft = (kst + B*ksx)/sqrt(1.0 - B*B); // the final, new photon energy
    kfx = (1.0/sqrt(1.0 - B*B))*((ksx + B*kst)*cos(phi_e)*sin(th_e)) - ksy*sin(phi_e) - ksz*cos(th_e)*cos(phi_e);
    kfy = (1.0/sqrt(1.0 - B*B))*((ksx + B*kst)*sin(phi_e)*sin(th_e)) + ksy*cos(phi_e) - ksz*cos(th_e)*sin(phi_e);
    kfz = (1.0/sqrt(1.0 - B*B))*(ksx + B*kst)*cos(th_e) + ksz*sin(th_e);

    // translate back into angles (and pass them up)
    *th_f = acos(kfz/sqrt(kfx*kfx + kfy*kfy + kfz*kfz));
    if ((kfx > 0.0) && (kfy > 0.0)) {
        // QI
        *phi_f = atan(kfy/kfx);
    }
    else if ((kfx < 0.0) && (kfy > 0.0)) {
        // QII
        *phi_f = M_PI - atan(kfy/fabs(kfx));
    }
    else if ((kfx < 0.0) && (kfy < 0.0)) {
        // QIII
        *phi_f = M_PI + atan(kfy/kfx);
    }
    else {
        // QIV
        *phi_f = 2.0*M_PI - atan(fabs(kfy)/kfx);
    }

    // pass up energy
    *E_f = kft;
}

void simulate_trajectory(double E_0, double T, double *E_record, double *mu_record) {
//  double *gam_grid, *gam_cdf;
    double dist, z, th, phi, E, E_f, th_f, phi_f;
    int i;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    if (T > 1.0e7) {
        gam_grid = NULL;
        while (gam_grid == NULL)
            gam_grid = malloc(GAM_N * sizeof(*gam_grid));
        gam_cdf = NULL;
        while (gam_cdf == NULL)
            gam_cdf = malloc(GAM_N * sizeof(*gam_cdf));
        make_gam_cdf(T, gam_grid, gam_cdf);
    }
    */

    z   = 0.0;
    th  = acos((double)rand()/(double)RAND_MAX);
    phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
    E   = E_0;

    i = 0;
    while (i < TOTAL_PHOTONS) {
        dist = -log((double)rand()/(double)RAND_MAX);

        z += sig_adjust(E) * dist * cos(th);

        if (z < 0.0) {
            z   = 0.0;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = E_0;
        } else if (z > TOTAL_THICKNESS) {
            E_record[i]  = E;
            mu_record[i] = cos(th);
            i++;
            printf("%d\n", i);
            z   = 0.0;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = E_0;
        } else {
//          simulate_scatter(E, th, phi, T, gam_grid, gam_cdf, &E_f, &th_f, &phi_f);
            simulate_scatter(E, th, phi, T, &E_f, &th_f, &phi_f);
            E   = E_f;
            th  = th_f;
            phi = phi_f;
        }
    }

    /*
    if (T > 1.0e7) {
        free(gam_grid);
        free(gam_cdf);
    }
    */
}

double draw_E_pl(double Gamma, double E_min, double E_max) {
    double u;

    u = (double)rand()/(double)RAND_MAX;

    return pow(u * pow(E_max, 1.0-Gamma) + (1.0 - u) * pow(E_min, 1.0-Gamma), 1.0/(1.0-Gamma));
}

void simulate_trajectory_pl(double Gamma, double E_min, double E_max, double T, double *E_record, double *mu_record) {
//  double *gam_grid, *gam_cdf;
    double dist, z, th, phi, E, E_f, th_f, phi_f;
    int i;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    if (T > 1.0e7) {
        gam_grid = NULL;
        while (gam_grid == NULL)
            gam_grid = malloc(GAM_N * sizeof(*gam_grid));
        gam_cdf = NULL;
        while (gam_cdf == NULL)
            gam_cdf = malloc(GAM_N * sizeof(*gam_cdf));
        make_gam_cdf(T, gam_grid, gam_cdf);
    }
    */

    z   = 0.0;
    th  = acos((double)rand()/(double)RAND_MAX);
    phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
    E   = draw_E_pl(Gamma, E_min, E_max);

    i = 0;
    while (i < TOTAL_PHOTONS) {
        dist = -log((double)rand()/(double)RAND_MAX);

        z += sig_adjust(E) * dist * cos(th);

        if (z < 0.0) {
            E_record[i]  = E;
            mu_record[i] = cos(th);
            i++;
            if (i % 1000 == 0) {
                printf("%.2f%%\n", 100.0 * (double)i/(double)(TOTAL_PHOTONS));
//              printf("sig_adjust(%e) = %e\n", E, sig_adjust(E));
            }
            z   = 0.0;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = draw_E_pl(Gamma, E_min, E_max);
        } else if (z > TOTAL_THICKNESS) {
            E_record[i]  = E;
            mu_record[i] = cos(th);
            i++;
            if (i % 1000 == 0) {
                printf("%.2f%%\n", 100.0 * (double)i/(double)(TOTAL_PHOTONS));
//              printf("sig_adjust(%e) = %e\n", E, sig_adjust(E));
            }
            z   = 0.0;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = draw_E_pl(Gamma, E_min, E_max);
        } else {
//          simulate_scatter(E, th, phi, T, gam_grid, gam_cdf, &E_f, &th_f, &phi_f);
            simulate_scatter(E, th, phi, T, &E_f, &th_f, &phi_f);
            E   = E_f;
            th  = th_f;
            phi = phi_f;
        }
    }

    /*
    if (T > 1.0e7) {
        free(gam_grid);
        free(gam_cdf);
    }
    */
}

void simulate_reflection(double E_0, double T, double *E_record, double *mu_record, double *scatter_record) {
    double dist, z, th, phi, E, E_f, th_f, phi_f;
    int i;
    long scatters;

    struct timeval timer_usec;

    // seed the random number generator with the number of microseconds past the second
    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    if (T > 1.0e7) {
        gam_grid = NULL;
        while (gam_grid == NULL)
            gam_grid = malloc(GAM_N * sizeof(*gam_grid));
        gam_cdf = NULL;
        while (gam_cdf == NULL)
            gam_cdf = malloc(GAM_N * sizeof(*gam_cdf));
        make_gam_cdf(T, gam_grid, gam_cdf);
    }
    */

    z   = 0.0001;
    th  = acos((double)rand()/(double)RAND_MAX);
    phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
    E   = E_0;
    scatters = 0;

    i = 0;
    while (i < TOTAL_PHOTONS) {
        dist = -log((double)rand()/(double)RAND_MAX);

        z += sig_adjust(E) * dist * cos(th);

        if (z < 0.0) {
            E_record[i]       = E;
            mu_record[i]      = cos(th);
            scatter_record[i] = (double) scatters;
            i++;
            printf("%d\n", i);
            z   = 0.0001;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = E_0;
            scatters = 0;
        } else if (z > TOTAL_THICKNESS) {
            E_record[i]       = E;
            mu_record[i]      = cos(th);
            scatter_record[i] = (double) -scatters;
            i++;
            printf("%d\n", i);
            z   = 0.0001;
            th  = acos((double)rand()/(double)RAND_MAX);
            phi = 2.0 * M_PI * (double)rand()/(double)RAND_MAX;
            E   = E_0;
            scatters = 0;
        } else {
            simulate_scatter(E, th, phi, T, &E_f, &th_f, &phi_f);
            E   = E_f;
            th  = th_f;
            phi = phi_f;
            scatters++;
        }
    }

    /*
    if (T > 1.0e7) {
        free(gam_grid);
        free(gam_cdf);
    }
    */
}

int new_z(double E, double *z_pos, double dtau, double *z, double *e, double *tau_grids, int D, int N, int *ndx) {
    double *tau_grid;
    double tau;

    int i, i_lb, i_ub;

    // find index i such that e[i] is nearest to E
    if (E < e[0]) {
        i = 0;
    } else if (E >= e[N-1]) {
        i = N-1;
    } else {
        i = (int)(log(E/e[0])/log(e[1]/e[0]));
        if (E - e[i] > e[i+1] - E)
            i++;
    }

    tau_grid = tau_grids + i*D;

    // find tau corresponding to z_pos
    i   = *ndx;
    tau = ((tau_grid[i+1] - tau_grid[i])/(z[i+1] - z[i]))*(*z_pos - z[i]) + tau_grid[i];

    // increment tau
    tau += dtau;

    // if the photon has escaped through the top, return 1
    if (tau <= tau_grid[0])
        return 1;
    // if through the bottom, return 2
    if (tau >= tau_grid[D-1])
        return 2;

    // otherwise, find new z_pos
    if (tau < tau_grid[D/2]) {
        i_lb = 0;
        i_ub = D/2;
    } else {
        i_lb = D/2;
        i_ub = D;
    }
    for (;;) {
        if (i_ub - i_lb == 1) {
            i = i_lb;
            break;
        } else {
            if (tau < tau_grid[(i_lb + i_ub)/2]) {
                i_ub = (i_lb + i_ub)/2;
            } else {
                i_lb = (i_lb + i_ub)/2;
            }
        }
    }
    *z_pos = ((z[i+1] - z[i])/(tau_grid[i+1] - tau_grid[i]))*(tau - tau_grid[i]) + z[i];

    // the photon is now in cell i; set ndx
    *ndx = i;

    return 0;
}

int response(int top_or_bot, double E_0, double *z, double *e, double *T, double *tau_grids, int D, int N, double *cdf) {
    double *gam_grids, *gam_cdfs, *count;
    double E, th, phi, z_pos, dtau;

    int i, ndx, code;

    long scatters;

    struct timeval timer_usec;

    gettimeofday(&timer_usec, NULL);
    srand(timer_usec.tv_usec);

    /*
    gam_grids = NULL;
    while (gam_grids == NULL)
        gam_grids = malloc((D-1) * GAM_N * sizeof(*gam_grids));
    gam_cdfs = NULL;
    while (gam_cdfs == NULL)
        gam_cdfs = malloc((D-1) * GAM_N * sizeof(*gam_cdfs));

    for (i = 0; i < D-1; i++) {
        if (T[i] > 1.0e7)
            make_gam_cdf(T[i], gam_grids + i*GAM_N, gam_cdfs + i*GAM_N);
    }
    */

    count = NULL;
    while (count == NULL)
        count = malloc(N * sizeof(*count));
    for (i = 0; i < N; i++)
        count[i] = 0.0;

    i        = 0;
    scatters = 0;
    while (i < (int)1e5) {
        E   = E_0;
        phi = 0.0;
        if (top_or_bot > 0) {
            z_pos = z[0];
            th    = acos((double)rand()/(double)RAND_MAX);
            ndx   = 0;
        } else {
            z_pos = z[D-1];
            th    = acos(-(double)rand()/(double)RAND_MAX);
            ndx   = D-2;
        }
        while (1) {
            // draw a random optical depth to traverse
            dtau = -log((double)rand()/(double)RAND_MAX);

            // where is the photon now
            code = new_z(E, &z_pos, cos(th)*dtau, z, e, tau_grids, D, N, &ndx);

            // has the photon escaped?
            if (code > 0)
                break;

            // the photon scatters, get new energy and direction
//          simulate_scatter(E, th, phi, T[ndx], gam_grids + ndx*GAM_N, gam_cdfs + ndx*GAM_N, &E, &th, &phi);
            single_scatter(T[ndx], E, &E, &th);
            th = acos(2.0*((double)rand()/(double)RAND_MAX) - 1.0);
//          printf("%e %e %e %e\n", E, T[ndx], th, phi);
            scatters++;
        }
        // if slab_type == 'upper' and the photon leaves through the bottom, don't count it
        if ((top_or_bot == 2) && (code == 2))
            continue;
        // if slab_type == 'lower' and the photon leaves through the top, don't count it either
        if ((top_or_bot == -2) && (code == 1))
            continue;
        i++;
        ndx = (int)(log(E/e[0])/log(e[1]/e[0]));
        if (ndx >= 0 && ndx < N)
            count[ndx] += 1.0;
    }

//  printf("mean scatters = %f\n", scatters/1.0e5);

    cdf[0] = 0.0;
    for (i = 1; i < N; i++)
        cdf[i] = cdf[i-1] + count[i-1];
    for (i = 0; i < N; i++) {
        cdf[i] /= 1e5;
        if (cdf[i] < 0.0)
            cdf[i] = 0.0;
        if (cdf[i] > 1.0)
            cdf[i] = 1.0;
    }
    cdf[0]   = 0.0;
    cdf[N-1] = 1.0;

//  free(gam_grids);
//  free(gam_cdfs);
    free(count);

    return 0;
}
