#include "fsolver.h"

// constructs "A" from Mihalas (1985), but only the nonzero diagonal elements as a single vector
int make_A(int d, double *mu, double *dtaub, double *dtauc, double *A, int D, int N) {
    int m, n, k;

    if ((d > 0) && (d < D-2)) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            A[k] = (mu[m] * mu[m])/(dtauc[d*(N-1) + n] * dtaub[d*(N-1) + n]);
        }
        return 0;
    }

    if (d == D-2) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            A[k] = mu[m]/dtaub[(D-2)*(N-1) + n];
        }
        return 0;
    }

    // we really shouldn't get here
    return -1;
}

// constructs "C" from Mihalas (1985), returns as a dense matrix
int make_C(int d, double *mu, double *dtaub, double *dtauc, double *C, int D, int N) {
    int i, j, k, m, n;

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            C[i*K + j] = 0.0;
        }
    }

    if (d == 0) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            C[k*K + k] = mu[m]/dtaub[N-1 + n];
        }
        return 0;
    }

    if ((d > 0) && (d < D-2)) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            C[k*K + k] = (mu[m] * mu[m])/(dtauc[d*(N-1) + n] * dtaub[(d+1)*(N-1) + n]);
        }
        return 0;
    }

    // we really shouldn't get here
    return -1;
}

// constructs "L" from Mihalas (1985)
int make_L(int d, double *mu, double *dtaub, double *dtauc, double *absorp, double *scatt, double *emis, double *inctop, double *incbot, double *L, int D, int N) {
    int m, n, k;

    if (d == 0) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            L[k] = inctop[k]/(1.0 + dtauc[n]/(2.0 * mu[m])) + (dtauc[n]/mu[m])*(1.0/(absorp[n] + scatt[n]))*emis[n];
        }
        return 0;
    }

    if ((d > 0) && (d < D-2)) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            L[k] = (1.0/(absorp[d*(N-1) + n] + scatt[d*(N-1) + n]))*emis[d*(N-1) + n];
        }
        return 0;
    }

    if (d == D-2) {
        for (k = 0; k < K; k++) {
            m = k % M;
            n = (k - m)/M;
            L[k] = incbot[k]/(1.0 + dtauc[(D-2)*(N-1) + n]/(2.0 * mu[m])) + (dtauc[(D-2)*(N-1) + n]/mu[m])*(1.0/(absorp[(D-2)*(N-1) + n] + scatt[(D-2)*(N-1) + n]))*emis[(D-2)*(N-1) + n];
        }
        return 0;
    }

    // we really shouldn't get here
    return -1;
}

int make_B(int d, double *mu, double *dtaub, double *dtauc, double *absorp, double *scatt, double *dmu, double *e_grid, double *tps, double *B, int D, int N) {
    int n_col, n_row, last_n_row, k, k_col, k_row, index, col_index, next_col_index, num_entries, mc, mr, nc, nr, i, j;
    int start_index, end_index;
    double tmp;

    int *cols          = NULL;
    int *rows          = NULL;
    double *phis       = NULL;
    double *phi_matrix = NULL;

    double phase[M][M];

    for (mr = 0; mr < M; mr++) {
        for (mc = 0; mc < M; mc++) {
            phase[mr][mc] = 0.5*(1.0 + pow(mu[mc] * mu[mr] + sqrt(1.0 - mu[mr] * mu[mr]) * sqrt(1.0 - mu[mc] * mu[mc]), 2));
        }
    }
    for (mc = 0; mc < M; mc++) {
        tmp = 0.0;
        for (mr = 0; mr < M; mr++) {
            tmp += phase[mr][mc] * dmu[mr];
        }
        for (mr = 0; mr < M; mr++) {
            phase[mr][mc] = phase[mr][mc]/tmp;
        }
    }

    if (d == 0) {
        start_index = 0;
        end_index   = (int)(tps[start_index] + 0.1);
    }
    else {
        start_index = 0;
        end_index   = (int)(tps[start_index] + 0.1);
        i = 0;
        while (i < d) {
            start_index = end_index;
            end_index   = (int)(tps[start_index] + 0.1);
            i++;
        }
    }
    num_entries = (end_index - start_index - 1)/3;

    while (cols == NULL)
        cols = malloc(num_entries * sizeof(*cols));
    while (rows == NULL)
        rows = malloc(num_entries * sizeof(*rows));
    while (phis == NULL)
        phis = malloc(num_entries * sizeof(*phis));

    for (i = start_index+1, j = 0; j < num_entries; i += 3, j++) {
        cols[j] = (int)(tps[i] + 0.1);
        rows[j] = (int)(tps[i+1] + 0.1);
        phis[j] = tps[i+2];
    }

    while (phi_matrix == NULL)
        phi_matrix = malloc((N-1)*(N-1) * sizeof(*phi_matrix));

    for (i = 0; i < N-1; i++) {
        for (j = 0; j < N-1; j++) {
            phi_matrix[i*(N-1) + j] = 0.0;
        }
    }

    if (0) {
        for (i = 0; i < N-1; i++) {
            if ((i == 0) || (i == N-2)) {
                phi_matrix[i*(N-1) + i] = 1.0;
            } else {
                phi_matrix[(i-1)*(N-1) + i] = 1.0/3.0;
                phi_matrix[i*(N-1) + i] = 1.0/3.0;
                phi_matrix[(i+1)*(N-1) + i] = 1.0/3.0;
            }
        }
    } else {
        for (i = 0; i < num_entries; i++) {
            if ((cols[i] >= 0) && (cols[i] < N-1) && (rows[i] >= 0) && (rows[i] < N-1)) {
                phi_matrix[rows[i]*(N-1) + cols[i]] = phis[i];
            } else {
                fprintf(stdout, ">>> set phi_matrix out of bounds: %d %d %d, %e\n", rows[i], cols[i], rows[i]*(N-1) + cols[i], phis[i]); fflush(stdout);
            }
            if (rows[i]*(N-1) + cols[i] >= (N-1)*(N-1)) {
                fprintf(stdout, ">>> set phi_matrix out of bounds (2): %d %d %d, %e\n", rows[i], cols[i], rows[i]*(N-1) + cols[i], phis[i]); fflush(stdout);
            }
        }
    }

    if (d == 0) {
        for (i = 0; i < K; i++) {
            mc = (i % M);
            nc = (i - mc)/M;
            mr = mc;
            nr = nc;
            B[i*K + i] = (mu[mr]/dtaub[N-1 + nr] + 1.0/(1.0 + dtauc[nr]/(2.0 * mu[mr])) + dtauc[nr]/mu[mr]);
        }
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                mc = (j % M);
                nc = (j - mc)/M;
                mr = (i % M);
                nr = (i - mr)/M;
                B[i*K + j] = B[i*K + j] - (dtauc[nr]/mu[mr]) * (scatt[nc]/(absorp[nr] + scatt[nr])) * phase[mr][mc] * dmu[mc] * phi_matrix[nr*(N-1) + nc];
            }
        }
    }

    if ((d > 0) && (d < D-2)) {
        for (i = 0; i < K; i++) {
            mc = (i % M);
            nc = (i - mc)/M;
            mr = mc;
            nr = nc;
            B[i*K + i] = (1.0 + ((mu[mr] * mu[mr])/dtauc[d*(N-1) + nr])*((1.0/dtaub[d*(N-1) + nr]) + (1.0/dtaub[(d+1)*(N-1) + nr])));
        }
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                mc = (j % M);
                nc = (j - mc)/M;
                mr = (i % M);
                nr = (i - mr)/M;
                B[i*K + j] = B[i*K + j] - (scatt[d*(N-1) + nc]/(absorp[d*(N-1) + nr] + scatt[d*(N-1) + nr])) * phase[mr][mc] * dmu[mc] * phi_matrix[nr*(N-1) + nc];
            }
        }
    }

    if (d == D-2) {
        for (i = 0; i < K; i++) {
            mc = (i % M);
            nc = (i - mc)/M;
            mr = mc;
            nr = nc;
            B[i*K + i] = (mu[mr]/dtaub[(D-2)*(N-1) + nr] + 1.0/(1.0 + dtauc[(D-2)*(N-1) + nr]/(2.0 * mu[mr])) + dtauc[(D-2)*(N-1) + nr]/mu[mr]);
        }
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                mc = (j % M);
                nc = (j - mc)/M;
                mr = (i % M);
                nr = (i - mr)/M;
                B[i*K + j] = B[i*K + j] - (dtauc[(D-2)*(N-1) + nr]/mu[mr]) * (scatt[(D-2)*(N-1) + nc]/(absorp[(D-2)*(N-1) + nr] + scatt[(D-2)*(N-1) + nr])) * phase[mr][mc] * dmu[mc] * phi_matrix[nr*(N-1) + nc];
            }
        }
    }

    free(cols);
    free(rows);
    free(phis);
    free(phi_matrix);

    return 0;
}

// performs A - B = C, where A, B, and C are dense matrices
int subtract(double *A, double *B, double *C, int N) {
    int i, j;

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            C[i*K + j] = A[i*K + j] - B[i*K + j];
        }
    }

    return 0;
}

// performs A * vec = res, where A is a dense matrix and vec and res are column vectors
int mult_mat_by_vec(double *A, double *vec, double *res, int N) {
    int i, j;

    for (i = 0; i < K; i++) {
        res[i] = 0.0;
    }

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            res[i] += A[i*K + j] * vec[j];
        }
    }

    return 0;
}

// performs diag_vec * A = B, where A and B are dense matrices and diag_vec is a diagonal matrix specified only by its nonzero elements
int mult_mat_by_diag(double *diag_vec, double *A, double *B, int N) {
    int i, j;

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            B[i*K + j] = diag_vec[i] * A[i*K + j];
        }
    }

    return 0;
}

// performs diag_vec * vec = res, where vec and res are column vectors and diag_vec is a diagonal matrix specified only by its nonzero elements
int mult_vec_by_diag(double *diag_vec, double *vec, double *res, int N) {
    int k;

    for (k = 0; k < K; k++) {
        res[k] = diag_vec[k] * vec[k];
    }

    return 0;
}

// performs vec1 + vec2 = res, where all are column vectors
int add_vecs(double *vec1, double *vec2, double *res, int N) {
    int k;

    for (k = 0; k < K; k++) {
        res[k] = vec1[k] + vec2[k];
    }

    return 0;
}

int copy_full(double *from, double *to, int N) {
    int i, j;

    for (i = 0; i < K; i++) {
        for (j = 0; j < K; j++) {
            to[i*K + j] = from[i*K + j];
        }
    }

    return 0;
}

int clean(double *arr, int dim) {
    int i;

    for (i = 0; i < dim; i++) {
        if (fabs(arr[i]) < LBND)
            arr[i] = 0.0;
    }

    return 0;
}

int fsolver(double *absorp, double *scatt, double *emis, double *dtaub, double *dtauc, double *inctop, double *incbot, double *mu, double *dmu, double *e_grid, double *tps, int D, int N, double *jf) {
    int d, k, i, j;
    lapack_int info = 0;

    double **B_mats = NULL;
    double **G_mats = NULL;

    double *A    = NULL;
    double *L    = NULL;
    double *v    = NULL;
    double *lhs  = NULL;
    double *rhs  = NULL;
    double *tmp1 = NULL;
    double *tmp2 = NULL;
    double *tmp3 = NULL;
    double *tmp4 = NULL;
    double *tmp5 = NULL;
    double *tmp6 = NULL;

    lapack_int *piv = NULL;

    while (B_mats == NULL)
        B_mats = malloc((D-1) * sizeof(double *));
    while (G_mats == NULL)
        G_mats = malloc((D-2) * sizeof(double *));

    while (A == NULL)
        A = malloc(K * sizeof(*A));
    while (L == NULL)
        L = malloc(K * sizeof(*L));

    while (v == NULL)
        v  = malloc((D-1) * K * sizeof(*v));

    while (lhs == NULL)
        lhs = malloc(K * K * sizeof(*lhs));
    while (rhs == NULL)
        rhs = malloc(K * K * sizeof(*rhs));

    while (tmp1 == NULL)
        tmp1 = malloc(K * K * sizeof(*tmp1));
    while (tmp2 == NULL)
        tmp2 = malloc(K * K * sizeof(*tmp2));

    while (tmp3 == NULL)
        tmp3 = malloc(K * sizeof(*tmp3));
    while (tmp4 == NULL)
        tmp4 = malloc(K * sizeof(*tmp4));
    while (tmp5 == NULL)
        tmp5 = malloc(K * sizeof(*tmp5));
    while (tmp6 == NULL)
        tmp6 = malloc(K * sizeof(*tmp6));

    while (piv == NULL)
        piv = malloc(K * sizeof(*piv));

    for (d = 0; d < D-1; d++) {
        B_mats[d] = NULL;
        while (B_mats[d] == NULL)
            B_mats[d] = malloc(K * K * sizeof(double));
        for (i = 0; i < K; i++) {
            for (j = 0; j < K; j++) {
                B_mats[d][i*K + j] = 0.0;
            }
        }
    }

    for (d = 0; d < D-2; d++) {
        G_mats[d] = NULL;
        while (G_mats[d] == NULL)
            G_mats[d] = malloc(K * K * sizeof(double));
    }

    for (d = 0; d < D-1; d++) {
        make_B(d, mu, dtaub, dtauc, absorp, scatt, dmu, e_grid, tps, B_mats[d], D, N);
    }

    for (d = 0; d < D-2; d++) {
        if (d == 0) {
            copy_full(B_mats[0], lhs, N);
            make_C(d, mu, dtaub, dtauc, rhs, D, N);
            info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int)K, (lapack_int)K, lhs, (lapack_int)K, piv, rhs, (lapack_int)K);
            clean(rhs, K * K);
            copy_full(rhs, G_mats[0], N);
        }
        else {
            make_A(d, mu, dtaub, dtauc, A, D, N);
            mult_mat_by_diag(A, G_mats[d-1], tmp1, N);
            subtract(B_mats[d], tmp1, tmp2, N);
            copy_full(tmp2, lhs, N);
            make_C(d, mu, dtaub, dtauc, rhs, D, N);
            info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int)K, (lapack_int)K, lhs, (lapack_int)K, piv, rhs, (lapack_int)K);
            clean(rhs, K * K);
            copy_full(rhs, G_mats[d], N);
        }
    }

    free(rhs);
    rhs = NULL;
    while (rhs == NULL)
        rhs = malloc(K * sizeof(*rhs));

    for (d = 0; d < D-1; d++) {
        if (d == 0) {
            copy_full(B_mats[0], lhs, N);
            make_L(d, mu, dtaub, dtauc, absorp, scatt, emis, inctop, incbot, rhs, D, N);
        }
        else {
            make_A(d, mu, dtaub, dtauc, A, D, N);
            mult_mat_by_diag(A, G_mats[d-1], tmp1, N);
            subtract(B_mats[d], tmp1, tmp2, N);
            copy_full(tmp2, lhs, N);
            make_L(d, mu, dtaub, dtauc, absorp, scatt, emis, inctop, incbot, L, D, N);
            for (k = 0; k < K; k++) {
                tmp3[k] = v[(d-1)*K + k];
            }
            mult_vec_by_diag(A, tmp3, tmp4, N);
            add_vecs(L, tmp4, rhs, N);
        }
        info = LAPACKE_dgesv(LAPACK_ROW_MAJOR, (lapack_int)K, (lapack_int)1, lhs, (lapack_int)K, piv, rhs, (lapack_int)1);
        clean(rhs, K);
        for (k = 0; k < K; k++) {
            v[d*K + k] = rhs[k];
        }
    }

    for (k = 0; k < K; k++) {
        jf[(D-2)*K + k] = v[(D-2)*K + k];
    }

    for (d = D-3; d > -1; d--) {
        for (k = 0; k < K; k++) {
            tmp3[k] = jf[(d+1)*K + k];
        }
        mult_mat_by_vec(G_mats[d], tmp3, tmp4, N);
        for (k = 0; k < K; k++) {
            tmp5[k] = v[d*K + k];
        }
        add_vecs(tmp4, tmp5, tmp6, N);
        for (k = 0; k < K; k++) {
            jf[d*K + k] = tmp6[k];
        }
    }

    for (d = 0; d < D-1; d++) {
        free(B_mats[d]);
    }
    free(B_mats);

    for (d = 0; d < D-2; d++) {
        free(G_mats[d]);
    }
    free(G_mats);

    free(A);
    free(L);
    free(v);
    free(lhs);
    free(rhs);
    free(tmp1);
    free(tmp2);
    free(tmp3);
    free(tmp4);
    free(tmp5);
    free(tmp6);
    free(piv);

    return 0;
}

