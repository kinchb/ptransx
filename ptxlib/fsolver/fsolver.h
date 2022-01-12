#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <lapacke.h>

#define M 8
#define K (M*(N-1))
#define LBND 1.0e-20

int make_A(int d, double *mu, double *dtaub, double *dtauc, double *A, int D, int N);
int make_C(int d, double *mu, double *dtaub, double *dtauc, double *C, int D, int N);
int make_L(int d, double *mu, double *dtaub, double *dtauc, double *absorp, double *scatt, double *emis, double *inctop, double *incbot, double *L, int D, int N);
int make_B(int d, double *mu, double *dtaub, double *dtauc, double *absorp, double *scatt, double *dmu, double *e_grid, double *tps, double *B, int D, int N);
int subtract(double *A, double *B, double *C, int N);
int mult_mat_by_vec(double *A, double *vec, double *res, int N);
int mult_mat_by_diag(double *diag_vec, double *A, double *B, int N);
int mult_vec_by_diag(double *diag_vec, double *vec, double *res, int N);
int add_vecs(double *vec1, double *vec2, double *res, int N);
int copy_full(double *from, double *to, int N);
int clean(double *arr, int dim);
int fsolver(double *absorp, double *scatt, double *emis, double *dtaub, double *dtauc, double *inctop, double *incbot, double *mu, double *dmu, double *e_grid, double *tps, int D, int N, double *jf);
