#include <Python.h>
#include <numpy/arrayobject.h>
#include "mc_compy.h"

// docstrings
static char mc_compy_docstring[] =
    "";

// available functions
static PyObject *mc_compy_mc_compy(PyObject *self, PyObject *args);
static PyObject *mc_compy_sa_ratio(PyObject *self, PyObject *args);
static PyObject *mc_compy_sa_sigma(PyObject *self, PyObject *args);
static PyObject *mc_compy_single_scatter(PyObject *self, PyObject *args);

// module specification
static PyMethodDef module_methods[] = {
    {"mc_compy", (PyCFunction) mc_compy_mc_compy, METH_VARARGS, mc_compy_docstring},
    {"sa_ratio", (PyCFunction) mc_compy_sa_ratio, METH_VARARGS, mc_compy_docstring},
    {"sa_sigma", (PyCFunction) mc_compy_sa_sigma, METH_VARARGS, mc_compy_docstring},
    {"single_scatter", (PyCFunction) mc_compy_single_scatter, METH_VARARGS, mc_compy_docstring},
    {NULL}
};

static struct PyModuleDef _mc_compy =
{
    PyModuleDef_HEAD_INIT,
    "_mc_compy",
    "usage: TODO\n",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__mc_compy(void)
{
    return PyModule_Create(&_mc_compy);
}

// the main workhorse for table generation
// input:
// n, the "n"-scatter redistribution function to calculate (i.e., 1-scatter, 2-scatter, etc.)
// T, the electron temperature (in Kelvin)
// E, the pre-scatter photon energy (in eV)
// output:
// mc_ratio, the statistically-determined n-scatter mean amplification ratio
// sa_ratio, the semi-analytically-determined 1-scatter mean amplification ratio (always gives 1-scatter value, regardless of n)
// sa_sigma, the semi-analytically-determined thermally-averaged fluid rest frame electron scattering opacity (relative to Thomson)
// cdf,      the statistically-deteremined, angle-averaged CDF for a photon to n-scatter from E_0 to eta * E_0, for the specified eta grid
// constants which control execution:
// MC_N,      the number of photons to simulate n-scattering when building up redistribution function
// N,         the number of grid points in eta to use when calculating E_0 -> eta * E_0 probability
// ETA_BND_L, the smallest resolved amplification factor eta
// ETA_BND_U, the largest resolved amplification
// GAM_M_N,   number of gamma grid-points used to determine gam_max
// GAM_N_SA,  number of gamma grid-points used for integration over MJ for semi-analytic mean amplification calculation
// MU_I_N,    number of supplied Gauss-Legendre quadrature points used to perform the fluid frame integral over angle for SA mean amplification calc
static PyObject *mc_compy_mc_compy(PyObject *self, PyObject *args) {
    int n;
    double T, E;
    PyArrayObject *mu_i_obj, *w_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "iddOO", &n, &T, &E, &mu_i_obj, &w_obj))
        return NULL;

    double *mu_i, *w;
    mu_i = (double *)mu_i_obj->data;
    w    = (double *)w_obj->data;

    // create pointers/variables to/for the output data from mc_compy
    double *E_record = NULL;
    while (E_record == NULL)
        E_record = malloc(MC_N * sizeof(*E_record));
    double *mu_record = NULL;
    while (mu_record == NULL)
        mu_record = malloc(MC_N * sizeof(*mu_record));

    double mc_ratio = mc_scatters(n, T, E, E_record, mu_record);
    double sa_ratio = sa_calc_ratio(T, E, mu_i, w);
    double sa_sigma = sa_calc_sigma(T, E, mu_i, w);

    double count[N];
    for (int i = 0; i < N; i++) {
        count[i] = 0.0;
    }
    for (int i = 0; i < MC_N; i++) {
        long ndx = (long)(log((E_record[i]/E)/ETA_BND_L)/(log(ETA_BND_U/ETA_BND_L)/(N-1)));
        if (ndx >= 0 && ndx < N) {
            count[ndx] += 1.0;
        }
    }

    double *cdf = NULL;
    while (cdf == NULL)
        cdf = malloc(N * sizeof(*cdf));

    cdf[0] = 0.0;
    for (int i = 1; i < N; i++)
        cdf[i] = cdf[i-1] + count[i-1];
    for (int i = 0; i < N; i++) {
        cdf[i] /= MC_N;
        if (cdf[i] < 0.0)
            cdf[i] = 0.0;
        if (cdf[i] > 1.0)
            cdf[i] = 1.0;
    }
    cdf[0]   = 0.0;
    cdf[N-1] = 1.0;

    PyObject *result_3 = PyTuple_New(N);
    for (int i = 0; i < N; i++) {
        PyTuple_SetItem(result_3, i, PyFloat_FromDouble(cdf[i]));
    }

    PyObject *result = PyTuple_New(4);
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(mc_ratio));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(sa_ratio));
    PyTuple_SetItem(result, 2, PyFloat_FromDouble(sa_sigma));
    PyTuple_SetItem(result, 3, result_3);

    free(E_record);
    free(mu_record);
    free(cdf);

    return result;
}

// takes nearly the same input as the mc_compy (minus the "n" as this only does a 1-scatter calculation),
// but only returns the mean amplification ratio component
static PyObject *mc_compy_sa_ratio(PyObject *self, PyObject *args) {
    double T, E;
    PyArrayObject *mu_i_obj, *w_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "ddOO", &T, &E, &mu_i_obj, &w_obj))
        return NULL;

    double *mu_i, *w;

    mu_i = (double *)mu_i_obj->data;
    w    = (double *)w_obj->data;

    double sa_ratio = sa_calc_ratio(T, E, mu_i, w);

    // build output
    PyObject *result = PyFloat_FromDouble(sa_ratio);

    return result;
}

// takes nearly the same input as the mc_compy (minus the "n"),
// but only returns the total cross section component
static PyObject *mc_compy_sa_sigma(PyObject *self, PyObject *args) {
    double T, E;
    PyArrayObject *mu_i_obj, *w_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "ddOO", &T, &E, &mu_i_obj, &w_obj))
        return NULL;

    double *mu_i, *w;

    mu_i = (double *)mu_i_obj->data;
    w    = (double *)w_obj->data;

    double sa_sigma = sa_calc_sigma(T, E, mu_i, w);

    // build output
    PyObject *result = PyFloat_FromDouble(sa_sigma);

    return result;
}

// simulate a single properly-sampled scattering event
static PyObject *mc_compy_single_scatter(PyObject *self, PyObject *args) {
    double T, E;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "dd", &T, &E))
        return NULL;

    double E_f = 0.0;
    double mu  = 0.0;

    single_scatter(T, E, &E_f, &mu);

    // build output tuple
    PyObject *result = PyTuple_New(2);
    PyTuple_SetItem(result, 0, PyFloat_FromDouble(E_f));
    PyTuple_SetItem(result, 1, PyFloat_FromDouble(mu));

    return result;
}
