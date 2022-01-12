#include <Python.h>
#include <numpy/arrayobject.h>
#include "fsolver.h"

// docstrings
static char fsolver_docstring[] =
    "";

// available functions
static PyObject *fsolver_fsolver(PyObject *self, PyObject *args);

// module specification
static PyMethodDef module_methods[] = {
    {"fsolver", (PyCFunction) fsolver_fsolver, METH_VARARGS, fsolver_docstring},
    {NULL}
};

static struct PyModuleDef _fsolver =
{
    PyModuleDef_HEAD_INIT,
    "_fsolver",
    "usage: TODO\n",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__fsolver(void)
{
    Py_Initialize();
    import_array();

    return PyModule_Create(&_fsolver);
}

static PyObject *fsolver_fsolver(PyObject *self, PyObject *args) {
    PyObject *absorp_obj, *scatt_obj, *emis_obj, *dtaub_obj, *dtauc_obj, *inctop_obj, *incbot_obj, *mu_obj, *dmu_obj, *e_grid_obj, *tps_obj, *D_obj, *N_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "OOOOOOOOOOOOO", &absorp_obj, &scatt_obj, &emis_obj, &dtaub_obj, &dtauc_obj, &inctop_obj, &incbot_obj, &mu_obj, &dmu_obj, &e_grid_obj, &tps_obj, &D_obj, &N_obj))
        return NULL;

    // interpret the input objects as numpy arrays
    PyObject *absorp_arr = PyArray_FROM_OTF(absorp_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *scatt_arr  = PyArray_FROM_OTF(scatt_obj,  NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *emis_arr   = PyArray_FROM_OTF(emis_obj,   NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *dtaub_arr  = PyArray_FROM_OTF(dtaub_obj,  NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *dtauc_arr  = PyArray_FROM_OTF(dtauc_obj,  NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *inctop_arr = PyArray_FROM_OTF(inctop_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *incbot_arr = PyArray_FROM_OTF(incbot_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mu_arr     = PyArray_FROM_OTF(mu_obj,     NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *dmu_arr    = PyArray_FROM_OTF(dmu_obj,    NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *e_grid_arr = PyArray_FROM_OTF(e_grid_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *tps_arr    = PyArray_FROM_OTF(tps_obj,    NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *D_arr      = PyArray_FROM_OTF(D_obj,      NPY_INT32,  NPY_IN_ARRAY);
    PyObject *N_arr      = PyArray_FROM_OTF(N_obj,      NPY_INT32,  NPY_IN_ARRAY);

    // if that didn't work, throw an exception
    if (absorp_arr == NULL || scatt_arr == NULL || emis_arr == NULL || dtaub_arr == NULL || dtauc_arr == NULL || inctop_arr == NULL || incbot_arr == NULL || mu_arr == NULL || dmu_arr == NULL || e_grid_arr == NULL || tps_arr == NULL || D_arr == NULL || N_arr == NULL) {
        Py_XDECREF(absorp_arr);
        Py_XDECREF(scatt_arr);
        Py_XDECREF(emis_arr);
        Py_XDECREF(dtaub_arr);
        Py_XDECREF(dtauc_arr);
        Py_XDECREF(inctop_arr);
        Py_XDECREF(incbot_arr);
        Py_XDECREF(mu_arr);
        Py_XDECREF(dmu_arr);
        Py_XDECREF(e_grid_arr);
        Py_XDECREF(tps_arr);
        Py_XDECREF(D_arr);
        Py_XDECREF(N_arr);
        return NULL;
    }

    // get pointers to the data as C types
    double *absorp = (double*)PyArray_DATA(absorp_arr);
    double *scatt  = (double*)PyArray_DATA(scatt_arr);
    double *emis   = (double*)PyArray_DATA(emis_arr);
    double *dtaub  = (double*)PyArray_DATA(dtaub_arr);
    double *dtauc  = (double*)PyArray_DATA(dtauc_arr);
    double *inctop = (double*)PyArray_DATA(inctop_arr);
    double *incbot = (double*)PyArray_DATA(incbot_arr);
    double *mu     = (double*)PyArray_DATA(mu_arr);
    double *dmu    = (double*)PyArray_DATA(dmu_arr);
    double *e_grid = (double*)PyArray_DATA(e_grid_arr);
    double *tps    = (double*)PyArray_DATA(tps_arr);
    int    D       = *((int*)PyArray_DATA(D_arr));
    int    N       = *((int*)PyArray_DATA(N_arr));

    // create a pointer to the output data from fsolver
    double *jf = NULL;

    while (jf == NULL)
        jf = malloc((D-1) * K * sizeof(*jf));

    int success = fsolver(absorp, scatt, emis, dtaub, dtauc, inctop, incbot, mu, dmu, e_grid, tps, D, N, jf);

    // clean up
    Py_DECREF(absorp);
    Py_DECREF(scatt);
    Py_DECREF(emis);
    Py_DECREF(dtaub);
    Py_DECREF(dtauc);
    Py_DECREF(inctop);
    Py_DECREF(incbot);
    Py_DECREF(mu);
    Py_DECREF(dmu);
    Py_DECREF(e_grid);
    Py_DECREF(tps);

    Py_DECREF(absorp_arr);
    Py_DECREF(scatt_arr);
    Py_DECREF(emis_arr);
    Py_DECREF(dtaub_arr);
    Py_DECREF(dtauc_arr);
    Py_DECREF(inctop_arr);
    Py_DECREF(incbot_arr);
    Py_DECREF(mu_arr);
    Py_DECREF(dmu_arr);
    Py_DECREF(e_grid_arr);
    Py_DECREF(tps_arr);
    Py_DECREF(D_arr);
    Py_DECREF(N_arr);

    // check for nonsense results
    if (success < 0 || jf == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "fsolver returned nonsense");
        return NULL;
    }

    // build output numpy array
    npy_intp dims[2] = {D-1, K};
    PyObject *jf_arr = PyArray_SimpleNewFromData(2, dims, NPY_DOUBLE, jf);
    PyArray_ENABLEFLAGS((PyArrayObject *)jf_arr, NPY_ARRAY_OWNDATA);

    return jf_arr;
}
