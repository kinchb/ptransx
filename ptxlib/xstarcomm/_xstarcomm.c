#include <Python.h>
#include <numpy/arrayobject.h>
#include "xstarcomm.h"

// docstrings
static char xstarcomm_docstring[] =
    "";

// available functions
static PyObject *xstarcomm_xstarcomm(PyObject *self, PyObject *args);
static PyObject *xstarcomm_gaunt(PyObject *self, PyObject *args);

// module specification
static PyMethodDef module_methods[] = {
    {"xstarcomm", (PyCFunction) xstarcomm_xstarcomm, METH_VARARGS, xstarcomm_docstring},
    {"fbg", (PyCFunction) xstarcomm_gaunt, METH_VARARGS, xstarcomm_docstring},
    {NULL}
};

static struct PyModuleDef _xstarcomm =
{
    PyModuleDef_HEAD_INIT,
    "_xstarcomm",
    "usage: TODO\n",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__xstarcomm(void)
{
    Py_Initialize();
    import_array();

    return PyModule_Create(&_xstarcomm);
}

static PyObject *xstarcomm_xstarcomm(PyObject *self, PyObject *args) {
    int niter, N;
    double density, heat, temp, Fe_abund;

    PyObject *energies_obj, *mi_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "OOddddii", &energies_obj, &mi_obj, &density, &heat, &temp, &Fe_abund, &niter, &N))
        return NULL;

    // interpret the input objects as numpy arrays
    PyObject *energies_arr = PyArray_FROM_OTF(energies_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *mi_arr       = PyArray_FROM_OTF(mi_obj,       NPY_DOUBLE, NPY_IN_ARRAY);

    // if that didn't work, throw an exception
    if (energies_arr == NULL || mi_arr == NULL) {
        Py_XDECREF(energies_arr);
        Py_XDECREF(mi_arr);
        return NULL;
    }

    // get pointers to the data as C types
    double *energies = (double *)PyArray_DATA(energies_arr);
    double *mi       = (double *)PyArray_DATA(mi_arr);

    // create pointers/variables to/for the output data from xstarcomm
    double *emis      = NULL;
    double *c_emis    = NULL;
    double *absorp    = NULL;
    double *heat_vals = NULL;
    double xee        = 0.0;
    double new_temp   = 0.0;

    while (emis == NULL)
        emis      = malloc((N-1) * sizeof(*emis));
    while (c_emis == NULL)
        c_emis    = malloc((N-1) * sizeof(*c_emis));
    while (absorp == NULL)
        absorp    = malloc((N-1) * sizeof(*absorp));
    while (heat_vals == NULL)
        heat_vals = malloc(9 * sizeof(*heat_vals));

    int success = xstarcomm(energies, mi, density, Fe_abund, heat, temp, niter, N, emis, c_emis, absorp, heat_vals, &xee, &new_temp);

    // check for nonsense results
    if (success != 0 || emis == NULL || c_emis == NULL || absorp == NULL || heat_vals == NULL || xee == 0.0 || new_temp == 0.0) {
        Py_DECREF(energies_arr);
        Py_DECREF(mi_arr);
        free(emis);
        free(c_emis);
        free(absorp);
        free(heat_vals);
        PyErr_SetString(PyExc_RuntimeError, "xstarcomm returned nonsense");
        return NULL;
    }

    // clean up
    Py_DECREF(energies_arr);
    Py_DECREF(mi_arr);

    // build output numpy arrays
    npy_intp dims[1]  = {N-1};
    npy_intp dims2[1] = {9};
    PyObject *emis_arr      = PyArray_SimpleNewFromData(1, dims,  NPY_DOUBLE, emis);
    PyObject *c_emis_arr    = PyArray_SimpleNewFromData(1, dims,  NPY_DOUBLE, c_emis);
    PyObject *absorp_arr    = PyArray_SimpleNewFromData(1, dims,  NPY_DOUBLE, absorp);
    PyObject *heat_vals_arr = PyArray_SimpleNewFromData(1, dims2, NPY_DOUBLE, heat_vals);

    // make sure these things got made
    if (emis_arr == NULL || c_emis_arr == NULL || absorp_arr == NULL || heat_vals_arr == NULL) {
        Py_XDECREF(emis_arr);
        Py_XDECREF(c_emis_arr);
        Py_XDECREF(absorp_arr);
        Py_XDECREF(heat_vals_arr);
        return NULL;
    }

    PyArray_ENABLEFLAGS((PyArrayObject *)emis_arr,      NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)c_emis_arr,    NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)absorp_arr,    NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)heat_vals_arr, NPY_ARRAY_OWNDATA);

    // build output tuple
    PyObject *result = Py_BuildValue("OOOOdd", emis_arr, c_emis_arr, absorp_arr, heat_vals_arr, xee, new_temp);

    Py_DECREF(emis_arr);
    Py_DECREF(c_emis_arr);
    Py_DECREF(absorp_arr);
    Py_DECREF(heat_vals_arr);

    return result;
}

static PyObject *xstarcomm_gaunt(PyObject *self, PyObject *args) {
    double u, gam, gaunt;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "dd", &u, &gam))
        return NULL;

    fbg2_(&u, &gam, &gaunt);

    PyObject *result = Py_BuildValue("d", gaunt);

    return result;
}

