#include <Python.h>
#include <numpy/arrayobject.h>
#include "response.h"

// docstrings
static char response_docstring[] =
    "";

// available functions
static PyObject *response_response(PyObject *self, PyObject *args);

// module specification
static PyMethodDef module_methods[] = {
    {"response", (PyCFunction) response_response, METH_VARARGS, response_docstring},
    {NULL}
};

static struct PyModuleDef _response =
{
    PyModuleDef_HEAD_INIT,
    "_response",
    "usage: TODO\n",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__response(void)
{
    Py_Initialize();
    import_array();

    return PyModule_Create(&_response);
}

static PyObject *response_response(PyObject *self, PyObject *args) {
    PyObject *top_or_bot_obj, *E_0_obj, *z_obj, *e_obj, *T_obj, *tau_grids_obj, *D_obj, *N_obj;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "OOOOOOOO", &top_or_bot_obj, &E_0_obj, &z_obj, &e_obj, &T_obj, &tau_grids_obj, &D_obj, &N_obj))
        return NULL;

    // interpret the input objects as numpy arrays
    PyObject *top_or_bot_arr = PyArray_FROM_OTF(top_or_bot_obj, NPY_INT32,  NPY_IN_ARRAY);
    PyObject *E_0_arr        = PyArray_FROM_OTF(E_0_obj,        NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *z_arr          = PyArray_FROM_OTF(z_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *e_arr          = PyArray_FROM_OTF(e_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *T_arr          = PyArray_FROM_OTF(T_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *tau_grids_arr  = PyArray_FROM_OTF(tau_grids_obj,  NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *D_arr          = PyArray_FROM_OTF(D_obj,          NPY_INT32,  NPY_IN_ARRAY);
    PyObject *N_arr          = PyArray_FROM_OTF(N_obj,          NPY_INT32,  NPY_IN_ARRAY);

    // if that didn't work, throw an exception
    if (top_or_bot_arr == NULL || E_0_arr == NULL || z_arr == NULL || e_arr == NULL || T_arr == NULL || tau_grids_arr == NULL || D_arr == NULL || N_arr == NULL) {
        Py_XDECREF(top_or_bot_arr);
        Py_XDECREF(E_0_arr);
        Py_XDECREF(z_arr);
        Py_XDECREF(e_arr);
        Py_XDECREF(T_arr);
        Py_XDECREF(tau_grids_arr);
        Py_XDECREF(D_arr);
        Py_XDECREF(N_arr);
        return NULL;
    }

    // get values of/pointers to the data as C types
    int    top_or_bot = *((int*)PyArray_DATA(top_or_bot_arr));
    double E_0        = *((double*)PyArray_DATA(E_0_arr));
    double *z         = (double*)PyArray_DATA(z_arr);
    double *e         = (double*)PyArray_DATA(e_arr);
    double *T         = (double*)PyArray_DATA(T_arr);
    double *tau_grids = (double*)PyArray_DATA(tau_grids_arr);
    int    D          = *((int*)PyArray_DATA(D_arr));
    int    N          = *((int*)PyArray_DATA(N_arr));

    // create pointers/variables to/for the output data from response
    double *cdf = NULL;
    while (cdf == NULL)
        cdf = malloc(N * sizeof(*cdf));

    int success = response(top_or_bot, E_0, z, e, T, tau_grids, D, N, cdf);

    // build output numpy array
    npy_intp dims[1]  = {N};
    PyObject *cdf_arr = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, cdf);

    // clean up
    Py_XDECREF(top_or_bot_arr);
    Py_XDECREF(E_0_arr);
    Py_XDECREF(z_arr);
    Py_XDECREF(e_arr);
    Py_XDECREF(T_arr);
    Py_XDECREF(tau_grids_arr);
    Py_XDECREF(D_arr);
    Py_XDECREF(N_arr);

    // build output and return
    return cdf_arr;
}
