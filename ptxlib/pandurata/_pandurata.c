#include <Python.h>
#include <numpy/arrayobject.h>
#include "panhead.h"

// docstrings
static char pandurata_docstring[] =
    "";

// available functions
static PyObject *pandurata_pandurata(PyObject *self, PyObject *args);

// module specification
static PyMethodDef module_methods[] = {
    {"pandurata", (PyCFunction) pandurata_pandurata, METH_VARARGS, pandurata_docstring},
    {NULL}
};

static struct PyModuleDef _pandurata =
{
    PyModuleDef_HEAD_INIT,
    "_pandurata",
    "usage: TODO\n",
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit__pandurata(void)
{
    Py_Initialize();
    import_array();

    return PyModule_Create(&_pandurata);
}

static PyObject *pandurata_pandurata(PyObject *self, PyObject *args) {
    int rank;
    PyObject *rr_obj, *tt_obj, *pp_obj, *rho_ijk_obj, *ut_ijk_obj, *ur_ijk_obj, *uz_ijk_obj, *up_ijk_obj, *diskbody_ik_obj, *Tdisk_ik_obj, *emtop_ik_obj, *embot_ik_obj;
    PyObject *T_list_obj, *E_list_obj, *ratio_obj, *cdf_obj;
    int d_Nphi, d_Nr;
    PyObject *phi_list_obj, *r_list_obj, *slab_exists_obj, *e_coarse_obj, *r_frac_top_disk_obj, *r_frac_bot_disk_obj, *refl_profs_top_disk_obj, *refl_profs_bot_disk_obj;
    PyObject *T_ijk_obj;
    PyObject *xstar_fluxtot_top_obj, *xstar_fluxtot_bot_obj;
    double thresh;

    // parse input tuple
    if (!PyArg_ParseTuple(args, "iOOOOOOOOOOOOOOOOiiOOOOOOOOOOOd",
                          &rank,
                          &rr_obj, &tt_obj, &pp_obj, &rho_ijk_obj, &ut_ijk_obj, &ur_ijk_obj, &uz_ijk_obj, &up_ijk_obj, &diskbody_ik_obj, &Tdisk_ik_obj, &emtop_ik_obj, &embot_ik_obj,
                          &T_list_obj, &E_list_obj, &ratio_obj, &cdf_obj,
                          &d_Nphi, &d_Nr,
                          &phi_list_obj, &r_list_obj, &slab_exists_obj, &e_coarse_obj, &r_frac_top_disk_obj, &r_frac_bot_disk_obj, &refl_profs_top_disk_obj, &refl_profs_bot_disk_obj,
                          &T_ijk_obj,
                          &xstar_fluxtot_top_obj, &xstar_fluxtot_bot_obj,
                          &thresh))
        return NULL;

    // interpret the input objects as numpy arrays
    PyObject *rr_arr          = PyArray_FROM_OTF(rr_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *tt_arr          = PyArray_FROM_OTF(tt_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *pp_arr          = PyArray_FROM_OTF(pp_obj,          NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *rho_ijk_arr     = PyArray_FROM_OTF(rho_ijk_obj,     NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ut_ijk_arr      = PyArray_FROM_OTF(ut_ijk_obj,      NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ur_ijk_arr      = PyArray_FROM_OTF(ur_ijk_obj,      NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *uz_ijk_arr      = PyArray_FROM_OTF(uz_ijk_obj,      NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *up_ijk_arr      = PyArray_FROM_OTF(up_ijk_obj,      NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *diskbody_ik_arr = PyArray_FROM_OTF(diskbody_ik_obj, NPY_INT64,  NPY_IN_ARRAY);
    PyObject *Tdisk_ik_arr    = PyArray_FROM_OTF(Tdisk_ik_obj,    NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *emtop_ik_arr    = PyArray_FROM_OTF(emtop_ik_obj,    NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *embot_ik_arr    = PyArray_FROM_OTF(embot_ik_obj,    NPY_DOUBLE, NPY_IN_ARRAY);

    PyObject *T_list_arr = PyArray_FROM_OTF(T_list_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *E_list_arr = PyArray_FROM_OTF(E_list_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *ratio_arr  = PyArray_FROM_OTF(ratio_obj,  NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *cdf_arr    = PyArray_FROM_OTF(cdf_obj,    NPY_DOUBLE, NPY_IN_ARRAY);

    PyObject *phi_list_arr            = PyArray_FROM_OTF(phi_list_obj,            NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *r_list_arr              = PyArray_FROM_OTF(r_list_obj,              NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *slab_exists_arr         = PyArray_FROM_OTF(slab_exists_obj,         NPY_INT64,  NPY_IN_ARRAY);
    PyObject *e_coarse_arr            = PyArray_FROM_OTF(e_coarse_obj,            NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *r_frac_top_disk_arr     = PyArray_FROM_OTF(r_frac_top_disk_obj,     NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *r_frac_bot_disk_arr     = PyArray_FROM_OTF(r_frac_bot_disk_obj,     NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *refl_profs_top_disk_arr = PyArray_FROM_OTF(refl_profs_top_disk_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *refl_profs_bot_disk_arr = PyArray_FROM_OTF(refl_profs_bot_disk_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    PyObject *T_ijk_arr = PyArray_FROM_OTF(T_ijk_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    PyObject *xstar_fluxtot_top_arr = PyArray_FROM_OTF(xstar_fluxtot_top_obj, NPY_DOUBLE, NPY_IN_ARRAY);
    PyObject *xstar_fluxtot_bot_arr = PyArray_FROM_OTF(xstar_fluxtot_bot_obj, NPY_DOUBLE, NPY_IN_ARRAY);

    // get pointers to the data as C types
    double *rr          = (double *)PyArray_DATA(rr_arr);
    double *tt          = (double *)PyArray_DATA(tt_arr);
    double *pp          = (double *)PyArray_DATA(pp_arr);
    double *rho_ijk     = (double *)PyArray_DATA(rho_ijk_arr);
    double *ut_ijk      = (double *)PyArray_DATA(ut_ijk_arr);
    double *ur_ijk      = (double *)PyArray_DATA(ur_ijk_arr);
    double *uz_ijk      = (double *)PyArray_DATA(uz_ijk_arr);
    double *up_ijk      = (double *)PyArray_DATA(up_ijk_arr);
    long   *diskbody_ik = (long  *)PyArray_DATA(diskbody_ik_arr);
    double *Tdisk_ik    = (double *)PyArray_DATA(Tdisk_ik_arr);
    double *emtop_ik    = (double *)PyArray_DATA(emtop_ik_arr);
    double *embot_ik    = (double *)PyArray_DATA(embot_ik_arr);

    double *T_list = (double *)PyArray_DATA(T_list_arr);
    double *E_list = (double *)PyArray_DATA(E_list_arr);
    double *ratio  = (double *)PyArray_DATA(ratio_arr);
    double *cdf    = (double *)PyArray_DATA(cdf_arr);

    double *phi_list            = (double *)PyArray_DATA(phi_list_arr);
    double *r_list              = (double *)PyArray_DATA(r_list_arr);
    long   *slab_exists         = (long   *)PyArray_DATA(slab_exists_arr);
    double *e_coarse            = (double *)PyArray_DATA(e_coarse_arr);
    double *r_frac_top_disk     = (double *)PyArray_DATA(r_frac_top_disk_arr);
    double *r_frac_bot_disk     = (double *)PyArray_DATA(r_frac_bot_disk_arr);
    double *refl_profs_top_disk = (double *)PyArray_DATA(refl_profs_top_disk_arr);
    double *refl_profs_bot_disk = (double *)PyArray_DATA(refl_profs_bot_disk_arr);

    double *T_ijk = (double *)PyArray_DATA(T_ijk_arr);

    double *xstar_fluxtot_top = (double *)PyArray_DATA(xstar_fluxtot_top_arr);
    double *xstar_fluxtot_bot = (double *)PyArray_DATA(xstar_fluxtot_bot_arr);

    double *reftop_ik = emtop_ik;
    double *refbot_ik = embot_ik;

    // create output arrays for pandurata results
    double *Ispec      = NULL;
    double *Rspecp_top = NULL;
    double *Rspecp_bot = NULL;
    double *corpow_ijk = NULL;

    while (Ispec == NULL)
        Ispec      = malloc((Ne+1)*(Nth_obs+1) * sizeof(*Ispec));
    while (Rspecp_top == NULL)
        Rspecp_top = malloc((Ne+1)*(Nph+1)*(Nr+1) * sizeof(*Rspecp_top));
    while (Rspecp_bot == NULL)
        Rspecp_bot = malloc((Ne+1)*(Nph+1)*(Nr+1) * sizeof(*Rspecp_bot));
    while (corpow_ijk == NULL)
        corpow_ijk = malloc((Nr+1)*(Nth+1)*(Nph+1) * sizeof(*corpow_ijk));

    int success = pandurata(rank,
                            rr, tt, pp, rho_ijk, ut_ijk, ur_ijk, uz_ijk, up_ijk, diskbody_ik, Tdisk_ik, emtop_ik, embot_ik, reftop_ik, refbot_ik,
                            T_list, E_list, ratio, cdf,
                            d_Nphi, d_Nr,
                            phi_list, r_list, slab_exists, e_coarse, r_frac_top_disk, r_frac_bot_disk, refl_profs_top_disk, refl_profs_bot_disk,
                            T_ijk,
                            xstar_fluxtot_top, xstar_fluxtot_bot,
                            Ispec, Rspecp_top, Rspecp_bot, corpow_ijk,
                            thresh);

    // check for failure
    if (success != 0) {
        PyErr_SetString(PyExc_RuntimeError, "pandurata failed");
        return NULL;
    }

    // clean up
    Py_DECREF(rr_arr);
    Py_DECREF(tt_arr);
    Py_DECREF(pp_arr);
    Py_DECREF(rho_ijk_arr);
    Py_DECREF(ut_ijk_arr);
    Py_DECREF(ur_ijk_arr);
    Py_DECREF(uz_ijk_arr);
    Py_DECREF(up_ijk_arr);
    Py_DECREF(diskbody_ik_arr);
    Py_DECREF(Tdisk_ik_arr);
    Py_DECREF(emtop_ik_arr);
    Py_DECREF(embot_ik_arr);

    Py_DECREF(T_list_arr);
    Py_DECREF(E_list_arr);
    Py_DECREF(ratio_arr);
    Py_DECREF(cdf_arr);

    Py_DECREF(phi_list_arr);
    Py_DECREF(r_list_arr);
    Py_DECREF(slab_exists_arr);
    Py_DECREF(e_coarse_arr);
    Py_DECREF(r_frac_top_disk_arr);
    Py_DECREF(r_frac_bot_disk_arr);
    Py_DECREF(refl_profs_top_disk_arr);
    Py_DECREF(refl_profs_bot_disk_arr);

    Py_DECREF(T_ijk_arr);

    Py_DECREF(xstar_fluxtot_top_arr);
    Py_DECREF(xstar_fluxtot_bot_arr);

    // build output numpy arrays
    npy_intp dims1[2] = {Ne+1, Nth_obs+1};
    npy_intp dims2[3] = {Ne+1, Nph+1, Nr+1};
    npy_intp dims3[3] = {Nr+1, Nth+1, Nph+1};
    PyObject *Ispec_arr      = PyArray_SimpleNewFromData(2, dims1, NPY_DOUBLE, Ispec);
    PyObject *Rspecp_top_arr = PyArray_SimpleNewFromData(3, dims2, NPY_DOUBLE, Rspecp_top);
    PyObject *Rspecp_bot_arr = PyArray_SimpleNewFromData(3, dims2, NPY_DOUBLE, Rspecp_bot);
    PyObject *corpow_ijk_arr = PyArray_SimpleNewFromData(3, dims3, NPY_DOUBLE, corpow_ijk);

    // make sure these got made
    if (Ispec_arr == NULL || Rspecp_top_arr == NULL || Rspecp_bot_arr == NULL || corpow_ijk_arr == NULL) {
        PyErr_SetString(PyExc_RuntimeError, "pandurata failed: making output arrays");
        return NULL;
    }

    PyArray_ENABLEFLAGS((PyArrayObject *)Ispec_arr,      NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)Rspecp_top_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)Rspecp_bot_arr, NPY_ARRAY_OWNDATA);
    PyArray_ENABLEFLAGS((PyArrayObject *)corpow_ijk_arr, NPY_ARRAY_OWNDATA);

    // build output tuple
    PyObject *result = Py_BuildValue("OOOO", Ispec_arr, Rspecp_top_arr, Rspecp_bot_arr, corpow_ijk_arr);

    Py_DECREF(Ispec_arr);
    Py_DECREF(Rspecp_top_arr);
    Py_DECREF(Rspecp_bot_arr);
    Py_DECREF(corpow_ijk_arr);

    return result;
}
