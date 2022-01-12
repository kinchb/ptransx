#!/bin/bash

rm -rf build
rm _mc_compy.*.so

module purge
module load intel/18.0    # NOTE: Change these to whatever modules on your target HPC platform load
module load intelmpi/2018 #       a C compiler, an MPI implementation, and a Python 3 interpreter.
module load python/3.6    #
pip install --user numpy
pip install --user scipy
pip install --user h5py
export MPICC=$(whereis mpicc | cut -d ' ' -f2) # user-level installation of mpi4py expects MPICC to point to location of mpicc
pip install --user mpi4py

CFLAGS="-O3 -std=c99" python3 setup.py build_ext --inplace
