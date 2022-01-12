import os

import sys
sys.path.insert(0, os.getcwd())

import numpy as np

import h5py

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

import _mc_compy

use_log_grid = True

if use_log_grid:
    N = 8000
else:
    N = 5000

if use_log_grid:
    eta_bound_l = 1.0e-4
    eta_bound_u = 1.0e4
else:
    eta_bound_l = 0.0
    eta_bound_u = 10.0

mu_i, w = np.polynomial.legendre.leggauss(128)

table_fname = 'compton_data.h5'

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

if (rank == 0):
    T_list = np.logspace(4, 12, 801, endpoint = True)
    E_list = np.logspace(0, 7,  701, endpoint = True)

    if (table_fname in os.listdir('.')):
        f = h5py.File(table_fname, 'a')
        restart = True
    else:
        f = h5py.File(table_fname, 'w')
        restart = False

    if not restart:
        f.create_dataset('T_list', (len(T_list),), dtype = 'f')
        f['/T_list'][...] = T_list
        f.create_dataset('E_list', (len(E_list),), dtype = 'f')
        f['/E_list'][...] = E_list
        f.create_dataset('eta_grid', (N,), dtype = 'f')
        if use_log_grid:
            f['/eta_grid'][...] = np.logspace(np.log10(eta_bound_l), np.log10(eta_bound_u), N, endpoint = True)
        else:
            f['/eta_grid'][...] = np.linspace(eta_bound_l, eta_bound_u, N, endpoint = True)
        f.create_dataset('mc_ratio', (len(T_list), len(E_list)), dtype = 'f')
        f.create_dataset('sa_ratio', (len(T_list), len(E_list)), dtype = 'f')
        f.create_dataset('sigma',    (len(T_list), len(E_list)), dtype = 'f')

    todo = []
    for i in range(0, len(T_list)):
        for j in range(0, len(E_list)):
            if restart:
                if (f['/sigma'][i, j] == 0.0):
                    todo.append([i, j, T_list[i], E_list[j]])
            else:
                todo.append([i, j, T_list[i], E_list[j]])
    print('len(todo) = ' + repr(len(todo)), flush = True)

    f.close()

    job_index = 0

    # send initial set of jobs
    for i in range(1, size):
        if (job_index == len(todo)):
            break
        COMM_WORLD.send([job_index, todo[job_index][2], todo[job_index][3]], dest = i, tag = 0)
        job_index += 1

    finished = 0

    # cycle through processors, receiving results and reassigning work
    while (finished != len(todo)):
        results = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
        finished += 1
        f = h5py.File(table_fname, 'a')
        i = todo[results[0]][0]
        j = todo[results[0]][1]
        k = results[5]
        print(repr(finished) + '/' + repr(len(todo)), flush = True)
        f.create_dataset('cdf_' + repr(i) + '_' + repr(j), (N,), dtype = 'f')
        f['/cdf_' + repr(i) + '_' + repr(j)][...] = results[1]
        f['/mc_ratio'][i,j, ...] = results[2]
        f['/sa_ratio'][i,j, ...] = results[3]
        f['/sigma'][i,j, ...]    = results[4]
        f.close()
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index][2], todo[job_index][3]], dest = k, tag = 0)
            job_index += 1

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send(-1, dest = i, tag = 0)
else:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = 0, tag = 0)

        # break out if sent done signal
        if (tmp == -1):
            break

        # get results
        mc_ratio, sa_ratio, sigma, cdf = _mc_compy.mc_compy(1, tmp[1], tmp[2], mu_i, w)

        cdf = np.array(cdf)

        # return results
        COMM_WORLD.send([tmp[0], cdf, mc_ratio, sa_ratio, sigma, rank], dest = 0, tag = 0)
