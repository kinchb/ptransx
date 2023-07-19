# command line execution:
# mpirun -n <numper of processors> python regrid.py <input directory> <output directory> <'panptx' input file> <Fe abund> >rgr.out 2>rgr.err

# ----------+---------- start: setup and configuration ----------+----------

import os

machine = ''
if ('bw.txt' in os.listdir('.')):
    machine = 'bw'
if ('marcc.txt' in os.listdir('.')):
    machine = 'marcc'
if ('lanl.txt' in os.listdir('.')):
    machine = 'lanl'
if ('protagoras.txt' in os.listdir('.')):
    machine = 'protagoras'

import sys

cwd = os.getcwd()

if (machine == 'bw'):
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    os.environ['LHEA_DATA'] = '/mnt/c/scratch/sciteam/kinch/xstar_lib'

if (machine == 'marcc'):
    sys.path.insert(0, '/home-1/bkinch1@jhu.edu/.local/lib/python2.7/site-packages')
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    os.environ['LHEA_DATA'] = '/home-1/bkinch1@jhu.edu/scratch/bkinch1/xstar_lib'

if (machine == 'lanl'):
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    os.environ['LHEA_DATA'] = cwd + '/ptxlib/xstarcomm'

if (machine == 'protagoras'):
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    os.environ['LHEA_DATA'] = cwd + '/ptxlib/xstarcomm'

skip = 4

import time

import numpy as np

from scipy.interpolate import interp1d

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

from _fsolver import fsolver

from _xstarcomm import xstarcomm

import ptxlib

from h5py_wrapper import GetFile

input_filename = sys.argv[3]

Fe_abund = float(sys.argv[4])

# ----------+----------  end: setup and configuration  ----------+----------

def recalc_slab(flnm):
    if ('lemis' + flnm[5:] in os.listdir('.')):
        return

    print('working on ' + flnm)

    start_time = time.time()

    # ----------+---------- start: actual calculation ----------+----------

    # --- --- --- regrid for separate line/continuum transfer --- --- ---

    slab_type = flnm.split('_')[1]
    phi_index = int(flnm.split('_')[2])
    r_index   = int(flnm.split('_')[3].split('.')[0])

    with GetFile(sys.argv[1] + '/' + flnm, open_type = 'r') as f:
        r        = float(f.read('r'))
        z        = f.read('z')
        zc       = f.read('zc')
        energies = f.read('energies')
        jf       = f.read('jf')
        dens     = f.read('dens')
        heat     = f.read('heat')
        temp     = f.read('temp')
        Tbnd     = float(f.read('Tbnd'))
        mu       = f.read('mu')
        dmu      = f.read('dmu')

    D = len(z)
    M = len(mu)
    N = len(energies)

    new_N = 801

    energies, inctop, incbot, emis, c_emis, absorp, xee = ptxlib.regrid(energies, jf, dens, Fe_abund, heat, temp, Tbnd, mu, dmu, new_N, D, M, N, slab_type)

    N = new_N

    # remove K-alpha line photons from continuum emissivity
    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5*(energies[n] + energies[n+1])

    e_lower = 6.3e3
    e_upper = 7.0e3

    l_emis = np.zeros((D-1, N-1))

    for d in range(0, D-1):
        for n in range(0, N-1):
            if ((ec[n] >= e_lower) and (ec[n] <= e_upper)):
                l_emis[d][n] = emis[d][n] - c_emis[d][n]
            else:
                c_emis[d][n] = emis[d][n]

    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * ptxlib.sigma(0.5*(energies[n+1] + energies[n]))

    # condense grids

    new_N = 401

    new_energies = np.logspace(1, 5, new_N, endpoint = True)

    new_ec = np.zeros(new_N-1)
    for n in range(0, new_N-1):
        new_ec[n] = 0.5*(new_energies[n+1] + new_energies[n])

    new_emis   = np.zeros((D-1, new_N-1))
    new_l_emis = np.zeros((D-1, new_N-1))
    new_c_emis = np.zeros((D-1, new_N-1))
    new_absorp = np.zeros((D-1, new_N-1))

    for d in range(0, D-1):
        new_emis[d]   = interp1d(ec, emis[d])(new_ec)
        new_l_emis[d] = interp1d(ec, l_emis[d])(new_ec)
        new_c_emis[d] = interp1d(ec, c_emis[d])(new_ec)
        new_absorp[d] = interp1d(ec, absorp[d])(new_ec)

    energies = new_energies
    ec       = new_ec
    emis     = new_emis
    l_emis   = new_l_emis
    c_emis   = new_c_emis
    absorp   = new_absorp

    N = new_N
    K = M*(new_N-1)

    inctop = np.zeros(K)
    incbot = np.zeros(K)

    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * ptxlib.sigma(0.5*(energies[n+1] + energies[n]))
    dtaub, dtauc = ptxlib.set_dtau(z, absorp, scatt, D, N)

    # write out for use in response.py

    ptxlib.write_output('opacs_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5', 0, D, M, N, r, 0, 0, z, zc, mu, dmu, energies, 0, 0, dens, xee, heat, temp, Tbnd, absorp, 0, scatt, emis, 0, c_emis, 0, 0, 0, 0, 0)

#   tps = maketps(energies.copy(), temp.copy(), D, N)
    tps = ptxlib.maketps_h5py(energies.copy(), temp.copy(), D, N)

    s_jf = fsolver(absorp.copy(), scatt.copy(), c_emis.copy(), dtaub.copy(), dtauc.copy(), inctop.copy(), incbot.copy(), mu.copy(), dmu.copy(), energies.copy(), tps.copy(), D, N)

    s_hf = np.zeros((D-2, K))
    for d in range(0, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            s_hf[d][k] = mu[m] * (s_jf[d+1][k] - s_jf[d][k])/dtaub[d+1][n]

    ptxlib.write_output('cemis_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5', 0, D, M, N, r, 0, 0, z, zc, mu, dmu, energies, 0, 0, dens, xee, heat, temp, Tbnd, absorp, 0, scatt, emis, 0, c_emis, 0, 0, s_jf, s_hf, 0)

    s_jf = fsolver(absorp.copy(), scatt.copy(), l_emis.copy(), dtaub.copy(), dtauc.copy(), inctop.copy(), incbot.copy(), mu.copy(), dmu.copy(), energies.copy(), tps.copy(), D, N)

    s_hf = np.zeros((D-2, K))
    for d in range(0, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            s_hf[d][k] = mu[m] * (s_jf[d+1][k] - s_jf[d][k])/dtaub[d+1][n]

    ptxlib.write_output('lemis_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5', 0, D, M, N, r, 0, 0, z, zc, mu, dmu, energies, 0, 0, dens, xee, heat, temp, Tbnd, absorp, 0, scatt, emis, 0, c_emis, 0, 0, s_jf, s_hf, 0)

    end_time = time.time()

    print('done with ' + flnm + ', time: ' + '{0:.2f}'.format((end_time - start_time)/60) + ' min')

    return

    # ----------+----------  end: actual calculation  ----------+----------

# ----------+---------- start: program execution ----------+----------

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

# ----------+---------- root process ----------+----------
if (rank == 0):
    file_list = []
    for flnm in os.listdir(sys.argv[1]):
        if (flnm[:5] == 'final'):
            file_list.append(flnm)

    print('detected: ')
    for flnm in file_list:
        print(flnm)

    print(repr(len(file_list)) + ' total')

    todo = range(0, len(file_list))

    # the ID numbers of the processors which will perform transfer solutions
    xfer_proc_ids = [rank + 1]
    while (xfer_proc_ids[-1] + skip < size):
        xfer_proc_ids.append(xfer_proc_ids[-1] + skip)

    # send initial set of jobs
    job_index = 0
    for i in xfer_proc_ids:
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index], file_list[todo[job_index]]], dest = i, tag = 0)
            job_index += 1

    # cycle through processors, receiving results and reassigning work
    finished = 0
    while (finished != len(todo)):
        results = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
        finished += 1
        # send next job to this processor
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index], file_list[todo[job_index]]], dest = results[2], tag = 0)
            job_index += 1

    # move all 'final' and input files into a specified directory
    os.system('mkdir ' + sys.argv[2])
    os.system('mv cemis_*.h5 ' + sys.argv[2])
    os.system('mv lemis_*.h5 ' + sys.argv[2])
    os.system('mv opacs_*.h5 ' + sys.argv[2])

    # run ptx2pan
    ptxlib.ptx2pan(input_filename, sys.argv[2], 'cemis', 'spec_out_c_hr.h5')
    ptxlib.ptx2pan(input_filename, sys.argv[2], 'lemis', 'spec_out_l_hr.h5')

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send([-1, -1], dest = i, tag = 0)

# ----------+---------- worker process ----------+----------
else:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        print('hello from ' + repr(rank))

        # break out if sent done signal
        if (tmp[0] == -1):
            break

        recalc_slab(tmp[2])

        reply = (0, 0, rank)

        # reply to root
        COMM_WORLD.send(reply, dest = 0, tag = 0)
