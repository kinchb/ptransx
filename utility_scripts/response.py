# command line execution:
# mpirun -n <numper of processors> python response.py <input directory> >rsp.out 2>rsp.err

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
    sys.path.insert(0, cwd + '/ptxlib/response')
    os.environ['LHEA_DATA'] = '/mnt/c/scratch/sciteam/kinch/xstar_lib'

if (machine == 'marcc'):
    sys.path.insert(0, '/home-1/bkinch1@jhu.edu/.local/lib/python2.7/site-packages')
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    sys.path.insert(0, cwd + '/ptxlib/response')
    os.environ['LHEA_DATA'] = '/home-1/bkinch1@jhu.edu/scratch/bkinch1/xstar_lib'

if (machine == 'lanl'):
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    sys.path.insert(0, cwd + '/ptxlib/response')
    os.environ['LHEA_DATA'] = cwd + '/ptxlib/xstarcomm'

if (machine == 'protagoras'):
    sys.path.insert(0, cwd + '/ptxlib')
    sys.path.insert(0, cwd + '/ptxlib/maketps')
    sys.path.insert(0, cwd + '/ptxlib/fsolver')
    sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
    sys.path.insert(0, cwd + '/ptxlib/response')
    os.environ['LHEA_DATA'] = cwd + '/ptxlib/xstarcomm'

import time

import numpy as np

from scipy.interpolate import interp1d

import h5py

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

from _response import response

from h5py_wrapper import GetFile

import ptxlib

num_ec = 33

eta_lb  = 1.0e-6
eta_ub  = 1.0e4
num_eta = 1000

# ----------+----------  end: setup and configuration  ----------+----------

def calc_profiles(flnm, E_0, eta, top_or_bot):
    slab_type = flnm.split('/')[-1].split('_')[1]

    with GetFile(sys.argv[2] + '/' + flnm, open_type = 'r') as f:
        z      = f.read('z')
        e      = f.read('energies')
        scatt  = f.read('scatt')
        absorp = f.read('absorp')
        temp   = f.read('temp')

    D = len(z)
    N = len(e)

    z = z - z[0]

    tmp_e = eta * E_0

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5 * (e[n] + e[n+1])

    tmp_scatt  = np.zeros((D-1, num_eta))
    tmp_absorp = np.zeros((D-1, num_eta))
    for d in range(0, D-1):
        tmp_scatt[d]  = interp1d(ec, scatt[d],  bounds_error = False, fill_value = (scatt[d][0],  scatt[d][-1]))(tmp_e)
        tmp_absorp[d] = interp1d(ec, absorp[d], bounds_error = False, fill_value = (absorp[d][0], absorp[d][-1]))(tmp_e)

    tau_grids = np.zeros((num_eta, D))
    for n in range(0, num_eta):
        for d in range(1, D):
            tau_grids[n][d] = tau_grids[n][d-1] + (tmp_scatt[d-1][n] + tmp_absorp[d-1][n])*(z[d] - z[d-1])
#           tau_grids[n][d] = tau_grids[n][d-1] + tmp_scatt[d-1][n]*(z[d] - z[d-1])

    if (top_or_bot == 'top'):
        if (slab_type == 'upper'):
#           print('rank = ' + repr(rank) + '; before response')
            cdf = response(2, E_0, z, tmp_e, temp, tau_grids, D, num_eta)
#           print('rank = ' + repr(rank) + '; after response')
        else:
#           print('rank = ' + repr(rank) + '; before response')
            cdf = response(1, E_0, z, tmp_e, temp, tau_grids, D, num_eta)
#           print('rank = ' + repr(rank) + '; after response')
    else:
        if (slab_type == 'lower'):
#           print('rank = ' + repr(rank) + '; before response')
            cdf = response(-2, E_0, z, tmp_e, temp, tau_grids, D, num_eta)
#           print('rank = ' + repr(rank) + '; after response')
        else:
#           print('rank = ' + repr(rank) + '; before response')
            cdf = response(-1, E_0, z, tmp_e, temp, tau_grids, D, num_eta)
#           print('rank = ' + repr(rank) + '; after response')

    return cdf

def calc_r_frac(flnm, top_or_bot):
    slab_type = flnm.split('/')[-1].split('_')[1]

    with GetFile(sys.argv[2] + '/' + flnm, open_type = 'r') as f:
        z        = f.read('z')
        mu       = f.read('mu')
        dmu      = f.read('dmu')
        energies = f.read('energies')
        absorp   = f.read('absorp')
        scatt    = f.read('scatt')

    D = len(z)
    M = len(mu)
    N = len(energies)
    K = M*(N-1)

    if (top_or_bot == 'top'):
        inctop = np.ones(K)
        incbot = np.zeros(K)
    else:
        inctop = np.zeros(K)
        incbot = np.ones(K)

    dtaub, dtauc = ptxlib.set_dtau(z, absorp, scatt, D, N)
    jf, hf = ptxlib.cfsolver(absorp, scatt, np.zeros((D-1, N-1)), dtaub, dtauc, inctop, incbot, mu, dmu, D, M, N)

    phs_out = np.zeros(N-1)
    for n in range(0, N-1):
        for m in range(0, M):
            phs_out[n] += (2*jf[0][n*M + m]  - inctop[n*M + m]) * dmu[m]
            phs_out[n] += (2*jf[-1][n*M + m] - incbot[n*M + m]) * dmu[m]

    phs_in = np.zeros(N-1)
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        if (top_or_bot == 'top'):
            phs_in[n] += inctop[k] * dmu[m]
        else:
            phs_in[n] += incbot[k] * dmu[m]

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5*(energies[n] + energies[n+1])

    f = interp1d(ec, phs_out/phs_in, fill_value = 'extrapolate')

    r_frac = f(energies)

    for i in range(0, len(r_frac)):
        if (r_frac[i] < 0.0):
            r_frac[i] = 0.0
        if (r_frac[i] > 1.0):
            r_frac[i] = 1.0

    return r_frac

def mc_for_slab(flnm):
    print('working on ' + flnm)

    start_time = time.time()

    slab_type = flnm.split('/')[-1].split('_')[1]

    with GetFile(sys.argv[2] + '/' + flnm, open_type = 'r') as f:
        r        = float(f.read('r'))
        energies = f.read('energies')

    num_ef = len(energies)
    e_fine = energies

    e_coarse = np.logspace(np.log10(energies[0]), np.log10(energies[-1]), num_ec, endpoint = True)

    eta = np.logspace(np.log10(eta_lb), np.log10(eta_ub), num_eta, endpoint = True)

    r_frac_top     = np.zeros(num_ef)
    r_frac_bot     = np.zeros(num_ef)
    refl_profs_top = np.zeros((num_ec, num_eta))
    refl_profs_bot = np.zeros((num_ec, num_eta))

    # for a "whole" slab, we must do two complete calculations, one for each side
    if (slab_type == 'whole'):
        r_frac_top = calc_r_frac(flnm, 'top')
        r_frac_bot = calc_r_frac(flnm, 'bot')
        # calculate the energy-dependent *shape* of the reflection+transmission profile on a relatively *coarse* grid (but with many photons)
        for i in range(0, num_ec):
            refl_profs_top[i] = calc_profiles(flnm, e_coarse[i], eta, 'top')
        for i in range(0, num_ec):
            refl_profs_bot[i] = calc_profiles(flnm, e_coarse[i], eta, 'bot')
    # otherwise we need only consider photons incident from the relevant side
    if (slab_type == 'upper'):
        r_frac_top = calc_r_frac(flnm, 'top')
        for i in range(0, num_ec):
            refl_profs_top[i] = calc_profiles(flnm, e_coarse[i], eta, 'top')
    if (slab_type == 'lower'):
        r_frac_bot = calc_r_frac(flnm, 'bot')
        for i in range(0, num_ec):
            refl_profs_bot[i] = calc_profiles(flnm, e_coarse[i], eta, 'bot')

    with GetFile(sys.argv[2] + '/mc' + flnm[5:], open_type = 'w') as f:
        f.write('dims',     np.array([num_ec, num_ef, num_eta], dtype = 'i'))
        f.write('r',        r)
        f.write('e_coarse', e_coarse)
        f.write('e_fine',   e_fine)
        f.write('eta',      eta)
        if (slab_type == 'whole'):
            f.write('r_frac_top',     r_frac_top)
            f.write('r_frac_bot',     r_frac_bot)
            f.write('refl_profs_top', refl_profs_top)
            f.write('refl_profs_bot', refl_profs_bot)
        if (slab_type == 'upper'):
            f.write('r_frac_top',     r_frac_top)
            f.write('refl_profs_top', refl_profs_top)
        if (slab_type == 'lower'):
            f.write('r_frac_bot',     r_frac_bot)
            f.write('refl_profs_bot', refl_profs_bot)

    end_time = time.time()

    print('finished with ' + flnm + ', time: ' + '{0:.2f}'.format((end_time - start_time)/60) + ' min')

# ----------+---------- start: program execution ----------+----------

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

# ----------+---------- root process ----------+----------
if (rank == 0):
    """
    if (sys.argv[2] != 'run_hr'):
        grand_iter = int(sys.argv[2].split('_')[-1])
        os.system('mkdir run_' + repr(grand_iter))
        os.system('mv final_*.h5 ./run_' + repr(grand_iter))
        os.system('mv scat_diskt.h5 ./run_' + repr(grand_iter))
        os.system('mv scat_diskb.h5 ./run_' + repr(grand_iter))
    """
    if (sys.argv[2] != 'run_hr'):
        # run ptx2pan
        grand_iter = int(sys.argv[2].split('_')[-1])
        ptxlib.ptx2pan(sys.argv[1], 'run_' + repr(grand_iter), 'final', 'spec_out_' + repr(grand_iter) + '.h5')

    file_list = []
    mc_list   = []
    for flnm in os.listdir(sys.argv[2]):
        if (((flnm[:5] == 'final') or (flnm[:5] == 'opacs')) and ('mc' + flnm[5:] not in os.listdir(sys.argv[2]))):
            file_list.append(flnm)
        if (((flnm[:5] == 'final') or (flnm[:5] == 'opacs')) and ('mc' + flnm[5:] in os.listdir(sys.argv[2]))):
            mc_list.append('mc' + flnm[5:])

    print('detected: ')
    for flnm in file_list:
        print(flnm)
    print(repr(len(file_list)) + ' total')
    print('mc files found = ' + repr(len(mc_list)))

    # is this the first time response is being run this round?
    if (len(mc_list) == 0):
        todo = range(0, len(file_list))

        # the ID numbers of the processors which will perform Monte Carlo solutions
        mc_proc_ids = [COMM_WORLD.Get_rank()+1]
        while (mc_proc_ids[-1] + 1 < COMM_WORLD.Get_size()):
            mc_proc_ids.append(mc_proc_ids[-1] + 1)

        # send initial set of jobs
        job_index = 0
        for i in mc_proc_ids:
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

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send([-1, -1], dest = i, tag = 0)

    # run rsp2pan on output
#   rsp2pan(sys.argv[2])

# ----------+---------- worker process ----------+----------
else:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        # break out if sent done signal
        if (tmp[0] == -1):
            break

        mc_for_slab(tmp[2])

        reply = (0, 0, rank)

        # reply to root
        COMM_WORLD.send(reply, dest = 0, tag = 0)
