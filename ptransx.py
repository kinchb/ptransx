# command line execution:
# mpirun -np <numper of processors> python ptransx.py <log file> <'panptx' input file> <Fe abund> <grand iteration> >ptx.out 2>ptx.err

# ----------+---------- start: setup and configuration ----------+----------

import os
import sys

cwd = os.getcwd()

sys.path.insert(0, cwd + '/ptxlib')
# sys.path.insert(0, cwd + '/ptxlib/maketps')
sys.path.insert(0, cwd + '/ptxlib/fsolver')
sys.path.insert(0, cwd + '/ptxlib/xstarcomm')
os.environ['LHEA_DATA'] = cwd + '/ptxlib/xstarlib'

# skip and jac_reserved_procs need to be configured on a per-machine basis depending on the
# number of cores and memory per node
skip = 2
jac_reserved_procs = 22

import time

import random

import numpy as np

from scipy.interpolate import interp1d

import h5py

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

# from _maketps import maketps

from _fsolver import fsolver

from _xstarcomm import xstarcomm
from _xstarcomm import fbg

from h5py_wrapper import GetFile

import ptxlib

input_filename = sys.argv[2]

# simulation parameters
with GetFile('./' + input_filename, open_type = 'r') as f:
    mass         = float(f.read('M'))     # central black hole mass (solar masses)
    a            = float(f.read('a'))     # black hole spin
#   slab_type_ik = f.read('slab_type_ik')

Fe_abund = float(sys.argv[3])           # the Fe abundance relative to solar

grand_iter = int(sys.argv[4])

# disk solution parameters:
T_init   = 1.0e7                         # initial temperature guess everywhere [K]
xee_init = ptxlib.get_xee_init(Fe_abund) # free electron fraction (electrons per nucleon)

M = 8   # number of angle groups
N = 141 # number of logarithmically spaced energy gridpoints (N-1 bins)

tau_photo = 1.0  # optical depth of upper/lower disk surface
dtau_tau  = 0.1  # logarithmic step size
dtau_max  = 0.25 # maximum step size
tau_max   = 10.0 # 1/2 maximum thickness of a 'whole' slab

e_min = 0 # energy grid boundaries, log10([eV])
e_max = 7

tol = 0.01 # desired tolerance

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

# the ID numbers of the processors which will receive and manage whole slabs
# 0 is of course taken, and 1 thru jac_reserved_procs+1 (inclusive) are taken
slab_mng_proc_ids = [jac_reserved_procs+2]
while (slab_mng_proc_ids[-1] + skip < size):
    slab_mng_proc_ids.append(slab_mng_proc_ids[-1] + skip)

# IDs of worker processors to tackle Jacobian columns
jac_worker_proc_ids = np.arange(2, jac_reserved_procs+2).tolist()

# ----------+----------  end: setup and configuration  ----------+----------

# this is *the* parallelized unit: many instances of calc_slab are spawned, run, and die;
# they know what to do because they are given the phi- and r-indices they are to process in the tuple "slab";
# they fail (relatively) gracefully by reporting an error code triplet; otherwise, they return (0, 0)
def calc_slab(slab):
    phi_index = slab[0]
    r_index   = slab[1]
    slab_type = slab[2]

    result_exists = 'final_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5' in os.listdir('.')
    if result_exists:
        return (slab, 'already_exists', -1)

    # ----------+---------- start: prepare input data ----------+----------

    tmp = ptxlib.get_vert_struct(input_filename, phi_index, r_index, mass, a, tau_photo, dtau_tau, dtau_max, tau_max, xee_init, Fe_abund)

    if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
        return (slab, 'get_vert_struct', tmp)

    to_cleave = tmp[0]
    if not to_cleave:
        r     = tmp[1]
        z     = tmp[2]
        zc    = tmp[3]
        dens  = tmp[4]
        heat  = tmp[5]
        D = len(z)
    else:
        r        = tmp[1]
        z_u      = tmp[2]
        zc_u     = tmp[3]
        dens_u   = tmp[4]
        heat_u   = tmp[5]
        z_l      = tmp[6]
        zc_l     = tmp[7]
        dens_l   = tmp[8]
        heat_l   = tmp[9]
        mid_heat = tmp[10]
        D = len(z_u)

    # we need at least 5 cells to work with
    if (D < 6):
        return (slab, 'too_few_cells', -1)

    # can't do at r < R_hor for obvious reasons
    if (r < 1.0 + np.sqrt(1.0 - a**2)):
        return (slab, 'too_far_in', -1)

    # not worth doing at r/M > 50
    if (r > 50):
        return (slab, 'too_far_out', -1)

    tmp = ptxlib.get_inc_fluxes(phi_index, r_index, e_min, e_max, N)

    if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
        return (slab, 'get_inc_fluxes', tmp)

    energies = tmp[0]
    flux_top = tmp[1]
    flux_bot = tmp[2]

    # ----------+----------  end: prepare input data  ----------+----------

    # only run slabs for which a corresponding output file does not exist in the current working directory
    if not to_cleave:
        tmp = perform_soln(phi_index, r_index, r, z, zc, dens, heat, 0.0, energies, flux_top, flux_bot, 'whole')
    else:
        if (slab_type == 'upper'):
            tmp = perform_soln(phi_index, r_index, r, z_u, zc_u, dens_u, heat_u, mid_heat, energies, flux_top, np.zeros(N-1), slab_type)
        if (slab_type == 'lower'):
            tmp = perform_soln(phi_index, r_index, r, z_l, zc_l, dens_l, heat_l, mid_heat, energies, np.zeros(N-1), flux_bot, slab_type)

    if (tmp < 0):
        return (slab, 'perform_soln', tmp)

    # this is the end of the processing for this slab, and we kick it back to the controlling processor
    return (0, 0)

def perform_soln(phi_index, r_index, r, z, zc, dens, heat, mid_heat, energies, flux_top, flux_bot, slab_type):
    # ----------+---------- start: actual calculation ----------+----------
    start_time = time.time()

    # let us review the relevant dimensions:
    # D       the number of cell edges in the slab
    # D-1     the number of cells, total
    # M       the number of angle groups
    # N       the number of energy bin *boundaries*
    # N-1     the number of energy bins

    D = len(z)
    K = M*(N-1) # the dimension of the "angle-frequency" space

    # the major players in our calculation:
#   mu      = np.linspace(0.5/M, 1.0-(0.5/M), M) # (cosine) angle bins, evenly spaced
#   dmu     = np.array(M*[1.0/M])                # the (uniform) width of each bin
    # --- --- ---
    mu, dmu = np.polynomial.legendre.leggauss(M)
    mu  = 0.5*mu + 0.5
    dmu = 0.5*dmu
    # --- --- ---
    inctop  = np.zeros(K)                        # the incident specific intensity at the top (photon number) [/s/cm^2/ster/eV]
    incbot  = np.zeros(K)                        # the incident specific intensity at the bottom (photon number) [/s/cm^2/ster/eV]
    jf      = np.zeros((D-1, K))                 # mean-intensity-like Feautrier variable, defined at cell centers (photon number) [/s/cm^2/ster/eV]
    hf      = np.zeros((D-2, K))                 # flux-like Feautrier variable, defined at cell boundaries, except first and last (photon number) [/s/cm^2/ster/eV]
    s_jf    = np.zeros((D-1, K))                 # "seed photon" versions of the above
    s_hf    = np.zeros((D-2, K))
    absorp  = np.zeros((D-1, N-1))               # absorption coefficient, defined at cell centers [/cm]
    scatt   = np.zeros((D-1, N-1))               # scattering coefficient, defined at cell centers [/cm]
    emis    = np.zeros((D-1, N-1))               # (presumably isotropic) total emissivity, defined at cell centers (photon number) [/s/cm^3/ster/eV]
    c_emis  = np.zeros((D-1, N-1))               # (presumably isotropic) continuum emissivity *only*, defined at cell centers (photon number) [/s/cm^3/ster/eV]
    dtaub   = np.zeros((D, N-1))                 # change in optical depth, defined at cell boundaries (thus "b")
    dtauc   = np.zeros((D-1, N-1))               # change in optical depth, defined at cell centers (thus "c")
    temp    = np.array((D-1) * [T_init])         # cell-center defined temperature [K], with an initial guess
    xee     = np.array((D-1) * [xee_init])       # free-electron fraction in each cell, with an initial guess

    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'r = ' + repr(r))
    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'D = ' + repr(D))

    prev_file_exists = False
    prev_filename = 'none'
    if ((grand_iter > 0) and os.path.isfile('run_' + repr(grand_iter-1) + '/final_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5')):
        prev_file_exists = True
        prev_filename = './run_' + repr(grand_iter-1) + '/final_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5'

    # now, we assume that the top and bottom fluxes are both isotropic in their respective half-spaces
    for k in range(0, K):
        m = k % M   # this is how we get the angle (m) and energy (n) indices out of their concatenated "k";
        n = int((k-m)/M) # it is super useful and we do it often
        if (slab_type == 'whole'):
            inctop[k] = flux_top[n]/np.pi
            incbot[k] = flux_bot[n]/np.pi
        elif (slab_type == 'upper'):
            inctop[k] = flux_top[n]/np.pi
        elif (slab_type == 'lower'):
            incbot[k] = flux_bot[n]/np.pi

    # set initial Tbnd (for cleaved slabs) and inctop/incbot; use previous iteration's data if possible
    if prev_file_exists:
        with GetFile(prev_filename, open_type = 'r') as f:
            Tbnd = float(f.read('Tbnd'))
    else:
        Tbnd = T_init
    if (slab_type == 'upper'):
        incbot = ptxlib.set_planck_I(energies, mu, Tbnd)
    if (slab_type == 'lower'):
        inctop = ptxlib.set_planck_I(energies, mu, Tbnd)

    # this ensures supplied flux equals the numerical integral (will by default on uniform mu grid)
    flux_bot = np.zeros(N-1)
    for n in range(0, N-1):
        for m in range(0, M):
            flux_bot[n] += 2*np.pi * incbot[n*M + m] * mu[m] * dmu[m]
    flux_top = np.zeros(N-1)
    for n in range(0, N-1):
        for m in range(0, M):
            flux_top[n] += 2*np.pi * inctop[n*M + m] * mu[m] * dmu[m]

    # perform coherent solution with free-free only first to get a starting guess at T vector; but use previous iteration's T if possible
    if prev_file_exists:
        with GetFile(prev_filename, open_type = 'r') as f:
            temp = f.read('temp')
        if (slab_type == 'whole'):
            T = temp.copy()[1:-1]
        else:
            T = temp.copy()[:-1]
            T[0] = Tbnd
    else:
        tmp = ptxlib.coherent_continuum_soln(z, zc, mu, dmu, energies, dens, xee, temp, Tbnd, heat, mid_heat, inctop, incbot, flux_top, flux_bot, slab_type, phi_index, r_index, D, M, N)
    #   if (tmp < 0):
    #       return tmp
        T = tmp

    if (slab_type == 'whole'):
        T = temp.copy()[1:-1]
    else:
        T = temp.copy()[:-1]
        T[0] = Tbnd

    # perform a full treatment, with compton scattering and xstar calls, to get the thermal and photoionization equilibrium solution
    tmp = ptxlib.full_compton_xstar_soln(z, zc, mu, dmu, energies, dens, xee, T, heat, mid_heat, inctop, incbot, flux_top, flux_bot, Fe_abund, slab_type, phi_index, r_index, D, M, N, prev_file_exists, prev_filename, rank)
#   tmp = ptxlib.tde_soln(z, zc, mu, dmu, energies, dens, xee, T, heat, mid_heat, inctop, incbot, flux_top, flux_bot, Fe_abund, slab_type, phi_index, r_index, D, M, N, prev_file_exists, prev_filename, rank)

#   if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
#       return tmp
    succ_code, T, xee, absorp, x_absorp, emis, x_emis, c_emis, jf, hf, s_jf, s_hf, y_flux, global_y = tmp

    if (succ_code == -1):
        print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'failed')
    if (succ_code == 0):
        print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'full success')
    if (succ_code == 1):
        print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'partial success')

    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'final global_y = ' + repr(global_y))

    if ((succ_code == 0) or (succ_code == 1)):
        temp = np.zeros(D-1)
        if (slab_type == 'whole'):
            temp[1:-1] = T
        else:
            temp[1:-1] = T[1:]
            Tbnd       = T[0]
        temp[0]  = temp[1]
        temp[-1] = temp[-2]

        scatt = np.zeros((D-1, N-1))
        for d in range(0, D-1):
            for n in range(0, N-1):
                scatt[d][n] = xee[d] * dens[d] * ptxlib.sigma(0.5*(energies[n+1] + energies[n]))

        end_time = time.time()

        duration = end_time - start_time

        # we're done, write the output file
        ptxlib.write_output('final_' + slab_type + '_' + repr(phi_index) + '_' + repr(r_index) + '.h5', succ_code, D, M, N, r, global_y, y_flux, z, zc, mu, dmu, energies, flux_top, flux_bot, dens, xee, heat, temp, Tbnd, absorp, x_absorp, scatt, emis, x_emis, c_emis, jf, hf, s_jf, s_hf, duration)

        return 0
    else:
        return -1
    # ----------+----------  end: actual calculation  ----------+----------

# ----------+---------- start: program execution ----------+----------

print('hello from ' + repr(rank))

# ----------+---------- root process ----------+----------
if (rank == 0):
    # start the clock
    start_time = time.time()

    # open the log file
    logname = sys.argv[1]
    log     = open(logname, 'w')

    with GetFile('./' + input_filename, open_type = 'r') as f:
        slab_type_ik = f.read('slab_type_ik')

#   if ((grand_iter == 0) and ('run_0' in os.listdir('.'))):
    if ('run_' + repr(grand_iter) in os.listdir('.')):
        # kill *all* processors
        for i in range(1, size):
            COMM_WORLD.send([-1, -1], dest = i, tag = 0)

        # presumably, we're done; clean up
        log = open(logname, 'a')
        log.write('run_' + repr(grand_iter) + ' found: already Finished.\n')

        # stop the clock, record time elapsed
        end_time = time.time()
        elapsed  = end_time - start_time
        log.write('Time elapsed: ' + '{0:.2f}'.format(elapsed/3600) + 'h.\n')

        log.close()

        sys.exit()

    if (grand_iter == 0):
        todo = []
        for i in range(0, len(slab_type_ik[0])):
            for j in range(0, len(slab_type_ik)):
                if (slab_type_ik[j][i] == 0):
                    continue
                if (slab_type_ik[j][i] == 1):
                    todo.append([j, i, 'whole'])
                if (slab_type_ik[j][i] == 2):
                    todo.append([j, i, 'lower'])
                    todo.append([j, i, 'upper'])
        todo = todo[::-1]
    else:
        prev_dir  = 'run_' + repr(grand_iter-1)
        todo = []
        for i in range(0, len(slab_type_ik)):
            for j in range(0, len(slab_type_ik[0])):
                if (slab_type_ik[i][j] == 0):
                    continue
                if (slab_type_ik[i][j] == 1):
                    prev_flnm = 'final_whole_' + repr(i) + '_' + repr(j) + '.h5'
                    if (prev_flnm in os.listdir(prev_dir)):
                        with GetFile('./' + prev_dir + '/' + prev_flnm, open_type = 'r') as f:
                            duration = float(f.read('duration'))
                        todo.append([duration, i, j, 'whole'])
                    else:
                        todo.append([1.0e9, i, j, 'whole'])
                if (slab_type_ik[i][j] == 2):
                    prev_flnm = 'final_upper_' + repr(i) + '_' + repr(j) + '.h5'
                    if (prev_flnm in os.listdir(prev_dir)):
                        with GetFile('./' + prev_dir + '/' + prev_flnm, open_type = 'r') as f:
                            duration = float(f.read('duration'))
                        todo.append([duration, i, j, 'upper'])
                    else:
                        todo.append([1.0e9, i, j, 'upper'])
                    prev_flnm = 'final_lower_' + repr(i) + '_' + repr(j) + '.h5'
                    if (prev_flnm in os.listdir(prev_dir)):
                        with GetFile('./' + prev_dir + '/' + prev_flnm, open_type = 'r') as f:
                            duration = float(f.read('duration'))
                        todo.append([duration, i, j, 'lower'])
                    else:
                        todo.append([1.0e9, i, j, 'lower'])
        todo.sort(key = lambda x: x[0])
        todo = todo[::-1]
        for el in todo:
            del el[0]

    print('len(todo) = ' + repr(len(todo)))

    todo = []
    todo.append([32, 150, 'whole'])

    # send initial set of jobs
    job_index = 0
    for i in slab_mng_proc_ids:
        if (job_index < len(todo)):
            log.write('Starting job at ' + repr(todo[job_index]) + ' on ' + repr(i) + '.\n')
            COMM_WORLD.send([job_index, todo[job_index]], dest = i, tag = 0)
            job_index += 1
    log.close()

    # cycle through processors, receiving results and reassigning work
    finished = 0
    while (finished != len(todo)):
        results = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
        finished += 1
        log = open(logname, 'a')
        # if "results" is an error code, it will have length 4
        if (len(results) == 4):
            log.write('Job at ' + repr(todo[results[0]]) + ' on ' + repr(results[2]) + ' issued: ' + results[1] + '.\n')
        # otherwise, a successful execution must've taken place
        else:
            log.write('Job at ' + repr(todo[results[0]]) + ' on ' + repr(results[2]) + ' finished with ' + repr(results[1]) + '.\n')
        # send next job to this processor
        if (job_index < len(todo)):
            log.write('Starting job at ' + repr(todo[job_index]) + ' on ' + repr(results[2]) + '.\n')
            COMM_WORLD.send([job_index, todo[job_index]], dest = results[2], tag = 0)
            job_index += 1
        log.close()

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send([-1, -1], dest = i, tag = 0)

    os.system('mkdir run_' + repr(grand_iter))
    os.system('mv final_*.h5 ./run_' + repr(grand_iter))
    os.system('mv scat_diskt.h5 ./run_' + repr(grand_iter))
    os.system('mv scat_diskb.h5 ./run_' + repr(grand_iter))

    # presumably, we're done; clean up
    log = open(logname, 'a')
    log.write('Finished.\n')

    # stop the clock, record time elapsed
    end_time = time.time()
    elapsed  = end_time - start_time
    log.write('Time elapsed: ' + '{0:.2f}'.format(elapsed/3600) + 'h.\n')

    log.close()
# ----------+---------- worker process ----------+----------
elif rank in slab_mng_proc_ids:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = 0, tag = 0)

        # break out if sent done signal
        if (not type(tmp[0]) is np.ndarray) and (tmp[0] == -1):
            print('rank ' + repr(rank) + ' quit!')
            break

        # perform calculation
        output = calc_slab(tmp[1])

        # if calc_slab returns an error code, output will have length 3;
        # else, 2; assemble reply tuple in either case
        if (len(output) == 3):
            reply = (tmp[0], output[1] + ' ' + repr(output[2]), rank, 0)
        else:
            reply = (tmp[0], output[1], rank)

        # reply to root
        COMM_WORLD.send(reply, dest = 0, tag = 0)
elif (rank == 1):
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        # break out if sent done signal
        if (not type(tmp[0]) is np.ndarray) and (tmp[0] == -1):
            print('rank ' + repr(rank) + ' quit!')
            break

        slab_mng_rank = tmp[0]
        data          = tmp[1]
        T             = data[0]

        print('Jac calc manager, received data from ' + repr(slab_mng_rank))

        Jac = np.zeros((len(T), len(T)))

        job_index = 0
        for i in jac_worker_proc_ids:
            if (job_index < len(T)):
                COMM_WORLD.send([job_index, data], dest = i, tag = 0)
                job_index += 1

        finished = 0
        while (finished != len(T)):
            col, dydT, jac_worker_rank = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 1)
            finished += 1
            print('Jac calc manager, received col data from ' + repr(jac_worker_rank))
            for i in range(0, len(T)):
                Jac[i][col] = dydT[i]
            if (job_index < len(T)):
                COMM_WORLD.send([job_index, data], dest = jac_worker_rank, tag = 0)
                job_index += 1

        COMM_WORLD.send(Jac, dest = slab_mng_rank, tag = 0)
else:
    while True:
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        # break out if sent done signal
        if (tmp[0] == -1):
            print('rank ' + repr(rank) + ' quit!')
            break

        col, data = tmp

        T, y_flux, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, x_absorp, x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N = data

        T_0       = T.copy()
        y_flux_0  = y_flux.copy()
        T         = T_0.copy()
        T[col]   += 1.0e-3 * T_0[col]

        temp = np.zeros(D-1)
        if (slab_type != 'whole'):
            temp[1:-1] = T[1:]
        else:
            temp[1:-1] = T
        temp[0]  = temp[1]
        temp[-1] = temp[-2]

        absorp = x_absorp + ptxlib.set_brems_absorp(energies, dens, xee, temp, D, N)
        emis   = x_emis   + ptxlib.set_brems_emis(energies, dens, xee, temp, D, N)

        y_flux, global_y, jf, hf = ptxlib.y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
        dydT = (y_flux - y_flux_0)/(1.0e-3 * T_0[col])

        COMM_WORLD.send([col, dydT, rank], dest = 1, tag = 1)
