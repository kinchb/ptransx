import sys
sys.path.append('/turquoise/users/kinch/spin_analysis')

import harm_utils

import os

import time

import numpy as np

import h5py

from scipy.interpolate import RegularGridInterpolator

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.0
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

if (a == 0.0):
    gdump_fname = './KDHARM0.gdump.a0.h5'
if (a == 0.5):
    gdump_fname = './KDHARM0.gdump.a05.h5'
if (a == 0.9):
    gdump_fname = './KDHARM0.gdump.a09.h5'
if (a == 0.99):
    gdump_fname = './KDHARM0.gdump.a099.h5'

dump_list, radflux_list = harm_utils.gather_files()

r, th, phi, dV_cgs = harm_utils.calc_dV_cgs(dump_list[0])

output_fname = os.getcwd().split('/')[-1] + '.h5'

def uconp2ucon_dynamic(data_file, ucon0, ucon1, ucon2, ucon3):
    dx_dxp10 = np.array(data_file['/dx_dxp10'])
    dx_dxp11 = np.array(data_file['/dx_dxp11'])
    dx_dxp12 = np.array(data_file['/dx_dxp12'])
    dx_dxp13 = np.array(data_file['/dx_dxp13'])
    dx_dxp20 = np.array(data_file['/dx_dxp20'])
    dx_dxp21 = np.array(data_file['/dx_dxp21'])
    dx_dxp22 = np.array(data_file['/dx_dxp22'])
    dx_dxp23 = np.array(data_file['/dx_dxp23'])
    dx_dxp30 = np.array(data_file['/dx_dxp30'])
    dx_dxp31 = np.array(data_file['/dx_dxp31'])
    dx_dxp32 = np.array(data_file['/dx_dxp32'])
    dx_dxp33 = np.array(data_file['/dx_dxp33'])

    ucontt = ucon0
    uconrr = dx_dxp10 * ucon0 + dx_dxp11 * ucon1 + dx_dxp12 * ucon2 + dx_dxp13 * ucon3
    uconth = dx_dxp20 * ucon0 + dx_dxp21 * ucon1 + dx_dxp22 * ucon2 + dx_dxp23 * ucon3
    uconph = dx_dxp30 * ucon0 + dx_dxp31 * ucon1 + dx_dxp32 * ucon2 + dx_dxp33 * ucon3

    return ucontt, uconrr, uconth, uconph

def ks2bl_con2(data_file, uconttKS, uconrrKS, uconthKS, uconphKS):
    spin = a

    rr = np.array(data_file['/x1'])
    th = np.array(data_file['/x2'])
    ph = np.array(data_file['/x3'])

    delta = rr*rr - 2.0*rr + spin*spin

    dtBL_drKS   = -2.0*rr/delta
    dphiBL_drKS = -spin/delta

    ucon0BL = uconttKS + dtBL_drKS*uconrrKS
    ucon1BL = uconrrKS
    ucon2BL = uconthKS
    ucon3BL = uconphKS + dphiBL_drKS*uconrrKS

    return ucon0BL, ucon1BL, ucon2BL, ucon3BL

def calc_ucon_bl(data_file, metric):
    v1 = np.array(data_file['/v1'])
    v2 = np.array(data_file['/v2'])
    v3 = np.array(data_file['/v3'])

    vsq = metric['gcov11']*v1*v1 + metric['gcov22']*v2*v2 + metric['gcov33']*v3*v3 + 2.0*(v1*(metric['gcov12']*v2 + metric['gcov13']*v3) + metric['gcov23']*v2*v3)

    gamma = np.sqrt(1.0 + vsq)
    alpha = metric['alpha']
    beta1 = metric['beta1']
    beta2 = metric['beta2']
    beta3 = metric['beta3']

    ucon0 = gamma/alpha
    ucon1 = v1 - ucon0 * beta1
    ucon2 = v2 - ucon0 * beta2
    ucon3 = v3 - ucon0 * beta3

    ucon0_ks, ucon1_ks, ucon2_ks, ucon3_ks = uconp2ucon_dynamic(data_file, ucon0, ucon1, ucon2, ucon3)

    ucon0_bl, ucon1_bl, ucon2_bl, ucon3_bl = ks2bl_con2(data_file, ucon0_ks, ucon1_ks, ucon2_ks, ucon3_ks)

    return ucon0_bl, ucon1_bl, ucon2_bl, ucon3_bl

rank = COMM_WORLD.Get_rank()
size = COMM_WORLD.Get_size()

filelist_1 = dump_list
filelist_2 = radflux_list

if ((rank == 0) or (rank % 2 == 1)):
    start = time.time()

    metric = harm_utils.get_metric(gdump_fname)

    end = time.time()
    print('at rank = ' + repr(rank) + ', setting metric took ' + '{0:.2f}'.format(end - start) + ' s')

    new_r     = np.logspace(np.log10(r[0]), np.log10(r[-1]), 192, endpoint = True)
    new_r[0]  = r[0]
    new_r[-1] = r[-1]
    new_th    = np.linspace(th[0], th[-1], len(th), endpoint = True)

    interp_points = []
    for i in range(0, len(new_r)):
        for j in range(0, len(new_th)):
            for k in range(0, len(phi)):
                interp_points.append([new_r[i], new_th[j], phi[k]])
    interp_points = np.array(interp_points)

# ----------+---------- parallel code ----------+----------

print 'hello from ' + repr(rank)
sys.stdout.flush()
sys.stderr.flush()

# ----------+---------- root process ----------+----------
if (rank == 0):
    output_file = h5py.File(output_fname, 'w')

    output_file['/M']    = np.array([M])
    output_file['/a']    = np.array([a])
    output_file['/mdot'] = np.array([mdot])
    output_file['/eta']  = np.array([eta])
    output_file['/r']    = new_r
    output_file['/th']   = new_th
    output_file['/phi']  = phi

    output_file.close()

    todo = []
    for i in range(0, len(filelist_1)):
        todo.append(i)

    times = np.zeros(len(todo))

    job_index = 0
    for i in range(1, size, 2):
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index]], dest = i, tag = 0)
            job_index += 1

    finished = 0
    while (finished != len(todo)):
        result    = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)
        finished += 1
        print repr(finished) + '/' + repr(len(todo))
        # record result
        i            = result[0]
        times[i]     = result[1]
        new_rho      = result[2]
        new_uu       = result[3]
        new_urad     = result[4]
        new_T_C      = result[5]
        new_coolfunc = result[6]
        new_ucon0    = result[7]
        new_ucon1    = result[8]
        new_ucon2    = result[9]
        new_ucon3    = result[10]
        start = time.time()
        while True:
            try:
                output_file = h5py.File(output_fname, 'r+')
            except:
                continue
            break
        output_file['/rho_'      + repr(i)] = new_rho
        output_file['/uu_'       + repr(i)] = new_uu
        output_file['/urad_'     + repr(i)] = new_urad
        output_file['/T_C_'      + repr(i)] = new_T_C
        output_file['/coolfunc_' + repr(i)] = new_coolfunc
        output_file['/ucon0_'    + repr(i)] = new_ucon0
        output_file['/ucon1_'    + repr(i)] = new_ucon1
        output_file['/ucon2_'    + repr(i)] = new_ucon2
        output_file['/ucon3_'    + repr(i)] = new_ucon3
        output_file.close()
        end = time.time()
        print('at rank = ' + repr(rank) + ', writing output took ' + '{0:.2f}'.format(end - start) + ' s')
        # send next job to this processor
        if (job_index < len(todo)):
            COMM_WORLD.send([job_index, todo[job_index]], dest = result[11], tag = 0)
            job_index += 1

    while True:
        try:
            output_file = h5py.File(output_fname, 'r+')
        except:
            continue
        break

    output_file['/times'] = times

    output_file.close()

    # kill *all* processors
    for i in range(1, size):
        COMM_WORLD.send([-1, -1], dest = i, tag = 0)

# ----------+---------- worker process ----------+----------
else:
    while True:
        # receive job
        tmp = COMM_WORLD.recv(source = MPI.ANY_SOURCE, tag = 0)

        # break out if sent done signal
        if (tmp[0] == -1):
            break

        harm_file_1 = h5py.File(filelist_1[tmp[0]], 'r')
        harm_file_2 = h5py.File(filelist_2[tmp[0]], 'r')

        t        = float(np.array(harm_file_1['/Header/Grid/t'])[0])
        rho      = np.array(harm_file_1['/rho'])
        uu       = np.array(harm_file_1['/uu'])
        urad     = np.array(harm_file_2['/urad'])
        T_C      = np.array(harm_file_2['/T_C'])
        coolfunc = np.array(harm_file_2['/coolfunc'])
        v1       = np.array(harm_file_1['/v1'])
        v2       = np.array(harm_file_1['/v2'])
        v3       = np.array(harm_file_1['/v3'])

        harm_file_2.close()

        start = time.time()

        ucon0, ucon1, ucon2, ucon3 = calc_ucon_bl(harm_file_1, metric)

        harm_file_1.close()

        end = time.time()
        print('at rank = ' + repr(rank) + ', calc_ucon_bl took ' + '{0:.2f}'.format(end - start) + ' s')

        start = time.time()

        new_rho     = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), rho)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_rho[i, j, k] = interp_vals[l]
                    l += 1

        new_uu      = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), uu)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_uu[i, j, k] = interp_vals[l]
                    l += 1

        new_urad    = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), urad)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_urad[i, j, k] = interp_vals[l]
                    l += 1

        new_T_C     = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), T_C)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_T_C[i, j, k] = interp_vals[l]
                    l += 1

        new_coolfunc = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj   = RegularGridInterpolator((r, th, phi), coolfunc)
        interp_vals  = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_coolfunc[i, j, k] = interp_vals[l]
                    l += 1

        new_ucon0   = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), ucon0)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_ucon0[i, j, k] = interp_vals[l]
                    l += 1

        new_ucon1   = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), ucon1)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_ucon1[i, j, k] = interp_vals[l]
                    l += 1

        new_ucon2   = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), ucon2)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_ucon2[i, j, k] = interp_vals[l]
                    l += 1

        new_ucon3   = np.zeros((len(new_r), len(new_th), len(phi)))
        interp_obj  = RegularGridInterpolator((r, th, phi), ucon3)
        interp_vals = interp_obj(interp_points)
        l = 0
        for i in range(0, len(new_r)):
            for j in range(0, len(new_th)):
                for k in range(0, len(phi)):
                    new_ucon3[i, j, k] = interp_vals[l]
                    l += 1

        end = time.time()
        print('at rank = ' + repr(rank) + ', interpolation took ' + '{0:.2f}'.format(end - start) + ' s')

        reply = (tmp[0], t, new_rho, new_uu, new_urad, new_T_C, new_coolfunc, new_ucon0, new_ucon1, new_ucon2, new_ucon3, rank)

        COMM_WORLD.send(reply, dest = 0, tag = 0)
