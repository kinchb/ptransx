#
# PTXLIB
# this file contains the implementation of many of the functions called in ptransx.py
#

import os

import sys

import time

import random

import numpy as np

from scipy.interpolate import interp1d
from scipy.interpolate import interp2d

from scipy.optimize import root

import h5py

from mpi4py import MPI
from mpi4py.MPI import COMM_WORLD

# from _maketps import maketps

from _fsolver import fsolver

from _xstarcomm import xstarcomm
from _xstarcomm import fbg

import matplotlib.pyplot as plt

from h5py_wrapper import GetFile

def spec_trapz(y, x):
    N = len(x)

    result = 0.0
    i = 0
    while (i < N-1):
        result += y[i] * (x[i+1] - x[i])
        i += 1

    return result

def get_xee_init(Fe_abund):
    Z_list = np.linspace(1, 30, num = 30, endpoint = True)

    mass_frac_list = [0.704,
                      0.28,
                      6.1E-11,
                      1.58E-10,
                      3.78E-09,
                      0.00278,
                      0.000814,
                      0.00756,
                      0.000000418,
                      0.00169,
                      0.0000343,
                      0.000645,
                      0.0000556,
                      0.000696,
                      0.0000061,
                      0.000479,
                      0.00000783,
                      0.0000701,
                      0.0000036,
                      0.0000641,
                      4.64E-08,
                      0.0000035,
                      0.000000356,
                      0.000017,
                      0.00000942,
                      0.00123,
                      0.00000342,
                      0.0000729,
                      0.00000072,
                      0.00000182]

    amu_list = [1.0079,
                4.0026,
                6.941,
                9.0122,
                10.811,
                12.0107,
                14.0067,
                15.9994,
                18.9984,
                20.1797,
                22.9897,
                24.305,
                26.9815,
                28.0855,
                30.9738,
                32.065,
                35.453,
                39.948,
                39.0983,
                40.078,
                44.9559,
                47.867,
                50.9415,
                51.9961,
                54.938,
                55.845,
                58.9332,
                58.6934,
                63.546,
                65.39]

    mass_frac_list = np.array(mass_frac_list)
    amu_list       = np.array(amu_list)

    mass_frac_list[25] = Fe_abund * mass_frac_list[25]
    mass_frac_list     = mass_frac_list/mass_frac_list.sum()

    num_frac_list = mass_frac_list/amu_list
    num_frac_list = num_frac_list/num_frac_list.sum()

    xee = (num_frac_list * Z_list).sum()

    return xee

def get_density(input_filename, phi_index, r_index, Fe_abund):
    with GetFile('./' + input_filename, open_type = 'r') as f:
        r   = f.read('r')
        th  = f.read('th')
        rho = f.read('rho')

    dens = np.zeros(len(th))
    for i in range(0, len(th)):
        dens[i] = rho[r_index][i][phi_index]

    # compute the conversion factor from g/cm^3 to number of nuclei per cm^3;
    # assume solar abundances for all elements but Fe, using data from Grevesse & Sauval 1998

    mass_frac_list = [0.704,
                      0.28,
                      6.1E-11,
                      1.58E-10,
                      3.78E-09,
                      0.00278,
                      0.000814,
                      0.00756,
                      0.000000418,
                      0.00169,
                      0.0000343,
                      0.000645,
                      0.0000556,
                      0.000696,
                      0.0000061,
                      0.000479,
                      0.00000783,
                      0.0000701,
                      0.0000036,
                      0.0000641,
                      4.64E-08,
                      0.0000035,
                      0.000000356,
                      0.000017,
                      0.00000942,
                      0.00123,
                      0.00000342,
                      0.0000729,
                      0.00000072,
                      0.00000182]

    amu_list = [1.0079,
                4.0026,
                6.941,
                9.0122,
                10.811,
                12.0107,
                14.0067,
                15.9994,
                18.9984,
                20.1797,
                22.9897,
                24.305,
                26.9815,
                28.0855,
                30.9738,
                32.065,
                35.453,
                39.948,
                39.0983,
                40.078,
                44.9559,
                47.867,
                50.9415,
                51.9961,
                54.938,
                55.845,
                58.9332,
                58.6934,
                63.546,
                65.39]

    mass_frac_list = np.array(mass_frac_list)
    amu_list       = np.array(amu_list)

    mass_frac_list[25] = Fe_abund * mass_frac_list[25]
    mass_frac_list     = mass_frac_list/mass_frac_list.sum()

    conv = (mass_frac_list/amu_list).sum()/1.661e-24

    dens = conv * np.array(dens)

    return (r[r_index], th, dens)

def get_heating(input_filename, phi_index, r_index):
    with GetFile('./' + input_filename, open_type = 'r') as f:
        th       = f.read('th')
        coolfunc = f.read('coolfunc')

    heat = np.zeros(len(th))
    for i in range(0, len(th)):
        heat[i] = coolfunc[r_index][i][phi_index]

    return heat

def get_vert_struct(input_filename, phi_index, r_index, mass, a, tau_photo, dtau_tau, dtau_max, tau_max, xee_init, Fe_abund):
    tmp = get_density(input_filename, phi_index, r_index, Fe_abund)

    if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
        return tmp

    r    = tmp[0]
    th   = tmp[1]
    dens = tmp[2]

    heat = get_heating(input_filename, phi_index, r_index)

    m_bh = mass * 2.0e33
    G    = 6.6726e-8
    c    = 3.0e10

    with GetFile('./' + input_filename, open_type = 'r') as f:
        emtop_ik = f.read('emtop_ik')
        embot_ik = f.read('embot_ik')

    th_top = emtop_ik[phi_index][r_index]
    th_bot = embot_ik[phi_index][r_index]

    tmp_th   = np.linspace(th_top, th_bot, 1e5, endpoint = True)
    tmp_dens = interp1d(th, dens, kind = 'linear')(tmp_th)
    tmp_heat = interp1d(th, heat, kind = 'linear')(tmp_th)

    th   = tmp_th
    dens = tmp_dens
    heat = tmp_heat

    s = np.zeros(len(th))

    for i in range(1, len(s)):
        s[i] = s[i-1] + ((G*m_bh)/(c*c)) * np.sqrt(r**2 + a**2 * np.cos(0.5*(th[i] + th[i-1]))**2) * (th[i] - th[i-1])

    tau_top    = np.zeros(len(s))
    tau_top[0] = tau_photo

    for i in range(1, len(s)):
        tau_top[i] = tau_top[i-1] + xee_init * 6.6524e-25 * 0.5 * (dens[i] + dens[i-1]) * (s[i] - s[i-1])

    tau_bot     = np.zeros(len(s))
    tau_bot[-1] = tau_photo

    for i in range(len(s)-2, -1, -1):
        tau_bot[i] = tau_bot[i+1] + xee_init * 6.6524e-25 * 0.5 * (dens[i] + dens[i+1]) * (s[i+1] - s[i])

    eq_index = 0
    curr_min = abs(tau_top[0] - tau_bot[0])
    for i in range(0, len(s)):
        if (abs(tau_top[i] - tau_bot[i]) < curr_min):
            curr_min = abs(tau_top[i] - tau_bot[i])
            eq_index = i

    mid_tau = min(tau_top[eq_index], tau_bot[eq_index])
    if (mid_tau < tau_max):
        to_cleave = False
    else:
        to_cleave = True

    tau_steps = []

    tau_steps.append(tau_photo)

    if not to_cleave:
        next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        while ((next_step < dtau_max) and (tau_steps[len(tau_steps)-1] + next_step < mid_tau)):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
            next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        next_step = dtau_max
        while (tau_steps[len(tau_steps)-1] + next_step < mid_tau):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
        if (len(tau_steps) < 2):
            return [-1, -1]
        tau_steps.insert(1, tau_photo + 0.0001)
    else:
        next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        while ((next_step < dtau_max) and (tau_steps[len(tau_steps)-1] + next_step < tau_max)):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
            next_step = tau_steps[len(tau_steps)-1] * dtau_tau
        next_step = dtau_max
        while (tau_steps[len(tau_steps)-1] + next_step < tau_max):
            tau_steps.append(tau_steps[len(tau_steps)-1] + next_step)
        if (tau_steps[len(tau_steps)-1] < tau_max - 0.0001):
            tau_steps.append(tau_max - 0.0001)
        tau_steps.append(tau_max)
        tau_steps.insert(1, tau_photo + 0.0001)

    tau_steps = np.array(tau_steps)

    tau_top_indices = []
    tau_top_indices.append(0)
    for i in range(1, len(tau_steps)):
        index = tau_top_indices[i-1] + 1
        curr_min = abs(tau_top[index] - tau_steps[i])
        for j in range(index+1, len(tau_top)):
            if (abs(tau_top[j] - tau_steps[i]) < curr_min):
                curr_min = abs(tau_top[j] - tau_steps[i])
                index = j
        tau_top_indices.append(index)

    tau_bot_indices = []
    tau_bot_indices.append(len(tau_bot)-1)
    for i in range(1, len(tau_steps)):
        index = tau_bot_indices[i-1] - 1
        curr_min = abs(tau_bot[index] - tau_steps[i])
        for j in range(index-1, -1, -1):
            if (abs(tau_bot[j] - tau_steps[i]) < curr_min):
                curr_min = abs(tau_bot[j] - tau_steps[i])
                index = j
        tau_bot_indices.append(index)

    tau_bot_indices = tau_bot_indices[::-1]

    if (tau_top_indices[-1] >= tau_bot_indices[0]):
        tau_top_indices = tau_top_indices[:-1]
        tau_bot_indices = tau_bot_indices[1:]

    if ((not to_cleave) and (tau_top_indices[-1] < eq_index) and (tau_bot_indices[0] > eq_index)):
        tau_top_indices.append(eq_index)

    # check that the combined indices are sorted and without duplicates
    if (not np.all(np.array(tau_top_indices + tau_bot_indices) == np.sort(np.array(tau_top_indices + tau_bot_indices)))):
        print('array poorly sorted')
        return [-1, -1]
    if (len(np.unique(np.array(tau_top_indices + tau_bot_indices))) != len(np.array(tau_top_indices + tau_bot_indices))):
        print('array contains duplicates')
        return [-1, -1]

    if not to_cleave:
        combined = tau_top_indices + tau_bot_indices

        z = s[combined]

        zc = np.zeros(len(z)-1)

        for i in range(0, len(zc)):
            zc[i] = 0.5 * (z[i] + z[i+1])

        dens = np.zeros(len(zc))

        for i in range(0, len(dens)):
            dens[i] = (tau_top[combined[i+1]] - tau_top[combined[i]])/(xee_init * 6.6524e-25 * (z[i+1] - z[i]))

        int_heat = np.zeros(len(s))

        for i in range(1, len(s)):
            int_heat[i] = int_heat[i-1] + 0.5 * (heat[i] + heat[i-1]) * (s[i] - s[i-1])

        heat = np.zeros(len(zc))

        for i in range(0, len(heat)):
            heat[i] = (int_heat[combined[i+1]] - int_heat[combined[i]])/(z[i+1] - z[i])

        dens[0]  = dens[1]
        dens[-1] = dens[-2]
        heat[0]  = heat[1]
        heat[-1] = heat[-2]

        return (to_cleave, r, z, zc, dens, heat)
    else:
        z_u = s[tau_top_indices]
        z_l = s[tau_bot_indices]

        zc_u = np.zeros(len(z_u)-1)
        zc_l = np.zeros(len(z_l)-1)

        for i in range(0, len(zc_u)):
            zc_u[i] = 0.5 * (z_u[i] + z_u[i+1])
        for i in range(0, len(zc_l)):
            zc_l[i] = 0.5 * (z_l[i] + z_l[i+1])

        dens_u = np.zeros(len(zc_u))
        dens_l = np.zeros(len(zc_l))

        for i in range(0, len(dens_u)):
            dens_u[i] = (tau_top[tau_top_indices[i+1]] - tau_top[tau_top_indices[i]])/(xee_init * 6.6524e-25 * (z_u[i+1] - z_u[i]))
        for i in range(0, len(dens_l)):
            dens_l[i] = (tau_bot[tau_bot_indices[i]] - tau_bot[tau_bot_indices[i+1]])/(xee_init * 6.6524e-25 * (z_l[i+1] - z_l[i]))

        int_heat = np.zeros(len(s))

        for i in range(1, len(s)):
            int_heat[i] = int_heat[i-1] + 0.5 * (heat[i] + heat[i-1]) * (s[i] - s[i-1])

        heat_u = np.zeros(len(zc_u))
        heat_l = np.zeros(len(zc_l))

        for i in range(0, len(heat_u)):
            heat_u[i] = (int_heat[tau_top_indices[i+1]] - int_heat[tau_top_indices[i]])/(z_u[i+1] - z_u[i])
        for i in range(0, len(heat_l)):
            heat_l[i] = (int_heat[tau_bot_indices[i+1]] - int_heat[tau_bot_indices[i]])/(z_l[i+1] - z_l[i])

        # each half of the full vertical slab gets half the in-between heating... thus the 0.5
        mid_heat = 0.5*(int_heat[tau_bot_indices[0]] - int_heat[tau_top_indices[-1]])

        dens_u[0]  = dens_u[1]
        dens_u[-1] = dens_u[-2]
        heat_u[0]  = heat_u[1]
        heat_u[-1] = heat_u[-2]
        dens_l[0]  = dens_l[1]
        dens_l[-1] = dens_l[-2]
        heat_l[0]  = heat_l[1]
        heat_l[-1] = heat_l[-2]

        return (to_cleave, r, z_u, zc_u, dens_u, heat_u, z_l, zc_l, dens_l, heat_l, mid_heat)

def get_spectra(phi_index, r_index, top_or_bot):
    if (top_or_bot == 'top'):
        flnm = './scat_diskt.h5'
    else:
        flnm = './scat_diskb.h5'

    with GetFile(flnm, open_type = 'r') as f:
        spec = f.read('data')

    Ne     = len(spec)
    Nphi   = len(spec[0])
    Nr     = len(spec[0][0])

    e_grid = np.logspace(0.0, 7.0, Ne, endpoint = True)

    spec2 = np.zeros((Nphi, Nr, Ne))
    for i in range(0, Nphi):
        for j in range(0, Nr):
            for k in range(0, Ne):
                spec2[i][j][k] = spec[k][i][j]

    # convert erg/s/cm^2/Hz to erg/s/cm^2/erg
    spec2 = 1.5091902e26 * spec2

    return (e_grid, spec2[phi_index][r_index])

def get_inc_fluxes(phi_index, r_index, e_min, e_max, N):
    tmp = get_spectra(phi_index, r_index, 'top')

    if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
        return tmp

    e_grid   = tmp[0]
    flux_top = tmp[1]

    tmp = get_spectra(phi_index, r_index, 'bot')

    if (not type(tmp[0]) is np.ndarray) and (tmp[0] < 0):
        return tmp

    e_grid   = tmp[0]
    flux_bot = tmp[1]

    energies = np.logspace(e_min, e_max, N, endpoint = True)

    f_top = interp1d(e_grid, flux_top, kind = 'linear')
    f_bot = interp1d(e_grid, flux_bot, kind = 'linear')

    flux_top = np.zeros(N-1)
    flux_bot = np.zeros(N-1)

    for n in range(0, N-1):
        flux_top[n] = (f_top(energies[n+1]) + f_top(energies[n]))/(energies[n+1] + energies[n])
        flux_bot[n] = (f_bot(energies[n+1]) + f_bot(energies[n]))/(energies[n+1] + energies[n])

    return (energies, flux_top, flux_bot)

def set_brems_emis(energies, dens, xee, temp, D, N):
    emis = np.zeros((D-1, N-1))

    for d in range(0, D-1):
        if (temp[d] == 0.0):
            continue
        for n in range(0, N-1):
            enrg = 0.5*(energies[n+1] + energies[n])
            gaunt = 1.0
            if (enrg/(8.617e-5 * temp[d]) < 100.0):
                gaunt = fbg(enrg/(8.617e-5 * temp[d]), 0.158/(temp[d]/1.0e6))
            if (np.isnan(gaunt) or np.isinf(gaunt)):
                gaunt = 1.0
            if (gaunt < 0.5):
                gaunt = 0.5
            if (gaunt > 5.0):
                gaunt = 5.0
            emis[d][n] = (1.0/(4.0*np.pi)) * (1.0/enrg) * 1.026e-11 * gaunt * xee[d] * (dens[d])**2 * (temp[d])**(-0.5) * np.exp(-enrg/(8.617e-5 * temp[d]))
            if (np.isnan(emis[d][n]) or np.isinf(emis[d][n])):
                emis[d][n] = 0.0

    return emis

def set_brems_absorp(energies, dens, xee, temp, D, N):
    absorp = np.zeros((D-1, N-1))

    coef = 2.0/((4.1357e-15)**3 * (3.0e10)**2)

    for d in range(0, D-1):
        if (temp[d] == 0.0):
            continue
        for n in range(0, N-1):
            enrg = 0.5*(energies[n+1] + energies[n])
            gaunt = 1.0
            if (enrg/(8.617e-5 * temp[d]) < 100.0):
                gaunt = fbg(enrg/(8.617e-5 * temp[d]), 0.158/(temp[d]/1.0e6))
            if (np.isnan(gaunt) or np.isinf(gaunt)):
                gaunt = 1.0
            if (gaunt < 0.5):
                gaunt = 0.5
            if (gaunt > 5.0):
                gaunt = 5.0
#           absorp[d][n] = ((1.0/(4.0*np.pi)) * 1.026e-11 * gaunt * xee[d] * (dens[d])**2 * (temp[d])**(-0.5) * np.exp(-enrg/(8.617e-5 * temp[d])) * (np.exp(enrg/(8.617e-5 * temp[d])) - 1.0))/(coef * enrg**3)
            absorp[d][n] = ((1.0/(4.0*np.pi)) * 1.026e-11 * gaunt * xee[d] * (dens[d])**2 * (temp[d])**(-0.5) * (1.0 - np.exp(-enrg/(8.617e-5 * temp[d]))))/(coef * enrg**3)
            if (np.isnan(absorp[d][n]) or np.isinf(absorp[d][n])):
                absorp[d][n] = 0.0

    return absorp

def sigma(e):
    if (e < 100.0):
        eps = e/511.0e3
        return np.pi*(2.81794032e-13)**2 * (8.0/3.0 - (16.0/3.0)*eps + (208.0/15.0)*eps**2 - (522.0/15.0)*eps**3)
    else:
        eps = e/511.0e3
        return np.pi*(2.81794032e-13)**2 * ((2*eps*(2 + eps*(1 + eps)*(8 + eps)))/(1 + 2*eps)**2 + (-2 + (-2 + eps)*eps)*np.log(1 + 2*eps))/eps**3

def set_planck_I(energies, mu, T):
    M = len(mu)
    N = len(energies)
    K = M*(N-1)

    I = np.zeros(K)

    coef = 2.0/((4.1357e-15)**3 * (3.0e10)**2)
    kT   = 8.6173e-5 * T
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        E = 0.5*(energies[n+1] + energies[n])
        I[k] = coef * (E**2/(np.exp(E/kT)-1.0))
        if (np.isnan(I[k]) or np.isinf(I[k])):
            I[k] = 0.0

    return I

# Mihalas (1985) 2.18-2.19
def set_dtau(z, absorp, scatt, D, N):
    dtaub = np.zeros((D, N-1))
    dtauc = np.zeros((D-1, N-1))

    for d in range(1, D-1):
        for n in range(0, N-1):
            dtaub[d, n] = 0.5*((absorp[d-1, n]+scatt[d-1, n])*(z[d]-z[d-1]) + (absorp[d, n]+scatt[d, n])*(z[d+1]-z[d]))
    for n in range(0, N-1):
        dtaub[0, n]   = 0.5*(absorp[0, n]+scatt[0, n])*(z[1]-z[0])
        dtaub[D-1, n] = 0.5*(absorp[D-2, n]+scatt[D-2, n])*(z[D-1]-z[D-2])
    for d in range(0, D-1):
        for n in range(0, N-1):
            dtauc[d, n] = 0.5*(dtaub[d, n] + dtaub[d+1, n])

    return (dtaub, dtauc)

# Mihalas (1985) 2.17/2.20/2.30
def set_A(mu, dtaub, dtauc, D, M, N):
    K = M*(N-1)

    A = np.zeros((D-1, N-1, M, M))

    for d in range(1, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            A[d, n, m, m] = mu[m]**2/(dtauc[d, n]*dtaub[d, n])
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        A[D-2, n, m, m] = mu[m]/dtaub[D-2, n]

    return A

# Mihalas (1985) 2.17/2.20/2.30
def set_B(mu, dtaub, dtauc, absorp, scatt, dmu, D, M, N):
    K = M*(N-1)

    phase = np.zeros((M, M))

    for m in range(0, M):
        for mm in range(0, M):
            phase[m, mm] = 0.5 * (1.0 + (mu[m] * mu[mm] + np.sqrt(1.0 - mu[m] * mu[m]) * np.sqrt(1.0 - mu[mm] * mu[mm]))**2)
    for m in range(0, M):
        tmp = 0.0
        for mm in range(0, M):
            tmp += phase[mm, m] * dmu[mm]
        for mm in range(0, M):
            phase[mm, m] = phase[mm, m]/tmp

    B = np.zeros((D-1, N-1, M, M))

    for d in range(1, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            B[d, n, m, m] = 1.0 + (mu[m]**2/dtauc[d, n])*(1/dtaub[d, n]+1/dtaub[d+1, n])
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            for mm in range(0, M):
                B[d, n, m, mm] -= (scatt[d, n]/(absorp[d, n]+scatt[d, n]))*dmu[mm] * phase[m, mm]
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        B[0, n, m, m]   = mu[m]/dtaub[1, n] + 1/(1+dtauc[0, n]/(2*mu[m])) + dtauc[0, n]/mu[m]
        B[D-2, n, m, m] = mu[m]/dtaub[D-2, n] + 1/(1+dtauc[D-2, n]/(2*mu[m])) + dtauc[D-2, n]/mu[m]
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        for mm in range(0, M):
            B[0, n, m, mm]   -= (dtauc[0, n]/mu[m])*(scatt[0, n]/(absorp[0, n]+scatt[0, n]))*dmu[mm] * phase[m, mm]
            B[D-2, n, m, mm] -= (dtauc[D-2, n]/mu[m])*(scatt[D-2, n]/(absorp[D-2, n]+scatt[D-2, n]))*dmu[mm] * phase[m, mm]

    return B

# Mihalas (1985) 2.17/2.20/2.30
def set_C(mu, dtaub, dtauc, D, M, N):
    K = M*(N-1)

    C = np.zeros((D-1, N-1, M, M))

    for d in range(1, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            C[d, n, m, m] = mu[m]**2/(dtauc[d, n]*dtaub[d+1, n])
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        C[0, n, m, m] = mu[m]/dtaub[1, n]

    return C

# Mihalas (1985) 2.17/2.20/2.30
def set_L(mu, dtaub, dtauc, absorp, scatt, emis, inctop, incbot, D, M, K):
    L = np.zeros((D-1, K))

    for d in range(1, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            L[d, k] = (1/(absorp[d, n]+scatt[d, n]))*emis[d, n]
    for k in range(0, K):
        m = k % M
        n = int((k-m)/M)
        L[0, k]   = inctop[k]/(1+dtauc[0, n]/(2*mu[m])) + (dtauc[0, n]/mu[m])*(1/(absorp[0, n]+scatt[0, n]))*emis[0, n]
        L[D-2, k] = incbot[k]/(1+dtauc[D-2, n]/(2*mu[m])) + (dtauc[D-2, n]/mu[m])*(1/(absorp[D-2, n]+scatt[D-2, n]))*emis[D-2, n]

    return L

def cfsolver(absorp, scatt, emis, dtaub, dtauc, inctop, incbot, mu, dmu, D, M, N):
    K = M*(N-1)

    G  = np.zeros((D-1, N-1, M, M))
    v  = np.zeros((D-1, K))
    jf = np.zeros((D-1, K))

    # fill out vectors and matrices; see Mihalas (1985) and the above for how these arrays are populated
    A = set_A(mu, dtaub, dtauc, D, M, N)
    B = set_B(mu, dtaub, dtauc, absorp, scatt, dmu, D, M, N)
    C = set_C(mu, dtaub, dtauc, D, M, N)
    L = set_L(mu, dtaub, dtauc, absorp, scatt, emis, inctop, incbot, D, M, K)

    # determine jf, essentially by solving 2.29-2.33 in Mihalas (1985);
    # this is just a straightforward implementation of that, though with the explicit assumption of coherent scattering
    for d in range(0, D-2):
        for n in range(0, N-1):
            tmp1    = np.dot(A[d, n], G[d-1, n])
            tmp2    = np.linalg.inv(B[d, n] - tmp1)
            G[d, n] = np.dot(tmp2, C[d, n])

    for d in range(0, D-1):
        for n in range(0, N-1):
            tmp1 = np.dot(A[d, n], G[d-1, n])
            tmp2 = np.linalg.inv(B[d, n] - tmp1)
            tmp3 = L[d, n*M:(n+1)*M] + np.dot(A[d, n], v[d-1, n*M:(n+1)*M])
            v[d, n*M:(n+1)*M] = np.dot(tmp2, tmp3)

    jf[D-2] = v[D-2]
    for d in range(D-3, -1, -1):
        for n in range(0, N-1):
            jf[d, n*M:(n+1)*M] = np.dot(G[d, n], jf[d+1, n*M:(n+1)*M]) + v[d, n*M:(n+1)*M]

    hf = np.zeros((D-2, K))
    for d in range(0, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            hf[d, k] = mu[m] * (jf[d+1, k] - jf[d, k])/dtaub[d+1, n]

    return (jf, hf)

def compute_ys(z, zc, mu, dmu, energies, heat, dtaub, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N):
    K = M*(N-1)

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5 * (energies[n+1] + energies[n])

    jfb = np.zeros((K, D-2))

    for k in range(0, K):
        jfb[k] = interp1d(zc, jf.T[k])(z[1:-1])

    jfb = jfb.T

    flux_up = np.zeros((D, N-1))
    flux_dn = np.zeros((D, N-1))

    for d in range(1, D-1):
        for n in range(0, N-1):
            for m in range(0, M):
                flux_up[d][n] += 2*np.pi * (jfb[d-1][n*M + m] + hf[d-1][n*M + m]) * mu[m] * dmu[m] * 1.602e-12 * ec[n]
                flux_dn[d][n] += 2*np.pi * (jfb[d-1][n*M + m] - hf[d-1][n*M + m]) * mu[m] * dmu[m] * 1.602e-12 * ec[n]

    for d in range(1, D-1):
        for n in range(0, N-1):
            flux_up[d][n] = 2*np.pi * ((jfb[d-1][n*M:n*M + M] + hf[d-1][n*M:n*M + M]) * mu * dmu).sum()
            flux_dn[d][n] = 2*np.pi * ((jfb[d-1][n*M:n*M + M] - hf[d-1][n*M:n*M + M]) * mu * dmu).sum()

    div_flux = np.zeros((D-1, N-1))
    for d in range(1, D-2):
        for n in range(0, N-1):
            div_flux[d][n] = (mu**2 * (-dtaub[d][n] * jf[d+1][n*M:n*M + M] + (dtaub[d][n] + dtaub[d+1][n]) * jf[d][n*M:n*M + M] - dtaub[d+1][n] * jf[d-1][n*M:n*M + M]) * dmu).sum()/(dtaub[d][n] * dtaub[d+1][n])

    int_heat = np.zeros(D-1)
    for d in range(0, D-1):
        int_heat[d] = heat[d] * (z[d+1] - z[d])

    int_div_flux = np.zeros(D-1)
    for d in range(1, D-2):
        int_div_flux[d] = 1.602e-12 * 4*np.pi * np.trapz(div_flux[d] * ec, ec)

    """
    flux_dn[0]  = flux_top * 1.602e-12 * ec
    flux_up[-1] = flux_bot * 1.602e-12 * ec
    """

    flux_dn[0]  = flux_top
    flux_up[-1] = flux_bot

    """
    for n in range(0, N-1):
        for m in range(0, M):
            flux_up[0][n]  += 4*np.pi * hf[0][n*M + m] * mu[m] * dmu[m] * 1.602e-12 * ec[n]
            flux_dn[-1][n] -= 4*np.pi * hf[-1][n*M + m] * mu[m] * dmu[m] * 1.602e-12 * ec[n]
        flux_up[0][n]  += flux_top[n] * 1.602e-12 * ec[n]
        flux_dn[-1][n] += flux_bot[n] * 1.602e-12 * ec[n]
    """

    for n in range(0, N-1):
        flux_up[0][n]  = 4*np.pi * (hf[0][n*M:n*M + M] * mu * dmu).sum()
        flux_dn[-1][n] = -4*np.pi * (hf[-1][n*M:n*M + M] * mu * dmu).sum()
    flux_up[0]  += flux_top
    flux_dn[-1] += flux_bot

    y      = np.zeros(D-1)
    y_flux = np.zeros(D-1)
    for d in range(0, D-1):
        """
        y[d]      = 2*(spec_trapz(flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1], energies) - int_heat[d])/(spec_trapz(flux_up[d] + flux_dn[d+1] + flux_dn[d] + flux_up[d+1], energies) + int_heat[d])
        y_flux[d] = spec_trapz(flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1], energies) - int_heat[d]
        """
        y[d]      = 2*(1.602e-12 * np.trapz((flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1]) * ec, ec) - int_heat[d])/(1.602e-12 * np.trapz((flux_up[d] + flux_dn[d+1] + flux_dn[d] + flux_up[d+1]) * ec, ec) + int_heat[d])
#       y_flux[d] = 1.602e-12 * np.trapz((flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1]) * ec, ec) - int_heat[d]
        y_flux[d] = int_div_flux[d] - int_heat[d]
#       print d
#       print 1.602e-12 * np.trapz((flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1]) * ec, ec)
#       print int_div_flux[d]

    if (slab_type == 'upper'):
        """
        y[0]      = 2*(spec_trapz(flux_up[D-1] - flux_dn[D-1], energies) - mid_heat)/(spec_trapz(flux_up[D-1] + flux_dn[D-1], energies) + mid_heat)
        y_flux[0] = spec_trapz(flux_up[D-1] - flux_dn[D-1], energies) - mid_heat
        """
        y[0]      = 2*(1.602e-12 * np.trapz((flux_up[D-1] - flux_dn[D-1]) * ec, ec) - mid_heat)/(1.602e-12 * np.trapz((flux_up[D-1] + flux_dn[D-1]) * ec, ec) + mid_heat)
        y_flux[0] = 1.602e-12 * np.trapz((flux_up[D-1] - flux_dn[D-1]) * ec, ec) - mid_heat
    if (slab_type == 'lower'):
        """
        y[0]      = 2*(spec_trapz(flux_dn[0] - flux_up[0], energies) - mid_heat)/(spec_trapz(flux_dn[0] + flux_up[0], energies) + mid_heat)
        y_flux[0] = spec_trapz(flux_dn[0] - flux_up[0], energies) - mid_heat
        """
        y[0]      = 2*(1.602e-12 * np.trapz((flux_dn[0] - flux_up[0]) * ec, ec) - mid_heat)/(1.602e-12 * np.trapz((flux_dn[0] + flux_up[0]) * ec, ec) + mid_heat)
        y_flux[0] = 1.602e-12 * np.trapz((flux_dn[0] - flux_up[0]) * ec, ec) - mid_heat

    if (slab_type == 'whole'):
        """
        flux_in  = spec_trapz((flux_top + flux_bot) * 1.602e-12 * ec, energies)
        flux_out = spec_trapz(flux_up[0] + flux_dn[D-1], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat)))/(flux_out + (flux_in + sum(int_heat)))
        """
        flux_in  = 1.602e-12 * np.trapz((flux_top + flux_bot) * ec, ec)
        flux_out = 1.602e-12 * np.trapz((flux_up[0] + flux_dn[D-1]) * ec, ec)
        global_y = 2*(flux_out - (flux_in + sum(int_heat)))/(flux_out + (flux_in + sum(int_heat)))
        return (y[1:-1], y_flux[1:-1], global_y)
    if (slab_type == 'upper'):
        """
        flux_in  = spec_trapz(flux_top * 1.602e-12 * ec, energies)
        flux_out = spec_trapz(flux_up[0], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))
        """
        flux_in  = 1.602e-12 * np.trapz(flux_top * ec, ec)
        flux_out = 1.602e-12 * np.trapz(flux_up[0] * ec, ec)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))
        return (y[:-1], y_flux[:-1], global_y)
    if (slab_type == 'lower'):
        """
        flux_in  = spec_trapz(flux_bot * 1.602e-12 * ec, energies)
        flux_out = spec_trapz(flux_dn[D-1], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))
        """
        flux_in  = 1.602e-12 * np.trapz(flux_bot * ec, ec)
        flux_out = 1.602e-12 * np.trapz(flux_dn[D-1] * ec, ec)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))

        return (y[:-1], y_flux[:-1], global_y)

def compute_ys_high_prec(z, zc, mu, dmu, energies, heat, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N):
    import mpmath as mp

    mp.mp.dps = 100

    K = M*(N-1)

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5 * (energies[n+1] + energies[n])

    ec = mp.mpf('1.0') * ec

    jfb = np.zeros((K, D-2))

    for k in range(0, K):
        jfb[k] = interp1d(zc, jf.T[k])(z[1:-1])

    jfb = jfb.T

    jfb = mp.mpf('1.0') * jfb

    flux_up = np.zeros((D, N-1))
    flux_dn = np.zeros((D, N-1))

    flux_up = mp.mpf('1.0') * flux_up
    flux_dn = mp.mpf('1.0') * flux_dn

    """
    for d in range(1, D-1):
        for n in range(0, N-1):
            for m in range(0, M):
                flux_up[d][n] += 2*np.pi * (jfb[d-1][n*M + m] + hf[d-1][n*M + m]) * mu[m] * dmu[m] * 1.602e-12 * ec[n]
                flux_dn[d][n] += 2*np.pi * (jfb[d-1][n*M + m] - hf[d-1][n*M + m]) * mu[m] * dmu[m] * 1.602e-12 * ec[n]
    """

    mu  = mp.mpf('1.0') * mu
    dmu = mp.mpf('1.0') * dmu

    for d in range(1, D-1):
        for n in range(0, N-1):
            flux_up[d][n] = 2*mp.pi * ((jfb[d-1][n*M:n*M + M] + hf[d-1][n*M:n*M + M]) * mu * dmu).sum()
            flux_dn[d][n] = 2*mp.pi * ((jfb[d-1][n*M:n*M + M] - hf[d-1][n*M:n*M + M]) * mu * dmu).sum()

    """
    div_flux = np.zeros((D-1, N-1)):
    for d in range(1, D-2):
        div_flux[d][n] = mu**2 * (dtaub[d] * (jf[d+1][n*M:n*M + M] - jf[d][n*M:n*M + M])
    """

    """
    flux_dn[0]  = flux_top * 1.602e-12 * ec
    flux_up[-1] = flux_bot * 1.602e-12 * ec
    """

    flux_dn[0]  = mp.mpf('1.0') * flux_top
    flux_up[-1] = mp.mpf('1.0') * flux_bot

    """
    for n in range(0, N-1):
        for m in range(0, M):
            flux_up[0][n]  += 4*np.pi * hf[0][n*M + m] * mu[m] * dmu[m] * 1.602e-12 * ec[n]
            flux_dn[-1][n] -= 4*np.pi * hf[-1][n*M + m] * mu[m] * dmu[m] * 1.602e-12 * ec[n]
        flux_up[0][n]  += flux_top[n] * 1.602e-12 * ec[n]
        flux_dn[-1][n] += flux_bot[n] * 1.602e-12 * ec[n]
    """

    for n in range(0, N-1):
        flux_up[0][n]  = 4*mp.pi * (hf[0][n*M:n*M + M] * mu * dmu).sum()
        flux_dn[-1][n] = -4*mp.pi * (hf[-1][n*M:n*M + M] * mu * dmu).sum()
    flux_up[0]  += flux_top
    flux_dn[-1] += flux_bot

    int_heat = mp.mpf('1.0') * np.zeros(D-1)

    heat = mp.mpf('1.0') * heat
    z    = mp.mpf('1.0') * z

    for d in range(0, D-1):
        int_heat[d] = heat[d] * (z[d+1] - z[d])

    y      = mp.mpf('1.0') * np.zeros(D-1)
    y_flux = mp.mpf('1.0') * np.zeros(D-1)
    for d in range(0, D-1):
        """
        y[d]      = 2*(spec_trapz(flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1], energies) - int_heat[d])/(spec_trapz(flux_up[d] + flux_dn[d+1] + flux_dn[d] + flux_up[d+1], energies) + int_heat[d])
        y_flux[d] = spec_trapz(flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1], energies) - int_heat[d]
        """
        y[d]      = 2*(mp.mpf('1.602e-12') * np.trapz((flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1]) * ec, ec) - int_heat[d])/(mp.mpf('1.602e-12') * np.trapz((flux_up[d] + flux_dn[d+1] + flux_dn[d] + flux_up[d+1]) * ec, ec) + int_heat[d])
        y_flux[d] = mp.mpf('1.602e-12') * np.trapz((flux_up[d] + flux_dn[d+1] - flux_dn[d] - flux_up[d+1]) * ec, ec) - int_heat[d]

    if (slab_type == 'upper'):
        y[0]      = 2*(spec_trapz(flux_up[D-1] - flux_dn[D-1], energies) - mid_heat)/(spec_trapz(flux_up[D-1] + flux_dn[D-1], energies) + mid_heat)
        y_flux[0] = spec_trapz(flux_up[D-1] - flux_dn[D-1], energies) - mid_heat
    if (slab_type == 'lower'):
        """
        y[0]      = 2*(spec_trapz(flux_dn[0] - flux_up[0], energies) - mid_heat)/(spec_trapz(flux_dn[0] + flux_up[0], energies) + mid_heat)
        y_flux[0] = spec_trapz(flux_dn[0] - flux_up[0], energies) - mid_heat
        """
        y[0]      = 2*(mp.mpf('1.602e-12') * np.trapz((flux_dn[0] - flux_up[0]) * ec, ec) - mid_heat)/(mp.mpf('1.602e-12') * np.trapz((flux_dn[0] + flux_up[0]) * ec, ec) + mid_heat)
        y_flux[0] = mp.mpf('1.602e-12') * np.trapz((flux_dn[0] - flux_up[0]) * ec, ec) - mid_heat

    if (slab_type == 'whole'):
        flux_in  = spec_trapz((flux_top + flux_bot) * mp.mpf('1.602e-12') * ec, energies)
        flux_out = spec_trapz(flux_up[0] + flux_dn[D-1], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat)))/(flux_out + (flux_in + sum(int_heat)))
        return (y[1:-1], y_flux[1:-1], global_y)
    if (slab_type == 'upper'):
        flux_in  = spec_trapz(flux_top * mp.mpf('1.602e-12') * ec, energies)
        flux_out = spec_trapz(flux_up[0], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))
        return (y[:-1], y_flux[:-1], global_y)
    if (slab_type == 'lower'):
        """
        flux_in  = spec_trapz(flux_bot * mp.mpf('1.602e-12') * ec, energies)
        flux_out = spec_trapz(flux_dn[D-1], energies)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))
        """
        flux_in  = mp.mpf('1.602e-12') * np.trapz(flux_bot * ec, ec)
        flux_out = mp.mpf('1.602e-12') * np.trapz(flux_dn[D-1] * ec, ec)
        global_y = 2*(flux_out - (flux_in + sum(int_heat) + mid_heat))/(flux_out + (flux_in + sum(int_heat) + mid_heat))

        return (np.array(y[:-1], dtype = float), np.array(y_flux[:-1], dtype = float), float(global_y))

def logtable_ndx(value_list, value):
    if (value <= value_list[0]):
        ndx   = 0
        value = value_list[0]
        return (ndx, value)
    elif (value >= value_list[-1]):
        ndx   = len(value_list)-2
        value = value_list[-1]
        return (ndx, value)
    else:
        ndx = int(np.log(value/value_list[0])/np.log(value_list[1]/value_list[0]))
        if ((value > value_list[ndx]) and (value <= value_list[ndx+1])):
            return (ndx, value)
        else:
            for ndx in range(0, len(value_list)-1):
                if ((value > value_list[ndx]) and (value <= value_list[ndx+1])):
                    return (ndx, value)
            print('in logtable_ndx: something went horribly wrong')
            print('in logtable_ndx: value = ' + repr(value) + ', value_list[0] = ' + repr(value_list[0]) + ', value_list[-1] = ' + repr(value_list[-1]))

def construct_cdf(T_list, E_list, eta_grid, T, E, f):
    i, T = logtable_ndx(T_list, T)
    j, E = logtable_ndx(E_list, E)

    x  = T
    x1 = T_list[i]
    x2 = T_list[i+1]
    y  = E
    y1 = E_list[j]
    y2 = E_list[j+1]

    q11 = f.read('cdf_' + repr(i)   + '_' + repr(j))
    q12 = f.read('cdf_' + repr(i)   + '_' + repr(j+1))
    q21 = f.read('cdf_' + repr(i+1) + '_' + repr(j))
    q22 = f.read('cdf_' + repr(i+1) + '_' + repr(j+1))

    cdf = (1.0/((x2 - x1)*(y2 - y1)))*(q11*(x2 - x)*(y2 - y) + q21*(x - x1)*(y2 - y) + q12*(x2 - x)*(y - y1) + q22*(x - x1)*(y - y1))

    e_grid = E * eta_grid

    cdf = interp1d(e_grid, cdf, bounds_error = False, fill_value = (0.0, 1.0))

    return cdf

def interp_sa_ratio(T_list, E_list, T, E, sa_ratio):
    i, T = logtable_ndx(T_list, T)
    j, E = logtable_ndx(E_list, E)

    x  = T
    x1 = T_list[i]
    x2 = T_list[i+1]
    y  = E
    y1 = E_list[j]
    y2 = E_list[j+1]

    q11 = sa_ratio[i,   j]
    q12 = sa_ratio[i,   j+1]
    q21 = sa_ratio[i+1, j]
    q22 = sa_ratio[i+1, j+1]

    return (1.0/((x2 - x1)*(y2 - y1)))*(q11*(x2 - x)*(y2 - y) + q21*(x - x1)*(y2 - y) + q12*(x2 - x)*(y - y1) + q22*(x - x1)*(y - y1))

def maketps_h5py(energies, temp, D, N):
    MIN_PROB = 0.5e-7

    compton_data_flnm = '/lustre/scratch4/turquoise/.mdt3/kinch/compy/compton_data.h5'

    with GetFile(compton_data_flnm, open_type = 'r') as f:
        T_list   = f.read('T_list')
        E_list   = f.read('E_list')
        eta_grid = f.read('eta_grid')
        sa_ratio = f.read('sa_ratio')

        ec = np.zeros(N-1)
        for n in range(0, N-1):
            ec[n] = 0.5 * (energies[n] + energies[n+1])

        tps = []

        for d in range(0, D-1):
            if (d > 0):
                tps[start_point] = len(tps)
            start_point = len(tps)
            tps.append(0.0)
            for n_col in range(0, N-1):
                e_sample = ec[n_col]
                cdf = construct_cdf(T_list, E_list, eta_grid, temp[d], e_sample, f)
                tmp_phis = np.zeros(N-1)
                # scattering into bins below n_col
                if (n_col > 0):
                    n_curr = n_col-1
                    while True:
                        prob = cdf(energies[n_curr+1]) - cdf(energies[n_curr])
                        if (prob >= MIN_PROB):
                            tmp_phis[n_curr] = prob
                        if ((prob < MIN_PROB) or (n_curr == 0)):
                            break
                        n_curr -= 1
                # scattering into n_col
                prob = cdf(energies[n_col+1]) - cdf(energies[n_col])
                tmp_phis[n_col] = prob
                # scattering into bins above n_col
                if (n_col < N-2):
                    n_curr = n_col+1
                    while True:
                        prob = cdf(energies[n_curr+1]) - cdf(energies[n_curr])
                        if (prob >= MIN_PROB):
                            tmp_phis[n_curr] = prob
                        if ((prob < MIN_PROB) or (n_curr == N-2)):
                            break
                        n_curr += 1
                # adjust tmp_phis so that its mean matches the target semi-analytic mean amplification ratio
                target = e_sample * interp_sa_ratio(T_list, E_list, temp[d], e_sample, sa_ratio)
                mean = (tmp_phis * ec).sum()
                for n in range(0, N-1):
                    tmp_phis[n] = tmp_phis[n]*(target/mean)*((energies[n_col+1] - energies[n_col])/(energies[n+1] - energies[n]))
                for n in range(0, N-1):
                    if (tmp_phis[n] > 0.0):
                        tps.append(float(n_col))
                        tps.append(float(n))
                        tps.append(tmp_phis[n])
        tps[start_point] = len(tps)

    return np.array(tps)

def y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N):
    K = M*(N-1)

    temp = np.zeros(D-1)
    if (slab_type != 'whole'):
        Tbnd = T[0]
        temp[1:-1] = T[1:]
        if (slab_type == 'upper'):
            incbot = set_planck_I(energies, mu, Tbnd)
        if (slab_type == 'lower'):
            inctop = set_planck_I(energies, mu, Tbnd)
    else:
        temp[1:-1] = T
    temp[0]  = temp[1]
    temp[-1] = temp[-2]

    # now, perform Compton solution at this T (+ Tbnd) with supplied source function
    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    dtaub, dtauc = set_dtau(z, absorp, scatt, D, N)

    start = time.time()

#   tps = maketps(energies.copy(), temp.copy(), D, N)
    tps = maketps_h5py(energies.copy(), temp.copy(), D, N)

    end = time.time()

#   print 'tps took ' + repr(end - start)

    start = time.time()

    jf = fsolver(absorp.copy(), scatt.copy(), emis.copy(), dtaub.copy(), dtauc.copy(), inctop.copy(), incbot.copy(), mu.copy(), dmu.copy(), energies.copy(), tps.copy(), D, N)

    end = time.time()

#   print 'fsolver took ' + repr(end - start)

    hf = np.zeros((D-2, K))
    for d in range(0, D-2):
        for k in range(0, K):
            m = k % M
            n = int((k-m)/M)
            hf[d][k] = mu[m] * (jf[d+1][k] - jf[d][k])/dtaub[d+1][n]

    # set flux_top or flux_bot if slab_type is not whole
    if (slab_type == 'upper'):
        flux_bot = np.zeros(N-1)
        for n in range(0, N-1):
            for m in range(0, M):
                flux_bot[n] += 2*np.pi * incbot[n*M + m] * mu[m] * dmu[m]
    if (slab_type == 'lower'):
        flux_top = np.zeros(N-1)
        for n in range(0, N-1):
            for m in range(0, M):
                flux_top[n] += 2*np.pi * inctop[n*M + m] * mu[m] * dmu[m]

    # get corresponding y vector
    y, y_flux, global_y = compute_ys(z, zc, mu, dmu, energies, heat, dtaub, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N)

    return (y_flux, global_y, jf, hf)

def new_T_from_Jac(T_0, y_flux_0, Jac):
    # if we're here, then presumably Jac cannot be inverted;
    # check first for rows or cols which are all zeros
    zero_col_indxs = []
    zero_row_indxs = []
    for i in range(0, len(T_0)):
        if np.all(Jac[i,::] == np.zeros(len(T_0))):
            zero_row_indxs.append(i)
        if np.all(Jac[::,i] == np.zeros(len(T_0))):
            zero_col_indxs.append(i)

    # check for positive definite (> 0) diagonal entries
    neg_diag_indxs = []
    for i in range(0, len(T_0)):
        if (Jac[i][i] <= 0.0):
            neg_diag_indxs.append(i)

    bad_indxs = list(set().union(zero_col_indxs, zero_row_indxs, neg_diag_indxs))

    diag_only = False

    if (len(bad_indxs) > 0):
        indx_map = []
        for i in range(0, len(T_0)):
            if (i not in bad_indxs):
                indx_map.append(i)
        if (len(indx_map) > 1):
            abbr_T_0      = np.zeros(len(indx_map))
            abbr_y_flux_0 = np.zeros(len(indx_map))
            abbr_Jac      = np.zeros((len(indx_map), len(indx_map)))
            for i in range(0, len(indx_map)):
                abbr_T_0[i]      = T_0[indx_map[i]]
                abbr_y_flux_0[i] = y_flux_0[indx_map[i]]
            for i in range(0, len(indx_map)):
                for j in range(0, len(indx_map)):
                    abbr_Jac[i][j] = Jac[indx_map[i]][indx_map[j]]
            try:
                abbr_dT = np.linalg.solve(abbr_Jac, -abbr_y_flux_0)
            except LinAlgError:
                diag_only = True
            T = T_0.copy()
            for i in range(0, len(indx_map)):
                T[indx_map[i]] += abbr_dT[i]
        else:
            diag_only = True
    else:
        try:
            dT = np.linalg.solve(Jac, -y_flux_0)
        except LinAlgError:
            diag_only = True
            print('detected singular matrix: ')
            print('y_flux_0 = ' + repr(y_flux_0))
            print('T_0 = ' + repr(T_0))
        T  = T_0.copy()
        T += dT

    if diag_only:
        T = T_0.copy()
        for i in range(0, len(T_0)):
            if (Jac[i][i] > 0.0):
                T[i] += -y_flux_0[i]/Jac[i][i]

    if np.all(T == T_0):
        for i in range(0, len(T)):
            if (np.random.rand() < 0.5):
                T[i] += 0.1 * T[i]
            else:
                T[i] -= 0.1 * T[i]

    return T

def coherent_continuum_soln(z, zc, mu, dmu, energies, dens, xee, temp, Tbnd, heat, mid_heat, inctop, incbot, flux_top, flux_bot, slab_type, phi_index, r_index, D, M, N):
    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    absorp       = set_brems_absorp(energies, dens, xee, temp, D, N)
    emis         = set_brems_emis(energies, dens, xee, temp, D, N)
    dtaub, dtauc = set_dtau(z, absorp, scatt, D, N)

    jf, hf = cfsolver(absorp, scatt, emis, dtaub, dtauc, inctop, incbot, mu, dmu, D, M, N)

    y, y_flux, global_y = compute_ys(z, zc, mu, dmu, energies, heat, dtaub, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N)

    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'coherent continuum soln, start, global_y = ' + repr(global_y))

    if (slab_type == 'whole'):
        T = temp.copy()[1:-1]
    else:
        T = temp.copy()[:-1]
        T[0] = Tbnd

    frac = 1.0e-5
    k = 0
    while True:
        T_0      = T.copy()
        y_flux_0 = y_flux.copy()
        Jac      = np.zeros((len(T_0), len(T_0)))
        for i in range(0, len(T_0)):
            T = T_0.copy()
            T[i] += frac * T_0[i]
            temp = np.zeros(D-1)
            if (slab_type != 'whole'):
                Tbnd = T[0]
                temp[1:-1] = T[1:]
                if (slab_type == 'upper'):
                    incbot = set_planck_I(energies, mu, Tbnd)
                if (slab_type == 'lower'):
                    inctop = set_planck_I(energies, mu, Tbnd)
            else:
                temp[1:-1] = T
            temp[0]  = temp[1]
            temp[-1] = temp[-2]
            absorp       = set_brems_absorp(energies, dens, xee, temp, D, N)
            emis         = set_brems_emis(energies, dens, xee, temp, D, N)
            dtaub, dtauc = set_dtau(z, absorp, scatt, D, N)
            jf, hf = cfsolver(absorp, scatt, emis, dtaub, dtauc, inctop, incbot, mu, dmu, D, M, N)
            if (slab_type == 'upper'):
                flux_bot = np.zeros(N-1)
                for n in range(0, N-1):
                    for m in range(0, M):
                        flux_bot[n] += 2*np.pi * incbot[n*M + m] * mu[m] * dmu[m]
            if (slab_type == 'lower'):
                flux_top = np.zeros(N-1)
                for n in range(0, N-1):
                    for m in range(0, M):
                        flux_top[n] += 2*np.pi * inctop[n*M + m] * mu[m] * dmu[m]
            y, y_flux, global_y = compute_ys(z, zc, mu, dmu, energies, heat, dtaub, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N)
            dydT = (y_flux - y_flux_0)/(frac * T_0[i])
            """
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'y_flux = ' + repr(y_flux)
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'y_flux_0 = ' + repr(y_flux_0)
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'T_0[' + repr(i) + '] = ' + repr(T_0[i])
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'dydT = ' + repr(dydT)
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'emis = ' + repr(emis)
            print '(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' +  'absorp = ' + repr(absorp)
            """
            for j in range(0, len(Jac[0])):
                Jac[j][i] = dydT[j]
        T = np.dot(np.linalg.inv(Jac), -y_flux_0) + T_0
        for i in range(0, len(T)):
            if (T[i] < 0.1*T_0[i]):
                T[i] = 0.1*T_0[i]
            if (T[i] > 10.0*T_0[i]):
                T[i] = 10.*T_0[i]
            if (T[i] < 1.0e3):
                T[i] = 1.0e3
            if (T[i] > 1.0e9):
                T[i] = 1.0e9
        temp = np.zeros(D-1)
        if (slab_type != 'whole'):
            Tbnd = T[0]
            temp[1:-1] = T[1:]
            if (slab_type == 'upper'):
                incbot = set_planck_I(energies, mu, Tbnd)
            if (slab_type == 'lower'):
                inctop = set_planck_I(energies, mu, Tbnd)
        else:
            temp[1:-1] = T
        temp[0]  = temp[1]
        temp[-1] = temp[-2]
        absorp       = set_brems_absorp(energies, dens, xee, temp, D, N)
        emis         = set_brems_emis(energies, dens, xee, temp, D, N)
        dtaub, dtauc = set_dtau(z, absorp, scatt, D, N)
        jf, hf = cfsolver(absorp, scatt, emis, dtaub, dtauc, inctop, incbot, mu, dmu, D, M, N)
        if (slab_type == 'upper'):
            flux_bot = np.zeros(N-1)
            for n in range(0, N-1):
                for m in range(0, M):
                    flux_bot[n] += 2*np.pi * incbot[n*M + m] * mu[m] * dmu[m]
        if (slab_type == 'lower'):
            flux_top = np.zeros(N-1)
            for n in range(0, N-1):
                for m in range(0, M):
                    flux_top[n] += 2*np.pi * inctop[n*M + m] * mu[m] * dmu[m]
        y, y_flux, global_y = compute_ys(z, zc, mu, dmu, energies, heat, dtaub, jf, hf, flux_top, flux_bot, mid_heat, slab_type, D, M, N)
        print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'coherent continuum soln, k = ' + repr(k) + ', global_y = ' + repr(global_y))
        k += 1
        if ((abs(global_y) < 0.001) or (k > 5)):
            break

    return T

def full_compton_xstar_soln(z, zc, mu, dmu, energies, dens, xee, T, heat, mid_heat, inctop, incbot, flux_top, flux_bot, Fe_abund, slab_type, phi_index, r_index, D, M, N, prev_file_exists, prev_file, rank):
    K = M*(N-1)

    temp = np.zeros(D-1)
    if (slab_type != 'whole'):
        Tbnd = T[0]
        temp[1:-1] = T[1:]
        if (slab_type == 'upper'):
            incbot = set_planck_I(energies, mu, Tbnd)
        if (slab_type == 'lower'):
            inctop = set_planck_I(energies, mu, Tbnd)
    else:
        temp[1:-1] = T
    temp[0]  = temp[1]
    temp[-1] = temp[-2]

    xstar_called = False
    all_elements = False
    if prev_file_exists:
#       succ_code = int(np.array(prev_file['/succ_code']))
#       if (succ_code == 0):
#           all_elements = True
        all_elements = True
        with GetFile(prev_file, open_type = 'r') as f:
            xee      = f.read('xee')
            x_absorp = f.read('x_absorp')
            x_emis   = f.read('x_emis')
            c_emis   = f.read('c_emis')
        xstar_called = True
    else:
        x_absorp = 0.0 * set_brems_absorp(energies, dens, xee, temp, D, N)
        x_emis   = 0.0 * set_brems_emis(energies, dens, xee, temp, D, N)
        c_emis   = x_emis.copy()
    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
    emis   = set_brems_emis(energies, dens, xee, temp, D, N)

    if xstar_called:
        y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp + x_absorp, emis + x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
    else:
        y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)

    s_jf = 0.0 * jf.copy()
    s_hf = 0.0 * hf.copy()

    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'after Compton xfer, global_y = ' + repr(global_y))

    outgoing_fluxes_top = []
    outgoing_fluxes_bot = []

    # succ_code == 0  : full success
    # succ_code == -1 : full failure
    # succ_code == +1 : partial success; solution found with only H, He, and Fe

    k = 0
    l = 0
    Jac_calls = 0
    finished  = False
    succ_code = 0
    while (not finished):
        if ((Jac_calls > 20) or (l > 20)):
            print('taking a long time')
            global_y = 0.049
#           finished = True
#           if (succ_code == 0):
#               succ_code = -1
#               break
        if (abs(global_y) >= 0.05):
            T_0      = T.copy()
            y_flux_0 = y_flux.copy()
            start_ndx = 0
            if (slab_type != 'whole'):
                start_ndx = 1
            for i in range(start_ndx, len(T_0)):
#               if (np.random.random() < 0.5):
                if (i % 2 == 0):
                    T[i] += 1.0e-3 * T_0[i]
                else:
                    T[i] -= 1.0e-3 * T_0[i]
            temp = np.zeros(D-1)
            if (slab_type != 'whole'):
                Tbnd = T[0]
                temp[1:-1] = T[1:]
                if (slab_type == 'upper'):
                    incbot = set_planck_I(energies, mu, Tbnd)
                if (slab_type == 'lower'):
                    inctop = set_planck_I(energies, mu, Tbnd)
            else:
                temp[1:-1] = T
            temp[0]  = temp[1]
            temp[-1] = temp[-2]
            if (not xstar_called):
                absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
                emis   = set_brems_emis(energies, dens, xee, temp, D, N)
            else:
                absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
            y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
            dydT = np.zeros(len(T))
            if (slab_type != 'whole'):
                dydT[1:] = (y_flux[1:] - y_flux_0[1:])/(T[1:] - T_0[1:])
            else:
                dydT = (y_flux - y_flux_0)/(T - T_0)
            if (slab_type == 'whole'):
                for i in range(0, len(T)):
                    if (dydT[i] > 0.0):
                        T[i] = -y_flux_0[i]/dydT[i] + T_0[i]
                for i in range(0, len(T)):
                    if (T[i] < 0.1*T_0[i]):
                        T[i] = 0.1*T_0[i]
                    if (T[i] > 10.0*T_0[i]):
                        T[i] = 10.0*T_0[i]
                    if (T[i] < 1.0e4):
                        T[i] = 1.0e4
                    if (T[i] > 1.0e12):
                        T[i] = 1.0e12
                temp = np.zeros(D-1)
                temp[1:-1] = T
                temp[0]  = temp[1]
                temp[-1] = temp[-2]
                if (not xstar_called):
                    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = set_brems_emis(energies, dens, xee, temp, D, N)
                else:
                    absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
                y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'full compton xstar soln, k = ' + repr(k) + ', global_y   = ' + repr(global_y))
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'full compton xstar soln, k = ' + repr(k) + ', quad_error = ' + repr(np.sqrt((y_flux**2).sum())))
            else:
                T = T_0.copy()
                T[0] += 1.0e-3 * T_0[i]
                temp = np.zeros(D-1)
                Tbnd = T[0]
                temp[1:-1] = T[1:]
                if (slab_type == 'upper'):
                    incbot = set_planck_I(energies, mu, Tbnd)
                if (slab_type == 'lower'):
                    inctop = set_planck_I(energies, mu, Tbnd)
                temp[0]  = temp[1]
                temp[-1] = temp[-2]
                if (not xstar_called):
                    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = set_brems_emis(energies, dens, xee, temp, D, N)
                else:
                    absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
                y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
                dydTbnd = np.zeros(len(T))
                dydTbnd = (y_flux - y_flux_0)/(T[0] - T_0[0])
                """
                T_map         = []
                T_0_abbr      = []
                y_flux_0_abbr = []
                T_map.append(0)
                T_0_abbr.append(T_0[0])
                y_flux_0_abbr.append(y_flux_0[0])
                for i in range(1, len(T)):
                    if (dydT[i] > 0.0):
                        T_map.append(i)
                        T_0_abbr.append(T_0[i])
                        y_flux_0_abbr.append(y_flux_0[i])
                T_0_abbr      = np.array(T_0_abbr)
                y_flux_0_abbr = np.array(y_flux_0_abbr)
#               print 'T_map: '
#               print T_map
                if (len(T_map) > 1):
                    Jac = np.zeros((len(T_map), len(T_map)))
                    for i in range(1, len(T_map)):
                        Jac[i][i] = dydT[T_map[i]]
                    for i in range(0, len(T_map)):
                        Jac[i][0] = dydTbnd[i]
                    T_abbr = np.dot(np.linalg.inv(Jac), -y_flux_0_abbr) + T_0_abbr
                    for i in range(0, len(T_map)):
                        T[T_map[i]] = T_abbr[i]
                else:
                    T[0] = -y_flux_0[0]/dydTbnd[0] + T_0[0]
                """
                Jac = np.zeros((len(T_0), len(T_0)))
                for i in range(1, len(T_0)):
                    Jac[i][i] = dydT[i]
                for i in range(0, len(T_0)):
                    Jac[i][0] = dydTbnd[i]
                T = new_T_from_Jac(T_0, y_flux_0, Jac)
                for i in range(0, len(T)):
                    if (T[i] < 0.1*T_0[i]):
                        T[i] = 0.1*T_0[i]
                    if (T[i] > 10.0*T_0[i]):
                        T[i] = 10.0*T_0[i]
                    if (T[i] < 1.0e4):
                        T[i] = 1.0e4
                    if (T[i] > 1.0e12):
                        T[i] = 1.0e12
                temp = np.zeros(D-1)
                Tbnd = T[0]
                temp[1:-1] = T[1:]
                if (slab_type == 'upper'):
                    incbot = set_planck_I(energies, mu, Tbnd)
                if (slab_type == 'lower'):
                    inctop = set_planck_I(energies, mu, Tbnd)
                temp[0]  = temp[1]
                temp[-1] = temp[-2]
                if (not xstar_called):
                    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = set_brems_emis(energies, dens, xee, temp, D, N)
                else:
                    absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
                y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'full compton xstar soln, k = ' + repr(k) + ', global_y   = ' + repr(global_y))
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'full compton xstar soln, k = ' + repr(k) + ', quad_error = ' + repr(np.sqrt((y_flux**2).sum())))
            k += 1
            if (k > 10):
                k = 0
                Jac_calls += 1
                T_0      = T.copy()
                y_flux_0 = y_flux.copy()
                if (not xstar_called):
                    data = [T, y_flux, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, np.zeros((D-1, N-1)), np.zeros((D-1, N-1)), inctop, incbot, flux_top, flux_bot, slab_type, D, M, N]
                else:
                    data = [T, y_flux, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, x_absorp, x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N]
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'at rank = ' + repr(rank) + ', about to send off for Jac processing, call = ' + repr(Jac_calls))
                COMM_WORLD.send([rank, data], dest = 1, tag = 0)
                Jac = COMM_WORLD.recv(source = 1, tag = 0)
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'at rank = ' + repr(rank) + ', received Jac')
#               T = np.dot(np.linalg.inv(Jac), -y_flux_0) + T_0
                T = new_T_from_Jac(T_0, y_flux_0, Jac)
                for i in range(0, len(T)):
                    if (T[i] < 0.1*T_0[i]):
                        T[i] = 0.1*T_0[i]
                    if (T[i] > 10.0*T_0[i]):
                        T[i] = 10.0*T_0[i]
                    if (T[i] < 1.0e4):
                        T[i] = 1.0e4
                    if (T[i] > 1.0e12):
                        T[i] = 1.0e12
                temp = np.zeros(D-1)
                if (slab_type != 'whole'):
                    Tbnd = T[0]
                    temp[1:-1] = T[1:]
                    if (slab_type == 'upper'):
                        incbot = set_planck_I(energies, mu, Tbnd)
                    if (slab_type == 'lower'):
                        inctop = set_planck_I(energies, mu, Tbnd)
                else:
                    temp[1:-1] = T
                temp[0]  = temp[1]
                temp[-1] = temp[-2]
                if (not xstar_called):
                    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = set_brems_emis(energies, dens, xee, temp, D, N)
                else:
                    absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                    emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
                y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'post Jac full compton xstar soln, global_y   = ' + repr(global_y))
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'post Jac full compton xstar soln, quad_error = ' + repr(np.sqrt((y_flux**2).sum())))
        if (abs(global_y) < 0.05):
            still_small = True
            while still_small:
                still_small = False
                k = 0
                Jac_calls = 0
                xstar_start = time.time()
                xee, x_absorp, x_emis, c_emis = call_xstar(energies, dmu, dens, temp, jf, Fe_abund, all_elements, D, M, N)
                xstar_called = True
                l += 1
                xstar_end = time.time()
#               print 'xstar calls took: ' + repr((xstar_end - xstar_start)/60.0)
                absorp = x_absorp + set_brems_absorp(energies, dens, xee, temp, D, N)
                emis   = x_emis   + set_brems_emis(energies, dens, xee, temp, D, N)
                y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)
                print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'after xstar update, l = ' + repr(l) + ', global_y = ' + repr(global_y))
                if (abs(global_y) < 0.05):
                    still_small = True
                if (slab_type == 'whole'):
                    tmp_y_flux, tmp_global_y, s_jf, s_hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, np.zeros(K), np.zeros(K), np.zeros(N-1), np.zeros(N-1), slab_type, D, M, N)
                if (slab_type == 'upper'):
                    tmp_y_flux, tmp_global_y, s_jf, s_hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, np.zeros(K), incbot, np.zeros(N-1), flux_bot, slab_type, D, M, N)
                if (slab_type == 'lower'):
                    tmp_y_flux, tmp_global_y, s_jf, s_hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, np.zeros(K), flux_top, np.zeros(N-1), slab_type, D, M, N)
                flux_out_top, flux_out_bot = compute_outgoing_seed_fluxes(energies, mu, dmu, s_jf, M, N)
                outgoing_fluxes_top.append(flux_out_top)
                outgoing_fluxes_bot.append(flux_out_bot)
                if (len(outgoing_fluxes_top) > 1):
                    e_lower_lim = 0.5e3
                    e_upper_lim = 30.0e3
                    for n in range(0, len(energies)):
                        if (energies[n] >= e_lower_lim):
                            l_ndx = n
                            break
                    for n in range(0, len(energies)):
                        if (energies[n] >= e_upper_lim):
                            u_ndx = n
                            break
                    max_diff_top = max(2*abs(outgoing_fluxes_top[-1][l_ndx:u_ndx] - outgoing_fluxes_top[-2][l_ndx:u_ndx])/(outgoing_fluxes_top[-1][l_ndx:u_ndx] + outgoing_fluxes_top[-2][l_ndx:u_ndx]))
                    max_diff_bot = max(2*abs(outgoing_fluxes_bot[-1][l_ndx:u_ndx] - outgoing_fluxes_bot[-2][l_ndx:u_ndx])/(outgoing_fluxes_bot[-1][l_ndx:u_ndx] + outgoing_fluxes_bot[-2][l_ndx:u_ndx]))
                    print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'outgoing fluxes max_diffs = ' + repr(max_diff_top) + ', ' + repr(max_diff_bot))
                    fluxes_converged = False
                    if ((slab_type == 'whole') and (max_diff_top < 0.05) and (max_diff_bot < 0.05)):
                        fluxes_converged = True
                    if ((slab_type == 'upper') and (max_diff_top < 0.05)):
                        fluxes_converged = True
                    if ((slab_type == 'lower') and (max_diff_bot < 0.05)):
                        fluxes_converged = True
                    if ((global_y < 0.05) and fluxes_converged):
                        if (not all_elements):
                            all_elements = True
                            l = 0
                            Jac_calls = 0
                            mid_output = (succ_code, T.copy(), xee, absorp.copy(), x_absorp.copy(), emis.copy(), x_emis.copy(), c_emis.copy(), jf.copy(), hf.copy(), s_jf.copy(), s_hf.copy(), y_flux.copy(), global_y)
                            succ_code = 1
                            print('(' + repr(phi_index) + ', ' + repr(r_index) + ', ' + slab_type + '):    ' + 'switching on the rest of the elements')
                            # instead of switching on the rest of the elements, just break here
                            finished = True
                            break
                        else:
                            finished = True
                            break

    if ((succ_code == 0) or (succ_code == -1)):
        return (succ_code, T, xee, absorp, x_absorp, emis, x_emis, c_emis, jf, hf, s_jf, s_hf, y_flux, global_y)
    else:
        return mid_output

def vector_error_func(T_np1, args):
    T, dt, jf_of_record, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, inctop, incbot, flux_top, flux_bot, slab_type, all_elements, Fe_abund, D, M, N = args

    print('in vector_error_func, T_np1 ' + repr(T_np1))

    temp = np.zeros(D-1)
    if (slab_type != 'whole'):
        Tbnd = T_np1[0]
        temp[1:-1] = T_np1[1:]
        if (slab_type == 'upper'):
            incbot = set_planck_I(energies, mu, Tbnd)
        if (slab_type == 'lower'):
            inctop = set_planck_I(energies, mu, Tbnd)
    else:
        temp[1:-1] = T_np1
    temp[0]  = temp[1]
    temp[-1] = temp[-2]

    xee, x_absorp, x_emis, c_emis = call_xstar(energies, dmu, dens, temp, jf_of_record, Fe_abund, all_elements, D, M, N)

    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
    emis   = set_brems_emis(energies, dens, xee, temp, D, N)

    y_flux, global_y, jf, hf = y_from_T_comp(T_np1, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp + x_absorp, emis + x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)

    error = np.zeros(len(T))
    for i in range(0, len(error)):
        if (slab_type == 'whole'):
            d = i + 1
        else:
            d = i
        if (d > 0):
            error[i] = (3.0/2.0) * (1.0 + xee[d]) * dens[d] * (1.380649e-16) * (z[d+1] - z[d]) * (T_np1[i] - T[i]) + y_flux[i] * dt

    print('in vector_error_func, error = ' + repr(error))

    return error

def expl_time_step(args):
    T, dt, jf_of_record, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, inctop, incbot, flux_top, flux_bot, slab_type, all_elements, Fe_abund, D, M, N = args

    print('in expl_time_step, T ' + repr(T))

    temp = np.zeros(D-1)
    if (slab_type != 'whole'):
        Tbnd = T[0]
        temp[1:-1] = T[1:]
        if (slab_type == 'upper'):
            incbot = set_planck_I(energies, mu, Tbnd)
        if (slab_type == 'lower'):
            inctop = set_planck_I(energies, mu, Tbnd)
    else:
        temp[1:-1] = T
    temp[0]  = temp[1]
    temp[-1] = temp[-2]

    xee, x_absorp, x_emis, c_emis = call_xstar(energies, dmu, dens, temp, jf_of_record, Fe_abund, all_elements, D, M, N)

    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
    emis   = set_brems_emis(energies, dens, xee, temp, D, N)

    y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp + x_absorp, emis + x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)

    T_np1 = np.zeros(len(T))
    for i in range(0, len(T_np1)):
        if (slab_type == 'whole'):
            d = i + 1
        else:
            d = i
        if (d > 0):
            cv = (3.0/2.0) * (1.0 + xee[d]) * dens[d] * (1.380649e-16) * (z[d+1] - z[d])
            T_np1[i] = ((-y_flux[i] * dt) + cv * T[i])/cv

    print('in expl_time_step, T_np1 ' + repr(T_np1))

    return T_np1

def tde_soln(z, zc, mu, dmu, energies, dens, xee, T, heat, mid_heat, inctop, incbot, flux_top, flux_bot, Fe_abund, slab_type, phi_index, r_index, D, M, N, prev_file_exists, prev_file, rank):
    K = M*(N-1)

    all_elements = False

    temp = np.zeros(D-1)
    if (slab_type != 'whole'):
        Tbnd = T[0]
        temp[1:-1] = T[1:]
        if (slab_type == 'upper'):
            incbot = set_planck_I(energies, mu, Tbnd)
        if (slab_type == 'lower'):
            inctop = set_planck_I(energies, mu, Tbnd)
    else:
        temp[1:-1] = T
    temp[0]  = temp[1]
    temp[-1] = temp[-2]

    scatt = np.zeros((D-1, N-1))
    for d in range(0, D-1):
        for n in range(0, N-1):
            scatt[d][n] = xee[d] * dens[d] * sigma(0.5*(energies[n+1] + energies[n]))
    absorp = set_brems_absorp(energies, dens, xee, temp, D, N)
    emis   = set_brems_emis(energies, dens, xee, temp, D, N)

    y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp, emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)

    xee, x_absorp, x_emis, c_emis = call_xstar(energies, dmu, dens, temp, jf, Fe_abund, all_elements, D, M, N)

    y_flux, global_y, jf, hf = y_from_T_comp(T, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, absorp + x_absorp, emis + x_emis, inctop, incbot, flux_top, flux_bot, slab_type, D, M, N)

    jf_of_record = jf.copy()

    time_step = 0
    dt = 0.01
    while True:
        args = [T, dt, jf_of_record, z, zc, mu, dmu, energies, dens, xee, heat, mid_heat, inctop, incbot, flux_top, flux_bot, slab_type, all_elements, Fe_abund, D, M, N]
#       T_np1 = root(vector_error_func, T, args = args, method = 'df-sane', jac = None)
        T_np1 = expl_time_step(args)
        print('after time step ' + repr(time_step) + ', T_np1/T = ' + repr(T_np1/T))
        T = T_np1.copy()
        jf_of_record = jf.copy()
        time_step += 1

# here is where XSTAR is actually called: this is our interface to it
def call_xstar(energies, dmu, dens, temp, jf, Fe_abund, all_elements, D, M, N):
    xee       = np.zeros(D-1)
    absorp    = np.zeros((D-1, N-1))
    emis      = np.zeros((D-1, N-1))
    c_emis    = np.zeros((D-1, N-1))

    for d in range(1, D-2):
        mi = np.zeros(N-1)
        for n in range(0, N-1):
            for m in range(0, M):
                mi[n] += 4.0 * np.pi * jf[d][n*M + m] * dmu[m] * 0.5 * (energies[n] + energies[n+1])

        if (not all_elements):
            tmp = xstarcomm(energies.copy(), mi.copy(), dens[d], 0.0, temp[d]/1.0e4, -Fe_abund, -99, N)
        else:
            tmp = xstarcomm(energies.copy(), mi.copy(), dens[d], 0.0, temp[d]/1.0e4, Fe_abund, -99, N)

        emis[d]   = tmp[0]
        c_emis[d] = tmp[1]
        absorp[d] = tmp[2]
        xee[d]    = tmp[4]

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5 * (energies[n+1] + energies[n])
    for d in range(1, D-2):
        f = interp1d(ec[1:], absorp[d][1:], kind = 'linear', fill_value = 'extrapolate')
        absorp[d][0] = f(ec[0])

    xee[0]  = xee[1]
    xee[-1] = xee[-2]

    absorp[0]  = absorp[1]
    absorp[-1] = absorp[-2]

    emis[0]  = emis[1]
    emis[-1] = emis[-2]

    c_emis[0]  = c_emis[1]
    c_emis[-1] = c_emis[-2]

    return (xee, absorp, emis, c_emis)

def compute_outgoing_seed_fluxes(energies, mu, dmu, jf, M, N):
    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5*(energies[n+1] + energies[n])

    flux_up = np.zeros(N-1)
    flux_dn = np.zeros(N-1)
    for n in range(0, N-1):
        for m in range(0, M):
            flux_up[n] += 4*np.pi * jf[0][n*M + m]  * mu[m] * dmu[m] * 1.602e-12 * ec[n]
            flux_dn[n] += 4*np.pi * jf[-1][n*M + m] * mu[m] * dmu[m] * 1.602e-12 * ec[n]

    return (flux_up, flux_dn)

def regrid(energies, jf, dens, Fe_abund, heat, temp, Tbnd, mu, dmu, new_N, D, M, N, slab_type):
    new_energies = np.logspace(np.log10(energies[0]), np.log10(energies[-1]), num = new_N, endpoint = True)

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5*(energies[n+1] + energies[n])

    new_ec = np.zeros(new_N-1)
    for n in range(0, new_N-1):
        new_ec[n] = 0.5*(new_energies[n+1] + new_energies[n])

    new_inctop = np.zeros(M*(new_N-1))
    new_incbot = np.zeros(M*(new_N-1))

    if (slab_type == 'upper'):
        new_incbot = set_planck_I(new_energies, mu, Tbnd)
    if (slab_type == 'lower'):
        new_inctop = set_planck_I(new_energies, mu, Tbnd)

    new_emis   = np.zeros((D-1, new_N-1))
    new_c_emis = np.zeros((D-1, new_N-1))
    new_absorp = np.zeros((D-1, new_N-1))
    new_xee    = np.zeros(D-1)

    for d in range(0, D-1):
        mi = np.zeros(N-1)
        for n in range(0, N-1):
            for m in range(0, M):
                mi[n] += 4.0 * np.pi * jf[d][n*M + m] * dmu[m] * 0.5 * (energies[n] + energies[n+1])
        f = interp1d(ec, mi, fill_value = 'extrapolate')
        new_mi = f(new_ec)

        tmp = xstarcomm(new_energies.copy(), new_mi.copy(), dens[d], heat[d], temp[d]/1.0e4, Fe_abund, -99, new_N)

        new_emis[d]   = tmp[0]
        new_c_emis[d] = tmp[1]
        new_absorp[d] = tmp[2]
        new_xee[d]    = tmp[4]

    # correction to the first energy bin
    for d in range(0, D-1):
        f = interp1d(new_ec[1:], new_absorp[d][1:], kind = 'linear', fill_value = 'extrapolate')
        new_absorp[d][0] = f(new_ec[0])

    new_emis   += set_brems_emis(new_energies, dens, new_xee, temp, D, new_N)
    new_c_emis += set_brems_emis(new_energies, dens, new_xee, temp, D, new_N)
    new_absorp += set_brems_absorp(new_energies, dens, new_xee, temp, D, new_N)

    return (new_energies, new_inctop, new_incbot, new_emis, new_c_emis, new_absorp, new_xee)

def write_output(outfile, succ_code, D, M, N, r, global_y, y_flux, z, zc, mu, dmu, energies, flux_top, flux_bot, dens, xee, heat, temp, Tbnd, absorp, x_absorp, scatt, emis, x_emis, c_emis, jf, hf, s_jf, s_hf, duration):
    with GetFile(outfile, open_type = 'w') as f:
        f.write('succ_code', succ_code)
        f.write('dims',      np.array([D, M, N], dtype = 'i'))
        f.write('r',         r)
        f.write('global_y',  global_y)
        f.write('y_flux',    y_flux)
        f.write('z',         z)
        f.write('zc',        zc)
        f.write('mu',        mu)
        f.write('dmu',       dmu)
        f.write('energies',  energies)
        f.write('flux_top',  flux_top)
        f.write('flux_bot',  flux_bot)
        f.write('dens',      dens)
        f.write('xee',       xee)
        f.write('heat',      heat)
        f.write('temp',      temp)
        f.write('Tbnd',      Tbnd)
        f.write('absorp',    absorp)
        f.write('x_absorp',  x_absorp)
        f.write('scatt',     scatt)
        f.write('emis',      emis)
        f.write('x_emis',    x_emis)
        f.write('c_emis',    c_emis)
        f.write('jf',        jf)
        f.write('hf',        hf)
        f.write('s_jf',      s_jf)
        f.write('s_hf',      s_hf)
        f.write('duration',  duration)

def set_tps_gaussian(energies, temp):
    kB  = 1.381e-16  # erg/K
    mc2 = 8.18712e-7 # ergs

    N = len(energies)
    D = len(temp) + 1

    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5 * (energies[n+1] + energies[n])

    tps = []

    for d in range(0, D-1):
        theta = (kB * temp[d])/mc2
        tps.append(0)
        start_point = len(tps)-1

        count = 0
        for n_in in range(0, N-1):
            e_in  = 0.5 * (energies[n_in+1] + energies[n_in])
            eps0  = e_in/(511.0e3)
            Ecomp = e_in*(1 + 4*theta - eps0)
            sigma = e_in*np.sqrt(2*theta + 0.4*(eps0**2))
            tmp   = (1.0/(sigma * np.sqrt(np.pi))) * np.exp(-(ec - Ecomp)**2/sigma**2)
            norm = 0.0
            for n_out in range(0, N-1):
                norm += tmp[n_out] * (energies[n_out+1] - energies[n_out])
            n_start = 0
            while True:
                if (tmp[n_start] * (energies[n_start+1] - energies[n_start]) > 0.5e-6):
                    break
                n_start += 1
            tmp2 = np.zeros(N-1)
            n = n_start
            while (n < N-1):
                if (tmp[n] * (energies[n+1] - energies[n]) > 0.5e-6):
                    tmp2[n] = tmp[n]
                    n += 1
                else:
                    break
            for n_out in range(0, N-1):
                if (tmp2[n_out] > 0.0):
                    tps.append(n_in)
                    tps.append(n_out)
                    tps.append((tmp2[n_out]/norm) * (energies[n_in+1] - energies[n_in]))
                    count += 1

        tps[start_point] = len(tps)

    return np.array(tps)

# --- --- --- code from ptx2pan.py --- --- ---

def get_r(flnm):
    with GetFile(flnm, open_type = 'r') as f:
        r = float(f.read('r'))
    return r

def fluxes_out(flnm, top_or_bot):
    with GetFile(flnm, open_type = 'r') as f:
        r        = float(f.read('r'))
        mu       = f.read('mu')
        dmu      = f.read('dmu')
        energies = f.read('energies')
        s_jf     = f.read('s_jf')
        s_hf     = f.read('s_hf')
    D = len(s_hf)
    M = len(mu)
    N = len(energies)
    ec = np.zeros(N-1)
    for n in range(0, N-1):
        ec[n] = 0.5*(energies[n+1] + energies[n])
    flux_out_top = np.zeros(N-1)
    flux_out_bot = np.zeros(N-1)
    for n in range(0, N-1):
        for m in range(0, M):
            flux_out_top[n] += 4*np.pi * s_jf[0][n*M + m] * mu[m] * dmu[m] * ec[n] * 1.602e-12
            flux_out_bot[n] += 4*np.pi * s_jf[-1][n*M + m] * mu[m] * dmu[m] * ec[n] * 1.602e-12
    ft = interp1d(ec, flux_out_top, fill_value = 'extrapolate')
    fb = interp1d(ec, flux_out_bot, fill_value = 'extrapolate')
    if (top_or_bot == 'top'):
        return ft(energies)
    else:
        return fb(energies)

    """
    else:
        new_energies = np.logspace(1, 5, 401, endpoint = True)
        if 'cemis' in flnm:
            ft = interp1d(ec, flux_out_top, fill_value = 'extrapolate')
            fb = interp1d(ec, flux_out_bot, fill_value = 'extrapolate')
            if (top_or_bot == 'top'):
                return ft(new_energies)
            else:
                return fb(new_energies)
        else:
            ft = interp1d(ec, flux_out_top, bounds_error = False, fill_value = 0.0)
            fb = interp1d(ec, flux_out_bot, bounds_error = False, fill_value = 0.0)
            if (top_or_bot == 'top'):
                return ft(new_energies)
            else:
                return fb(new_energies)
    """

def ptx2pan(input_filename, dir, prefix, out_flnm):
    for top_or_bot in ['top', 'bot']:
        phi_index_list = []
        r_index_list   = []
        phi_list = []
        r_list   = []
        for flnm in os.listdir(dir):
            if (flnm[:len(prefix)] == prefix):
                phi_index = int(flnm.split('_')[2])
                if phi_index not in phi_index_list:
                    phi_index_list.append(phi_index)
                    phi_list.append((phi_index/64.0) * 0.5 * np.pi)
                r_index = int(flnm.split('_')[3].split('.')[0])
                if r_index not in r_index_list:
                    r_index_list.append(r_index)
                    r_list.append(get_r('./' + dir + '/' + flnm))
        phi_index_list.sort()
        r_index_list.sort()
        phi_list.sort()
        r_list.sort()
        Nphi = len(phi_index_list)
        Nr   = len(r_index_list)
        Nnu  = 0
        for i in range(0, Nphi):
            for j in range(0, Nr):
                flnm = prefix + '_whole_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
                if flnm in os.listdir(dir):
                    tmp = fluxes_out('./' + dir + '/' + flnm, top_or_bot)
                    if (Nnu == 0):
                        Nnu = len(tmp)
                        e_spec_out = np.zeros((Nphi, Nr, Nnu))
                    e_spec_out[i][j] = tmp
                    continue
                flnm = prefix + '_upper_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
                if ((top_or_bot == 'top') and flnm in os.listdir(dir)):
                    tmp = fluxes_out('./' + dir + '/' + flnm, top_or_bot)
                    if (Nnu == 0):
                        Nnu = len(tmp)
                        e_spec_out = np.zeros((Nphi, Nr, Nnu))
                    e_spec_out[i][j] = tmp
                    continue
                flnm = prefix + '_lower_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
                if ((top_or_bot == 'bot') and flnm in os.listdir(dir)):
                    tmp = fluxes_out('./' + dir + '/' + flnm, top_or_bot)
                    if (Nnu == 0):
                        Nnu = len(tmp)
                        e_spec_out = np.zeros((Nphi, Nr, Nnu))
                    e_spec_out[i][j] = tmp
                    continue
        # cull any negative values or NaNs or infs (just in case)
        for i in range(0, Nphi):
            for j in range(0, Nr):
                for k in range(0, Nnu):
                    if ((e_spec_out[i][j][k] < 0.0) or np.isnan(e_spec_out[i][j][k]) or np.isinf(e_spec_out[i][j][k])):
                        e_spec_out[i][j][k] = 0.0

        # --- --- ---

        # read in grid file to match
        with GetFile('./' + input_filename, open_type = 'r') as f:
            m_r   = f.read('r')
            m_phi = f.read('phi')

        m_Nr   = len(m_r)
        m_Nphi = len(m_phi)

        m_e_spec_out = np.zeros((m_Nphi, m_Nr, Nnu))

        for i in range(0, m_Nphi):
            nearest_i = 0
            for k in range(1, Nphi):
                if (abs(phi_list[k] - m_phi[i]) < abs(phi_list[nearest_i] - m_phi[i])):
                    nearest_i = k
            for j in range(0, m_Nr):
                nearest_j = 0
                # if trying to match to cells within the very thin (not computed by PTRANSX) region, just move on
                if (m_r[j] < r_list[0]):
                    continue
                for k in range(1, Nr):
                    if (abs(r_list[k] - m_r[j]) < abs(r_list[nearest_j] - m_r[j])):
                        nearest_j = k
                m_e_spec_out[i][j] = e_spec_out[nearest_i][nearest_j]

        # --- --- ---

        # convert from erg/s/cm^2/eV to erg/s/cm^2/Hz/ster (assuming half-space isotropy)
        m_e_spec_out *= 4.135668e-15/np.pi

        # reformat for pandurata
        trans_spec = np.zeros((Nnu, m_Nphi, m_Nr))

        for i in range(0, Nnu):
            for j in range(0, m_Nphi):
                for k in range(0, m_Nr):
                    trans_spec[i][j][k] = m_e_spec_out[j][k][i]

        if (top_or_bot == 'top'):
            open_type = 'w'
        else:
            open_type = 'r+'
        with GetFile(out_flnm, open_type = open_type) as f:
            f.write('seed_spec_' + top_or_bot, trans_spec)

#       for i in range(0, m_Nr):
#           for j in range(0, m_Nphi):
#               for k in range(0, Nnu):
#                   e_outfile.write('{0:.6E}'.format(m_e_spec_out[j][i][k]) + '\n')

# --- from rsp2pan ---

def collect_data(flnm):
    with GetFile(flnm, open_type = 'r') as f:
        e_coarse = f.read('e_coarse')
        e_fine   = f.read('e_fine')
        eta_grid = f.read('eta')

        slab_type = flnm.split('/')[-1].split('_')[1]

        if (slab_type == 'whole'):
            r_frac_top     = f.read('r_frac_top')
            r_frac_bot     = f.read('r_frac_bot')
            refl_profs_top = f.read('refl_profs_top')
            refl_profs_bot = f.read('refl_profs_bot')
        else:
            r_frac_top     = f.read('r_frac_top')
            refl_profs_top = f.read('refl_profs_top')

    if (slab_type != 'whole'):
        with GetFile(flnm.replace('upper', 'lower'), open_type = 'r') as f:
            r_frac_bot     = f.read('r_frac_bot')
            refl_profs_bot = f.read('refl_profs_bot')
    
    Ne = len(e_fine)-1
    A  = 10.0**((1.0/Ne)*np.log10(e_fine[-1]/e_fine[0]))
    new_eta_grid = np.zeros(2*Ne+1)
    for i in range(0, 2*Ne+1):
        new_eta_grid[i] = A**(i-Ne)
    bnd_eta_grid = np.zeros(2*Ne+2)
    for i in range(1, 2*Ne+1):
        bnd_eta_grid[i] = 0.5*(new_eta_grid[i] + new_eta_grid[i-1])
    bnd_eta_grid[0]  = new_eta_grid[0] - 0.5*(new_eta_grid[1] - new_eta_grid[0])
    bnd_eta_grid[-1] = new_eta_grid[-1] + 0.5*(new_eta_grid[-1] - new_eta_grid[-2])

    new_refl_profs_top = np.zeros((len(e_coarse), 2*Ne+2))
    new_refl_profs_bot = np.zeros((len(e_coarse), 2*Ne+2))

    for i in range(0, len(e_coarse)):
        f = interp1d(eta_grid, refl_profs_top[i], bounds_error = False, fill_value = (0.0, 1.0))
        new_refl_profs_top[i] = f(bnd_eta_grid)
        f = interp1d(eta_grid, refl_profs_bot[i], bounds_error = False, fill_value = (0.0, 1.0))
        new_refl_profs_bot[i] = f(bnd_eta_grid)

    return r_frac_top, r_frac_bot, new_refl_profs_top, new_refl_profs_bot

def rsp2pan(dir_name):
    dir = dir_name

    phi_index_list = []
    r_index_list   = []

    phi_list = []
    r_list   = []

    for flnm in os.listdir(dir):
        if (flnm[:2] == 'mc'):
            phi_index = int(flnm.split('_')[2])
            if phi_index not in phi_index_list:
                phi_index_list.append(phi_index)
                phi_list.append((phi_index/64.0) * 0.5 * np.pi)
            r_index = int(flnm.split('_')[3].split('.')[0])
            if r_index not in r_index_list:
                r_index_list.append(r_index)
                r_list.append(get_r('./' + dir + '/' + flnm))

    phi_index_list.sort()
    r_index_list.sort()

    phi_list.sort()
    r_list.sort()

    Nphi = len(phi_index_list)
    Nr   = len(r_index_list)

    # get energy grids, presumably the same across all files

    for flnm in os.listdir(dir):
        if (flnm[:2] == 'mc'):
            with GetFile('./' + dir + '/' + flnm) as f:
                e_coarse = f.read('e_coarse')
                e_fine   = f.read('e_fine')
                break

    # not all pairs of (phi_index, r_index) will have slabs;
    # record existence in slab_exists: 1 indicates existence

    slab_exists = np.zeros((Nphi, Nr), dtype = 'i')

    for i in range(0, Nphi):
        for j in range(0, Nr):
            # an extant slab will either have a "whole" file or an "upper" *and* "lower"
            flnm = 'mc_whole_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            if flnm in os.listdir(dir):
                slab_exists[i][j] = 1
                r_frac_top, r_frac_bot, refl_profs_top, refl_profs_bot = collect_data('./' + dir + '/' + flnm)
            flnm1 = 'mc_lower_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            flnm2 = 'mc_upper_' + repr(phi_index_list[i]) + '_' + repr(r_index_list[j]) + '.h5'
            if (flnm1 in os.listdir(dir)) and (flnm2 in os.listdir(dir)):
                slab_exists[i][j] = 1
                r_frac_top, r_frac_bot, refl_profs_top, refl_profs_bot = collect_data('./' + dir + '/' + flnm2)
            # only make disk-covering data structures the one time
            if (slab_exists.sum() == 1):
                r_frac_top_disk     = np.zeros((Nphi, Nr, len(e_fine)))
                r_frac_bot_disk     = np.zeros((Nphi, Nr, len(e_fine)))
                refl_profs_top_disk = np.zeros((Nphi, Nr, len(e_coarse), 2*(len(e_fine)-1)+2))
                refl_profs_bot_disk = np.zeros((Nphi, Nr, len(e_coarse), 2*(len(e_fine)-1)+2))
            if (slab_exists[i][j] > 0):
                r_frac_top_disk[i][j]     = r_frac_top
                r_frac_bot_disk[i][j]     = r_frac_bot
                refl_profs_top_disk[i][j] = refl_profs_top
                refl_profs_bot_disk[i][j] = refl_profs_bot

    print(Nphi, Nr, slab_exists.sum())

    out_flnm = dir.split('/')[-1].replace('run', 'greens') + '.h5'

    with GetFile(out_flnm, open_type = 'w') as f:        
        f.write('Nphi',                Nphi)
        f.write('Nr',                  Nr)
        f.write('slab_exists',         slab_exists)
        f.write('phi_list',            np.array(phi_list))
        f.write('r_list',              np.array(r_list))
        f.write('e_coarse',            e_coarse)
        f.write('e_fine',              e_fine)
        f.write('r_frac_top_disk',     r_frac_top_disk)
        f.write('r_frac_bot_disk',     r_frac_bot_disk)
        f.write('refl_profs_top_disk', refl_profs_top_disk)
        f.write('refl_profs_bot_disk', refl_profs_bot_disk)
