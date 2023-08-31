import os

import numpy as np

from scipy.interpolate import interp1d

import h5py

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc

import imageio

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 15)

sys.path.append('/Users/kinch/spin_analysis/util_scripts')
import harm_utils

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.99
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
m_i    = 1.67262171e-24
m_e    = 9.1093826e-28
# --- --- --- --- --- --- --- --- --- ---

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

dir = '../data/h3d_10_a099_01_gcc_later/dumps'

harm_data_1 = h5py.File(dir + '/KDHARM0.000520.h5', 'r')
harm_data_2 = h5py.File(dir + '/../KDHARM0.RADFLUX.016800.h5', 'r')
gdump_fname = dir + '/../KDHARM0.gdump.a09.h5'

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

harm_data_1.close()
harm_data_1 = h5py.File(dir + '/../KDHARM0.000840.h5', 'r')

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

for i in range(0, len(x1)):
    r[i]   = x1[i][0][0]
for i in range(0, len(x2[0])):
    th[i]  = x2[0][i][0]
for i in range(0, len(x3[0][0])):
    phi[i] = x3[0][0][i]

max_r = 11.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i
r = r[:max_r_index]

def pad(r, th, data):
    new_r    = np.zeros(len(r)+1)
    new_data = np.zeros((len(r)+1, len(th)))

    for i in range(1, len(r)+1):
        new_r[i] = r[i-1]
        for j in range(0, len(th)):
            new_data[i][j] = data[i-1][j]

    return (new_r, new_data)

def get_mean_photo_thetas(rho):
    if os.path.isfile(dir + '/photo_locs.npz'):
        r_abbr = np.load(dir + '/photo_locs.npz')['r_abbr']
        emtop  = np.load(dir + '/photo_locs.npz')['emtop']
        embot  = np.load(dir + '/photo_locs.npz')['embot']
        if (max_r < r_abbr[-1]):
            for i in range(0, len(r_abbr)):
                if (r_abbr[i] > max_r):
                    break
            max_r_index = i
            r_abbr = r_abbr[:max_r_index]
            emtop  = emtop[:max_r_index]
            embot  = embot[:max_r_index]
        return (r_abbr, emtop, embot)

    # convert rho and coolfunc to cgs units
    rho *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

    # determine tau = 1 surface location and related quantities
    th_b = np.zeros(len(th)+1)
    for i in range(1, len(th_b)-1):
        th_b[i] = 0.5*(th[i-1] + th[i])
    th_b[0]  = 0.0
    th_b[-1] = np.pi

    tau_from_top = np.zeros((len(r), len(th_b), len(phi)))
    tau_from_bot = np.zeros((len(r), len(th_b), len(phi)))

    for i in range(0, len(r)):
        for j in range(1, len(th_b)):
            for k in range(0, len(phi)):
                tau_from_top[i][j][k] = tau_from_top[i][j-1][k] + kappa * rho[i][j-1][k] * ((G*m_bh)/(c*c)) * np.sqrt(r[i]**2 + a**2 * np.cos(th[j-1])**2) * (th_b[j] - th_b[j-1])

    for i in range(0, len(r)):
        for j in range(len(th_b)-2, -1, -1):
            for k in range(0, len(phi)):
                tau_from_bot[i][j][k] = tau_from_bot[i][j+1][k] + kappa * rho[i][j][k] * ((G*m_bh)/(c*c)) * np.sqrt(r[i]**2 + a**2 * np.cos(th[j])**2) * (th_b[j+1] - th[j])

    emtop_ik = (np.pi/2) * np.ones((len(r), len(phi)))
    embot_ik = (np.pi/2) * np.ones((len(r), len(phi)))

    for i in range(0, len(r)-1):
        for k in range(0, len(phi)):
            if ((tau_from_top[i,::,k].max() > 1.0) and (tau_from_bot[i,::,k].max() > 1.0)):
                emtop_ik[i][k] = interp1d(tau_from_top[i,::,k], th_b)(1.0)
                embot_ik[i][k] = interp1d(tau_from_bot[i,::,k][::-1], th_b[::-1])(1.0)
                if (emtop_ik[i][k] > embot_ik[i][k]):
                    emtop_ik[i][k] = np.pi/2
                    embot_ik[i][k] = np.pi/2
    for k in range(0, len(phi)):
        emtop_ik[len(r)-1][k] = emtop_ik[len(r)-2][k]
        embot_ik[len(r)-1][k] = embot_ik[len(r)-2][k]

    emtop = np.mean(emtop_ik, axis = 1)
    embot = np.mean(embot_ik, axis = 1)

    r_abbr = r.copy()

    start_ndx = 0
    while ((emtop[start_ndx] == np.pi/2) or (embot[start_ndx] == np.pi/2)):
        start_ndx += 1

    np.savez(dir + '/photo_locs.npz', r_abbr = r_abbr[start_ndx:], emtop = emtop[start_ndx:], embot = embot[start_ndx:])

    return (r_abbr[start_ndx:], emtop[start_ndx:], embot[start_ndx:])

rho = np.array(harm_data_2['rho'])

rho = rho[:max_r_index,::,::]

r_abbr, emtop, embot = get_mean_photo_thetas(rho.copy())

# --- --- ---
"""
rho        *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)
rho         = np.mean(rho, axis = 2)
slice_data  = rho
vmin  = -10.0
vmax  = -5.0
title = r'\textsc{harm3d}'
label = r'$\log \rho \left[\mathrm{g}\ \mathrm{cm}^{-3} \right]$'
flnm  = 'harm_rho.png'
"""
# --- --- ---
"""
bsq         = np.array(harm_data_1['bsq'])
bsq         = np.mean(bsq, axis = 2)
slice_data  = bsq
vmin  = -10.0
vmax  = -5.0
title = r'\textsc{harm3d}'
label = r'$\log B^2$'
flnm  = 'harm_bsq.png'
"""
# --- --- ---
bsq         = np.array(harm_data_1['bsq'])
rho         = np.array(harm_data_1['rho'])
uu          = np.array(harm_data_1['uu'])
bsq_o_rho   = np.mean(bsq/rho, axis = 2)
bsq_o_uu    = np.mean(bsq/uu,  axis = 2)
slice_data  = bsq_o_rho
vmin  = -4.0
vmax  = 6.0
title = r'\textsc{harm3d}: $a = 0.99$'
label = r'$\log B^2/\rho$'
flnm = 'harm_bsq_o_rho_a99_zoomed.png'
# flnm   = 'harm_bsq_o_uu_a99.png'
# --- --- ---
"""
cool        = np.array(harm_data_2['coolfunc'])
r, cool     = pad(r, th, phi, cool)
cool       *= (4 * np.pi * c**7 * mdot)/(kappa * G**2 * m_bh**2 * 3.0e-4 * eta)
# cool        = np.ma.masked_where(cool == 0.0, cool)
cool        = np.mean(cool, axis = 2)
slice_data  = cool
vmin  = 12.0
vmax  = 16.0
title = r'\textsc{harm3d}'
label = r'$\log \mathcal{L} \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-3} \right]$'
flnm  = 'harm_cool_a99_uc.png'
"""
# --- --- ---
"""
cool        = np.array(harm_data_2['coolfunc'])
uu          = np.array(harm_data_1['uu'])
rho         = np.array(harm_data_2['rho'])
T_e         = (1.0/(1.0 + 1.2)) * ((5.0/3.0) - 1.0) * (m_i/m_e) * (uu/rho)
T_e         = np.ma.masked_where(cool == 0.0, T_e)
T_e         = np.mean(T_e, axis = 2)
slice_data  = T_e
vmin  = -2.0
vmax  = 2.0
title = r'\textsc{harm3d}: $a = 0.99, \dot{m} = 0.01, t = -1M$'
label = r'$\log \Theta_e$'
flnm  = 'harm_Te_a99_t0.png'
"""
# --- --- ---
"""
urad        = np.array(harm_data_2['urad'])
urad       *= (4 * np.pi * c**4 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)
urad        = np.mean(urad, axis = 2)
slice_data  = urad
vmin  = 9.0
vmax  = 13.0
title = r'\textsc{harm3d}'
label = r'$\log u_\mathrm{rad} \left[\mathrm{erg}\ \mathrm{cm}^{-3} \right]$'
flnm  = 'harm_urad.png'
"""
# --- --- ---
"""
gamma          = np.array(harm_data_1['gamma'])
gamma          = np.mean(gamma, axis = 2)
slice_data     = gamma - 1.0
vmin  = -3.0
vmax  = 1.0
title = r'\textsc{harm3d}: $a = 0.99, \dot{m} = 0.01$'
label = r'$\log \left(\gamma - 1\right)$'
flnm  = 'harm_gamma_a99_nocool_unbound.png'
"""
# --- --- ---
"""
metric         = harm_utils.get_metric(gdump_fname)
u0, u1, u2, u3 = harm_utils.calc_ucon_bl(harm_data_1, metric)
beta           = np.mean(u1/u0, axis = 2)
gamma_pos      = np.ma.masked_where(beta < 0.0, 1.0/np.sqrt(1.0 - beta**2))
slice_data     = gamma_pos - 1.0
vmin  = -3.0
vmax  = 1.0
title = r'\textsc{harm3d}'
label = r'$\log \left(\gamma_r - 1\right)$'
flnm  = 'harm_ucon1.png'
"""
# --- --- ---

r, slice_data = pad(r, th, slice_data)

rad, theta = np.meshgrid(r, th)

fig = plt.figure(figsize = (4.5, 8))
ax  = Axes3D(fig)

rad, theta = np.meshgrid(r, th)

ax = plt.subplot(projection = 'polar')
ax.set_rmax(50.0)
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rlabel_position(80)
ax.set_xticks(np.pi/180.0 * np.array([0, 30, 60, 90, 120, 150, 180]))
ax.set_yticks(np.array([0, 10, 20, 30, 40, 50]))

cmap = plt.get_cmap('jet')
# cmap = plt.get_cmap('RdBu')

plt.title(title)
plt.plot(emtop, r_abbr, 'w-')
plt.plot(embot, r_abbr, 'w-')

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = vmin, vmax = vmax)
# plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap)#, vmin = vmin, vmax = vmax)

plt.colorbar(label = label, fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform = ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()

f.savefig(flnm, bbox_inches = 'tight', dpi = 600)

plt.show()
