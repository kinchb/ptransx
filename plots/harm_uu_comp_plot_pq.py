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

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.0
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

def calc_incor(rho, r, th, phi):
    rho *= (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

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

    emtop_ik    = (np.pi/2) * np.ones((len(r), len(phi)))
    embot_ik    = (np.pi/2) * np.ones((len(r), len(phi)))
    diskbody_ik = np.zeros((len(r), len(phi)), dtype = int)

    for i in range(0, len(r)-1):
        for k in range(0, len(phi)):
            if ((tau_from_top[i,::,k].max() > 1.0) and (tau_from_bot[i,::,k].max() > 1.0)):
                emtop_ik[i][k] = interp1d(tau_from_top[i,::,k], th_b)(1.0)
                embot_ik[i][k] = interp1d(tau_from_bot[i,::,k][::-1], th_b[::-1])(1.0)
                if (emtop_ik[i][k] > embot_ik[i][k]):
                    emtop_ik[i][k] = np.pi/2
                    embot_ik[i][k] = np.pi/2
                else:
                    diskbody_ik[i][k] = 2
    for k in range(0, len(phi)):
        emtop_ik[len(r)-1][k]    = emtop_ik[len(r)-2][k]
        embot_ik[len(r)-1][k]    = embot_ik[len(r)-2][k]
        diskbody_ik[len(r)-1][k] = 2

    incor = np.ones((len(r), len(th), len(phi)), dtype = int)

    for i in range(0, len(r)):
        for j in range(0, len(th)):
            for k in range(0, len(phi)):
                if (diskbody_ik[i][k] != 0):
                    if ((th[j] > emtop_ik[i][k]) and (th[j] < embot_ik[i][k])):
                        incor[i][j][k] = 0

    return incor

dir = '../harm_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

for i in range(0, len(x1)):
    r[i]   = x1[i][0][0]
for i in range(0, len(x2[0])):
    th[i]  = x2[0][i][0]
for i in range(0, len(x3[0][0])):
    phi[i] = x3[0][0][i]

max_r = 51.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i
r = r[:max_r_index]

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')

# harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

# urad = np.array(harm_data_2['urad'])

uu = np.array(harm_data_1['uu'])

rho = np.array(harm_data_1['rho'])

# uu = np.ma.masked_where(urad == 0.0, uu)

# incor = calc_incor(rho, r, th, phi)

# uu = uu * incor

# uu = np.mean(uu[:max_r_index,::,::], axis = 2)
# rho = rho[:max_r_index,::,32]

slice_data_new = uu.copy()

dir = '../harm_data/old_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.000500.h5', 'r')

# harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.010000.h5', 'r')

# urad = np.array(harm_data_2['urad'])

uu = np.array(harm_data_1['uu'])

# uu = np.ma.masked_where(urad == 0.0, uu)

rho = np.array(harm_data_1['rho'])

# uu = np.ma.masked_where(urad == 0.0, uu)

# incor = calc_incor(rho, r, th, phi)

# uu = uu * incor

slice_data_old = uu.copy()

dV_cgs = np.load('../harm_data/dV_cgs.npz')['dV_cgs']

uconv = c**2 * (4 * np.pi * c**2 * mdot)/(kappa * G * m_bh * 3.0e-4 * eta)

change_in_utot = uconv * ((slice_data_new * dV_cgs).sum() - (slice_data_old * dV_cgs).sum())

net_heat_rate = change_in_utot/(1000.0 * 4.9e-5)

# slice_data = np.mean((slice_data_new[:max_r_index,::,::] - slice_data_old[:max_r_index,::,::])/slice_data_old[:max_r_index,::,::], axis = 2)
slice_data = np.mean(slice_data_new[:max_r_index,::,::]/slice_data_old[:max_r_index,::,::], axis = 2)

rad, theta = np.meshgrid(r, th)

fig = plt.figure(figsize = (4.5, 8))
ax  = Axes3D(fig)

rad, theta = np.meshgrid(r, th)

ax = plt.subplot(projection = 'polar')
ax.set_theta_zero_location('N')
ax.set_theta_direction(-1)
ax.set_thetamin(0)
ax.set_thetamax(180)
ax.set_rlabel_position(80)
ax.set_xticks(np.pi/180.0 * np.array([0, 30, 60, 90, 120, 150, 180]))
ax.set_yticks(np.array([0, 10, 20, 30, 40, 50]))

cmap = plt.get_cmap('RdBu')

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -1.5, vmax = 1.5)
plt.colorbar(label = r'$(p_{\mathrm{gas},1000M} - p_{\mathrm{gas},-1M})/p_{\mathrm{gas},-1M}$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_p_gas_comp.png', bbox_inches = 'tight', dpi = 600)

plt.show()
