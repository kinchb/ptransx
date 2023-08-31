import numpy as np

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

dir = '../harm_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

harm_data_1.close()

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

harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

rho = np.array(harm_data_2['rho'])

# rho = np.mean(rho[:max_r_index,::,::], axis = 2)
# rho = rho[:max_r_index,::,32]

slice_data_new = rho.copy()

dir = '../harm_data/old_data'

harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.010000.h5', 'r')

rho = np.array(harm_data_2['rho'])

# rho = np.mean(rho[:max_r_index,::,::], axis = 2)
# rho = rho[:max_r_index,::,32]

slice_data_old = rho.copy()

slice_data = np.mean((slice_data_new[:max_r_index,::,::] - slice_data_old[:max_r_index,::,::])/slice_data_old[:max_r_index,::,::], axis = 2)

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

r_abbr = np.load('photo_locs.npz')['r_abbr']
emtop  = np.load('photo_locs.npz')['emtop']
embot  = np.load('photo_locs.npz')['embot']
plt.plot(emtop, r_abbr, 'k-')
plt.plot(embot, r_abbr, 'k-')

plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap, vmin = -1.0, vmax = 1.0)
plt.colorbar(label = r'$( \rho_{1000M} - \rho_{-1M} )/\rho_{-1M}$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_rho_comp.png', bbox_inches = 'tight', dpi = 600)

plt.show()
