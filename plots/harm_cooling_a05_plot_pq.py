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
a    = 0.5
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

dir = '../harm_data/a05_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.002581.h5', 'r')

print(repr(np.array(harm_data_1['/Header/Grid/t'])))

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

harm_data_1 = h5py.File(dir + '/KDHARM0.002581.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.012081.h5', 'r')

cool = np.array(harm_data_2['coolfunc'])

slice_data = np.mean(cool[:max_r_index,::,::], axis = 2)
# slice_data = cool[:max_r_index,::,32]

slice_data = np.ma.masked_where(slice_data == 0.0, slice_data)

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

cmap = plt.get_cmap('jet')

lum_conv = (4 * np.pi * c**7 * mdot)/(kappa * G**2 * m_bh**2 * 3.0e-4 * eta)

slice_data *= lum_conv

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = 10.0, vmax = 15.0)
plt.colorbar(label = r'$\log \mathcal{L} \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-3} \right]$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

# f = plt.gcf()
# f.savefig('harm_cool_a05.png', bbox_inches = 'tight', dpi = 600)

plt.show()
