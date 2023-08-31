import numpy as np
import numpy.ma as ma

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

dir = '/mnt/archive/h3d_10_a0_01_1o5000_beta/dumps'
a   = 0.0

harm_data_1 = h5py.File(dir + '/KDHARM0.000700.h5', 'r')

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

harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.010200.h5', 'r')

rho = np.array(harm_data_2['rho'])

# rho = np.mean(rho[:max_r_index,::,::], axis = 2)
rho = rho[:max_r_index,::,::]

dlnrho = np.zeros((len(rho), len(rho[0]), len(rho[0][0])))

for i in range(0, max_r_index):
    for k in range(0, len(phi)):
        for j in range(1, len(th)-1):
            dlnrho[i, j, k] = abs((np.log(rho[i, j-1, k]) - np.log(rho[i, j+1, k]))/2)
        dlnrho[i, 0, k]  = abs(np.log(rho[i, 1, k]) - np.log(rho[i, 0, k]))
        dlnrho[i, -1, k] = abs(np.log(rho[i, -1, k]) - np.log(rho[i, -2, k]))

coolfunc_corona = np.array(harm_data_2['coolfunc_corona'])[:max_r_index,::,::]

# dlnrho = ma.masked_where(coolfunc_corona == 0.0, dlnrho)

print dlnrho[::,::,32].mean()

# dlnrho = np.mean(dlnrho, axis = 2)
dlnrho = dlnrho[::,::,32]

slice_data = dlnrho

# --- --- ---

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

plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap, vmin = 0.0, vmax = 2.0)
# plt.pcolormesh(theta, rad, slice_data.T, cmap = cmap, vmin = -4.0, vmax = 5.0)
plt.colorbar(label = r'$| \Delta \ln \rho |$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.title(r'$\phi = \pi/4, t = 200 M$')

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_Dlnrho_t200.png', bbox_inches = 'tight', dpi = 600)

plt.show()

