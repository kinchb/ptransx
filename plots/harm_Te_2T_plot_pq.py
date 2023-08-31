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

m_e = 9.1093826e-28
m_i = 1.67262171e-24

dir = '../harm_data/2T_data'
a   = 0.0

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

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

T_e  = np.array(harm_data_2['Te_eV'])/1000.0

rho  = np.array(harm_data_1['rho'])
uu   = np.array(harm_data_1['uu'])

T_i = 938.3e3 * ((2.0/3.0) * (uu/rho) - 1.21 * (m_e/m_i) * (T_e/511.0))

slice_data = T_e[:max_r_index,::,::]

slice_data = np.ma.masked_where(np.array(harm_data_2['urad'])[:max_r_index,::,::] == 0.0, slice_data)

slice_data = np.mean(slice_data, axis = 2)

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

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -1.0, vmax = 2.0)
plt.colorbar(label = r'$\log T_e \left[\mathrm{keV}\right]$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('harm_Te_2T.png', bbox_inches = 'tight', dpi = 600)

plt.show()
