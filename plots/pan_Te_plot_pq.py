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

a = 0.0

flnm = '../panptx_data/RADFLUX_panptx_1000.h5'

panptx_data = h5py.File(flnm, 'r')

r   = np.array(panptx_data['/r'])
th  = np.array(panptx_data['/th'])
phi = np.array(panptx_data['/phi'])

panptx_data.close()

max_r = 51.0
for i in range(0, len(r)):
    if (r[i] > max_r):
        break
max_r_index = i
r = r[:max_r_index]

pan_data = h5py.File('../panptx_data/te.1000.0.h5', 'r')

T_e = 8.617e-8 * np.array(pan_data['T_e'])

slice_data = np.mean(T_e[:max_r_index,::,::], axis = 2)

slice_data = np.ma.masked_where(slice_data == 0.0 , slice_data)

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

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = 0.0, vmax = 3.0)
plt.colorbar(label = r'$\log T_e \left[\mathrm{keV}\right]$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('pan_Te.png', bbox_inches = 'tight', dpi = 600)

plt.show()
