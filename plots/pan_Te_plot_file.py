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

flnm = '../data/h3d_10_a09_10_gcc/radfluxes/RADFLUX_panptx_1000.h5'

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

def pad(r, th, phi, data):
    new_r    = np.zeros(len(r)+1)
    new_data = np.zeros((len(r)+1, len(th), len(phi)))

    for i in range(1, len(r)+1):
        new_r[i] = r[i-1]
        for j in range(0, len(th)):
            for k in range(0, len(phi)):
                new_data[i][j][k] = data[i-1][j][k]

    return (new_r, new_data)

"""
slices = []
for i in range(20, 1001, 20):
    flnm = '../data/h3d_10_a09_10_gcc/tes/te.' + repr(i) + '.0.h5'
    pan_data = h5py.File(flnm, 'r')
    T_e = (8.617e-8 * np.array(pan_data['T_e']))/511.0
    slice_data = np.ma.masked_where(T_e[:max_r_index,::,::] == 0.0 , T_e[:max_r_index,::,::])
    slices.append(slice_data)
    pan_data.close()
slices = np.array(slices)
"""

flnm = '../data/h3d_10_a09_10_gcc/tes/te.1000.0.h5'
pan_data = h5py.File(flnm, 'r')
T_e = (8.617e-8 * np.array(pan_data['T_e']))/511.0
pan_data.close()

flnm = '../data/h3d_10_a09_10_gcc/radfluxes/RADFLUX_panptx_1000.h5'
radflux_data = h5py.File(flnm, 'r')
cool = np.array(radflux_data['coolfunc'])
slice_data = np.ma.masked_where(T_e[:max_r_index,::,::] == 0.0 , T_e[:max_r_index,::,::])
slice_data = np.ma.masked_where(cool[:max_r_index,::,::] == 0.0, slice_data)

# slice_data = np.mean(slices, axis = 0)

r, slice_data = pad(r, th, phi, slice_data)

# slice_data = np.mean(slice_data, axis = 2)
slice_data = slice_data[::,::,32]

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

plt.title(r'\textsc{pandurata} (single $\phi$ slice)')

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -3.0, vmax = 3.0)
plt.colorbar(label = r'$\log \Theta_e$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

f = plt.gcf()
f.savefig('pan_Te_slice.png', bbox_inches = 'tight', dpi = 600)

plt.show()
