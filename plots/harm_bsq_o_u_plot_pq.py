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

dir = '../harm_data'
a   = 0.0

harm_data_1 = h5py.File(dir + '/KDHARM0.000500.h5', 'r')

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

dir = '../harm_data/old_data'

harm_data_1 = h5py.File(dir + '/KDHARM0.000500.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.010000.h5', 'r')

# T_e  = np.array(harm_data_2['Te_eV'])/1000.0

bsq  = np.array(harm_data_1['bsq'])
rho  = np.array(harm_data_1['rho'])
uu   = np.array(harm_data_1['uu'])

bsq_o_rho = bsq/rho
bsq_o_uu  = bsq/uu

# slice_data = np.mean(bsq_o_rho[:max_r_index,::,::], axis = 2)
slice_data = np.mean(bsq_o_uu[:max_r_index,::,::], axis = 2)
# slice_data = bsq_o_rho[:max_r_index,::,32]
# slice_data = bsq_o_uu[:max_r_index,::,32]

# cool = np.mean(cool[:max_r_index,::,::], axis = 2)
# slice_data = np.ma.masked_where(cool == 0.0, slice_data)

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

r_abbr = np.load('photo_locs_old.npz')['r_abbr']
emtop  = np.load('photo_locs_old.npz')['emtop']
embot  = np.load('photo_locs_old.npz')['embot']
plt.plot(emtop, r_abbr, 'w-')
plt.plot(embot, r_abbr, 'w-')

plt.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -1.0, vmax = 5.0)
plt.colorbar(label = r'$\log \left( B^2/ u \right)$', fraction = 0.046, pad = -0.15)

circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5)
ax.add_artist(circle)

plt.tight_layout()

plt.title(r'$t = -1M$')

f = plt.gcf()
f.savefig('harm_bsq_o_u_old.png', bbox_inches = 'tight', dpi = 600)

plt.show()
