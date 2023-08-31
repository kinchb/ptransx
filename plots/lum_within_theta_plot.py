import os

import numpy as np

import h5py

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 15)

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.9
mdot = 0.1
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

dir = '../data/h3d_10_a09_10_gcc/dumps'

harm_data_1 = h5py.File(dir + '/KDHARM0.001500.h5', 'r')
harm_data_2 = h5py.File(dir + '/KDHARM0.RADFLUX.011000.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

r      = np.load(dir + '/../dV_cgs.npz')['r']
th     = np.load(dir + '/../dV_cgs.npz')['th']
phi    = np.load(dir + '/../dV_cgs.npz')['phi']
dV_cgs = np.load(dir + '/../dV_cgs.npz')['dV_cgs']

cool = np.array(harm_data_2['coolfunc_corona'])

r_h = 1.0 + np.sqrt(1.0 - a**2)

cool = np.where(x1 < r_h, 0.0, cool)

L_edd = (1.26e38) * M

lum_conv = (4 * np.pi * c**7 * mdot)/(kappa * G**2 * m_bh**2 * 3.0e-4 * eta)

th = th[int(len(th)/2):] - th[int(len(th)/2)]

L_within_theta_1 = np.zeros(len(th))

for i in range(0, len(th)):
    tmp_cool = np.where(np.logical_or(x2 > np.pi/2 + th[i], x2 < np.pi/2 - th[i]), 0.0, cool)
    tmp_cool = np.where(np.logical_or(x1 < 2.0, x1 > 10.0), 0.0, tmp_cool)
    L_within_theta_1[i] = 4.0 * lum_conv * (tmp_cool * dV_cgs).sum()

L_within_theta_2 = np.zeros(len(th))

for i in range(0, len(th)):
    tmp_cool = np.where(np.logical_or(x2 > np.pi/2 + th[i], x2 < np.pi/2 - th[i]), 0.0, cool)
    tmp_cool = np.where(np.logical_or(x1 < 10.0, x1 > 20.0), 0.0, tmp_cool)
    L_within_theta_2[i] = 4.0 * lum_conv * (tmp_cool * dV_cgs).sum()

L_within_theta_3 = np.zeros(len(th))

for i in range(0, len(th)):
    tmp_cool = np.where(np.logical_or(x2 > np.pi/2 + th[i], x2 < np.pi/2 - th[i]), 0.0, cool)
    tmp_cool = np.where(np.logical_or(x1 < 20.0, x1 > 40.0), 0.0, tmp_cool)
    L_within_theta_3[i] = 4.0 * lum_conv * (tmp_cool * dV_cgs).sum()

L_within_theta_4 = np.zeros(len(th))

for i in range(0, len(th)):
    tmp_cool = np.where(np.logical_or(x2 > np.pi/2 + th[i], x2 < np.pi/2 - th[i]), 0.0, cool)
    tmp_cool = np.where(np.logical_or(x1 < 40.0, x1 > 70.0), 0.0, tmp_cool)
    L_within_theta_4[i] = 4.0 * lum_conv * (tmp_cool * dV_cgs).sum()

L_tot = L_within_theta_1[-1] + L_within_theta_2[-1] + L_within_theta_3[-1] + L_within_theta_4[-1]

plt.plot(th * (180.0/np.pi), L_within_theta_1/L_tot, 'r-', label = r'$2 < r/M \leq 10$')
plt.plot(th * (180.0/np.pi), L_within_theta_2/L_tot, 'g-', label = r'$10 < r/M \leq 20$')
plt.plot(th * (180.0/np.pi), L_within_theta_3/L_tot, 'b-', label = r'$20 < r/M \leq 40$')
plt.plot(th * (180.0/np.pi), L_within_theta_4/L_tot, 'm-', label = r'$40 < r/M \leq 70$')

plt.xlim([0, 90])
# plt.ylim([0, 1.0])

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((0, 15, 30, 45, 60, 75, 90))
ax.set_xticklabels(('$0^\circ$', '$15^\circ$', '$30^\circ$', '$45^\circ$', '$60^\circ$', '$75^\circ$', '$90^\circ$'))

plt.xlabel(r'$\theta\ (\mathrm{from\ midplane})$')
plt.ylabel(r'$L_\mathrm{IC}(\vartheta < \theta)/L_\mathrm{IC,tot}$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('L_within_theta.pdf', bbox_inches = 'tight')

plt.show()
