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

L_within_r = np.zeros(len(r))

for i in range(1, len(L_within_r)):
    print(i)
    L_within_r[i] = (4.0 * lum_conv * (cool[:i,::,::] * dV_cgs[:i,::,::]).sum())/L_edd

plt.plot(r, L_within_r/L_within_r[-1], 'k-')

for i in range(0, len(r)):
    if (L_within_r[i] > 0.5*L_within_r[-1]):
        break

plt.plot([1.0, r[i]], [L_within_r[i]/L_within_r[-1], L_within_r[i]/L_within_r[-1]], 'r--')
plt.plot([r[i], r[i]], [0.0, L_within_r[i]/L_within_r[-1]], 'r--')

plt.xlim([r_h, 70])
plt.ylim([0.0, 1.0])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('1', '2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$R/M$')
plt.ylabel(r'$L(r < R)/L_\mathrm{tot}$')

plt.tight_layout()

f = plt.gcf()
f.savefig('L_within_r.pdf', bbox_inches = 'tight')

plt.show()
