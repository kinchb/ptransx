import os

import numpy as np

import h5py

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

th             = np.load('../harm_data/lum_within_theta_2-10.npz')['th']
L_within_theta = np.load('../harm_data/lum_within_theta_2-10.npz')['L_within_theta']

# plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'g-', label = r'$2$--$10 M$')
plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'g-', label = r'$2 < r/M \leq 10$')

th             = np.load('../harm_data/lum_within_theta_10-20.npz')['th']
L_within_theta = np.load('../harm_data/lum_within_theta_10-20.npz')['L_within_theta']

# plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'b-', label = r'$10$--$20 M$')
plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'b-', label = r'$10 < r/M \leq 20$')

th             = np.load('../harm_data/lum_within_theta_20-40.npz')['th']
L_within_theta = np.load('../harm_data/lum_within_theta_20-40.npz')['L_within_theta']

# plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'c-', label = r'$20$--$40 M$')
plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'c-', label = r'$20 < r/M \leq 40$')

th             = np.load('../harm_data/lum_within_theta_40-65.npz')['th']
L_within_theta = np.load('../harm_data/lum_within_theta_40-65.npz')['L_within_theta']

# plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'm-', label = r'$40$--$65 M$')
plt.plot(th * (180.0/np.pi), L_within_theta/L_within_theta[-1], 'm-', label = r'$40 < r/M \leq 65$')

plt.xlim([0, 90])
plt.ylim([0, 1.0])

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
