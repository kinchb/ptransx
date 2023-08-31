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

flnm = '../harm_data/lum_within_r_a0.0.npz'

r         = np.load(flnm)['r']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

L_within_r = cor_lum

dLdr = np.zeros(len(r))

for i in range(1, len(dLdr)):
    dLdr[i] = r[i] * (L_within_r[i] - L_within_r[i-1])/(r[i] - r[i-1])

plt.plot(r, dLdr, 'k-', label = r'$\mathrm{IC}$')

flnm = '../harm_data/old_data/lum_within_r_a0.0.npz'

r         = np.load(flnm)['r']
cor_lum   = np.load(flnm)['cor_lum']

L_within_r = cor_lum

dLdr = np.zeros(len(r))

for i in range(1, len(dLdr)):
    dLdr[i] = r[i] * (L_within_r[i] - L_within_r[i-1])/(r[i] - r[i-1])

plt.plot(r, dLdr, 'r-', label = r'$\mathrm{target}$-$\mathrm{temperature}$')

plt.xlim([2, 70])
plt.ylim([0.0, 0.007])

plt.semilogx()

plt.legend(loc = 'best', frameon = False)

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$r/M$')
plt.ylabel(r'$\partial (L_\mathrm{cor}/L_\mathrm{Edd}) / \partial \log r$')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdr_comp.pdf', bbox_inches = 'tight')

plt.show()
