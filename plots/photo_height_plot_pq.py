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

flnm = '../harm_data/photo_height_old.npz'

t       = np.load(flnm)['times']
r       = np.load(flnm)['r']
indices = np.load(flnm)['indices'].tolist()
height  = np.load(flnm)['height']

height_0 = height[indices].copy()
height_0 = np.mean(height_0, axis = 0)

flnm = '../harm_data/photo_height.npz'

t       = np.load(flnm)['times']
r       = np.load(flnm)['r']
indices = np.load(flnm)['indices'].tolist()
height  = np.load(flnm)['height']

height_1 = height[indices[:int(len(indices)/2)]].copy()
height_1 = np.mean(height_1, axis = 0)

height_2 = height[indices[int(len(indices)/2):]].copy()
height_2 = np.mean(height_2, axis = 0)

plt.plot(r, height_0/r, 'k-', label = r'$-1000$--$0 M$')
plt.plot(r, height_1/r, 'r-', label = r'$0$--$1000 M$')
plt.plot(r, height_2/r, 'b-', label = r'$1000$--$2000 M$')

plt.xlim([2, 70])
plt.ylim([0.0, 0.175])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$r/M$')
plt.ylabel(r'$H_\mathrm{phot}/r$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('photo_height.pdf', bbox_inches = 'tight')

plt.show()
