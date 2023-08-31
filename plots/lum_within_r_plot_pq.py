import numpy as np

from scipy.interpolate import interp1d

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

r_half = interp1d(total_lum/total_lum[-1], r)(0.5)

print(r_half)

plt.plot(r, total_lum, 'k-', label = r'$\mathrm{total}$')
plt.plot(r, cor_lum,   'r-', label = r'$\mathrm{corona}$')
plt.plot(r, disk_lum/(r**2),  'b-', label = r'$\mathrm{disk}$')

plt.plot([0.0, r_half], [0.5 * total_lum[-1], 0.5 * total_lum[-1]], 'k--')
plt.plot([r_half, r_half], [0.0, 0.5 * total_lum[-1]], 'k--')

plt.xlim([2, 70])
plt.ylim([0.0, 0.021])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$R/M$')
plt.ylabel(r'$L(r < R)/L_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'upper left')

plt.tight_layout()

# f = plt.gcf()
# f.savefig('lum_within_r.pdf', bbox_inches = 'tight')

plt.show()
