import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

def deriv(x):
    dx = np.zeros(len(x))

    for i in range(0, len(x)):
        if (i == 0):
            dx[i] = (-3*x[0] + 4*x[1] - x[2])/2
        elif (i == len(x)-1):
            dx[i] = (3*x[len(x)-1] - 4*x[len(x)-2] + x[len(x)-3])/2
        else:
            dx[i] = (x[i+1] - x[i-1])/2

    return dx

flnm = '../harm_data/lum_within_r_a0.0.npz'

# t         = np.load(flnm)['t']
r         = np.load(flnm)['r']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

# total_lum = np.mean(disk_lum, axis = 0)
# disk_lum  = np.mean(disk_lum, axis = 0)

plt.plot(r, r*(deriv(total_lum)/deriv(r)), 'k-')

plt.xlim([2, 70])
# plt.ylim([0.0, 1.0])

plt.loglog()

plt.show()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$r/M$')
# plt.ylabel(r'$L(r < R)/L_\mathrm{total}$')

plt.legend(frameon = False, loc = 'upper left')

plt.tight_layout()

# f = plt.gcf()
# f.savefig('lum_within_r.pdf', bbox_inches = 'tight')

plt.show()
