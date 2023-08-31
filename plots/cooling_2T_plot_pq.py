import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

flnm = '../harm_data/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

plt.plot([0.0, 0.0], [-1.0, 1.0], 'k--')

t_adj = 10000.0

# plt.plot(t - t_adj, total_lum, 'k-', label = r'$\mathrm{total}$')
plt.plot(t - t_adj, cor_lum,   'r-', label = r'$\mathrm{corona}$')
# plt.plot(t - t_adj, disk_lum,  'b-', label = r'$\mathrm{disk}$')

flnm = '../harm_data/lum_data_2T_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

# plt.plot(t - t_adj, total_lum, 'r-') #, label = r'$\mathrm{total}$')
plt.plot(t - t_adj, cor_lum,   'r-') #, label = r'$\mathrm{corona}$')
# plt.plot(t - t_adj, disk_lum,  'b-') #, label = r'$\mathrm{disk}$')

# plt.semilogy()

plt.xlim([0.0, 1000.0])
# plt.ylim([4.0e-3, 3.0e-1])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L/L_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'upper right')

plt.tight_layout()

f = plt.gcf()
f.savefig('cooling_2T.pdf', bbox_inches = 'tight')

plt.show()
