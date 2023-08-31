import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

flnm = '/mnt/archive/h3d_10_a0_01_c1.5/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_1 = t.copy()
total_lum_1 = total_lum.copy()

t_adj = 10000.0

plt.plot(t - t_adj, cor_lum, 'r-', label = r'$B^2/\rho > 1.5$')

flnm = '/mnt/archive/h3d_10_a0_01_redo/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_2 = t.copy()
total_lum_2 = total_lum.copy()

t_adj = 10000.0

plt.plot(t - t_adj, cor_lum, 'k-', label = r'$B^2/\rho > 1$')

flnm = '/mnt/archive/h3d_10_a0_01_c0.5/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, cor_lum, 'b-', label = r'$B^2/\rho > 0.5$')

# plt.semilogy()

plt.xlim([0, 200])
plt.ylim([0.016, 0.03])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L_\mathrm{IC}/L_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'upper right')

plt.tight_layout()

f = plt.gcf()
f.savefig('cooling_cutoff.pdf', bbox_inches = 'tight')

plt.show()

