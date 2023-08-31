import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

flnm = '/mnt/archive/h3d_10_a0_01_higher_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'r-', label = r'$1/\beta_\mathrm{min} = 50$')

# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_redo/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'g-', label = r'$1/\beta_\mathrm{min} = 100$')

# --- --- ---
"""
flnm = '/mnt/archive/h3d_10_a0_01_lower_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'b-', label = r'$1/\beta_\mathrm{min} = 150$')
"""
# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_lowest_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'b-', label = r'$1/\beta_\mathrm{min} = 300$')

# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_even_lower_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'c-', label = r'$1/\beta_\mathrm{min} = 500$')

# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_super_low_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'm-', label = r'$1/\beta_\mathrm{min} = 1000$')

# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_1o2000_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'y-', label = r'$1/\beta_\mathrm{min} = 2000$')

# --- --- ---

flnm = '/mnt/archive/h3d_10_a0_01_1o5000_beta/lum_data_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k-', label = r'$1/\beta_\mathrm{min} = 5000$')

# --- --- ---

plt.semilogy()

plt.xlim([0, 200])
plt.ylim([0.02, 0.1])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L/L_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'upper right', fontsize = 14)

plt.tight_layout()

f = plt.gcf()
f.savefig('beta_cutoff.pdf', bbox_inches = 'tight')

plt.show()

