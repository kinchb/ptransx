import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

def make_Mdot(M, mdot, a, flnm):
    if (a == 0.0):
        eta = 0.0572
    if (a == 0.5):
        eta = 0.0821
    if (a == 0.9):
        eta = 0.1558
    if (a == 0.99):
        eta = 0.2640

    Mdot_edd = (1.26e38 * M)/(eta * c**2)

    t    = np.load(flnm)['t']
    r    = np.load(flnm)['r']
    Mdot = np.load(flnm)['Mdot']

    t -= t[0]

    t_1000_ndx = 0
    i = 0
    while (i < len(t)):
        if (t[i] > 1000.0):
            break
        else:
            i += 1
    t_1000_ndx = i+1

    sigma = np.std(Mdot[:t_1000_ndx,::],  axis = 0)
    Mdot  = np.mean(Mdot[:t_1000_ndx,::], axis = 0)

    rh = 1.0 + np.sqrt(1.0 - a**2)

    rh_ndx = 0
    while True:
        if (r[rh_ndx] > rh):
            break
        rh_ndx += 1

    print('===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ')
    print('M = ' + '{:.2f}'.format(M) + ', a = ' + '{:.2f}'.format(a) + ', mdot = ' + '{:.2f}'.format(mdot))
    print(repr(Mdot[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma[rh_ndx]/Mdot_edd))

    return (r, Mdot/Mdot_edd)

def make_plot(ax, r, Mdot, mdot, a):
    ax.plot(r, Mdot, 'k-')
    ax.plot([-1.0, 100.0], [mdot, mdot], 'k:')
    ax.set_xscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks((1, 2, 6, 10, 20, 30, 40, 50, 60, 70))
    ax.set_xticklabels(('1', '2', '6', '10', '20', '', '', '', '', '70'))
    ax.set_xlim([1.0, 70])
    ax.set_xlabel(r'$r/M$')
    ax.set_ylabel(r'$\dot{M}(r)/\dot{M}_\mathrm{Edd}$')

    if (mdot == 0.01):
        ax.set_ylim([-0.05, 0.05])
    if (mdot == 0.1):
        ax.set_ylim([-0.5, 0.5])

    mdot_label = r'$\dot{m} = ' + r'{:g}'.format(mdot) + r'$'
    spin_label = r'$a = '       + r'{:g}'.format(a)    + r'$'
    ax.text(0.1, 0.1, mdot_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
    ax.text(0.1, 0.2, spin_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/spin_analysis/data/'

flnm = [[root + 'h3d_10_a0_01_long_run/mdot_data_a0.0.npz', root + 'h3d_10_a05_01_gcc/mdot_data_a0.5.npz', root + 'h3d_10_a09_01_gcc/mdot_data_a0.9.npz'],
        [root + 'h3d_10_a0_10_intel/mdot_data_a0.0.npz',    root + 'h3d_10_a05_10_gcc/mdot_data_a0.5.npz', root + 'h3d_10_a09_10_gcc/mdot_data_a0.9.npz']]

fig, axs = plt.subplots(len(mdot), len(spin), sharex = 'col', sharey = 'row', figsize = (8, 6))

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        r, Mdot = make_Mdot(10.0, mdot[i], spin[j], flnm[i][j])
        make_plot(axs[i][j], r, Mdot, mdot[i], spin[j])

for ax in axs.flat:
    ax.label_outer()

plt.tight_layout()

plt.show()
