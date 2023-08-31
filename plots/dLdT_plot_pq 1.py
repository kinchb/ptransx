import numpy as np

import h5py

import matplotlib.pyplot as plt
custom_preamble = {
    "text.usetex": True,
    "text.latex.preamble": [
        r"\usepackage{amsmath}", # for the align enivironment
        ],
    }
plt.rcParams.update(custom_preamble)

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

M     = 10.0
L_edd = (1.26e38) * M

def get_dLdT(dLdT_flnm):
    T_grids = np.load(dLdT_flnm)['T_grids']
    pdfs    = np.load(dLdT_flnm)['cool_pdfs']

    pdf = np.mean(pdfs, axis = 0)
#   pdf = pdfs[-1]

    return (T_grids[0], pdf/L_edd)

def make_plot(ax, T_grid, pdf, mdot, a):
    ax.plot(T_grid, pdf, 'k-')

    ax.set_xscale('log')

    ax.set_xlim([1.0e3, 1.0e8])
#   ax.set_ylim([1.0e32, 5.0e36])

    ax.set_xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
    ax.set_ylabel(r'$\partial \left(L_\mathrm{IC}/L_\mathrm{Edd}\right) / \partial \log T_e$')

    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks((1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7, 1.0e8))
    ax.set_xticklabels((r'$10^3$', r'$10^4$', r'$10^5$', r'$10^6$', r'$10^7$', r'$10^8$'))

    if (mdot == 0.01):
        ax.set_ylim([0.0, 0.075])
    if (mdot == 0.1):
        ax.set_ylim([0.0, 0.2])

    mdot_label  = r'\dot{m} &= ' + r'{:g}'.format(mdot)    + r'\\[-1ex]'
    spin_label  = r'a &= '       + r'{:g}'.format(a)       + r'\\'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    ax.text(0.1, 0.8, text,  horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

#   gamma_label = r'$\Gamma = '  + r'{:.2f}'.format(gamma) + r'$'
#   ax.text(0.1, 0.9, mdot_label,  horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.8, spin_label,  horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.5, 0.9, gamma_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/spin_analysis/data/'

flnm = [[root + 'h3d_10_a0_01_long_run/dLdT_pan_data_a0.0.npz', root + 'h3d_10_a05_01_gcc/dLdT_pan_data_a0.5.npz', root + 'h3d_10_a09_01_gcc/dLdT_pan_data_a0.9.npz'],
        [root + 'h3d_10_a0_10_intel/dLdT_pan_data_a0.0.npz',    root + 'h3d_10_a05_10_gcc/dLdT_pan_data_a0.5.npz', root + 'h3d_10_a09_10_gcc/dLdT_pan_data_a0.9.npz']]

fig, axs = plt.subplots(len(mdot), len(spin), sharex = True, sharey = 'row', figsize = (9, 6))

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        T_grid, pdf = get_dLdT(flnm[i][j])
        make_plot(axs[i][j], T_grid, pdf, mdot[i], spin[j])

for ax in axs.flat:
    ax.label_outer()

plt.tight_layout()

# lgd = fig.legend(handles = handle, loc = 'upper center', ncol = 4)

fig.savefig('dLdT_grid.pdf', bbox_inches = 'tight')

plt.show()
