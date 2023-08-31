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
rc('font', family = 'serif', size = 10)

e_min = 3.0
e_max = 79.0

def get_spec(M, flnm, ndx):
    L_edd = (1.26e38) * M

    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    angle = np.arccos(np.linspace(1, -1, len(spec), endpoint = True)) * (180.0/np.pi)

    print(len(spec))

    print(angle[ndx])

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return e_grid, 4.0*np.pi*(spec[ndx] * 2.417989e14 * (1.0/(2.0 * np.pi)) * (2.0/len(spec)))/L_edd

def pl_fit(e, f, e_min, e_max):
    for i in range(0, len(e)):
        if (e[i] > e_min):
            break
    min_ndx = i - 1

    for i in range(0, len(e)):
        if (e[i] > e_max):
            break
    max_ndx = i + 1

    e = e[min_ndx:max_ndx]
    f = f[min_ndx:max_ndx]

    gamma = -np.polyfit(np.log10(e), np.log10(f/e), 1)[0]

    p = np.poly1d(np.polyfit(np.log10(e), np.log10(f), 1))

    return (e, 10.0**p(np.log10(e)), gamma)

spec_ndx_list = [0, 3, 10, 15, 20]
angle_list    = [0, 30, 60, 75, 90]

flags = ['c-', 'g-', 'b-', 'y-', 'm-']

def make_plot(ax, mass, mdot, a, flnm):
    for i in range(0, len(spec_ndx_list)):
        nu, Lnu = get_spec(10.0, flnm, spec_ndx_list[i])
        e, f, gamma = pl_fit(nu, Lnu, e_min * 1.0e3, e_max * 1.0e3)
        if (i == 0):
            ax.plot(nu/1.0e3, nu * Lnu, flags[i], label = r'$i = \phantom{0}' + r'{:g}'.format(angle_list[i]) + r'^\circ : \Gamma = ' + r'{:.2f}'.format(gamma) + r'$')
        else:
            ax.plot(nu/1.0e3, nu * Lnu, flags[i], label = r'$i = ' + r'{:g}'.format(angle_list[i]) + r'^\circ : \Gamma = ' + r'{:.2f}'.format(gamma) + r'$')
        ax.plot(e/1.0e3, e * f, 'r--')

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
    ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))
    ax.set_xlim([0.1, 1.0e3])
    ax.set_ylim([1.0e-7, 1.0e-3])
    ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$', fontsize = 12)
    # ax.set_ylabel(r'$\varepsilon I_\varepsilon\ \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\ \mathrm{sr}^{-1}\right]$')
    ax.set_ylabel(r'$4\pi \varepsilon I_\varepsilon / L_\mathrm{Edd}$', fontsize = 12)
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.98, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.7), frameon = False)

    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(10)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot)  + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)     + r'\\'
    # text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    if (mdot == 0.01):
        ax.text(0.6, 0.75, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes, fontsize = 12)
    else:
        ax.text(0.2, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes, fontsize = 12)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

flnm = dir.copy()

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        flnm[i][j] = dir[i][j] + 'scat_spec.1000.5.h5'

fig, axs = plt.subplots(nrows = len(spin), ncols = len(mdot), sharex = True, sharey = True, figsize = (9, 9))

for i in range(0, len(spin)):
    for j in range(0, len(mdot)):
        make_plot(axs[i][j], 10.0, mdot[j], spin[i], flnm[j][i])

for ax in axs.flat:
    ax.label_outer()

plt.tight_layout()

plt.savefig('inc_grid.pdf', bbox_inches = 'tight')

plt.show()
