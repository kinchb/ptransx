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

e_min = 3.0
e_max = 79.0

def get_spec(M, flnm, ndx):
    L_edd = (1.26e38) * M

    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    angle = np.arccos(np.linspace(1, -1, len(spec), endpoint = True)) * (180.0/np.pi)

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

def make_plot(ax, angle, gamma, mdot, a):
    ax.plot(np.cos(angle), gamma, 'k-')

#   ax.set_xscale('log')
#   ax.set_yscale('log')

#   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#   ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
#   ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))

#   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#   ax.set_xticks((1.0, 5.0, 10.0, 15.0, 20.0))
#   ax.set_xticklabels((r'$3$', r'$10$', r'$20$', r'$40$', r'$79$'))

    ax.set_xlim([-1.0, 1.0])
    if (mdot == 0.01):
        ax.set_ylim([1.5, 3.0])
    else:
        ax.set_ylim([2.5, 4.0])
#   ax.set_ylim([0.0, 7.0e15])

    ax.set_xlabel(r'$\cos i$')
    ax.set_ylabel(r'$\Gamma$')
#   ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
#   ax.set_ylabel(r'$\mathrm{total/power\ law\ fit}$')
#   ax.set_ylabel(r'$\mathrm{total/continuum}$')
#   ax.set_ylabel(r'$\mathrm{line\ flux} \varepsilon I^\mathrm{line}_\varepsilon\ \left[\mathrm{erg}\ \mathrm{s}^{-1}\ \mathrm{cm}^{-2}\ \mathrm{sr}^{-1}\right]$')

#   if (mdot == 0.01):
#       ax.set_ylim([0.001, 0.1])
#   if (mdot == 0.1):
#       ax.set_ylim([0.01, 1.0])

    mass = 10.0
    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(mass)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot)    + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)       + r'\\'
#   text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    ax.text(0.7, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

#   gamma_label = r'$\Gamma = '  + r'{:.2f}'.format(gamma) + r'$'
#   ax.text(0.5, 0.9, gamma_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

flnm = dir.copy()

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        flnm[i][j] = dir[i][j] + 'scat_spec.1000.5.h5'

fig, axs = plt.subplots(len(mdot), len(spin), sharex = True, sharey = False, figsize = (12, 6))

angle = np.arccos(np.linspace(1, -1, 41, endpoint = True))

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        gamma = np.zeros(len(angle))
        for k in range(0, len(angle)):
            nu, Lnu = get_spec(10.0, flnm[i][j], k)
            e, f, gam = pl_fit(nu, Lnu, e_min * 1.0e3, e_max * 1.0e3)
            gamma[k] = gam
        make_plot(axs[i][j], angle, gamma, mdot[i], spin[j])

for ax in axs.flat:
    ax.label_outer()

plt.tight_layout()

# lgd = fig.legend(handles = handle, loc = 'upper center', ncol = 4)
fig.savefig('gammas.pdf', bbox_inches = 'tight')

plt.show()
