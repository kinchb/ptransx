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
e_max = 12.0

def get_ews(M, mdot, a, cont_flnm):
    try:
        infile = h5py.File(cont_flnm, 'r')
    except:
        e_grid = np.logspace(1, 5, 401, endpoint = True)
        return (e_grid, np.ones(len(e_grid)), 0.0)
    cont_spec = np.array(infile['/data']).T
    infile.close()
    line_flnm = cont_flnm.replace('_c_', '_l_')
    try:
        infile = h5py.File(line_flnm, 'r')
    except:
        e_grid = np.logspace(1, 5, 401, endpoint = True)
        return (e_grid, np.ones(len(e_grid)), 0.0)
    line_spec = np.array(infile['/data']).T
    infile.close()

    e_grid = np.logspace(1, 5, 401, endpoint = True)

    angle = np.arccos(np.linspace(1, -1, len(line_spec), endpoint = True)) #* (180.0/np.pi)

    e_min_ndx = 0
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_min*1.0e3):
            e_min_ndx = i-1
            break
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_max*1.0e3):
            e_max_ndx = i+1
            break

    ews = np.zeros(len(angle))
    for ndx in range(0, len(ews)):
        ew_integrand = (line_spec[ndx]/cont_spec[ndx])[e_min_ndx:e_max_ndx]
        ew = np.trapz(ew_integrand, e_grid[e_min_ndx:e_max_ndx])
        ews[ndx] = ew

    return (angle, ews)
#   return (e_grid, line_spec[ndx], ew)

def make_plot(ax, angle, ews, mdot, a):
    ax.plot(np.cos(angle), ews, 'k-')

#   ax.set_xscale('log')
#   ax.set_yscale('log')

#   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#   ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
#   ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))

#   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#   ax.set_xticks((1.0, 5.0, 10.0, 15.0, 20.0))
#   ax.set_xticklabels((r'$3$', r'$10$', r'$20$', r'$40$', r'$79$'))

    ax.set_xlim([-1.0, 1.0])
    ax.set_ylim([0.0, 600.0])
#   ax.set_ylim([0.0, 7.0e15])

    ax.set_xlabel(r'$\cos i$')
    ax.set_ylabel(r'$\mathrm{Fe\ K}\alpha\ \mathrm{EW}\ \left[\mathrm{eV}\right]$')
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

root = '/Users/kinch/panptx/panptx_data/survey/'

abund = [0.5, 1.0, 3.0]

dir = [root + 'h3d_10_a0_01_loFe/', root + 'h3d_10_a0_01/', root + 'h3d_10_a0_01_hiFe/']

flnm = dir.copy()

for i in range(0, len(abund)):
    flnm[i] = dir[i] + 'scat_spec_c_hr.h5'

flags = [r'r-', r'g-', r'b-']

for i in range(0, len(abund)):
    angle, ews = get_ews(10.0, 0.01, 0.0, flnm[i])
    plt.plot(np.cos(angle), ews, flags[i], label = r'$[\mathrm{Fe}]/[\mathrm{Fe}]_\odot = ' + r'{:g}'.format(abund[i]) + r'$')
plt.legend(loc = 'upper left', frameon = False)
ax = plt.gca()
ax.set_xlim([-1.0, 1.0])
ax.set_ylim([0.0, 500.0])
ax.set_xlabel(r'$\cos i$')
ax.set_ylabel(r'$\mathrm{Fe\ K}\alpha\ \mathrm{EW}\ \left[\mathrm{eV}\right]$')
plt.tight_layout()
plt.savefig('ews_abund.pdf', bbox_inches = 'tight')
plt.show()
