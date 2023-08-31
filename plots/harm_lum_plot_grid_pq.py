# --- --- --- preamble --- --- ---

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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import rc
from matplotlib import gridspec

import imageio

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 10)

import os

from scipy.interpolate import interp1d

# --- --- ---

def get_lum(flnm):
    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return np.trapz(lum * 2.417989e14, e_grid)

def make_plot(ax, dir, M, mdot, a):
    # useful constants (cgs)
    m_bh   = M * 2.0e33
    G      = 6.6726e-8
    c      = 3.0e10
    kappa  = 0.4

    if (a == 0.0):
        eta = 0.0572
    if (a == 0.5):
        eta = 0.0821
    if (a == 0.9):
        eta = 0.1558
    if (a == 0.99):
        eta = 0.2640

    flnm = dir + 'lum_data.npz'

    t         = np.load(flnm)['t']
    total_lum = np.load(flnm)['total_lum']
    cor_lum   = np.load(flnm)['cor_lum']
    disk_lum  = np.load(flnm)['disk_lum']

#   plt.plot([0.0, 0.0], [-1.0, 1.0], 'k--')

    L_edd = (1.26e38) * M
    t_rt = []
    L_rt = []
    for i in range(100, 1001, 100):
        t_rt.append(i)
        L_rt.append(get_lum(dir + 'scat_spec.' + repr(i) + '.0.h5')/L_edd)
    t_rt = np.array(t_rt)
    L_rt = np.array(L_rt)

    print(L_edd * L_rt.mean())

    ax.plot([0.0, 1000.0], [mdot, mdot], 'k:')
    ax.plot(t - t[0], total_lum, 'k-', label = r'$\mathrm{total}$')
    ax.plot(t_rt, L_rt, 'm-', label = r'$\mathrm{ray}$-$\mathrm{traced}$')
    ax.plot(t - t[0], cor_lum,   'r-', label = r'$\mathrm{corona}$')
    ax.plot(t - t[0], disk_lum,  'b-', label = r'$\mathrm{disk}$')

    disk_fraction = disk_lum[100:1000].mean()/total_lum[100:1000].mean()

#   print(disk_fraction)

    mean_total_lum = L_edd * total_lum[100:1000].mean()

#   print(mean_total_lum)

    ax.set_yscale('log')

    ax.set_xlim([0, 1000])
    ax.set_ylim([1.0e-3, 1.0])

    ax.set_ylabel(r'$L/L_\mathrm{Edd}$')
    ax.set_xlabel(r'$t/M$')

    ax.set_yticks([1.0e-3, 1.0e-2, 1.0e-1, 1.0])
    ax.set_xticks([0.0, 500.0, 1000.0])

    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(M)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot) + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)    + r'\\'
#   text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    if mdot == 0.01:
        ax.text(0.2, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
    else:
        ax.text(0.2, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

# --- --- ---

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

fig, axs = plt.subplots(len(mdot), len(spin), sharex = True, sharey = True, figsize = (6.5, 4))

plotted_already = False
for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        make_plot(axs[i][j], dir[i][j], 10.0, mdot[i], spin[j])

for ax in axs.flat:
    ax.label_outer()

ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.98, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.7), frameon = False)

plt.tight_layout()

fig.savefig('lum_grid.pdf', bbox_inches = 'tight')

plt.show()
