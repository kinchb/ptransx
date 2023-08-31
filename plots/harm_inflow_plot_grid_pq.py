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

def isco_radius(a):
    z1 = 1.0 + ((1.0 - a*a)**(1.0/3.0))*((1.0 + a)**(1.0/3.0) + (1.0 - a)**(1.0/3.0))
    z2 = np.sqrt(3.0*a*a + z1*z1)
    r  = 3.0 + z2 - np.sqrt((3.0 - z1)*(3.0 + z1 + 2.0*z2))
    print(a, r)
    return r

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

    Mdot_edd = (1.26e38 * M)/(eta * c**2)

    flnm = dir + 'mdot_data.npz'

    t    = np.load(flnm)['t']
    r    = np.load(flnm)['r']
    Mdot = np.load(flnm)['Mdot']

    t -= t[0]

    sigma = np.std(Mdot, axis = 0)
    Mdot  = np.mean(Mdot, axis = 0)

    r_isco = isco_radius(a)

    ax.plot(r/r_isco, Mdot/Mdot_edd, 'k-')

    print(Mdot[0])

    ax.plot(r/r_isco, ((Mdot + sigma)/Mdot_edd), 'k-', linewidth = 0.5)
    ax.plot(r/r_isco, ((Mdot - sigma)/Mdot_edd), 'k-', linewidth = 0.5)

    ax.plot([0.0, 100.0], [mdot, mdot], 'k:')

    ax.set_xticks([1.0, 2.0, 3.0, 6.0, 10.0, 20.0, 30.0, 50.0, 70.0])
    ax.set_xlim([0.2, 70.0])
    ax.set_xscale('log')

#   ax.set_xlabel(r'$r/M$')
    ax.set_xlabel(r'$r/r_\mathrm{ISCO}$')
    ax.set_ylabel(r'$\dot{M}(r)/\dot{M}_\mathrm{Edd}$')

    if (mdot == 0.01):
        ax.set_ylim([-0.04, 0.04])
    else:
        ax.set_ylim([-0.4, 0.4])

    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(M)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot) + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)    + r'\\'
#   text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    ax.text(0.2, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

# --- --- ---

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

fig, axs = plt.subplots(len(mdot), len(spin), sharex = True, sharey = False, figsize = (6.5, 4))

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

fig.savefig('inflow_grid.pdf', bbox_inches = 'tight')

plt.show()
