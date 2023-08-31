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

def make_plot(ax, dir, M, mdot, a, plotted_already):
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

    harm_data_1 = h5py.File(dir + 'KDHARM0.001500.h5', 'r')

    x1 = np.array(harm_data_1['x1'])
    x2 = np.array(harm_data_1['x2'])
    x3 = np.array(harm_data_1['x3'])

    r   = np.zeros(len(x1))
    th  = np.zeros(len(x2[0]))
    phi = np.zeros(len(x3[0][0]))

    for i in range(0, len(x1)):
        r[i]   = x1[i][0][0]
    for i in range(0, len(x2[0])):
        th[i]  = x2[0][i][0]
    for i in range(0, len(x3[0][0])):
        phi[i] = x3[0][0][i]

    max_r = 31.0
    for i in range(0, len(r)):
        if (r[i] > max_r):
            break
    max_r_index = i
    r = r[:max_r_index]

    rho = np.array(harm_data_1['rho'])
    rho = rho[:max_r_index,::,::]
    uu  = np.array(harm_data_1['uu'])
    uu  = uu[:max_r_index,::,::]

    m_e = 9.1093826e-28
    m_i = 1.67262171e-24
    T_e = 511.0 * (1.0/(1.0 + 1.21)) * ((5.0/3.0) - 1) * (m_i/m_e) * (uu/rho)

    r_abbr = np.load(dir + 'photo_locs.npz')['r_abbr']
    emtop  = np.load(dir + 'photo_locs.npz')['emtop']
    embot  = np.load(dir + 'photo_locs.npz')['embot']

    T_e = np.mean(T_e, axis = 2)

    slice_data = T_e

    rad, theta = np.meshgrid(r, th)

    ax.set_theta_zero_location('N')
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(180)
    if not plotted_already:
        ax.set_rlabel_position(80)
#       ax.set_xticks(np.pi/180.0 * np.array([0, 30, 60, 90, 120, 150, 180]))
        ax.set_xticks(np.pi/180.0 * np.array([0, 45, 90, 135, 180]))
        ax.set_yticks(np.array([0, 10, 20, 30, 40, 50]))
        plotted_already = True
    else:
        ax.set_xticks([])
        ax.set_yticks([])

    cmap = plt.get_cmap('jet')

#   im = ax.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = -10.0, vmax = -5.0, shading = 'gourand')
    im = ax.pcolormesh(theta, rad, np.log10(slice_data.T), cmap = cmap, vmin = 0.0, vmax = 5.0, shading = 'gourand')

    ax.plot(emtop, r_abbr, 'w-', linewidth = 1)
    ax.plot(embot, r_abbr, 'w-', linewidth = 1)

    circle = plt.Circle((0.0, 0.0), 1.0 + np.sqrt(1.0 - a**2), transform=ax.transData._b, color = 'black', fill = True, linewidth = 0.5, zorder = 1000)
    ax.add_artist(circle)

    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(M)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot) + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)    + r'\\'
#   text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
    ax.text(-0.15, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

    return im, plotted_already

# --- --- ---

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

nrow = len(mdot)
ncol = len(spin)

# fig = plt.figure(figsize=(3*ncol+1, 3*nrow+1))
fig = plt.figure(figsize = (6.5, 5.5))

gs = gridspec.GridSpec(nrow, ncol,
         wspace=0.0, hspace=0.0,
         top=1.-0.5/(nrow+1), bottom=0.5/(nrow+1),
         left=0.5/(ncol+1), right=1-0.5/(ncol+1)) 
#        left=0.25/(ncol+1), right=1-0.25/(ncol+1))

plotted_already = False
for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        ax = plt.subplot(gs[i, j], polar = True)
        im, plotted_already = make_plot(ax, dir[i][j], 10.0, mdot[i], spin[j], plotted_already)

axs = fig.axes

ticks = [0, 1, 2, 3, 4, 5]

fig.colorbar(im, ax = axs, label = r'$\log T_e \left[\mathrm{keV}\right]$', shrink = 0.6, pad = -0.05, ticks = ticks)

plt.savefig('harm_Te_grid.png', dpi = 600, bbox_inches = 'tight')

plt.show()
