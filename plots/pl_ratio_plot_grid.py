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

def get_line(M, mdot, a, spec_dir, cont_flnm, ndx):
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
    ratio = 1.0 + (line_spec[ndx]/cont_spec[ndx])

    angle = np.arccos(np.linspace(1, -1, len(line_spec), endpoint = True)) * (180.0/np.pi)
    print(angle[ndx])

    e_min_ndx = 0
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_min*1.0e3):
            e_min_ndx = i-1
            break
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_max*1.0e3):
            e_max_ndx = i+1
            break

    ew_integrand = (line_spec[ndx]/cont_spec[ndx])[e_min_ndx:e_max_ndx]

    ew = np.trapz(ew_integrand, e_grid[e_min_ndx:e_max_ndx])

    print(ew)

    print('===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ')
    print('M = ' + '{:.2f}'.format(M) + ', a = ' + '{:.2f}'.format(a) + ', mdot = ' + '{:.2f}'.format(mdot))

    return (e_grid, ratio, ew)
#   return (e_grid, line_spec[ndx], ew)

def get_line_by_pl(M, mdot, a, spec_dir, cont_flnm, ndx):
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
    total = line_spec[ndx] + cont_spec[ndx]
    ratio = 1.0 + (line_spec[ndx]/cont_spec[ndx])

    angle = np.arccos(np.linspace(1, -1, len(line_spec), endpoint = True)) * (180.0/np.pi)
    print(angle[ndx])

    e, f, gamma = pl_fit(e_grid, total, e_min * 1.0e3, e_max * 1.0e3)

#   plt.plot(e_grid, e_grid * total)
#   plt.loglog()
#   plt.show()

    e_min_ndx = 0
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_min*1.0e3):
            e_min_ndx = i-1
            break
    for i in range(0, len(e_grid)):
        if (e_grid[i] > e_max*1.0e3):
            e_max_ndx = i+1
            break

#   ratio = ratio[e_min_ndx:e_max_ndx]
    ratio = total[e_min_ndx:e_max_ndx]/f
#   ratio = (line_spec[ndx][e_min_ndx:e_max_ndx])#/(np.trapz(line_spec[ndx][e_min_ndx:e_max_ndx], e_grid[e_min_ndx:e_max_ndx]))

#   plt.plot(e_grid[e_min_ndx:e_max_ndx], ratio)
#   plt.show()

    ew_integrand = (line_spec[ndx]/cont_spec[ndx])[e_min_ndx:e_max_ndx]

    ew = np.trapz(ew_integrand, e_grid[e_min_ndx:e_max_ndx])

    print(ew)

    print('===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ')
    print('M = ' + '{:.2f}'.format(M) + ', a = ' + '{:.2f}'.format(a) + ', mdot = ' + '{:.2f}'.format(mdot))

    return (e, ratio, ew)

def make_plot(ax, nu, ratio, mdot, a, flag, i):
    ax.plot(nu/1.0e3, np.ones(len(nu)), 'k:')
    ax.plot(nu/1.0e3, ratio, flag, label = r'$i = ' + r'{:g}'.format(angle_list[i]) + r'^\circ$')

#   ax.set_xscale('log')
#   ax.set_yscale('log')

#   ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
#   ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
#   ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))

    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks((3.0, 10.0, 20.0, 40.0, 79.0))
    ax.set_xticklabels((r'$3$', r'$10$', r'$20$', r'$40$', r'$79$'))

    ax.set_xlim([e_min, e_max])
    ax.set_ylim([0.5, 1.5])
#   ax.set_ylim([0.0, 7.0e15])

    ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
    ax.set_ylabel(r'$\mathrm{total/power\ law\ fit}$')
#   ax.set_ylabel(r'$\mathrm{total/continuum}$')
#   ax.set_ylabel(r'$\mathrm{line\ flux}$')

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
    ax.text(0.6, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
#   ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

#   gamma_label = r'$\Gamma = '  + r'{:.2f}'.format(gamma) + r'$'
#   ax.text(0.5, 0.9, gamma_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/panptx/panptx_data/survey/'

dir  = [[root + 'h3d_10_a0_01/', root + 'h3d_10_a05_01/', root + 'h3d_10_a09_01/'],
        [root + 'h3d_10_a0_10/', root + 'h3d_10_a05_10/', root + 'h3d_10_a09_10/']]

spec_ndx_list = [0, 3, 10, 15, 20]
angle_list    = [0, 30, 60, 75, 90]

flags = ['c-', 'g-', 'b-', 'y-', 'm-']

flnm = dir.copy()

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        flnm[i][j] = dir[i][j] + 'scat_spec_c_hr.h5'

fig, axs = plt.subplots(len(mdot), len(spin), sharex = True, sharey = True, figsize = (12, 6))

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        for k in range(0, len(angle_list)):
            nu, ratio, ew = get_line_by_pl(10.0, mdot[i], spin[j], dir[i][j], flnm[i][j], spec_ndx_list[k])
            make_plot(axs[i][j], nu, ratio, mdot[i], spin[j], flags[k], k)

for ax in axs.flat:
    ax.label_outer()

ax = plt.gca()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.98, box.height])
# Put a legend to the right of the current axis
ax.legend(loc='center left', bbox_to_anchor=(1, 0.7), frameon = False)

plt.tight_layout()

# lgd = fig.legend(handles = handle, loc = 'upper center', ncol = 4)
fig.savefig('pl_ratio.pdf', bbox_inches = 'tight')

plt.show()
