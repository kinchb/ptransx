import numpy as np

import h5py

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

def get_lum(M, mdot, a, spec_dir, lum_flnm):
    L_edd = (1.26e38) * M

    t_rt   = []
    lum_rt = []
    for k in range(20, 1001, 20):
        t_rt.append(k)
        spec_flnm = spec_dir + 'scat_spec.' + repr(k) + '.0.h5'
        infile    = h5py.File(spec_flnm, 'r')
        spec      = np.array(infile['/data']).T
        infile.close()
        lum       = np.zeros(len(spec[0]))
        for i in range(0, len(spec[0])):
            for j in range(0, len(spec)):
                lum[i] += spec[j][i]
        e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)
        lum_rt.append(np.trapz(lum * 2.417989e14, e_grid)/L_edd)
    t_rt   = np.array(t_rt)
    lum_rt = np.array(lum_rt)

    t         = np.load(lum_flnm)['t']
    total_lum = np.load(lum_flnm)['total_lum']
    cor_lum   = np.load(lum_flnm)['cor_lum']
    disk_lum  = np.load(lum_flnm)['disk_lum']

    t -= t[0]

    print('===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ')
    print('M = ' + '{:.2f}'.format(M) + ', a = ' + '{:.2f}'.format(a) + ', mdot = ' + '{:.2f}'.format(mdot))
    print('L = ' + repr(np.mean(lum_rt)) + ', disk frac = ' + repr((disk_lum/total_lum).mean()))

    return (t_rt, lum_rt, t, total_lum, cor_lum, disk_lum)

def make_plot(ax, t_rt, lum_rt, t, total_lum, cor_lum, disk_lum, mdot, a):
    lbl2, = ax.plot(t,    total_lum, 'k-', label = 'total')
    lbl3, = ax.plot(t,    cor_lum,   'r-', label = 'corona')
    lbl4, = ax.plot(t,    disk_lum,  'b-', label = 'disk')
    lbl1, = ax.plot(t_rt, lum_rt,    'm-', label = 'ray-traced')

    ax.plot([-1.0, 1000.0], [mdot, mdot], 'k:')
    ax.set_yscale('log')
    ax.set_xlim([0.0, 1000.0])
    ax.set_xlabel(r'$t/M$')
    ax.set_ylabel(r'$L/L_\mathrm{Edd}$')

    if (mdot == 0.01):
        ax.set_ylim([0.001, 0.1])
    if (mdot == 0.1):
        ax.set_ylim([0.01, 1.0])

    mdot_label = r'$\dot{m} = ' + r'{:g}'.format(mdot) + r'$'
    spin_label = r'$a = '       + r'{:g}'.format(a)    + r'$'
    ax.text(0.1, 0.05, mdot_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
    ax.text(0.1, 0.15, spin_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

    return (lbl1, lbl2, lbl3, lbl4)

mdot = [0.01, 0.1]
spin = [0.0, 0.5, 0.9]

root = '/Users/kinch/spin_analysis/data/'

dir  = [[root + 'h3d_10_a0_01_long_run/scat_specs/', root + 'h3d_10_a05_01_gcc/scat_specs/', root + 'h3d_10_a09_01_gcc/scat_specs/'],
        [root + 'h3d_10_a0_10_intel/scat_specs/',    root + 'h3d_10_a05_10_gcc/scat_specs/', root + 'h3d_10_a09_10_gcc/scat_specs/']]

flnm = [[root + 'h3d_10_a0_01_long_run/lum_data_ucon0_a0.0.npz', root + 'h3d_10_a05_01_gcc/lum_data_ucon0_a0.5.npz', root + 'h3d_10_a09_01_gcc/lum_data_ucon0_a0.9.npz'],
        [root + 'h3d_10_a0_10_intel/lum_data_ucon0_a0.0.npz',    root + 'h3d_10_a05_10_gcc/lum_data_ucon0_a0.5.npz', root + 'h3d_10_a09_10_gcc/lum_data_ucon0_a0.9.npz']]

fig, axs = plt.subplots(len(mdot), len(spin), sharex = 'col', sharey = 'row', figsize = (9, 6))

for i in range(0, len(mdot)):
    for j in range(0, len(spin)):
        t_rt, lum_rt, t, total_lum, cor_lum, disk_lum = get_lum(10.0, mdot[i], spin[j], dir[i][j], flnm[i][j])
        handle = make_plot(axs[i][j], t_rt, lum_rt, t, total_lum, cor_lum, disk_lum, mdot[i], spin[j])

for ax in axs.flat:
    ax.label_outer()

lgd = fig.legend(handles = handle, loc = 'upper center', ncol = 4)

fig.savefig('cooling_rate_grid.pdf', bbox_extra_artists = (lgd,), bbox_inches = 'tight')

plt.show()
