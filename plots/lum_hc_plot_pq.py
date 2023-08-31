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

def get_spec(M, mdot, a, spec_dir, spec_flnm):
    L_edd = (1.26e38) * M

    """
    first_time = True
    den = 0
    for k in range(20, 1001, 20):
        den += 1
        spec_flnm = spec_dir + 'scat_spec.' + repr(k) + '.0.h5'
        infile    = h5py.File(spec_flnm, 'r')
        spec      = np.array(infile['/data']).T
        infile.close()
        if first_time:
            lum    = np.zeros(len(spec[0]))
            e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)
        for i in range(0, len(spec[0])):
            for j in range(0, len(spec)):
                lum[i] += spec[j][i]
    lum /= den
    """

    try:
        infile = h5py.File(spec_flnm, 'r')
    except:
        e_grid = np.logspace(0, 7, 141, endpoint = True)
        return (e_grid, np.zeros(len(e_grid)))
    spec   = np.array(infile['/data']).T
    infile.close()
    lum    = np.zeros(len(spec[0]))
    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)
    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    print('===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ===== ')
    print('M = ' + '{:.2f}'.format(M) + ', a = ' + '{:.2f}'.format(a) + ', mdot = ' + '{:.2f}'.format(mdot))

    print(np.trapz((lum * 2.417989e14)/L_edd, e_grid))

    return (e_grid, (lum * 2.417989e14)/L_edd)

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

def make_plot(ax, nu, Lnu, abund, mdot, a):
    ax.plot(nu/1.0e3, nu * Lnu, 'k-')
#   ax.plot(nu/1.0e3, Lnu, 'k-')

    e, f, gamma = pl_fit(nu, Lnu, e_min * 1.0e3, e_max * 1.0e3)

    ax.plot(e/1.0e3, e * f, 'r--')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
    ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))

    ax.set_xlim([0.1, 1.0e3])
    ax.set_ylim([1.0e-5, 1.0e-2])

    ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
    ax.set_ylabel(r'$\varepsilon L_\varepsilon / L_\mathrm{Edd}$')

#   if (mdot == 0.01):
#       ax.set_ylim([0.001, 0.1])
#   if (mdot == 0.1):
#       ax.set_ylim([0.01, 1.0])

    mass_label  = r'M/M_\odot &= ' + r'{:g}'.format(abund)    + r'\\[-1ex]'
    mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(mdot)    + r'\\[-1ex]'
    spin_label  = r'a &= '         + r'{:g}'.format(a)       + r'\\'
    text        = r'\begin{align*} ' + mass_label + mdot_label + spin_label + r'\end{align*}'
    ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

    gamma_label = r'$\Gamma = '  + r'{:.2f}'.format(gamma) + r'$'
    ax.text(0.5, 0.9, gamma_label, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

root = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a09_01/'

flnm = [root + 'scat_spec.1000.th3.h5', root + 'scat_spec.1000.th3.hc.h5']

flags = [r'r-', r'k-']

for i in range(0, 2):
    nu, Lnu = get_spec(10.0, 0.01, 0.9, 'tmp', flnm[i])
    plt.plot(nu/1.0e3, nu * Lnu, flags[i])#, label = r'$\Theta_e^\mathrm{cut} = ' + r'{:g}'.format(theta[i]) + r'$')
# plt.legend(loc = 'lower center', frameon = False)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
# ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))
ax.set_xlim([0.1, 2.0e3])
ax.set_ylim([1.0e-4, 1.0e-2])
ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
ax.set_ylabel(r'$\varepsilon L_\varepsilon / L_\mathrm{Edd}$')
mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(0.01)    + r'\\[-1ex]'
spin_label  = r'a &= '         + r'{:g}'.format(0.9)       + r'\\'
text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
plt.tight_layout()
plt.savefig('lum_count.pdf', bbox_inches = 'tight')
plt.show()

# --- --- ---
"""
root = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a0_10/'

theta = [1, 2, 3]

flnm = [root + 'scat_spec.1000.th1.h5', root + 'scat_spec.1000.th2.h5', root + 'scat_spec.1000.th3.h5']

flags = [r'r-', r'g-', r'b-']

for i in range(0, len(theta)):
    nu, Lnu = get_spec(10.0, 0.01, 0.9, 'tmp', flnm[i])
    plt.plot(nu/1.0e3, nu * Lnu, flags[i], label = r'$\Theta_e^\mathrm{cut} = ' + r'{:g}'.format(theta[i]) + r'$')
plt.legend(loc = 'lower center', frameon = False)
ax = plt.gca()
ax.set_xscale('log')
ax.set_yscale('log')
ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.set_xticks((1.0e-1, 1.0e0, 1.0e1, 1.0e2, 1.0e3))
# ax.set_xticklabels((r'$10^{-1}$', r'$10^0$', r'$10^1$', r'$10^2$', r'$10^3$'))
ax.set_xlim([0.1, 2.0e3])
ax.set_ylim([1.0e-7, 1.0e-1])
ax.set_xlabel(r'$\mathrm{energy}\ \left[\mathrm{keV}\right]$')
ax.set_ylabel(r'$\varepsilon L_\varepsilon / L_\mathrm{Edd}$')
mdot_label  = r'\dot{m} &= '   + r'{:g}'.format(0.1)    + r'\\[-1ex]'
spin_label  = r'a &= '         + r'{:g}'.format(0)       + r'\\'
text        = r'\begin{align*} ' + mdot_label + spin_label + r'\end{align*}'
ax.text(0.1, 0.2, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)
plt.tight_layout()
plt.savefig('theta_cut_plot_v2.pdf', bbox_inches = 'tight')
plt.show()
"""