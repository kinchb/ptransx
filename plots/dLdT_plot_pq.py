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

e_min = 1.0e-1
e_max = 1.0e5

M = 10.0

L_edd = (1.26e38) * M

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a0_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'r-', label = '$a = 0$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a05_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'g-', label = '$a = 0.5$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a09_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'b-', label = '$a = 0.9$')

plt.semilogx()

plt.xlim([e_min, e_max])
plt.ylim([0, 0.7])

plt.xlabel(r'$T_e\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\partial L_\mathrm{IC} / \partial \log T_e\ (\mathrm{normalized})$')

ax = plt.gca()
text = r'$\dot{m} = ' + r'{:g}'.format(0.01) + r'$'
ax.text(0.05, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

# plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_01.pdf', bbox_inches = 'tight')

plt.show()

# --- --- --- --- ---

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a0_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'r-', label = '$a = 0$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a05_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'g-', label = '$a = 0.5$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a09_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['cool_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'b-', label = '$a = 0.9$')

plt.semilogx()

plt.xlim([e_min, e_max])
plt.ylim([0, 0.7])

plt.xlabel(r'$T_e\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\partial L_\mathrm{IC} / \partial \log T_e\ (\mathrm{normalized})$')

ax = plt.gca()
text = r'$\dot{m} = ' + r'{:g}'.format(0.1) + r'$'
ax.text(0.05, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_10.pdf', bbox_inches = 'tight')

plt.show()

# --- --- ---

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a0_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'r-', label = '$a = 0$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a05_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'g-', label = '$a = 0.5$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a09_01/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.e03
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'b-', label = '$a = 0.9$')

plt.semilogx()

plt.xlim([e_min, e_max])
plt.ylim([0, 1.4])

plt.xlabel(r'$T_e\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\partial \mathcal{M}_\mathrm{cor} / \partial \log T_e\ (\mathrm{normalized})$')

ax = plt.gca()
text = r'$\dot{m} = ' + r'{:g}'.format(0.01) + r'$'
ax.text(0.05, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

# plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dMdT_01.pdf', bbox_inches = 'tight')

plt.show()

# --- --- ---

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a0_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'r-', label = '$a = 0$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a05_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'g-', label = '$a = 0.5$')

flnm = '/Users/kinch/panptx/panptx_data/survey/h3d_10_a09_10/dLdT_collect_cor_data.npz'

T_grids = np.load(flnm)['T_grids']/1.0e3
pdfs    = np.load(flnm)['rho_pdfs']
times   = np.load(flnm)['times']

L_tot = np.trapz(np.nan_to_num(np.mean(pdfs[1:], axis = 0)), np.log10(T_grids[0]))

plt.plot(T_grids[0], np.mean(pdfs[1:], axis = 0)/L_tot, 'b-', label = '$a = 0.9$')

plt.semilogx()

plt.xlim([e_min, e_max])
plt.ylim([0, 1.4])

plt.xlabel(r'$T_e\ \left[\mathrm{keV}\right]$')
plt.ylabel(r'$\partial \mathcal{M}_\mathrm{cor} / \partial \log T_e\ (\mathrm{normalized})$')

ax = plt.gca()
text = r'$\dot{m} = ' + r'{:g}'.format(0.1) + r'$'
ax.text(0.05, 0.8, text, horizontalalignment = 'left', verticalalignment = 'center', transform = ax.transAxes)

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dMdT_10.pdf', bbox_inches = 'tight')

plt.show()