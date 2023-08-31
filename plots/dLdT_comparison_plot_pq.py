import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

M = 10.0

L_edd = (1.26e38) * M

T_grids = np.load('../harm_scripts/dLdT_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_scripts/dLdT_data_a0.0.npz')['pdfs']
times   = np.load('../harm_scripts/dLdT_data_a0.0.npz')['times']

# plt.plot(T_grids[2], pdfs[2]/L_edd, 'r-', label = '$100M$')

# plt.plot(T_grids[10], pdfs[10]/L_edd, 'g-', label = '$500M$')

plt.plot(T_grids[20], pdfs[20]/L_edd, 'r-', label = r'\textsc{harm3d}')

# plt.plot(T_grids[40], pdfs[40]/L_edd, 'c-', label = '$2000M$')

print times[20]

T_grids = np.load('../harm_scripts/dLdT_pan_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_scripts/dLdT_pan_data_a0.0.npz')['pdfs']
times   = np.load('../harm_scripts/dLdT_pan_data_a0.0.npz')['times']

# sys.exit()

plt.plot(T_grids[4], pdfs[4]/L_edd, 'b-', label = r'\textsc{pandurata}')

# plt.plot(T_grids[-1], pdfs[-1]/L_edd, 'b--', label = '$1000M$')

print times[4]

plt.semilogx()

plt.xlim([1.0e3, 1.0e8])
plt.ylim([0, 0.005])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial \left(L_\mathrm{IC}/L_\mathrm{Edd}\right) / \partial \log T_e$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_comparison.pdf', bbox_inches = 'tight')

plt.show()

# --- --- --- --- ---
"""
M = 10.0

L_edd = (1.26e38) * M

T_grids = np.load('../harm_scripts/dLdT_past_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_scripts/dLdT_past_data_a0.0.npz')['pdfs']
times   = np.load('../harm_scripts/dLdT_past_data_a0.0.npz')['times']

# plt.plot(T_grids[0], pdfs[0]/L_edd, 'k-', label = '$-1000M$')
plt.plot(T_grids[-1], pdfs[-1]/np.trapz(pdfs[-1][1:], np.log10(T_grids[-1][1:])), 'k-', label = '$0M$')

print times[0]
print times[-1]

T_grids = np.load('../harm_scripts/dLdT_data_a0.0.npz')['T_grids']
pdfs    = np.load('../harm_scripts/dLdT_data_a0.0.npz')['pdfs']
times   = np.load('../harm_scripts/dLdT_data_a0.0.npz')['times']

plt.plot(T_grids[2], pdfs[2]/np.trapz(pdfs[2][1:], np.log10(T_grids[2][1:])), 'r-', label = '$100M$')

plt.plot(T_grids[10], pdfs[10]/np.trapz(pdfs[10][1:], np.log10(T_grids[10][1:])), 'g-', label = '$500M$')

plt.plot(T_grids[20], pdfs[20]/np.trapz(pdfs[20][1:], np.log10(T_grids[20][1:])), 'b-', label = '$1000M$')

plt.plot(T_grids[40], pdfs[40]/np.trapz(pdfs[40][1:], np.log10(T_grids[40][1:])), 'c-', label = '$2000M$')

plt.semilogx()

plt.xlim([1.0e3, 1.0e8])
plt.ylim([0, 0.8])

plt.xlabel(r'$T_e\ \left[\mathrm{eV}\right]$')
plt.ylabel(r'$\partial L_\mathrm{IC} / \partial \log T_e\ (\mathrm{normalized})$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('dLdT_norm.pdf', bbox_inches = 'tight')

plt.show()
"""

