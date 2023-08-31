import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

thetas   = np.load('../harm_data/Te_errors_polar_rho_wght_2-10.npz')['thetas']
errors_1 = np.load('../harm_data/Te_errors_polar_rho_wght_2-10.npz')['errors']
errors_2 = np.load('../harm_data/Te_errors_polar_rho_wght_10-20.npz')['errors']
errors_3 = np.load('../harm_data/Te_errors_polar_rho_wght_20-40.npz')['errors']
errors_4 = np.load('../harm_data/Te_errors_polar_rho_wght_40-65.npz')['errors']

errors_1 = np.ma.masked_where(errors_1 == 0.0, errors_1)
errors_2 = np.ma.masked_where(errors_2 == 0.0, errors_2)
errors_3 = np.ma.masked_where(errors_3 == 0.0, errors_3)
errors_4 = np.ma.masked_where(errors_4 == 0.0, errors_4)

plt.plot([-1.0, 181.0], [1.0, 1.0], 'k--')
plt.plot([-1.0, 181.0], [2.0, 2.0], 'r--')
plt.plot([-1.0, 181.0], [0.5, 0.5], 'r--')

"""
plt.plot(180 * (thetas/np.pi), errors_1, 'g-', label = r'$r = 2$--$10 M$')
plt.plot(180 * (thetas/np.pi), errors_2, 'b-', label = r'$r = 10$--$20 M$')
plt.plot(180 * (thetas/np.pi), errors_3, 'c-', label = r'$r = 20$--$40 M$')
plt.plot(180 * (thetas/np.pi), errors_4, 'm-', label = r'$r = 40$--$65 M$')
"""
plt.plot(180 * (thetas/np.pi), errors_1, 'g-', label = r'$2 < r/M \leq 10$')
plt.plot(180 * (thetas/np.pi), errors_2, 'b-', label = r'$10 < r/M \leq 20$')
plt.plot(180 * (thetas/np.pi), errors_3, 'c-', label = r'$20 < r/M \leq 40$')
plt.plot(180 * (thetas/np.pi), errors_4, 'm-', label = r'$40 < r/M \leq 65$')

plt.xlim([0.0, 180])
plt.ylim([0.4, 40.0])
# plt.ylim([0.4, 100.0])

plt.semilogy()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((0, 30, 60, 90, 120, 150, 180))
ax.set_xticklabels(('$0^\circ$', '$30^\circ$', '$60^\circ$', '$90^\circ$', '$120^\circ$', '$150^\circ$', '$180^\circ$'))

plt.xlabel(r'$\theta$')
plt.ylabel(r'$\langle T_e^{\textsc{harm}} / T_e^{\textsc{pandurata}} \rangle$')

plt.legend(frameon = False, loc = 'best')

ax.legend(bbox_to_anchor = (0.50, 0.60), frameon = False)

plt.tight_layout()

f = plt.gcf()
f.savefig('Te_errors_polar.pdf', bbox_inches = 'tight')

plt.show()
