import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

# simulation parameters
# --- --- --- --- --- --- --- --- --- ---
M    = 10.0
a    = 0.0
mdot = 0.01
# --- --- --- --- --- --- --- --- --- ---

# useful constants (cgs)
# --- --- --- --- --- --- --- --- --- ---
m_bh   = M * 2.0e33
G      = 6.6726e-8
c      = 3.0e10
kappa  = 0.4
# --- --- --- --- --- --- --- --- --- ---

L_edd = (1.26e38) * M

if (a == 0.0):
    eta = 0.0572
if (a == 0.5):
    eta = 0.0821
if (a == 0.9):
    eta = 0.1558
if (a == 0.99):
    eta = 0.2640

Mdot_edd = (1.26e38 * M)/(eta * c**2)

flnm = '../harm_data/mdot_data_past_a0.0.npz'

t    = np.load(flnm)['t']
r    = np.load(flnm)['r']
Mdot = np.load(flnm)['Mdot']

flnm = '../harm_data/old_data/heat_flux_data_a0.0.npz'

heat_flux = np.load(flnm)['heat_flux']
# heat_flux = np.load(flnm)['disk_flux']

flnm = '../harm_data/lum_data_past_a0.0.npz'

total_lum = np.load(flnm)['total_lum'][89:]
cor_lum   = np.load(flnm)['cor_lum'][89:]
disk_lum  = np.load(flnm)['disk_lum'][89:]

t -= t[0]
t -= -1000.0

print(cor_lum[-1]/disk_lum[-1])

Mdot      = np.mean(Mdot, axis = 0)
heat_flux_o1 = np.mean(heat_flux[:50,::], axis = 0)
heat_flux_o2 = np.mean(heat_flux[50:,::], axis = 0)
total_lum = np.mean(total_lum)
cor_lum   = np.mean(cor_lum)

rh_ndx = 0
while True:
    if (r[rh_ndx] > 2.0):
        break
    rh_ndx += 1

Mdot_past      = Mdot.copy()
heat_flux_past = heat_flux.copy()
total_lum_past = L_edd * total_lum
cor_lum_past   = L_edd * cor_lum

flnm = '../harm_data/mdot_data_a0.0.npz'

t    = np.load(flnm)['t']
r    = np.load(flnm)['r']
Mdot = np.load(flnm)['Mdot']

flnm = '../harm_data/heat_flux_data_a0.0.npz'

heat_flux = np.load(flnm)['heat_flux']
# heat_flux = np.load(flnm)['disk_flux']

flnm = '../harm_data/lum_data_a0.0.npz'

total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t -= t[0]

t_1000_ndx = 0
i = 0
while (i < len(t)):
    if (t[i] > 1000.0):
        break
    else:
        i += 1
t_1000_ndx = i+1

t_2000_ndx = 0
while (i < len(t)):
    if (t[i] > 2000.0):
        break
    else:
        i += 1
t_2000_ndx = i+1

print(cor_lum[2000]/disk_lum[2000])

Mdot_1      = np.mean(Mdot[:t_1000_ndx,::], axis = 0)
Mdot_2      = np.mean(Mdot[t_1000_ndx-1:t_2000_ndx,::], axis = 0)
heat_flux_1 = np.mean(heat_flux[:t_1000_ndx,::], axis = 0)
heat_flux_2 = np.mean(heat_flux[t_1000_ndx-1:t_2000_ndx,::], axis = 0)
total_lum_1 = L_edd * np.mean(total_lum[:t_1000_ndx])
total_lum_2 = L_edd * np.mean(total_lum[t_1000_ndx-1:t_2000_ndx])
cor_lum_1   = L_edd * np.mean(cor_lum[:t_1000_ndx])
cor_lum_2   = L_edd * np.mean(cor_lum[t_1000_ndx-1:t_2000_ndx])

plt.plot([0.0, 100.0], [0.0, 0.0], 'k--')

# plt.plot(r, heat_flux_past/L_edd, 'k-', label = r'$-1000$--$0 M$')
plt.plot(r, heat_flux_1/L_edd,    'r-', label = r'$0$--$1000 M$')
plt.plot(r, heat_flux_2/L_edd,    'b-', label = r'$1000$--$2000 M$')

plt.plot(r, heat_flux_o1/L_edd,    'r--', label = r'$0$--$1000 M$')
plt.plot(r, heat_flux_o2/L_edd,    'b--', label = r'$1000$--$2000 M$')


# plt.plot(r, (Mdot_edd/L_edd) * heat_flux_past/Mdot_past, 'k-', label = r'$-1000$--$0 M$')
# plt.plot(r, (Mdot_edd/L_edd) * heat_flux_1/Mdot_1,       'r-', label = r'$0$--$1000 M$')
# plt.plot(r, (Mdot_edd/L_edd) * heat_flux_2/Mdot_2,       'b-', label = r'$1000$--$2000 M$')

plt.xlim([2, 70])
plt.ylim([-0.0005, 0.003])
# plt.ylim([-0.05, 0.2])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))

plt.xlabel(r'$r/M$')
plt.ylabel(r'$\dot{U}(r)/L_\mathrm{Edd}$')
# plt.ylabel(r'$(\dot{U}(r)/L_\mathrm{Edd})/(\dot{M}(r)/\dot{M}_\mathrm{Edd})$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('heat_flux_disk.pdf', bbox_inches = 'tight')

plt.show()

"""
net_heat_rate_past = -3.902511454332309e+34
net_heat_rate_1    = -1.9524814993694006e+35

heat_rate_past = net_heat_rate_past + (heat_flux_past[rh_ndx] - heat_flux_past[-1] + total_lum_past)
heat_rate_1    = net_heat_rate_1    + (heat_flux_1[rh_ndx]    - heat_flux_1[-1]    + total_lum_1)

(heat_flux_past[rh_ndx] - heat_flux_past[-1])/(c**2 * Mdot_past[rh_ndx])
(heat_flux_1[rh_ndx]    - heat_flux_1[-1])/(c**2 * Mdot_1[rh_ndx])

(total_lum_past)/(c**2 * Mdot_past[rh_ndx])
(total_lum_1)/(c**2 * Mdot_1[rh_ndx])

mod_eta_past = (heat_flux_past[rh_ndx] + heat_flux_past[-1] + total_lum_past)/(c**2 * Mdot_past[rh_ndx])
mod_eta_1    = (heat_flux_1[rh_ndx]    + heat_flux_1[-1]    + total_lum_1)/(c**2 * Mdot_1[rh_ndx])
mod_eta_2    = (heat_flux_2[rh_ndx]    + total_lum_2)/(c**2 * Mdot_2[rh_ndx])
"""
