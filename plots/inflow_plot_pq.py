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

t -= t[0]
t -= -1000.0

sigma = np.std(Mdot, axis = 0)
Mdot  = np.mean(Mdot, axis = 0)

rh_ndx = 0
while True:
    if (r[rh_ndx] > 2.0):
        break
    rh_ndx += 1

print('old data: ' + repr(Mdot[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma[rh_ndx]/Mdot_edd))

plt.plot(r, Mdot/Mdot_edd, 'k-', label = r'$-1000$--$0 M$')
plt.plot(r, (Mdot + sigma)/Mdot_edd, 'k-', linewidth = 0.5)
plt.plot(r, (Mdot - sigma)/Mdot_edd, 'k-', linewidth = 0.5)

# --- --- ---

flnm = '../harm_data/mdot_data_a0.0.npz'

t    = np.load(flnm)['t']
r    = np.load(flnm)['r']
Mdot = np.load(flnm)['Mdot']

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

print(t[t_1000_ndx])
print(t[t_2000_ndx])

sigma_1 = np.std(Mdot[:t_1000_ndx,::], axis = 0)
sigma_2 = np.std(Mdot[t_1000_ndx-1:t_2000_ndx,::], axis = 0)

Mdot_1 = np.mean(Mdot[:t_1000_ndx,::], axis = 0)
Mdot_2 = np.mean(Mdot[t_1000_ndx-1:t_2000_ndx,::], axis = 0)

plt.plot([0.0, 100.0], [0.0, 0.0], 'k--')

print('IC (1--1000):    ' + repr(Mdot_1[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma_1[rh_ndx]/Mdot_edd))
print('IC (1000--2000): ' + repr(Mdot_2[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma_2[rh_ndx]/Mdot_edd))

plt.plot(r, Mdot_1/Mdot_edd, 'r-', label = r'$0$--$1000 M$')
plt.plot(r, (Mdot_1 + sigma)/Mdot_edd, 'r-', linewidth = 0.5)
plt.plot(r, (Mdot_1 - sigma)/Mdot_edd, 'r-', linewidth = 0.5)

plt.plot(r, Mdot_2/Mdot_edd, 'b-', label = r'$1000$--$2000 M$')
plt.plot(r, ((Mdot_2 + sigma)/Mdot_edd), 'b-', linewidth = 0.5)
plt.plot(r, ((Mdot_2 - sigma)/Mdot_edd), 'b-', linewidth = 0.5)

print(Mdot_1[rh_ndx])
print(Mdot_2[rh_ndx])

# --- --- ---

flnm = '../harm_data/old_data/mdot_data_a0.0.npz'

t    = np.load(flnm)['t']
r    = np.load(flnm)['r']
Mdot = np.load(flnm)['Mdot']

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
t_2000_ndx = i-1

print(t[t_1000_ndx])
print(t[t_2000_ndx])

sigma_1 = np.std(Mdot[:t_1000_ndx,::], axis = 0)
sigma_2 = np.std(Mdot[t_1000_ndx-1:t_2000_ndx,::], axis = 0)

Mdot_1 = np.mean(Mdot[:t_1000_ndx,::], axis = 0)
Mdot_2 = np.mean(Mdot[t_1000_ndx-1:t_2000_ndx,::], axis = 0)

plt.plot([0.0, 100.0], [0.0, 0.0], 'k--')

plt.plot(r, Mdot_1/Mdot_edd, 'r--') #, label = r'$0$--$1000 M$')
plt.plot(r, Mdot_2/Mdot_edd, 'b--') #, label = r'$1000$--$2000 M$')

print('TT (1--1000):    ' + repr(Mdot_1[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma_1[rh_ndx]/Mdot_edd))
print('TT (1000--2000): ' + repr(Mdot_2[rh_ndx]/Mdot_edd) + ' +/- ' + repr(sigma_2[rh_ndx]/Mdot_edd))

print(Mdot_1[rh_ndx])
print(Mdot_2[rh_ndx])

# --- --- ---

"""
plt.xlim([2, 70])
plt.ylim([-0.02, 0.04])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 8, 10, 15, 20, 30, 40, 50, 70))
ax.set_xticklabels(('2', '3', '4', '5', '6', '8', '10', '15', '20', '30', '40', '50', '70'))
"""

plt.xlim([2, 10])
plt.ylim([0.00, 0.04])

plt.semilogx()

ax = plt.gca()

ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax.set_xticks((2, 3, 4, 5, 6, 7, 8, 9, 10))
ax.set_xticklabels(('2', '3', '4', '5', '6', '7', '8', '9', '10'))

plt.xlabel(r'$r/M$')
plt.ylabel(r'$\dot{M}(r)/\dot{M}_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('inflow.pdf', bbox_inches = 'tight')

plt.show()
