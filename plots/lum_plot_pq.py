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

M = 10.0

def get_lum(flnm):
    infile = h5py.File(flnm, 'r')

    spec = np.array(infile['/data']).T

    infile.close()

    lum = np.zeros(len(spec[0]))

    for i in range(0, len(spec[0])):
        for j in range(0, len(spec)):
            lum[i] += spec[j][i]

    e_grid = np.logspace(0, 7, len(spec[0]), endpoint = True)

    return np.trapz(lum * 2.417989e14, e_grid)

t = []
L = []

times = np.load('../panptx_data/h3d_10_a0_01_old/old_a0_data_times.npz')['times']
for i in range(0, 51, 1):
#   t.append(-1000.0 + 20*float(i))
    t.append(times[i])
    L.append(get_lum('../panptx_data/h3d_10_a0_01_old/scat_spec.' + repr(i) + '.0.h5'))

t_rt_old = np.array(t)
L_rt_old = np.array(L)

print(np.mean(L_rt_old))

t = []
L = []

times = np.load('../panptx_data/h3d_10_a0_01/h3d_10_a0_01_times.npz')['times']
for i in range(50, 2001, 50):
    t.append(times[i])
    L.append(get_lum('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5'))

t_rt = np.array(t)
L_rt = np.array(L)

print(np.mean(L_rt[:int(len(t_rt)/2)]))
print(np.mean(L_rt[int(len(t_rt)/2):]))

L_edd = (1.26e38) * M

flnm = '../../harm3d/harm_data/lum_data_ucon0_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

disk_frac = (disk_lum/total_lum).copy()

plt.plot([0.0, 0.0], [-1.0, 1.0], 'k--')

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k-', label = r'$\mathrm{total}$')

plt.plot(t_rt_old - t_adj, L_rt_old/L_edd, 'm-')
plt.plot(t_rt - t_adj, L_rt/L_edd, 'm-', label = r'$\mathrm{ray}$-$\mathrm{traced}$')

plt.plot(t - t_adj, cor_lum,   'r-', label = r'$\mathrm{corona}$')
plt.plot(t - t_adj, disk_lum,  'b-', label = r'$\mathrm{disk}$')

flnm = '../../harm3d/harm_data/old_data/lum_data_ucon0_past_a0.0_1.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

disk_frac = (disk_lum/total_lum).copy()

print(disk_frac)

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k-')
plt.plot(t - t_adj, cor_lum,   'r-')
plt.plot(t - t_adj, disk_lum,  'b-')

flnm = '../../harm3d/harm_data/old_data/lum_data_ucon0_past_a0.0_2.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

disk_frac = (disk_lum/total_lum).copy()

print(disk_frac)

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k--')
plt.plot(t - t_adj, cor_lum,   'r--')
plt.plot(t - t_adj, disk_lum,  'b--')

flnm = '../../harm3d/harm_data/old_data/lum_data_ucon0_past_a0.0_3.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

disk_frac = (disk_lum/total_lum).copy()

print(disk_frac)

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k--')
plt.plot(t - t_adj, cor_lum,   'r--')
plt.plot(t - t_adj, disk_lum,  'b--')

print(L_edd * np.mean(total_lum))

# plt.semilogy()

plt.xlim([-1000, 2000])
# plt.ylim([5.0e-3, 1.0e-1])
plt.ylim([0.0, 0.02])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L/L_\mathrm{Edd}$')

plt.legend(bbox_to_anchor = (1.0, 1.0), frameon = False, loc = 'upper right', ncol = 2, handletextpad = 0.5, columnspacing = 1)

plt.tight_layout()

f = plt.gcf()
f.savefig('cooling_w_rt.pdf', bbox_inches = 'tight')

plt.show()

print(old_disk_frac.mean())
plt.plot(old_disk_frac)
plt.show()

print(disk_frac.mean())
plt.plot(disk_frac)
plt.show()
