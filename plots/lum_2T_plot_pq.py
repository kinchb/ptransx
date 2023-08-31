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

    ndx = 50
    print(e_grid[ndx]/1000.0)

    return np.trapz(lum[ndx:] * 2.417989e14, e_grid[ndx:])

t = []
L = []

times = np.load('../panptx_data/h3d_10_a0_01/h3d_10_a0_01_times.npz')['times']
for i in range(50, 1001, 50):
    t.append(times[i])
    L.append(get_lum('../panptx_data/h3d2T_10_a0_01/scat_spec.' + repr(i) + '.0.h5'))

t_rt_2T = np.array(t)
L_rt_2T = np.array(L)

t = []
L = []

times = np.load('../panptx_data/h3d_10_a0_01/h3d_10_a0_01_times.npz')['times']
for i in range(50, 2001, 50):
    t.append(times[i])
    L.append(get_lum('../panptx_data/h3d_10_a0_01/scat_spec.' + repr(i) + '.0.h5'))

t_rt_1T = np.array(t)
L_rt_1T = np.array(L)

print(np.mean(L_rt_2T))

L_edd = (1.26e38) * M

flnm = '../../harm3d/harm_data/lum_data_ucon0_2T_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

disk_frac = (disk_lum/total_lum).copy()

t_2T = t.copy()

t_adj = 10000.0

plt.plot(t - t_adj, total_lum, 'k-', label = r'$\mathrm{total\ (2T)}$')
# plt.plot(t_rt_2T - t_adj, L_rt_2T/L_edd, 'm-', label = r'$\mathrm{ray}$-$\mathrm{traced}\ (2T)$')
plt.plot(t - t_adj, cor_lum, 'r-', label = r'$\mathrm{corona\ (2T)}$')
# plt.plot(t - t_adj, disk_lum, 'b-', label = r'$\mathrm{disk}$')

cor_lum_2T  = cor_lum.copy()
disk_lum_2T = disk_lum.copy()

flnm = '../../harm3d/harm_data/lum_data_ucon0_a0.0.npz'

t         = np.load(flnm)['t']
total_lum = np.load(flnm)['total_lum']
cor_lum   = np.load(flnm)['cor_lum']
disk_lum  = np.load(flnm)['disk_lum']

t_1T = t.copy()

plt.plot(t - t_adj, total_lum, 'k--', label = r'$\mathrm{total\ (1T)}$')
# plt.plot(t_rt_1T - t_adj, L_rt_1T/L_edd, 'm-', label = r'$\mathrm{ray}$-$\mathrm{traced}\ (2T)$')
plt.plot(t - t_adj, cor_lum, 'r--', label = r'$\mathrm{corona\ (1T)}$')
# plt.plot(t - t_adj, disk_lum, 'b-', label = r'$\mathrm{disk\ (1T)}$')

cor_lum_1T  = cor_lum.copy()
disk_lum_1T = disk_lum.copy()
# plt.semilogy()

plt.xlim([0, 1000])
plt.ylim([0.0, 0.06])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L/L_\mathrm{Edd}$')

plt.legend(frameon = False, loc = 'upper right')

plt.tight_layout()

f = plt.gcf()
f.savefig('cooling_w_rt_2T.pdf', bbox_inches = 'tight')

plt.show()

def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

def make_pspec(data, start=0, end=-1):
    pspec = np.square(np.abs(np.fft.rfft(data[start:end])))
    time_step = 10.0 * 4.9e-6
    rate = 1.0/time_step
    freq = np.linspace(0.0, rate/2, len(pspec))
    return (freq, pspec)

freq, pspec = make_pspec(cor_lum_1T, start = 100, end = 1000)
# freq, pspec = make_pspec(cor_lum_1T/cor_lum_1T[100:1000].mean(), start = 100, end = 1000)
plt.plot(freq, pspec, 'r-', label = r'$\mathrm{1T}$')

freq, pspec = make_pspec(cor_lum_2T, start = 100, end = 1000)
# freq, pspec = make_pspec(cor_lum_2T/cor_lum_2T[100:1000].mean(), start = 100, end = 1000)
plt.plot(freq, pspec, 'k-', label = r'$\mathrm{2T}$')

plt.xlim([1.0e1, 1.1e4])

plt.xlabel(r'$\mathrm{frequency}\ [\mathrm{Hz}]$')
plt.ylabel(r'$L_\mathrm{IC}\ \mathrm{power\ spectrum}\ \left[\mathrm{Hz}^{-1}\right]$')

plt.loglog()

plt.legend(frameon = False, loc = 'best')

plt.tight_layout()

f = plt.gcf()
f.savefig('corona_power_spec.pdf', bbox_inches = 'tight')

plt.show()

N = 50

cor_lum_2T_ma = moving_average(cor_lum_2T, n = N)
cor_lum_1T_ma = moving_average(cor_lum_1T, n = N)

plt.plot([-1000, 2000], [1.0, 1.0], 'k--')
plt.plot(t_2T[N-1:] - t_adj, cor_lum_2T[N-1:]/cor_lum_2T_ma, 'k-', label = r'$\mathrm{2T}$')
plt.plot(t_1T[N-1:] - t_adj, cor_lum_1T[N-1:]/cor_lum_1T_ma, 'r-', label = r'$\mathrm{1T}$')

"""
N = 3

lum_2T_ma = moving_average(L_rt_2T, n = N)
lum_1T_ma = moving_average(L_rt_1T, n = N)

lum_1T_ma = np.mean(L_rt_1T)
lum_2T_ma = np.mean(L_rt_2T)

plt.plot([-1000, 2000], [1.0, 1.0], 'k--')
plt.plot(t_rt_2T[N-1:] - t_adj, L_rt_2T[N-1:]/lum_2T_ma, 'k-', label = r'$\mathrm{2T}$')
plt.plot(t_rt_1T[N-1:] - t_adj, L_rt_1T[N-1:]/lum_1T_ma, 'r-', label = r'$\mathrm{1T}$')
"""

cor_2T = (cor_lum_2T[N-1:]/cor_lum_2T_ma)[100:1000]
cor_1T = (cor_lum_1T[N-1:]/cor_lum_1T_ma)[100:1000]

print(cor_2T.std())
print(cor_1T.std())

plt.xlim([0, 1000])
plt.ylim([0.9, 1.1])

plt.xlabel(r'$t/M$')
plt.ylabel(r'$L_\mathrm{IC}/\langle L_\mathrm{IC} \rangle$')

plt.legend(frameon = False, loc = 'upper right')

plt.tight_layout()

f = plt.gcf()
f.savefig('cooling_variation.pdf', bbox_inches = 'tight')

plt.show()
