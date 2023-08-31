import numpy as np

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/pdflatex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

flnm = '../harm_data/Te_and_Ti_data_0.0.npz'

t   = np.load(flnm)['t']
T_e = np.load(flnm)['T_e_pw']
T_i = np.load(flnm)['T_i_pw']

plt.plot(t - t[0], T_e/T_i, 'k-')
# plt.semilogy()
plt.show()

from scipy.fft import fft

T_e_f = fft(T_e)

plt.plot(T_e_f)

plt.show()
