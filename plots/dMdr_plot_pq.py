import os

import numpy as np

import h5py

import matplotlib.pyplot as plt

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc

import sys
sys.path.append('/usr/bin/latex')
rc('text', usetex = True)
rc('font', family = 'serif', size = 16)

dir = '/mnt/archive/h3d_10_a0_01/dumps'

harm_data_1 = h5py.File(dir + '/KDHARM0.000500.h5', 'r')

x1 = np.array(harm_data_1['x1'])
x2 = np.array(harm_data_1['x2'])
x3 = np.array(harm_data_1['x3'])

harm_data_1.close()

r   = np.zeros(len(x1))
th  = np.zeros(len(x2[0]))
phi = np.zeros(len(x3[0][0]))

for i in range(0, len(x1)):
    r[i]   = x1[i][0][0]
for i in range(0, len(x2[0])):
    th[i]  = x2[0][i][0]
for i in range(0, len(x3[0][0])):
    phi[i] = x3[0][0][i]

filelist = []

i = 0
while (i < 10000):
    filename = dir + '/KDHARM0.RADFLUX.01' + '{:04}'.format(i) + '.h5'
    if os.path.isfile(filename):
        print filename
        filelist.append(filename)
    i += 1

rho = np.zeros((len(r), len(th), len(phi)))

j = 0
for i in range(0, len(filelist)):
    print repr(i) + '/' + repr(len(filelist))
    harm_data_2 = h5py.File(filelist[i], 'r')
    rho += np.array(harm_data_2['rho'])
    j += 1
    harm_data_2.close()

rho /= j

