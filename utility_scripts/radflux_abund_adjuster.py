import sys

import os

import h5py

import numpy as np

new_Fe_abund = float(sys.argv[1])

files = os.listdir('./')
for file in files:
    if 'RADFLUX_panptx' not in file:
        files.remove(file)
files.sort()

for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    print(i)
    f = h5py.File('./' + file, 'a')
    del f['/Fe_abund']
    f['/Fe_abund'] = np.array(new_Fe_abund)
    f.close()
