import os

import shutil

import sys

import numpy as np

import h5py

pan_time = '8:00:00'

nodes           = 5
ntasks_per_node = 36
ntasks          = nodes * ntasks_per_node
account         = 'w20_xraysims'
job_name        = 'pan'

module_list = ['python/2.7-anaconda-5.0.1', 'intel/19.0.4', 'mkl/2019.0.4', 'openmpi/2.1.2', 'hdf5-parallel/1.8.16']

prefix = 'mpirun -np ' + repr(ntasks) + ' python '

def make_pan_line(harm_file):
    run_id = (harm_file.split('_')[-1]).split('.')[0]

    line = prefix + 'nr_iterate_corona.py ' + harm_file + ' ' + run_id + ' ' + repr(0) + ' x.h5 normal >pan.out 2>pan.err\n'

    return line

def make_job_script(harm_file):
    run_id = (harm_file.split('_')[-1]).split('.')[0]

    script = '#!/bin/bash' + '\n\n'

    script += '#SBATCH --time='            + pan_time              + '\n'
    script += '#SBATCH --nodes='           + repr(nodes)           + '\n'
    script += '#SBATCH --ntasks='          + repr(ntasks)          + '\n'
    script += '#SBATCH --ntasks-per-node=' + repr(ntasks_per_node) + '\n'
    script += '#SBATCH --account='         + account               + '\n'

    script += '#SBATCH --job-name=' + job_name + '_' + run_id + '\n'

    # these lines may have to be changed depending on the machine
    script += '#SBATCH --qos=standard'             + '\n'
    script += '#SBATCH --partition=standard'       + '\n'
    script += '#SBATCH --mail-user=kinch@lanl.gov' + '\n'   # set to preferred email (please change)
    script += '#SBATCH --mail-type=BEGIN'          + '\n'
    script += '#SBATCH --mail-type=END'            + '\n'
    script += '#SBATCH --mail-type=FAIL'           + '\n'
    script += '#SBATCH --no-requeue'               + '\n'   # if node failure recovery procedure doesn't work, set to --no-requeue
    script += '#SBATCH --signal=23@60'             + '\n\n'
#   script += '#SBATCH --output=%j_panptx.out'     + '\n'
#   script += '#SBATCH --error=%j_panptx.err'      + '\n\n'

    for module in module_list:
        script += 'module load ' + module + '\n'
    script += '\n'

    script += make_pan_line(harm_file)

    sfile = open('pan_' + run_id + '.sh', 'w')

    sfile.write(script + '\n')

    sfile.close()

files = os.listdir('./')
for file in files:
    if 'RADFLUX_panptx' not in file:
        files.remove(file)
files.sort()

for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    print i
    os.mkdir('run_' + repr(i))
    os.chdir('run_' + repr(i))
    os.rename('../RADFLUX_panptx_' + repr(i) + '.h5', './RADFLUX_panptx_' + repr(i) + '.h5')
    os.system('cp -r /turquoise/users/kinch/panptx/* ./')
    os.system('./build_all.sh')
    make_job_script('RADFLUX_panptx_' + repr(i) + '.h5')
    os.chdir('../')

script = '#!/bin/bash' + '\n\n'

for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    script += 'cd run_' + repr(i) + '\n'
    script += 'sbatch pan_' + repr(i) + '.sh' + '\n'
    script += 'cd ..' + '\n'

file = open('submit_all_jobs.sh', 'w')
file.write(script)
file.close()
os.chmod('./submit_all_jobs.sh', 0777)
