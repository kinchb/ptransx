import os

import sys

import numpy as np

import h5py

harm_file   = sys.argv[1]
run_id      = (harm_file.split('_')[-1]).split('.')[0]
tag         = run_id
results_dir = 'panptx_' + tag

f = h5py.File('./' + harm_file, 'r')
Fe_abund = repr(float(np.array(f['/Fe_abund'])))
f.close()

pan_time = '8:00:00'
ptx_time = '16:00:00'
rsp_time = '6:00:00'
rgr_time = '8:00:00'

iterations = 5

nodes           = 5
ntasks_per_node = 36
ntasks          = nodes * ntasks_per_node
account         = 'w21_xraysims'
job_name        = 'ptx_' + tag

pan_module_list = ['python/3.6-anaconda-5.0.1', 'gcc/9.3.0', 'openmpi/2.1.2']
ptx_module_list = ['python/3.6-anaconda-5.0.1', 'intel/19.0.4', 'mkl/2019.0.4', 'openmpi/2.1.2']
rsp_module_list = ['python/3.6-anaconda-5.0.1', 'intel/19.0.4', 'openmpi/2.1.2']
rgr_module_list = ['python/3.6-anaconda-5.0.1', 'intel/19.0.4', 'mkl/2019.0.4', 'openmpi/2.1.2']

prefix = 'mpirun -np ' + repr(ntasks) + ' python '

def make_pan_line(itr):
    if isinstance(itr, int):
        if (itr == 0):
            line = prefix + 'nr_iterate_corona.py ' + harm_file + ' ' + run_id + ' ' + repr(itr) + ' x.h5 normal >pan.out 2>pan.err\n'
        else:
            line = prefix + 'nr_iterate_corona.py ' + harm_file + ' ' + run_id + ' ' + repr(itr) + ' spec_out_' + repr(itr-1) + '.h5 normal >pan.out 2>pan.err\n'
    else:
        line  = prefix + 'nr_iterate_corona.py ' + harm_file + ' ' + run_id + ' ' + repr(iterations+1) + ' spec_out_c_hr.h5 hires >panhr.out 2>panhr.err\n'
        line += 'mv scat_spec.' + run_id + '.' + repr(iterations+1) + '.h5 scat_spec_c_hr.h5\n'
        line += prefix + 'nr_iterate_corona.py ' + harm_file + ' ' + run_id + ' ' + repr(iterations+1) + ' spec_out_l_hr.h5 hires >panhr.out 2>panhr.err\n'
        line += 'mv scat_spec.' + run_id + '.' + repr(iterations+1) + '.h5 scat_spec_l_hr.h5\n'

    return line

def make_ptx_line(itr):
    line = prefix + 'ptransx.py log.out ' + harm_file + ' ' + Fe_abund + ' ' + repr(itr) + ' >ptx.out 2>ptx.err\n'

    return line

def make_rsp_line(itr):
    if isinstance(itr, int):
        line = prefix + 'response.py ' + harm_file + ' run_' + repr(itr) + ' >rsp.out 2>rsp.err\n'
        line += 'mv rsp.out rsp.err ./run_' + repr(itr) + '/\n'
    else:
        line = prefix + 'response.py ' + harm_file + ' run_' + itr + ' >rsp.out 2>rsp.err\n'
        line += 'mv rsp.out rsp.err ./run_' + itr + '/\n'

    return line

def make_rgr_line():
    line = prefix + 'regrid.py run_' + repr(iterations-1) + ' run_hr ' + harm_file + ' ' + Fe_abund + ' >rgr.out 2>rgr.err\n'

    return line

def make_job_script(type, itr = 0):
    script = '#!/bin/bash' + '\n\n'

    if (type == 'pan'):
        script += '#SBATCH --time=' + pan_time + '\n'
    if (type == 'ptx'):
        script += '#SBATCH --time=' + ptx_time + '\n'
    if (type == 'rsp'):
        script += '#SBATCH --time=' + rsp_time + '\n'
    if (type == 'rgr'):
        script += '#SBATCH --time=' + rgr_time + '\n'
    script += '#SBATCH --nodes='           + repr(nodes)           + '\n'
    script += '#SBATCH --ntasks='          + repr(ntasks)          + '\n'
    script += '#SBATCH --ntasks-per-node=' + repr(ntasks_per_node) + '\n'
    script += '#SBATCH --account='         + account               + '\n'
    if isinstance(itr, int):
        script += '#SBATCH --job-name=' + job_name + '_' + type + '_' + repr(itr) + '\n'
    else:
        script += '#SBATCH --job-name=' + job_name + '_' + type + '_' + itr + '\n'
    # these lines may have to be changed depending on the machine
    script += '#SBATCH --qos=standard'             + '\n'
    script += '#SBATCH --partition=standard'       + '\n'
    script += '#SBATCH --mail-user=kinch@lanl.gov' + '\n'   # set to preferred email (please change)
    script += '#SBATCH --mail-type=BEGIN'          + '\n'
    script += '#SBATCH --mail-type=END'            + '\n'
    script += '#SBATCH --mail-type=FAIL'           + '\n'
    script += '#SBATCH --no-requeue'               + '\n'   # if node failure recovery procedure doesn't work, set to --no-requeue
    script += '#SBATCH --signal=23@60'             + '\n\n'

    script += 'module purge\n'
    if (type == 'pan'):
        module_list = pan_module_list
    if (type == 'ptx'):
        module_list = ptx_module_list
    if (type =='rsp'):
        module_list = rsp_module_list
    if (type == 'rgr'):
        module_list = rgr_module_list
    for module in module_list:
        script += 'module load ' + module + '\n'
    script += '\n'

    script += 'source /turquoise/users/kinch/virtual_python/bin/activate\n\n'

    if (type == 'pan'):
        script += make_pan_line(itr)
    if (type == 'ptx'):
        script += make_ptx_line(itr)
    if (type == 'rsp'):
        script += make_rsp_line(itr)
    if (type == 'rgr'):
        script += make_rgr_line()

    if isinstance(itr, int):
        sfile = open(type + '_' + repr(itr) + '.sh', 'w')
    else:
        sfile = open(type + '_' + itr + '.sh', 'w')

    sfile.write(script + '\n')

    sfile.close()

for itr in range(0, iterations):
    make_job_script('pan', itr)
    make_job_script('ptx', itr)
    make_job_script('rsp', itr)
make_job_script('pan', itr+1)
make_job_script('rgr')
make_job_script('rsp', 'hr')
make_job_script('pan', 'hr')

script = '#!/bin/bash' + '\n\n'

script += 'jid0=$(sbatch pan_0.sh | grep -o -E \'[0-9]+\')' + '\n'
script += 'jid1=$(sbatch --dependency=afterok:$jid0 ptx_0.sh | grep -o -E \'[0-9]+\')' + '\n'
script += 'jid2=$(sbatch --dependency=afterany:$jid1 ptx_0.sh | grep -o -E \'[0-9]+\')' + '\n'
script += 'jid3=$(sbatch --dependency=afterany:$jid2 ptx_0.sh | grep -o -E \'[0-9]+\')' + '\n'
script += 'jid4=$(sbatch --dependency=afterany:$jid3 ptx_0.sh | grep -o -E \'[0-9]+\')' + '\n'
script += 'jid5=$(sbatch --dependency=afterok:$jid4 rsp_0.sh | grep -o -E \'[0-9]+\')' + '\n'
j = 6
for i in range(1, iterations):
    script += 'jid' + repr(j) + '=$(sbatch --dependency=afterany:$jid' + repr(j-1) + ' pan_' + repr(i) + '.sh | grep -o -E \'[0-9]+\')' + '\n'
    j += 1
    script += 'jid' + repr(j) + '=$(sbatch --dependency=afterok:$jid' + repr(j-1) + ' ptx_' + repr(i) + '.sh | grep -o -E \'[0-9]+\')' + '\n'
    j += 1
    script += 'jid' + repr(j) + '=$(sbatch --dependency=afterok:$jid' + repr(j-1) + ' rsp_' + repr(i) + '.sh | grep -o -E \'[0-9]+\')' + '\n'
    j += 1
script += 'jid' + repr(j) + '=$(sbatch --dependency=afterany:$jid' + repr(j-1) + ' pan_' + repr(i+1) + '.sh | grep -o -E \'[0-9]+\')' + '\n'
j += 1
script += 'jid' + repr(j) + '=$(sbatch --dependency=afterok:$jid' + repr(j-1) + ' rgr_0.sh | grep -o -E \'[0-9]+\')' + '\n'
j += 1
script += 'jid' + repr(j) + '=$(sbatch --dependency=afterany:$jid' + repr(j-1) + ' rsp_hr.sh | grep -o -E \'[0-9]+\')' + '\n'
j += 1
script += 'jid' + repr(j) + '=$(sbatch --dependency=afterany:$jid' + repr(j-1) + ' pan_hr.sh | grep -o -E \'[0-9]+\')' + '\n'

script += '\n'

file = open('submit_all_jobs.sh', 'w')
file.write(script)
file.close()
os.chmod('./submit_all_jobs.sh', 777)
