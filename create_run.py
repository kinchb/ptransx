import sys
import os
import shutil
import glob
import json
import time

def copy_ptransx(src_dir, target_dir):
    exclude = ['compy']
    os.chdir(target_dir)
    for entry in os.listdir(src_dir):
        if entry not in exclude:
            os.system('cp -r ' + src_dir + '/' + entry + ' ./')
    os.chdir('..')

def build_fsolver(config, rebuild = False):
    os.chdir('ptxlib/fsolver')
    if rebuild:
        os.system('rm -rf build')
        os.system('rm _fsolver.*.so')
    build_command = 'module purge\n'
    for key in config['module list']:
        build_command += 'module load ' + config['module list'][key] + '\n'
    build_command += config['fsolver build line'] + '\n'
    os.system(build_command)
    os.chdir('..')
    os.chdir('..')

def build_response(config, rebuild = False):
    os.chdir('ptxlib/response')
    if rebuild:
        os.system('rm -rf build')
        os.system('rm _response.*.so')
    build_command = 'module purge\n'
    for key in config['module list']:
        build_command += 'module load ' + config['module list'][key] + '\n'
    build_command += config['response build line'] + '\n'
    os.system(build_command)
    os.chdir('..')
    os.chdir('..')

def build_pandurata(config, rebuild = False):
    os.chdir('ptxlib/pandurata')
    if rebuild:
        os.system('rm -rf build')
        os.system('rm _pandurata.*.so')
        os.system('rm ./*.o')
    os.system('dos2unix panhead.h')
    build_command = 'module purge\n'
    for key in config['module list']:
        build_command += 'module load ' + config['module list'][key] + '\n'
    build_command += config['pan file build line'] + ' '
    for file in glob.glob('./*.c'):
        if file not in ['./_pandurata.c', './pandurata.c']:
            os.system(build_command + file + '\n')
    os.system(build_command + '\n' + config['pandurata build line'] + '\n')
    os.chdir('..')
    os.chdir('..')

def build_xstarcomm(config, rebuild = False):
    os.chdir('ptxlib/xstarcomm')
    if rebuild:
        os.system('rm -rf build')
        os.system('rm _xstarcomm.*.so')
        os.system('rm ./*.o')
    build_command = 'module purge\n'
    for key in config['module list']:
        build_command += 'module load ' + config['module list'][key] + '\n'
    os.system(build_command + config['xstarsub build line'] + '\n')
    os.system(build_command + config['xstarcomm build line'] + '\n')
    os.chdir('..')
    os.chdir('..')

def make_test(config):
    test = '#!/bin/bash\n\n'
    for key in config['module list']:
        test += 'module load ' + config['module list'][key] + '\n'
    test += '\n'
    test += 'python3 test.py\n'
    with open('test.sh', 'w') as f:
        f.write(test)
    os.system('chmod a+x test.sh')

# load the JSON config file
config_file_name = sys.argv[1]
with open(config_file_name, 'r') as f:
    config = json.load(f)

# make the run directory
run_dir_name = sys.argv[2]
os.mkdir(run_dir_name)

# copy the relevant parts of the PTRANSX repo into the new run directory
copy_ptransx(config['ptransx repo dir'], run_dir_name)

# copy the HARM3D snapshot file
try:
    harm3d_file_name = sys.argv[3]
    os.system('cp ' + harm3d_file_name + ' ./ ' + run_dir_name)
except:
    print("HARM3D snapshot file not found or provided.")

# move into the new run directory
os.chdir(run_dir_name)

# make a symlink from the indicated config file to "config.json"
os.symlink(sys.argv[1].split('/')[-1], 'config.json')

# build the custom python modules in ptxlib
build_fsolver(config,   rebuild = True)
build_response(config,  rebuild = True)
build_pandurata(config, rebuild = True)
build_xstarcomm(config, rebuild = True)

# make test script
make_test(config)
