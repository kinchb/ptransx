import os

files = os.listdir('./')
for file in files:
    if 'RADFLUX_panptx' not in file:
        files.remove(file)
files.sort()

new_files = []
for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    if (i % 100 == 0):
        print(i)
        new_files.append(file)
files = new_files

for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    print(i)
    os.mkdir('run_' + repr(i))
    os.chdir('run_' + repr(i))
    os.rename('../RADFLUX_panptx_' + repr(i) + '.h5', './RADFLUX_panptx_' + repr(i) + '.h5')
    os.system('cp -r /turquoise/users/kinch/panptx/* ./')
#   os.system('./build_all.sh')
    os.system('python panptx_full_run_generator.py RADFLUX_panptx_' + repr(i) + '.h5')
    os.chdir('../')

script = '#!/bin/bash' + '\n\n'

for file in files:
    i = int(file.split('_')[-1].split('.')[0])
    script += 'cd run_' + repr(i) + '\n'
    script += './submit_all_jobs.sh' + '\n'
    script += 'cd ..' + '\n'

file = open('submit_all_jobs.sh', 'w')
file.write(script)
file.close()
os.chmod('./submit_all_jobs.sh', 777)
