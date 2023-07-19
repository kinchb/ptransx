import os

target_dir_list =  ['h3d_3_a0_01',
                    'h3d_10_a0_01', 'h3d_10_a05_01', 'h3d_10_a09_01',
                    'h3d_10_a0_10', 'h3d_10_a05_10', 'h3d_10_a09_10',
                    'h3d_30_a0_01',
                    'h3d_10_a0_01_loFe', 'h3d_10_a0_01_hiFe']

target_base_dir = '/lustre/scratch4/turquoise/.mdt3/kinch'
# target_base_dir = '/mnt/archive'

"""
target_files = ['KDHARM0.001500.h5']
target_files = ['KDHARM0.RADFLUX.011000.h5']
target_files = ['mdot_data.npz', 'lum_data.npz', 'dLdT_collect_data.npz', 'dLdT_collect_cor_data.npz']
"""
target_files = ['run_1000/scat_diskt.h5', 'run_1000/scat_diskb.h5',
                'run_1000/spec_out_c_hr.h5', 'run_1000/spec_out_l_hr.h5', 'run_1000/scat_spec.1000.1.h5', 'run_1000/scat_spec.1000.5.h5', 'run_1000/scat_spec_c_hr.h5', 'run_1000/scat_spec_l_hr.h5']
"""
target_files = ['run_100/scat_spec.100.0.h5',
                'run_200/scat_spec.200.0.h5',
                'run_300/scat_spec.300.0.h5',
                'run_400/scat_spec.400.0.h5',
                'run_500/scat_spec.500.0.h5',
                'run_600/scat_spec.600.0.h5',
                'run_700/scat_spec.700.0.h5',
                'run_800/scat_spec.800.0.h5',
                'run_900/scat_spec.900.0.h5',
                'run_1000/scat_spec.1000.0.h5',
]
"""

local_base_dir = '/Users/kinch/panptx/panptx_data/survey'

for dir in target_dir_list:
    for target_file in target_files:
        line = 'scp wtrw:gr-fe:\'' + target_base_dir + '/' + dir + '/' + target_file + '\' ' + local_base_dir + '/' + dir + '/'
#       line = 'scp protagoras:\'' + target_base_dir + '/' + dir + '/' + target_file + '\' ' + local_base_dir + '/' + dir + '/'
        if dir not in os.listdir(local_base_dir):
            os.system('mkdir ' + local_base_dir + '/' + dir)
        os.system(line)
