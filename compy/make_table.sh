#!/bin/bash

#SBATCH --job-name=compy
#SBATCH --time=24:00:00
#SBATCH --partition=parallel
#SBATCH --nodes=3
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
#SBATCH --mail-type=end
#SBATCH --mail-user=brooks.e.kinch@gmail.com

#### load and unload modules you may need
module purge
module load intel/18.0
module load intelmpi/2018
module load python/3.6

#### execute code and write output file to compy.log
mpirun -n 72 python3 p_mc_comp.py > compy.log
