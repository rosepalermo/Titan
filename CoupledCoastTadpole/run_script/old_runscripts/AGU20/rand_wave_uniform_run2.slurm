#!/bin/bash  
#SBATCH --gres=gpu:1  
#SBATCH -N 1  
#SBATCH -n 16  
#SBATCH --time=03:00:00  
#SBATCH --exclusive  
#SBATCH -p sched_engaging_default
. /etc/profile.d/modules.sh
module load mit/matlab/2019a 
matlab -nodisplay -nodesktop -nosplash -r "run rand_wave_uniform_run2.m"