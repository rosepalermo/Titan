#!/bin/bash  
#SBATCH --gres=gpu:1  
#SBATCH -N 1  
#SBATCH -n 16  
#SBATCH --time=01:00:00  
#SBATCH --exclusive  
#SBATCH -p sched_engaging_default
. /etc/profile.d/modules.sh
module load mit/matlab/2019a 
matlab -nodisplay -nodesktop -nosplash -r "run river_wave_runs3.m"