#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH -t 4-00:00:00      
#SBATCH --mem=50GB
#SBATCH --output=job_%a.out
#SBATCH --error=job_%a.err


module unload matlab
module load matlab/R2018a
matlab -nodesktop -singleCompThread -r  "server_bloch_multi(${SLURM_ARRAY_TASK_ID},40)" < /dev/null
