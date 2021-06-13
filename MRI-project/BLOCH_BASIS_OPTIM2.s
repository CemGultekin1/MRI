#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH -t 8:00:00      
#SBATCH --mem=40GB
#SBATCH --output=job_bss2_%a.out
#SBATCH --error=job_bss2_%a.err


module purge
module load matlab/2020b
matlab -nodesktop -singleCompThread -r  "server_bloch_basis_optimization(${SLURM_ARRAY_TASK_ID},40,'')" < /dev/null