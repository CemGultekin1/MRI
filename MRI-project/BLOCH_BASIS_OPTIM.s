#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH -t 8:00:00      
#SBATCH --mem=40GB
#SBATCH --output=job_bss_%a.out
#SBATCH --error=job_bss_%a.err


module unload matlab
module load matlab/R2018a
matlab -nodesktop -singleCompThread -r  "server_bloch_basis_optimization(${SLURM_ARRAY_TASK_ID},[9,13],40)" < /dev/null
