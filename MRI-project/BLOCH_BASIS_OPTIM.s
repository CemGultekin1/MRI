#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=40
#SBATCH -t 8:00:00      
#SBATCH --mem=40GB
#SBATCH --output=job_bss_%a.out
#SBATCH --error=job_bss_%a.err
#SBATCH --array=25,38,39

module purge
module load matlab/2020b
if %SLURM_ARRAY_TASK_ID%<38 (
matlab -nodesktop -singleCompThread -r  "server_bloch_basis_optimization(${SLURM_ARRAY_TASK_ID},40,'m0s=0.1 tag=m0s0p1')" < /dev/null
) else (
matlab -nodesktop -singleCompThread -r  "server_bloch_basis_optimization(${SLURM_ARRAY_TASK_ID},40,'')"
) 
