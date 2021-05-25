#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 7:00:00      
#SBATCH --mem=5GB
#SBATCH --output=job_crb_%a.out
#SBATCH --error=job_crb_%a.err


module unload matlab
module load matlab/R2018a
matlab -nodesktop -singleCompThread -r  "server_bloch_crb_compute(${SLURM_ARRAY_TASK_ID},[9,13])" < /dev/null
