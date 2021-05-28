#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -t 4-00:00:00      
#SBATCH --mem=50GB
#SBATCH --output=job_stat_%a.out
#SBATCH --error=job_stat_%a.err


module purge
module load matlab/2020b
matlab -nodesktop -singleCompThread -r  "server_m0svalue_tests(${SLURM_ARRAY_TASK_ID},5)" < /dev/null
