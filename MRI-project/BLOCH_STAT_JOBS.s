#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5
#SBATCH -t 4-00:00:00      
#SBATCH --mem=50GB
#SBATCH --output=job_stat_%a.out
#SBATCH --error=job_stat_%a.err
#SBATCH --array=23-26,38-40

module purge
module load matlab/2020b
if %SLURM_ARRAY_TASK_ID%<38 (
matlab -nodesktop -singleCompThread -r  "server_bloch_multi_statistics(${SLURM_ARRAY_TASK_ID},5,'tag=m0s0p1')" < /dev/null
matlab -nodesktop -singleCompThread -r  "server_bloch_multi_statistics(${SLURM_ARRAY_TASK_ID},5,'')" < /dev/null
) else (
matlab -nodesktop -singleCompThread -r  "server_bloch_multi_statistics(${SLURM_ARRAY_TASK_ID},5,'')"< /dev/null
) 

