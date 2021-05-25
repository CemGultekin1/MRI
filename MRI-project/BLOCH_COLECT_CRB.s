#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -t 5:00:00      
#SBATCH --mem=5GB
#SBATCH --output=job_colcrb.out
#SBATCH --error=job_colcrb.err

module unload matlab
module load matlab/R2018a
matlab -nodesktop -singleCompThread -r  "collect_crb" 
