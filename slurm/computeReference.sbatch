#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=ce
#SBATCH --output=ce.out

module purge

cd /scratch/taa357/mlis/

module load matlab/2022b

matlab -nodisplay -nosplash -nodesktop -r "run('computeReference.m');exit;"
