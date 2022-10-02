#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=ce_%a
#SBATCH --output=ce_%a.out

module purge

cd /scratch/taa357/mlis/

module load matlab/2022b

matlab -batch "run('computeReference()');exit;"
