import sys


jobnum = int(sys.argv[1])



# Generate the SBATCH script with the parameters
script = """#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=1GB
#SBATCH --job-name=mlis-{:d}
#SBATCH --output=mlis-{:d}.out

module purge
module load anaconda3/2020.07
module load matlab/2022b
cd /scratch/taa357/mlis/

matlab -nodisplay -nosplash -nodesktop -r "run('computeBatchSamples({:d})');exit;" | tail -n +11

""".format(jobnum, jobnum, jobnum)




with open("runBatch.sbatch", "w") as sbatch_file:
    sbatch_file.write(script)
sbatch_file.close()

