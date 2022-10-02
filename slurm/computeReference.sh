module load anaconda3/2020.07

for n in {0..99}
do
    python3 makeSbatch.py $n
    sbatch runBatch.sbatch
    rm runBatch.sbatch
done