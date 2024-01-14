#!/bin/bash -l

# job name
SBATCH -J myfft
# account
SBATCH -A edu19.SF2568
# email notification
SBATCH --mail-type=BEGIN,END
# 10 minutes wall-clock time will be given to this job
SBATCH -t 00:10:00
# set tasks per node to 24 in order to disablr hyperthreading
SBATCH --ntasks-per-node=24

module add i-compilers intelmpi

N = 8388608 # 2**23 : number of points of the array
processors=(1 2 4 8 16 32 64 128 256)
for P in "${processors[@]}"; do
    echo "$P processors"
    # 5 experiments for each number of processors
    for i in {0..5}; do
        echo "$i"
        time mpirun -np $P ./fft $N
    done
done
