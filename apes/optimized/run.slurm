#!/bin/bash -l
#SBATCH -N 1 
#SBATCH -C knl
#SBATCH -q premium
#SBATCH -t 48:00:00

export OMP_NUM_THREADS=64
export OMP_PLACES=threads
export OMP_PROC_BIND=spread

echo running with 64 threads
./test.x

