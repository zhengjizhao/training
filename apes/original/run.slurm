#!/bin/bash -l
#SBATCH -N 1 
#SBATCH -C knl
#SBATCH -q regular
#SBATCH -t 3:00:00

export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=64

echo running with $OMP_NUM_THREADS threads
./test.x 



