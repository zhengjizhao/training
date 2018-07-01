#!/bin/bash -l
#SBATCH -N 1 
#SBATCH -C knl
#SBATCH -q premium
#SBATCH -t 24:00:00

export OMP_PLACES=threads
export OMP_PROC_BIND=spread
export OMP_NUM_THREADS=64

module load gcc/7.2.0
module load advisor/2018.integrated_roofline.up1

advixe-cl -collect survey --project-dir advisor_par --data-limit=0 -- ./test.x

advixe-cl -collect tripcounts -flop --project-dir advisor_par --data-limit=0 -- ./test.x

