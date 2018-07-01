#!/bin/bash -l
#SBATCH -N 1 
#SBATCH -C knl
#SBATCH -q premium
#SBATCH -t 24:00:00
#SBATCH --perf=vtune

export OMP_PLACES=threads
export OMP_PROC_BIND=spread

export OMP_NUM_THREADS=64

module load gcc/7.2.0
module load vtune
amplxe-cl -collect advanced-hotspots -knob analyze-openmp=true -r vtune_ah_par --search-dir=`pwd` --source-search-dir=`pwd` -- ./test.x

amplxe-cl -collect general-exploration -knob analyze-openmp=true -r vtune_ge_par --search-dir=`pwd` --source-search-dir=`pwd` -- ./test.x

amplxe-cl -collect memory-access -knob analyze-openmp=true -knob analyze-mem-objects=true -r vtune_ma_par --search-dir=`pwd` --source-search-dir=`pwd` -- ./test.x
