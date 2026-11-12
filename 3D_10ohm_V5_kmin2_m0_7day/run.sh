#!/bin/bash
#SBATCH --account=nn9872k
#SBATCH --job-name=TDGP3D
#SBATCH --partition=large
#SBATCH --time=7-00:00:00
#SBATCH --nodes=9            
#SBATCH --ntasks-per-node=6
#SBATCH --cpus-per-task=42       # cores for each frequency
##SBATCH --mem-per-cpu=2G
#SBATCH --output=%x-%j.out

module --quiet reset
module load NRIS/CPU
module load intel/2023b
module list

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PROC_BIND=close
export OMP_PLACES=cores
export OMP_DYNAMIC=FALSE
export MKL_NUM_THREADS=1
export MKL_DYNAMIC=FALSE

cd "$SLURM_SUBMIT_DIR"
srun --ntasks=49 --cpus-per-task=$SLURM_CPUS_PER_TASK ~/CODES_TDGP3D_m0/tdgp3d_mpi


