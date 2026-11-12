#!/bin/bash
# Slurm job for the parallel-tempered 3D TDGP sampler (tdgp3d_mpi) on Olivia.
# Submit from the run folder (tdgp.settings, MT3D_V6.dat and the data file):  sbatch run_tdgp3d.sh
#
# Ranks:   the program needs exactly 1 + N_chains*nproc_per_team ranks (8 chains x 6 = 48, + 1 = 49);
#          srun below reads both from tdgp.settings and starts that many, so --ntasks must be at
#          least that and may be more (e.g. to fill whole nodes): extra task slots stay idle.
#          Each chain's nproc_per_team ranks share its frequencies' x and y polarizations
#          (nproc_per_team = 2 x number of periods: one polarization per rank, the fastest).
# Threads: --cpus-per-task = OpenMP threads, i.e. cores, for each polarization solve. They share
#          the nnx columns of MT3D_V6.dat, so use a divisor of nnx (nnx = 42: 7, 14, 21 or 42).
#          A 256-core node holds 6 tasks of 42 cores, so 49 tasks take 9 nodes; Slurm picks them.
#          (The PT coordinator, rank 0, needs only one core; its other 41 stay idle.)

#SBATCH --account=nn9872k
#SBATCH --job-name=TDGP3D
#SBATCH --time=7-00:00:00
#SBATCH --ntasks=49
#SBATCH --cpus-per-task=42
#SBATCH --hint=nomultithread
#SBATCH --mem-per-cpu=1G
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

# Ranks the program needs, from tdgp.settings (a setting left out has the default 8 or 6):
NCH=$(awk -F= 'tolower($1) ~ /^ *n_chains *$/ {print $2+0}' tdgp.settings)
NPT=$(awk -F= 'tolower($1) ~ /^ *nproc_per_team *$/ {print $2+0}' tdgp.settings)
NRANKS=$((1 + ${NCH:-8}*${NPT:-6}))
echo "tdgp.settings: N_chains = ${NCH:-8}, nproc_per_team = ${NPT:-6} -> $NRANKS ranks; the job has $SLURM_NTASKS tasks"
if [ "$NRANKS" -gt "$SLURM_NTASKS" ]; then
   echo "FATAL: raise --ntasks to at least $NRANKS (or lower N_chains / nproc_per_team in tdgp.settings)"
   exit 1
fi

# srun (not mpirun): gives every rank its --cpus-per-task cores. Newer Slurm versions do not pass
# --cpus-per-task from #SBATCH on to srun by themselves, hence the explicit option.
srun --ntasks=$NRANKS --cpus-per-task=$SLURM_CPUS_PER_TASK ~/CODES_TDGP3D/tdgp3d_mpi
