#!/usr/bin/env python
filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """#!/bin/bash
#SBATCH -J 10_structure_optimization
#SBATCH -A m3337
#SBATCH -C cpu
#SBATCH -q regular
#SBATCH -t 13:00:00
#SBATCH --nodes 2
#SBATCH --ntasks-per-node=128
#SBATCH -c 2
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=gabriel.graf@duke.edu

module load PrgEnv-intel

source /opt/intel/oneapi/setvars.sh > /dev/null
ulimit -s unlimited

# Load the working intel tuning file
export I_MPI_TUNING_BIN=/opt/intel/oneapi/mpi/2021.6.0/etc/tuning_generic_shm-ofi_mlx_hcoll.dat

# Need library for srun to function
export I_MPI_PMI_LIBRARY=/usr/lib/shifter/mpich-2.2/dep/libpmi.so.0

srun -n ${SLURM_NTASKS} ~/fhi-aims.240920/build/aims.240920.scalapack.mpi.x > aims.out
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

