#!/usr/bin/env python
# Running this script will create a sample jobscript for SLURM 
# submission on Perlmutter in the working directory. Please make 
# sure all the SBATCH settings are as you want them.

filename = "jobscript.sh"

# Define the content of the job script
job_script_content = """#!/bin/bash
#SBATCH -J forgot_to_name
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
module load intel/2023.2.0
ulimit -s unlimited

# Load the working intel tuning file
export I_MPI_TUNING_BIN=/opt/intel/oneapi/mpi/2021.6.0/etc/tuning_generic_shm-ofi_mlx_hcoll.dat

srun -n ${SLURM_NTASKS} ~/FHIaims/build/aims.250131.scalapack.mpi.x > aims.out
"""

# Write the content to the file
with open(filename, "w") as file:
    file.write(job_script_content)

print(f"Job script written to {filename}")

